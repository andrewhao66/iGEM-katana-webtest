/* index.js — load Katana's own audit modules into Pyodide and run them.
 *
 * The point of this file is what it does NOT do: it does not reimplement any check. It
 * fetches core/ and the audit half of kagami/ as plain .py files, writes them into
 * Pyodide's virtual filesystem, and calls them. Two implementations of one check drift
 * apart, and that drift is the failure the whole project exists to prevent -- so the
 * browser runs the same lines the terminal does.
 *
 * MODULES is pinned by tests/test_web_importable.py, which imports every one of them
 * with subprocess blocked. If a module is added there it must be added here, and the
 * deploy script copies exactly this set.
 */

const MODULES = [
  "core/__init__.py",
  "core/hashing.py",
  "core/lock.py",
  "core/parts.py",
  "core/result.py",
  "kagami/kg_parse.py",
  "kagami/kg_refs.py",
  "kagami/kg_seedmatch.py",
  "kagami/kg_identify.py",
  "kagami/kg_audit.py",
  "kagami/kg_bridge.py",
];

// The reference set, and the one bundled genome. Fetched on demand: the genome is only
// read when a host is chosen, so somebody checking a plasmid for restriction sites
// never pays for it.
const DATA = [
  "kagami/refs/reference_parts.tsv",
  "kagami/refs/reference_parts.fasta",
];
const GENOME = "kagami/genomes/MG1655_ecoli_NC_000913.3.fna";

const HOSTS = {
  "mg1655-pos": { file: GENOME, reca: true },
  "mg1655-neg": { file: GENOME, reca: false },
};

// Any chassis other than the bundled one. The page ships one genome because the file is
// 4.6 MB and shipping 23 of them would be most of the download; `get_genome.py` in the
// bundle fetches the rest. Without this the browser could only ever answer the
// off-target question about MG1655 -- and a team whose chassis is Nissle 1917 would have
// been reading an answer about a different organism, with nothing on the page saying so.
// The genome keeps the person's own filename inside the engine's filesystem. Two files
// go into an audit -- the sequence and the genome -- and an error only ever names one of
// them, so a fixed internal path meant a genome the reader could not parse produced
// "_local_genome: the LOCUS line declares 100 bp", which does not tell a student which
// of their two files is the problem.
const LOCAL_GENOME_DIR = "/genome";
let localGenomePath = "";
// recA is a property of the strain, not of the file, so a genome supplied this way has
// no recA status attached. null is what kg_audit reads as "assume recA+", the worse of
// the two, matching what the desktop GUI does for a browsed file.
const UPLOAD_RECA = null;
let localGenomeLoaded = false;

const VENDOR_CAP = { Twist: 5000, IDT: 3000, GenScript: 10000 };

const el = (id) => document.getElementById(id);
let pyodide = null;
let genomeLoaded = false;

function say(text) {
  el("status").textContent = text;
}

function progress(fraction) {
  el("bar").style.width = Math.round(fraction * 100) + "%";
}

async function fetchInto(path) {
  const res = await fetch(path);
  if (!res.ok) throw new Error(path + " -> HTTP " + res.status);
  return res.arrayBuffer();
}

async function boot() {
  try {
    say("Loading Python… (about 10 MB, once — your browser caches it)");
    progress(0.05);
    pyodide = await loadPyodide();

    progress(0.35);
    say("Loading Katana's own audit modules…");
    // Directory layout inside Pyodide mirrors the repository, so `import kg_parse` and
    // `from core import hashing` resolve exactly as they do on disk.
    pyodide.FS.mkdirTree("/katana/core");
    pyodide.FS.mkdirTree("/katana/kagami/refs");
    pyodide.FS.mkdirTree("/katana/kagami/genomes");
    for (const m of MODULES) {
      const buf = await fetchInto(m);
      pyodide.FS.writeFile("/katana/" + m, new Uint8Array(buf));
    }

    progress(0.55);
    say("Loading the reference set… (about 3 MB compressed)");
    for (const d of DATA) {
      const buf = await fetchInto(d);
      pyodide.FS.writeFile("/katana/" + d, new Uint8Array(buf));
    }

    progress(0.85);
    await pyodide.runPythonAsync(`
import sys
sys.path.insert(0, "/katana")
sys.path.insert(0, "/katana/kagami")
import kg_parse, kg_refs, kg_identify, kg_audit
`);
    const n = await pyodide.runPythonAsync("len(kg_refs.REFERENCE_PARTS)");
    progress(1);
    say("Ready. " + n.toLocaleString() + " reference parts loaded. "
        + "Nothing you check here is uploaded.");
    el("drop").classList.remove("hide");
  } catch (err) {
    progress(0);
    say("The engine could not load: " + err.message
        + "\n\nThis page needs to be served over http, not opened as a file:// URL — "
        + "browsers block a page from fetching its own files that way. "
        + "From the bundle, run:  ./katana web");
  }
}

async function ensureGenome() {
  if (genomeLoaded) return;
  say("Loading the host genome… (about 1.4 MB compressed, once)");
  const buf = await fetchInto(GENOME);
  pyodide.FS.writeFile("/katana/" + GENOME, new Uint8Array(buf));
  genomeLoaded = true;
}

// Reads the chosen genome into the in-browser filesystem. Nothing is sent anywhere —
// Pyodide's FS is a block of memory in this tab, and it is gone when the tab closes.
async function ensureLocalGenome() {
  const f = el("genomefile").files[0];
  if (!f) return false;
  if (localGenomeLoaded === f.name + ":" + f.size) return true;
  say("Reading " + f.name + " (" + Math.round(f.size / 1e6) + " MB)…");
  const bytes = new Uint8Array(await f.arrayBuffer());
  pyodide.FS.mkdirTree(LOCAL_GENOME_DIR);
  localGenomePath = LOCAL_GENOME_DIR + "/" + f.name.replace(/[^\w.\-]/g, "_");
  pyodide.FS.writeFile(localGenomePath, bytes);
  localGenomeLoaded = f.name + ":" + f.size;
  return true;
}

// The engine's kind -> the pill's class. PASS is the only thing that may look like a
// pass; everything else, including anything this page does not recognise, does not.
const PILL = { FAIL: "bad", REVIEW: "cond", PASS: "ok" };

const STAT = { PASS: ["p", "✓"], FLAG: ["w", "!"], FAIL: ["f", "✕"],
               NOTE: ["n", "○"], SKIP: ["s", "–"] };
// Unrecognised statuses sort FIRST, not last. An old page meeting a newer engine's tier
// must surface it, not bury it below the passes where nobody scrolls -- the same reason
// SKIP exists at all: "we did not check" must never read as "checked and fine".
const ORDER = { FAIL: 0, FLAG: 1, NOTE: 2, SKIP: 3, PASS: 4 };
const UNKNOWN_FIRST = -1;

function esc(s) {
  return String(s == null ? "" : s).replace(/[&<>"]/g,
    (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
}

function render(report) {
  el("rname").textContent = report.name + " · " + report.length.toLocaleString()
    + " bp · " + report.topology;

  // The colour comes from the engine's kind, never from the headline. verdict_kind() is
  // documented as "the canonical verdict token for logic, colour and exit codes"; the
  // headline beside it is a human sentence and will be reworded. Sniffing the sentence to
  // decide the colour is the project's own central failure -- a name and a fact that can
  // drift apart -- and it drifted the optimistic way: the ternary chain this replaces
  // fell through to "ok", so a reworded FAIL headline rendered GREEN.
  //
  // An unrecognised kind, from an engine newer than this page, renders as review: we do
  // not know what it is, so a person should look at it. It must never render as pass.
  const pill = el("pill");
  pill.className = "vpill " + (PILL[report.kind] || "cond");
  pill.textContent = report.verdict;

  let html = '<table class="blocks"><tr><th>where</th><th>role</th>'
    + "<th>what it actually is</th><th>claimed</th></tr>";
  for (const b of report.blocks) {
    const ident = b.identity || b.note || "unidentified";
    const conf = b.pident == null ? ""
      : " <span style='color:var(--faint)'>" + b.pident + "% id, "
        + Math.round((b.coverage || 0) * 100) + "% cov</span>";
    const alts = b.alternatives && b.alternatives.length
      ? " <span style='color:var(--faint)'>(= " + b.alternatives.length
        + " other ref" + (b.alternatives.length === 1 ? "" : "s") + ")</span>"
      : "";
    html += "<tr><td>" + b.start + "–" + b.end
      + (b.strand === -1 ? " −" : " +") + "</td><td>" + esc(b.role || "-")
      + "</td><td>" + esc(ident) + conf + alts + "</td><td>"
      + esc(b.claim || "") + "</td></tr>";
  }
  el("blocks").innerHTML = html + "</table>";

  const sorted = report.findings.slice().sort(
    (a, b) => (ORDER[a.status] ?? UNKNOWN_FIRST) - (ORDER[b.status] ?? UNKNOWN_FIRST));
  el("findings").innerHTML = sorted.map((f) => {
    // An unrecognised status gets the flag's glyph, not the note's circle: an unknown
    // tier needs attention, and a neutral circle is a claim that it does not.
    const [cls, gly] = STAT[f.status] || ["w", "?"];
    const note = f.detail ? '<span class="note">' + esc(f.detail) + "</span>" : "";
    // A SKIP's `fix` is the ACTIONABLE line -- "Choose a host to run the >40 bp
    // recombination-substrate check" -- and suppressing it meant the CLI showed the one
    // thing a person could do about a not-run check and the web page did not. A SKIP is
    // the tier that most needs its next step shown, not least.
    const fix = (f.fix && f.status !== "PASS")
      ? '<span class="fix">→ ' + esc(f.fix) + "</span>" : "";
    const loc = f.loc ? '<span class="loc">' + esc(f.loc) + "</span>"
      : '<span class="loc">—</span>';
    return '<div class="rrow"><span class="st ' + cls + '">' + gly + "</span>"
      + '<span class="what"><b>' + esc(f.summary) + "</b>" + note + fix + "</span>"
      + loc + "</div>";
  }).join("");

  el("result").classList.remove("hide");
}

async function audit(file) {
  if (!pyodide) { say("The engine is still loading — try again in a moment."); return; }
  el("result").classList.add("hide");

  const hostKey = el("host").value;
  if (hostKey === "local") {
    // A not-run check must never read as a clean one, so refusing to start is the only
    // honest response to "another host" with no file behind it.
    if (!(await ensureLocalGenome())) {
      say("Choose a genome file for the host, or set Host back to one of the presets.");
      return;
    }
  } else if (hostKey) {
    await ensureGenome();
  }

  say("Reading " + file.name + "…");
  progress(0.2);
  const bytes = new Uint8Array(await file.arrayBuffer());
  pyodide.FS.mkdirTree("/input");
  const inPath = "/input/" + file.name.replace(/[^\w.\-]/g, "_");
  pyodide.FS.writeFile(inPath, bytes);

  say("Identifying parts against " + (await pyodide.runPythonAsync(
    "len(kg_refs.REFERENCE_PARTS)")).toLocaleString() + " references…");
  progress(0.5);

  const host = HOSTS[hostKey];
  pyodide.globals.set("_in_path", inPath);
  if (hostKey === "local") {
    pyodide.globals.set("_host_file", localGenomePath);
    pyodide.globals.set("_host_reca", UPLOAD_RECA);
  } else {
    pyodide.globals.set("_host_file", host ? "/katana/" + host.file : "");
    pyodide.globals.set("_host_reca", host ? host.reca : null);
  }
  pyodide.globals.set("_assembly", el("assembly").value || null);
  pyodide.globals.set("_vendor", el("vendor").value || null);
  pyodide.globals.set("_cap", VENDOR_CAP[el("vendor").value] || null);

  let json;
  try {
    json = await pyodide.runPythonAsync(`
import json, tempfile
import kg_parse, kg_identify, kg_audit

record = kg_parse.parse(_in_path)

host_seq = None
if _host_file:
    # Two files reach an audit and an error only names one. Saying which is the host
    # genome costs a line here and saves a student reading a LOCUS complaint about a
    # chromosome as though it were about their plasmid. Not caught and turned into a
    # SKIP on purpose: they asked for this check against this file, and quietly not
    # running it is the failure this project is about.
    try:
        host_seq = kg_parse.parse(_host_file).seq
    except Exception as _e:
        raise RuntimeError(
            "The host genome file could not be read, so nothing was audited. "
            "The sequence you dropped was not the problem -- this is the genome: %s"
            % _e)

status = {}
with tempfile.TemporaryDirectory() as wd:
    blocks = kg_identify.identify(record, wd, status=status)

findings = kg_audit.audit(record, blocks, vendor=_vendor, fragment_bp_max=_cap,
                          host_seq=host_seq, host_reca=_host_reca,
                          assembly=_assembly, identify_status=status)

json.dumps(dict(
    name=record.name,
    length=len(record.seq),
    topology=record.topology,
    joined=list(getattr(record, "assembled_from", None) or []),
    verdict=kg_audit.verdict(findings),
    kind=kg_audit.verdict_kind(findings),
    identification_ran=bool(status.get("ran", True)),
    blocks=[dict(start=b.start, end=b.end, strand=b.strand, claim=b.claim_label,
                 identity=b.ident_id, role=(b.ident_role or b.claim_role),
                 pident=b.pident, coverage=b.coverage, note=b.note,
                 alternatives=list(getattr(b, "alternatives", None) or []))
            for b in blocks],
    findings=[dict(category=f.category, status=f.status, summary=f.summary,
                   loc=f.loc, detail=f.detail, fix=f.fix) for f in findings],
))
`);
  } catch (err) {
    progress(0);
    // A file that is not a sequence is not an engine failure, and must not read like
    // one. kg_parse raises NotASequenceFile with a message that already names the file
    // and the next thing to do, so pass it through whole.
    const m = String(err && err.message || err);
    const notSeq = m.match(/NotASequenceFile:\s*([\s\S]*?)(?:\n\n|$)/);
    say(notSeq ? "Katana could not read that file.\n\n" + notSeq[1].trim()
               : "The audit could not run: " + m);
    return;
  }

  progress(1);
  const report = JSON.parse(json);
  const notes = [];
  if (report.joined.length > 1) {
    notes.push(report.joined.length + " DNA fragments were found and joined in file "
      + "order (" + report.joined.join(" bp, ") + " bp). If that is not the order they "
      + "belong in, the audit below is of a construct you did not build.");
  }
  if (!report.identification_ran) {
    notes.push("Identification did not fully run — see the finding below.");
  }
  say("Done." + (notes.length ? "\n\n" + notes.join("\n\n") : ""));
  render(report);
}

// ---- wiring -----------------------------------------------------------------
const drop = el("drop");
const input = el("file");

drop.addEventListener("click", () => input.click());
drop.addEventListener("keydown", (e) => {
  if (e.key === "Enter" || e.key === " ") { e.preventDefault(); input.click(); }
});
input.addEventListener("change", () => {
  if (input.files && input.files[0]) audit(input.files[0]);
});
["dragenter", "dragover"].forEach((ev) =>
  drop.addEventListener(ev, (e) => { e.preventDefault(); drop.classList.add("over"); }));
["dragleave", "drop"].forEach((ev) =>
  drop.addEventListener(ev, (e) => { e.preventDefault(); drop.classList.remove("over"); }));
drop.addEventListener("drop", (e) => {
  const f = e.dataTransfer && e.dataTransfer.files && e.dataTransfer.files[0];
  if (f) audit(f);
});

// The genome row is only shown for "another host", so the presets stay a one-click
// choice and the file picker does not sit there implying the bundled genome needs one.
el("host").addEventListener("change", () => {
  el("genomerow").classList.toggle("hide", el("host").value !== "local");
});

boot();
