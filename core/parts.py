"""core.parts — reading and writing sealed part files, and auditing a library.

One sequence reader, because there were four: katana_build.extract_gb_sequence,
katana_lock.parse_sequence, add_part.extract_gb_sequence (whose docstring said
"byte-for-byte the same logic as katana_build" -- a promise maintained by hand) and
add_part.read_fasta. They agreed, which is the point: four implementations of one
convention agree only until one of them is edited.

Standard library only. Loaded by Pyodide in the browser front end.
"""
import os
import re
from datetime import date

from . import hashing

_GENBANK = (".gb", ".gbk", ".genbank")


def extract_sequence(text, suffix=""):
    """The sequence a part file holds, uppercase, letters only.

    Dispatches on `suffix` when given, and otherwise on content: a GenBank ORIGIN block
    if there is one, else a FASTA body. Content-based fallback matters because a part
    file's extension is a naming convention, not a guarantee.
    """
    low = suffix.lower()
    if low in _GENBANK or (not low and "ORIGIN" in text):
        out, in_origin = [], False
        for line in text.splitlines():
            if line.startswith("ORIGIN"):
                in_origin = True
                continue
            if in_origin:
                if line.startswith("//"):
                    break
                out.append(re.sub(r"[^A-Za-z]", "", line))
        return "".join(out).upper()
    return "".join(re.sub(r"[^A-Za-z]", "", ln)
                   for ln in text.splitlines() if not ln.startswith(">")).upper()


def read_sequence(path):
    """The sequence a sealed part file on disk holds."""
    with open(path, "r", encoding="utf-8") as f:
        return extract_sequence(f.read(), os.path.splitext(path)[1])


def render_genbank(part_id, version, seq, source, klass):
    """A minimal, valid GenBank record for a part entering a library.

    GenBank rather than FASTA on purpose: the build engine reads part files through an
    ORIGIN-block reader, so a .fasta part would read as empty and be blocked.
    """
    today = date.today().strftime("%d-%b-%Y").upper()
    lines = [
        "LOCUS       %-24s%d bp    DNA     linear   SYN %s" % (part_id, len(seq), today),
        "DEFINITION  %s v%s, admitted to a Katana Parts Library." % (part_id, version),
        "ACCESSION   %s" % part_id,
        "VERSION     %s.%s" % (part_id, version),
        "KEYWORDS    .",
        "SOURCE      %s" % source,
        "COMMENT     Admitted by add_part.py. The manifest row in LOCK.tsv, not this file,",
        "            is the record of provenance; this file is the sequence it points at.",
        "            class: %s" % klass,
        "FEATURES             Location/Qualifiers",
        "     source          1..%d" % len(seq),
        '                     /note="%s"' % source,
        '                     /label="%s"' % part_id,
        "ORIGIN",
    ]
    low = seq.lower()
    for i in range(0, len(low), 60):
        chunk = low[i:i + 60]
        blocks = " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10))
        lines.append("%9d %s" % (i + 1, blocks))
    lines.append("//")
    return "\n".join(lines) + "\n"


def verify_library(lock_path):
    """Audit a whole library: every row's file, bytes, sequence, filename and row hash,
    plus orphan files and the root. Returns a list of problems; empty means sealed.

    Fail-closed: a file it cannot read is a problem, never a skip.
    """
    from . import lock as _lock

    d = os.path.dirname(os.path.abspath(lock_path))
    problems = []
    try:
        _header, rows = _lock.read(lock_path)
    except _lock.LockError as exc:
        return [str(exc)]

    seen = set()
    for row in rows:
        rid = "%s v%s" % (row.get("id", "?"), row.get("version", "?"))
        outfile = row.get("outfile", "")
        path = os.path.join(d, outfile.replace("\\", os.sep))
        if not os.path.exists(path):
            problems.append("%s: file MISSING (%s)" % (rid, outfile))
            continue
        seen.add(os.path.normpath(path))

        if hashing.file_sha256(path) != row.get("file_sha256"):
            problems.append("%s: file_sha256 MISMATCH (bytes changed)" % rid)
        try:
            if hashing.seq_sha256(read_sequence(path)) != row.get("seq_sha256"):
                problems.append("%s: seq_sha256 MISMATCH (sequence changed)" % rid)
        except Exception as exc:
            problems.append("%s: seq parse error: %s" % (rid, exc))

        fn12 = _lock.filename_sha12(outfile)
        if fn12 and fn12 != row.get("seq_sha256", "")[:12].lower():
            problems.append("%s: filename sha12 != seq_sha256" % rid)

        if hashing.row_sha256(row) != row.get("row_sha256"):
            problems.append("%s: row_sha256 MISMATCH (a trust field -- source, version, "
                            "class or outfile -- was edited)" % rid)

    # Orphans: a part file with no row is trust by dropping a file into the store.
    exts = (".gb", ".gbk", ".faa", ".fa", ".fasta")
    for root, dirs, files in os.walk(d):
        # Staging directories (names starting with "_", e.g. _incoming) hold raw
        # pre-seal fetches, which legitimately have no row and are NOT orphans.
        dirs[:] = [x for x in dirs if not x.startswith("_")]
        for fn in files:
            if fn.lower().endswith(exts):
                fp = os.path.normpath(os.path.join(root, fn))
                if fp not in seen:
                    problems.append("ORPHAN part file with no LOCK row: %s"
                                    % os.path.relpath(fp, d))

    ok, msg = _lock.verify_root(lock_path, os.path.join(d, "LOCK.root"))
    if not ok:
        problems.append(msg)

    return problems


def count_pending(lock_path):
    """Part files sitting in a staging directory, awaiting intake. Not orphans."""
    d = os.path.dirname(os.path.abspath(lock_path))
    exts = (".gb", ".gbk", ".faa", ".fa", ".fasta")
    pending = 0
    for root, dirs, _files in os.walk(d):
        for pd in [x for x in dirs if x.startswith("_")]:
            for _r, _ds, fs in os.walk(os.path.join(root, pd)):
                pending += sum(1 for f in fs if f.lower().endswith(exts))
        dirs[:] = [x for x in dirs if not x.startswith("_")]
    return pending
