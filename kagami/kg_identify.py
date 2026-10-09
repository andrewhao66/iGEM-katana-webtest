"""
kg_identify.py — decompose a construct into blocks and IDENTIFY each block by
what its sequence ACTUALLY is (via blastn against the public reference seed set),
independent of what the construct's own annotation claims.

This is the "identify" half of Kagami. The audit (kg_audit) then compares the
identity to the claim — that comparison, not the identification, is the point.
"""
import os

# shutil, subprocess and tempfile are imported INSIDE the --deep helpers below, not here.
# They exist only on the blastn path, and this module is loaded by Pyodide in the browser
# front end, which has no process tools. Pyodide does ship a `subprocess` that exists and
# raises when used, so a module-level import happens to survive today -- but relying on
# that is relying on a stub's politeness, and tests/test_web_importable.py blocks the
# modules outright to make the stricter claim.

from kg_parse import revcomp
import kg_refs
import kg_seedmatch


class Block:
    def __init__(self, start, end, strand):
        self.start = start            # 1-based inclusive, on the given (+) strand
        self.end = end
        self.strand = strand          # +1/-1 : orientation the identified part sits in
        # claim (from the construct's own annotation) — may be None
        self.claim_label = None
        self.claim_role = None
        # The orientation the ANNOTATION claimed, which is not the same thing as
        # self.strand above: that one is where the sequence actually matched. Keeping
        # both is what lets the audit notice they disagree. Before this field existed
        # there was no orientation check anywhere, so a part whose label was right and
        # whose strand was wrong went unreported -- and once identical-sequence synonyms
        # were accepted, a forward part annotated complement was told it "is correct".
        self.claim_strand = None
        # identity (what blastn says the sequence really is) — may be None
        self.ident_id = None
        self.ident_name = None
        self.ident_role = None
        self.ident_variant = None
        self.ident_registry = None
        self.pident = None            # % identity over the reference's own span
        self.core_pident = None       # % identity over the best-supported segment only
        self.ref_in_query = None      # did the whole reference fit inside this sequence?
        self.wraps_origin = False     # on a plasmid, does this part cross the origin?
        self.indel = 0                # bases inserted or deleted inside the part
        self.indel_net = 0            # signed net shift: what the reading frame feels
        self.indel_events = 0         # how many separate shifts
        self.coverage = None          # matched_len / reference_len
        self.ref_len = None
        self.matched_len = None
        self.provenance = None
        # other references that fit these bases equally well (identical-sequence re-deposits)
        self.alternatives = []
        # for unidentified regions
        self.note = None

    @property
    def length(self):
        return self.end - self.start + 1

    def label(self):
        if self.ident_id:
            return self.ident_id
        if self.claim_label:
            return self.claim_label
        return "unidentified"


def _have_blast():
    """Is NCBI BLAST+ on PATH? Only the opt-in --deep path asks."""
    try:
        import shutil
    except ImportError:
        return False
    return bool(shutil.which("blastn") and shutil.which("makeblastdb"))


def _write_ref_db(workdir):
    import subprocess
    fa = os.path.join(workdir, "refs.fasta")
    with open(fa, "w", encoding="utf-8") as fh:
        for p in kg_refs.REFERENCE_PARTS:
            fh.write(f">{p['id']}\n{kg_refs.normalise(p['seq'])}\n")
    subprocess.run(
        ["makeblastdb", "-in", fa, "-dbtype", "nucl", "-out",
         os.path.join(workdir, "refs")],
        check=True, capture_output=True, text=True,
    )
    return os.path.join(workdir, "refs")


def _blast(construct_seq, dbpath, workdir):
    import subprocess
    q = os.path.join(workdir, "query.fasta")
    with open(q, "w", encoding="utf-8") as fh:
        fh.write(f">construct\n{construct_seq}\n")
    cols = "qseqid sseqid pident length qstart qend sstart send evalue bitscore"
    out = subprocess.run(
        ["blastn", "-task", "blastn-short", "-query", q, "-db", dbpath,
         "-word_size", "7", "-dust", "no", "-soft_masking", "false",
         "-evalue", "1000", "-perc_identity", "80", "-strand", "both",
         # Default is 500. With ~18k references a genuine 18 bp hit need not make the top 500,
         # so a real part silently became "unidentified" as the set grew. The alignment was fine;
         # the reporting cap was not.
         "-max_target_seqs", "5000",
         "-outfmt", f"6 {cols}"],
        check=True, capture_output=True, text=True,
    )
    hits = []
    reflen = {p["id"]: len(kg_refs.normalise(p["seq"])) for p in kg_refs.REFERENCE_PARTS}
    for line in out.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 10:
            continue
        sid = f[1]
        pid = float(f[2]); length = int(f[3])
        qs, qe = int(f[4]), int(f[5])
        ss, se = int(f[6]), int(f[7])
        bit = float(f[9])
        strand = 1 if se >= ss else -1
        qstart, qend = min(qs, qe), max(qs, qe)
        cov = length / reflen[sid] if reflen.get(sid) else 0.0
        hits.append(dict(sid=sid, pident=pid, length=length, qstart=qstart,
                         qend=qend, strand=strand, cov=cov, bit=bit))
    return hits


# A reference matched this much of its length counts as matched END TO END. Not 1.0, because a
# single base of vendor scar or a one-off trim should not demote a real identification.
COMPLETE_COV = 0.97


def _tile(hits, min_cov=0.6):
    """Greedy non-overlapping selection, preferring references matched end to end.

    Ranking by bitscore alone is what let one composite swallow three exact parts: bitscore grows
    with alignment length, so an 863 bp hit at 99.4% identity beat a 720 bp hit at 100% and then
    blocked it by overlap. The construct decomposed into 2 blocks instead of 4 and the RBS-spacing
    and ORF checks silently stopped running.

    Completeness comes first instead. A reference matched end to end is an identification; one
    matched at 81% is a FRAGMENT of something larger, which is what a composite looks like when the
    real parts are present too. Labels cannot be trusted to mark devices - thousands of authors over
    twenty years, and many entries carry no useful type at all - but this is structural and needs no
    label.

    It only bites when hits compete for the same bases. A genuinely truncated part with no rival
    still wins its region and is still reported as truncated.
    """
    kept = [h for h in hits if h["cov"] >= min_cov]
    # Completeness, then identity, then how much we trust the reference's origin, then bitscore.
    # The tier matters because identical sequences are common: fifteen references carry the exact
    # B0015 sequence, and without a deterministic tie-break the reported name changed run to run.
    kept.sort(key=lambda h: (h["cov"] >= COMPLETE_COV, h["pident"],
                             kg_refs.tier(h["sid"]), h["bit"]), reverse=True)
    chosen = []
    for h in kept:
        overlap = False
        for c in chosen:
            lo = max(h["qstart"], c["qstart"])
            hi = min(h["qend"], c["qend"])
            ov = max(0, hi - lo + 1)
            shorter = min(h["qend"] - h["qstart"] + 1, c["qend"] - c["qstart"] + 1)
            if ov > 0.5 * shorter:
                overlap = True
                break
        if not overlap:
            chosen.append(h)

    # Say when other references fit the same bases equally well. Naming one of fifteen identical
    # candidates with no hint that the rest exist is the tool sounding more certain than it is.
    for c in chosen:
        rivals = [h["sid"] for h in kept
                  if h["sid"] != c["sid"]
                  and abs(h["qstart"] - c["qstart"]) <= 2 and abs(h["qend"] - c["qend"]) <= 2
                  and h["pident"] >= c["pident"] - 0.01 and h["cov"] >= COMPLETE_COV]
        c["alternatives"] = sorted(set(rivals))

    chosen.sort(key=lambda h: h["qstart"])
    return chosen


def _find_orfs(seq, min_aa=50):
    """Return (start,end,strand) of the longest clean ORF in each strand region,
    to label unidentified stretches that look like a CDS."""
    stops = {"TAA", "TAG", "TGA"}
    found = []
    for strand, s in ((1, seq), (-1, revcomp(seq))):
        n = len(s)
        for frame in range(3):
            i = frame
            while i < n - 2:
                if s[i:i+3] == "ATG":
                    j = i + 3
                    while j < n - 2:
                        if s[j:j+3] in stops:
                            aa = (j - i) // 3
                            if aa >= min_aa:
                                if strand == 1:
                                    found.append((i + 1, j + 3, 1, aa))
                                else:
                                    found.append((n - (j + 3) + 1, n - i, -1, aa))
                            i = j
                            break
                        j += 3
                i += 3
    found.sort(key=lambda t: t[3], reverse=True)
    return found


# The no-blastn fallback that used to live here -- _exact_index() and
# _exact_hits_from_features(), which checked an annotated record's own claims by exact
# comparison -- was deleted on 2026-10-07. kg_seedmatch supersedes it: it searches
# unannotated sequence as well as it verifies annotated claims, and needs no binary. An
# unreachable ALTERNATIVE identification path is exactly the thing that gets rewired
# later and then disagrees with the live one, which is this project's cardinal sin.


def _nested_claim(feat, sub):
    """Identify `sub` on its own and return what the claim at those coordinates really is.

    Used for a claimed feature that sits wholly inside an identified block. Returning None
    means the bases could not be identified, which is not the same as the claim being
    wrong and is left to the audit to phrase.
    """
    if len(sub) < 8:
        return None
    try:
        hits = kg_seedmatch.identify_hits(sub, kg_refs.REFERENCE_PARTS)
    except Exception:
        return None
    if not hits:
        return dict(claim=feat.label or feat.kind, role=feat.kind,
                    start=feat.start, end=feat.end, strand=feat.strand,
                    ident_id=None, pident=None, coverage=None, alternatives=[])
    best = max(hits, key=lambda h: (h["cov"], h["pident"]))
    # Everything that fits these bases equally well, so a correct claim naming a
    # re-deposit of the same sequence is not called a mislabel.
    alts = sorted({h["sid"] for h in hits
                   if abs(h["pident"] - best["pident"]) < 0.05
                   and abs(h["cov"] - best["cov"]) < 0.02})
    return dict(claim=feat.label or feat.kind, role=feat.kind,
                start=feat.start, end=feat.end, strand=feat.strand,
                ident_id=best["sid"], pident=round(best["pident"], 1),
                coverage=round(best["cov"], 2), alternatives=alts)


def identify(record, workdir, status=None, deep=False):
    """Return an ordered list[Block] covering the construct, with claim + identity.

    status: an optional dict the caller passes in. Filled with {"ran": True, "reason": ""}.
    Identification now ALWAYS runs: it is pure Python with no external dependency, so
    the old "BLAST+ is missing" path -- which let a construct carrying a planted
    mislabel report PASS, clean to order, exit 0 -- cannot occur. When the caller asked
    for --deep and BLAST+ was unavailable or failed, status["deep_failed"] carries a
    sentence saying so; the built-in identifier ran regardless.

    deep: opt in to blastn's gapped local alignment for distant homologs. Not automatic:
    an ambient dependency makes two machines disagree about the same file, which is the
    works-on-my-machine failure the engine's own --expect-root flag exists to prevent.
    """
    seq = record.seq
    blocks = []

    if status is not None:
        status.clear()
        status.update({"ran": True, "reason": ""})

    def _deep_failed(reason):
        if status is not None:
            status["deep_failed"] = reason

    id_hits = []
    if len(seq) < 8:
        # Not a failure: there is nothing to identify. Left as ran=True so a 4 bp input does not
        # produce an alarming "identification did not run" on top of its real findings.
        pass
    else:
        # The built-in identifier ALWAYS runs. --deep adds blastn's gapped alignment on
        # top; it does not replace this. An earlier version ran only blastn when --deep
        # succeeded, which made a built-in hit that blastn did not report DISAPPEAR when
        # a user asked for a deeper search -- the opposite of what the flag promises.
        # The record's own topology, not an assumption. A part sitting across a
        # plasmid's origin was read as two truncated pieces before this was passed.
        raw = kg_seedmatch.identify_hits(
            seq, kg_refs.REFERENCE_PARTS,
            circular=(getattr(record, "topology", "linear") == "circular"))
        if status is not None:
            _capped = {h["sid"]: h["loci_capped"] for h in raw if h.get("loci_capped")}
            if _capped:
                status["loci_capped"] = _capped
        if deep:
            if not _have_blast():
                _deep_failed("--deep was requested but NCBI BLAST+ is not installed "
                             "(blastn/makeblastdb are not on PATH). The built-in "
                             "identifier ran, so the audit is complete for everything "
                             "except distant homologs.")
            else:
                try:
                    db = _write_ref_db(workdir)
                    raw = raw + _blast(seq, db, workdir)
                except Exception as exc:
                    _deep_failed("--deep was requested but BLAST+ failed to run "
                                 "(%s: %s). The built-in identifier ran, so the audit "
                                 "is complete for everything except distant homologs."
                                 % (exc.__class__.__name__, exc))
        id_hits = _tile(raw)

    refs = kg_refs.by_id()

    def make_ident_block(h):
        b = Block(h["qstart"], h["qend"], h["strand"])
        r = refs[h["sid"]]
        b.ident_id = r["id"]; b.ident_name = r["name"]; b.ident_role = r["role"]
        b.ident_variant = r.get("variant"); b.ident_registry = r.get("registry")
        b.pident = round(h["pident"], 1); b.coverage = round(h["cov"], 2)
        # Both identity measures, and whether the reference even fit, so the audit
        # can describe a shortfall instead of guessing at its cause.
        b.core_pident = h.get("core_pident")
        b.wraps_origin = bool(h.get("wraps_origin"))
        # How many bases the alignment had to shift to fit. Non-zero means an
        # insertion or a deletion inside the part, which is a frameshift risk in a
        # CDS and is not the same thing as the part being short.
        b.indel = h.get("indel") or 0
        # The net shift decides the frame; the count of events decides how many
        # places there are to look.
        b.indel_net = h.get("indel_net") or 0
        b.indel_events = h.get("indel_events") or 0
        b.ref_in_query = h.get("ref_in_query")
        b.alternatives = h.get("alternatives") or []
        b.ref_len = len(kg_refs.normalise(r["seq"])); b.matched_len = h["length"]
        b.provenance = r.get("provenance")
        return b

    # Attach the construct's own claims (features) to identity blocks by overlap.
    used_features = set()

    for h in id_hits:
        b = make_ident_block(h)
        best_f, best_ov = None, 0
        for fi, feat in enumerate(record.features):
            if feat.kind in ("source",):
                continue
            lo = max(b.start, feat.start); hi = min(b.end, feat.end)
            ov = max(0, hi - lo + 1)
            if ov > best_ov:
                best_ov, best_f = ov, fi
        if best_f is not None and best_ov > 0.3 * b.length:
            feat = record.features[best_f]
            b.claim_label = feat.label or feat.kind
            b.claim_role = feat.kind
            b.claim_strand = feat.strand
            used_features.add(best_f)
        blocks.append(b)

    # Features with no identity block (e.g. a CDS not in the seed set): keep as
    # claim-only blocks so they are still audited.
    for fi, feat in enumerate(record.features):
        if fi in used_features or feat.kind in ("source", "misc"):
            continue
        # A feature fully inside an identified block does not become a block of its own:
        # the blocks table is "what this construct is made of", and two rows covering the
        # same bases makes it unreadable. But its CLAIM is still a claim, and skipping it
        # here -- three lines under a comment promising such features are still audited --
        # meant it vanished. A file annotating B0015 correctly AND its inner eighty bases
        # as "B0034", a flat falsehood, reported no identity findings at all and a verdict
        # of PASS.
        #
        # Greedy decomposition must not decide WHICH of the submitter's claims get
        # verified. So the claim is identified on its own coordinates and handed to the
        # audit through status, which checks it exactly as it checks a block's.
        covered = any(b.start <= feat.start and b.end >= feat.end for b in blocks)
        if covered:
            if status is not None:
                sub = seq[feat.start - 1:feat.end]
                nested = _nested_claim(feat, sub)
                if nested is not None:
                    status.setdefault("nested_claims", []).append(nested)
            continue
        b = Block(feat.start, feat.end, feat.strand)
        b.claim_label = feat.label or feat.kind
        b.claim_role = feat.kind
        b.claim_strand = feat.strand
        blocks.append(b)

    blocks.sort(key=lambda b: b.start)

    # The smallest leftover worth showing as its own block. Below this it is a scar or a couple
    # of spare bases; at this size it is a real element (an RBS is ~18 bp) and the reader wants it.
    MIN_REMAINDER = 6

    # Fill unidentified gaps; label CDS-like ones via ORF scan.
    #
    # The threshold here was 30, four lines below the constant above whose comment argues
    # for 6 -- two thresholds disagreeing, and the one in force was longer than an RBS. So
    # any gap under 30 bases vanished from the decomposition entirely: not an identified
    # block, not an unidentified row, not a number. Measured on Katana's own build of
    # pSense-Nit, 1000 of 1013 bases were accounted for and the table simply began at 14.
    # The missing 13 were the start of the construct's own promoter, not a scar. On the
    # demo, 13 bases went the same way across two gaps of 7 and 6.
    #
    # MIN_REMAINDER is the threshold now, since it is the one whose reasoning is written
    # down. Anything smaller than that still does not vanish: identify() records the
    # shortfall in `status` and the report states it, because "we did not account for
    # these" must not read as nothing at all.
    covered = [(b.start, b.end) for b in blocks]
    gaps = []
    cursor = 1
    for s, e in sorted(covered):
        if s - cursor >= MIN_REMAINDER:
            gaps.append((cursor, s - 1))
        cursor = max(cursor, e + 1)
    if len(seq) - cursor + 1 >= MIN_REMAINDER:
        gaps.append((cursor, len(seq)))

    def _unidentified(a, b_):
        blk = Block(a, b_, 1)
        blk.note = "unidentified region (no reference match)"
        return blk

    for gs, ge in gaps:
        sub = seq[gs - 1:ge]
        orfs = _find_orfs(sub)
        if orfs and orfs[0][3] * 3 >= 0.7 * len(sub):
            os_, oe_, ostr, aa = orfs[0]
            # The CDS block is the ORF ITSELF, not the whole gap it was found in. Anything left
            # over on either side becomes its own unidentified block rather than being absorbed.
            # Absorbing it used to run the CDS past its own stop codon and report a sound
            # construct as truncated - see the S4 RBS case in the commit message.
            a_start, a_end = gs + os_ - 1, gs + oe_ - 1
            if a_start - gs >= MIN_REMAINDER:
                blocks.append(_unidentified(gs, a_start - 1))
            b = Block(a_start, a_end, ostr)
            b.ident_role = "cds"
            b.note = f"unannotated CDS-like ORF ({aa} aa, {'+' if ostr==1 else '-'} strand)"
            blocks.append(b)
            if ge - a_end >= MIN_REMAINDER:
                blocks.append(_unidentified(a_end + 1, ge))
        else:
            blocks.append(_unidentified(gs, ge))

    blocks.sort(key=lambda b: b.start)

    # What the decomposition does NOT account for. A leftover smaller than

    # MIN_REMAINDER gets no row of its own -- a 2 bp scar between two parts is not

    # worth a line -- but it must not vanish either, so the count goes out through

    # status and the report states it. A table that looks like a complete

    # accounting of the sequence has to be one, or say how much it is short.

    if status is not None:

        _acc = sum(b.end - b.start + 1 for b in blocks)

        if _acc < len(seq):

            status["unaccounted_bp"] = len(seq) - _acc

            status["accounted_bp"] = _acc

    return blocks
