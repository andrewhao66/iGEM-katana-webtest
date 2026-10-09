"""
kg_audit.py — the AUDIT. This is Kagami's reason to exist: not "what parts are in
here" (annotation, already solved) but "is each part what it claims to be, and is
the construct actually buildable and clean" — with a reason and a fix per finding.

The checks reuse Katana's forward conventions verbatim so a Kagami PASS means the
same thing a katana-validate PASS does:
  - forbidden restriction / Type IIS site definitions  (forbid_sites_check.py)
  - >~40 bp host-homology recombination-substrate rule (host_homology_check.py)
  - junction rules: RBS→ATG spacing 5–9 nt, no spurious internal ATG/stop
"""
from kg_parse import revcomp
from kg_seedmatch import MAX_LOCI

# Mirrors kg_identify's MIN_REMAINDER: the shortest leftover that gets its own row.
_MIN_ROW = 6
import kg_refs

# --- enzyme sites, copied from parts-library/_tools/forbid_sites_check.py ---
TYPEIIS = {"BsmBI": "CGTCTC", "BsaI": "GGTCTC", "SapI": "GCTCTTC"}
SIXCUT = {"NcoI": "CCATGG", "EcoRI": "GAATTC", "XbaI": "TCTAGA",
          "SpeI": "ACTAGT", "PstI": "CTGCAG"}
# RFC[10] BioBrick-incompatible sites: presence breaks BioBrick assembly.
BIOBRICK_FORBIDDEN = ("EcoRI", "XbaI", "SpeI", "PstI")

STOPS = {"TAA", "TAG", "TGA"}

# Host off-target thresholds. All three are lengths of an EXACT match to the host
# chromosome, and each one marks a point where the question being asked changes.
#
#   <= 40            not a recombination substrate. PASS.
#   41 .. 499        a substrate, and whether that matters depends on recA. The recA-
#                    cloning strain downgrades it to a NOTE.
#   >= 500           gene scale. recA is no longer the question: a construct carrying
#                    half a kilobase of verbatim chromosome is either carrying an
#                    intended host-derived part or the wrong part, and an audit with no
#                    Design Spec in front of it cannot tell which. Both readings get
#                    reported, and the recA- downgrade stops applying -- a cloning
#                    strain makes a 60 bp homology inert, not a kilobase of chromosome
#                    correct.
#
# HOST_MATCH_CAP bounds the measurement, because the search is a seed-and-extend over
# a 4.6 Mb string and the advice does not change past gene scale. It is a reported
# FLOOR, not a measurement: past it the finding says "at least". The previous cap was
# 60 and was printed as though it were the length -- 3000 bp of verbatim chromosome
# read as "60 bp exact match", with "recode the stretch" as the advice.
HOST_SUBSTRATE_MIN = 40
HOST_GENE_SCALE = 500
HOST_MATCH_CAP = 1000

# Five tiers, and the distinction is the whole point of the 2026-09-15 rework:
#   FAIL  — a defect in the DNA itself (internal stop, empty). Blocks.
#   FLAG  — a real property of the DNA that is actionable in essentially any context
#           (mislabel, truncation, bad spacing, GC/homopolymer/repeat). Holds back PASS.
#   NOTE  — a real observation whose importance depends on the reader's assembly/host
#           context (a Type IIS site with no matching enzyme; host homology in a recA- strain).
#           Surfaced, never hidden, but does NOT hold back PASS.
#   SKIP  — a check that did not run (no host chosen). Shown as not-run, never counted as a
#           finding, because "we did not check" must never read as "we checked and it is fine".
#   PASS  — a positive result.
# The verdict word a reader sees is derived from these: any FAIL -> FAIL; else any FLAG ->
# "REVIEW — N to resolve"; else "PASS — N notes" (or a bare PASS). This is what lets a genuinely
# clean construct read PASS while its contextual notes stay visible.
PASS, FLAG, FAIL = "PASS", "FLAG", "FAIL"
NOTE, SKIP = "NOTE", "SKIP"


class Finding:
    def __init__(self, category, status, summary, loc="", detail="", fix=""):
        self.category = category
        self.status = status      # PASS / FLAG / FAIL
        self.summary = summary
        self.loc = loc
        self.detail = detail
        self.fix = fix


def _site_positions(seq, site):
    rc = revcomp(site)
    L = len(site)
    pos = []
    for i in range(len(seq) - L + 1):
        w = seq[i:i + L]
        if w == site or w == rc:
            pos.append(i + 1)
    return pos


def _gc(seq):
    if not seq:
        return 0.0
    g = sum(seq.count(b) for b in "GC")
    return 100.0 * g / len(seq)


def _longest_homopolymer(seq):
    best, run, prev = 1, 1, ""
    for c in seq:
        if c == prev:
            run += 1
            best = max(best, run)
        else:
            run, prev = 1, c
    return best


def _longest_direct_repeat(seq, kmin=20, cap=4000):
    """Cheap direct-repeat probe: does any k-mer occur twice? Binary search k.
    Capped for speed on large inputs (a coarse misassembly flag, not an aligner)."""
    s = seq[:cap] if len(seq) > cap else seq
    lo, hi, found = kmin, min(len(s) // 2, 200), 0
    while lo <= hi:
        k = (lo + hi) // 2
        seen = set()
        hit = False
        for i in range(len(s) - k + 1):
            km = s[i:i + k]
            if km in seen:
                hit = True
                break
            seen.add(km)
        if hit:
            found = k
            lo = k + 1
        else:
            hi = k - 1
    return found


def _next_atg(seq, from_idx, window=40):
    for i in range(from_idx, min(len(seq) - 2, from_idx + window)):
        if seq[i:i + 3] == "ATG":
            return i
    return -1


def _accuses_the_design(b):
    """Is a FLAG or FAIL about this block a statement about the DESIGN, or about a guess?

    A FLAG accuses. It is earned when somebody LABELLED this region -- then the finding is
    about what they said -- or when the identification is complete and clean enough that
    the block really is that part, whatever the file does or does not say about it.

    Otherwise the block is the identifier's own best guess at a stretch nobody named, and a
    finding about a partial guess is a fact about the guess. Those are still reported, as a
    NOTE, which never changes a PASS.

    Measured on a circular record, which is where this shows up with no ambiguity: the same
    molecule rotated to six different origins gave three different verdicts. The parts that
    were actually present -- sfGFP and B0015 -- were identified completely at every
    rotation, coverage 1.000, which is correct and rotation-invariant. What moved was which
    incidental composite the tiler picked up in the filler, each then carrying FLAGs and
    once an ORF FAIL: at one origin K1758101 was matched at 1049-1899 on a 1049-base record
    -- a block wrapping the whole molecule, 94% covered with eight indels -- and the ORF
    check read 851 bases across a junction nobody intended as a reading frame and found a
    stop in it. A plasmid has no canonical start, so a FAIL that depends on where someone
    chose to cut the file is noise with a severity attached.
    """
    if (b.claim_label or "").strip():
        return True
    cov = b.coverage
    return (cov is not None and cov >= 0.95
            and not (getattr(b, "indel", 0) or 0))


SD_CONSENSUS = "AGGAGG"


def _find_sd(block_seq):
    """Best Shine-Dalgarno core inside an RBS block: (start, end) 0-based, end exclusive.

    Katana seals an RBS as a functional unit, SD plus spacer, so the block's own end is NOT the SD's
    end and cannot be used to measure spacing. Score every 6-mer against AGGAGG; on a tie prefer the
    3'-most, since that is the one the ribosome reads. Returns None if nothing scores well enough to
    be called an SD, which is honest rather than guessing.
    """
    if len(block_seq) < len(SD_CONSENSUS):
        return None
    best, best_score = None, 0
    for i in range(len(block_seq) - len(SD_CONSENSUS) + 1):
        window = block_seq[i:i + len(SD_CONSENSUS)]
        score = sum(1 for a, b in zip(window, SD_CONSENSUS) if a == b)
        if score >= best_score:            # >= so later (3'-most) wins a tie
            best, best_score = i, score
    # 4 of 6 is the floor for calling something an SD. Below that it is just purine-ish noise.
    if best is None or best_score < 4:
        return None
    return best, best + len(SD_CONSENSUS)


def audit(record, blocks, vendor=None, fragment_bp_max=None, host_seq=None,
          host_reca=None, assembly=None, identify_status=None):
    """identify_status: the dict filled by kg_identify.identify(..., status=...), or None if the
    caller did not ask. When it says identification did NOT run, that is raised as a FLAG, which
    holds the verdict back to REVIEW.

    Why FLAG and not SKIP, since SKIP is what the host scan uses. SKIP is for a check the user
    OPTED OUT of: no host chosen means no host scan, and reporting that as a problem would be
    wrong. A missing BLAST+ is not an opt-out — the user asked for the audit and silently did not
    get its central half. It matches the Registry precedent instead: an unreachable Registry is
    already a FLAG ("could not reach" is not "checked and fine"), and identification is the more
    important of the two, so it cannot grade softer. Without this, a construct carrying a known
    mislabel read PASS — clean to order, exit 0.

    assembly: the assembly method context, one of {None, "BsaI", "BsmBI", "SapI",
    "BioBrick"}. It decides SEVERITY, not presence: a restriction/Type IIS site is only an
    actionable FLAG when the chosen method's enzyme would cut it; otherwise it is a NOTE.
    host_reca: recA status of the chosen host (True/False/None). A >40 bp host match is a FLAG
    in a recA+ (or unknown) host and only a NOTE in a recA- cloning strain, because that is the
    background in which it is actually a recombination substrate."""
    seq = record.seq
    findings = []

    # ---- POSITIVE INVARIANTS FIRST (Katana §4.1: never PASS by absence-of-bad) ----
    if not seq:
        findings.append(Finding("invariant", FAIL, "Empty sequence",
                                fix="Provide a non-empty FASTA/GenBank."))
        return findings

    # ---- DID THE IDENTIFICATION HALF ACTUALLY RUN? ----
    # As of 2026-10-07 it always does: identification is pure Python (kg_seedmatch) with no
    # external dependency, so the path that used to skip it -- NCBI BLAST+ absent, or present
    # and raising -- no longer exists. That path is what let a construct carrying a planted
    # mislabel report "PASS, clean to order" with exit 0.
    #
    # This branch is kept rather than deleted, and the distinction matters. It is now a guard
    # against a FUTURE identification path that cannot run, not a report on a missing binary:
    # a caller that hands us ran=False is telling us the audit is incomplete, and that must
    # never read as "checked and fine" -- the same rule the SKIP tier encodes. What it no longer
    # does is carry install instructions for a program this software does not need.
    _ident_ran = not (identify_status is not None and identify_status.get("ran") is False)
    if not _ident_ran:
        findings.append(Finding(
            "identification", FLAG,
            "Part identification did NOT run — this audit cannot catch a mislabel",
            detail=(identify_status.get("reason") or "the caller reported it did not run") +
                   ". Everything that does not depend on naming a part still ran (reading frame, "
                   "restriction sites, GC, repeats, size), and those results stand. But no block "
                   "was compared against a reference, so a part whose label disagrees with its "
                   "bases would NOT have been reported. Blocks shown as \"unidentified\" below "
                   "mean 'not checked', not 'checked and unmatched'.",
            fix="Re-run the audit. Identification needs nothing installed, so this should not "
                "happen; if it persists, report it."))
    # A requested --deep search that could not run is a check the user ASKED FOR and did
    # not get. Printing it to the terminal is not enough: the JSON and HTML reports and
    # the exit code would still read PASS, which is "we did not check" reading as
    # "checked and fine" -- the exact failure this project refuses, on the new opt-in
    # path. FLAG, matching the Registry precedent, not SKIP: SKIP is for a check the
    # user opted OUT of.
    _deep_failed = (identify_status or {}).get("deep_failed")
    if _deep_failed:
        findings.append(Finding(
            "deep-search", FLAG,
            "The deeper homology search you asked for did not run",
            detail=_deep_failed + " Everything the built-in identifier covers was still "
                                  "checked, and those results stand; what was not done "
                                  "is the gapped search for homologs below about 90% "
                                  "identity.",
            fix="Install NCBI BLAST+ if you need distant-homology search, or drop "
                "--deep: the audit is complete without it for parts that match a "
                "reference closely."))

    if not blocks:
        findings.append(Finding("invariant", FLAG,
                                "No blocks identified against the reference seed set",
                                detail=("Identification did not run, so nothing could match."
                                        if not _ident_ran else
                                        "Sequence parsed but nothing matched the public "
                                        "seed parts. Expand the reference set (fetch from "
                                        "the Registry) or supply an annotated GenBank."),
                                fix=("Re-run the audit; identification needs nothing installed."
                                     if not _ident_ran else
                                     "Add the relevant parts to the seed set via "
                                     "katana-parts-library intake, then re-run.")))

    # ---- THE DECOMPOSITION DOES NOT ACCOUNT FOR EVERY BASE ----
    # A leftover under MIN_REMAINDER gets no row of its own, which is right: a 2 bp scar
    # between two parts is not worth a line. But a decomposition that reads like a
    # complete accounting of the sequence has to be one, or say how much it is short --
    # and it used to say neither. On Katana's own build of pSense-Nit the table began at
    # base 14 and 13 bases appeared nowhere at all, not even as a number. A finding rather
    # than a header line, so the web page and --json carry it too.
    _short = (identify_status or {}).get("unaccounted_bp") or 0
    if _short:
        _acc = (identify_status or {}).get("accounted_bp") or 0
        findings.append(Finding(
            "decomposition", NOTE,
            f"{_short} base(s) are not in any block",
            detail=(f"The blocks above account for {_acc} of {len(record.seq)} bases. "
                    f"The remaining {_short} sit in fragments shorter than "
                    f"{_MIN_ROW} bp, between or beside the parts -- usually assembly "
                    f"scars or spare bases from a primer. They were looked at; they are "
                    f"too short to identify and too short to list as their own row."),
            fix=None))

    # ---- NOT EVERY OCCURRENCE WAS EXAMINED ----
    # The identifier extends at most MAX_LOCI diagonals per reference per strand. Beyond
    # that it used to drop the rest in silence: five tandem copies of a part got four
    # identity comparisons and the fifth got none, so a wrong label on it was never
    # checked. Raising the cap would not fix it -- any cap drops the next one -- so the
    # cap is reported instead. "We did not check" must never read as "checked and fine".
    # Only for references the construct ACTUALLY contains. With 18,538 references, some
    # short repetitive one seeds more than MAX_LOCI diagonals in almost any sequence --
    # measured on the demo: 355 references hit the cap and not one of them was identified
    # anywhere in the construct. Reporting those would put a FLAG on every clean audit,
    # which is the cry-wolf failure this finding exists to prevent.
    _identified = {b.ident_id for b in blocks if b.ident_id}
    _capped = {sid: k for sid, k in
               ((identify_status or {}).get("loci_capped") or {}).items()
               if sid in _identified}
    if _capped:
        _worst = sorted(_capped.items(), key=lambda kv: -kv[1])
        _names = ", ".join("%s (%d places)" % (sid, k) for sid, k in _worst[:4])
        findings.append(Finding(
            "identification-partial", FLAG,
            "Not every occurrence of a repeated part was examined",
            detail=(f"These references occur in more places than Katana examines: "
                    f"{_names}. It compares the best {MAX_LOCI} per reference per "
                    f"strand, so a wrong label on a further copy would not be caught "
                    f"here. This is a limit of the search, not a finding about your "
                    f"sequence."),
            fix="Check the repeated copies against each other by hand, or audit the "
                "repeated region on its own so each copy gets compared."))

    # ---- IDENTITY: claim vs sequence (the headline check) ----
    for b in blocks:
        if not b.claim_label:
            continue
        if not b.ident_id:
            # A block that HAS a label and no identity was dropped here in silence: no
            # FLAG, no SKIP, no finding of any kind. A construct with four identified
            # blocks and one labelled-but-unidentified block reported PASS without
            # mentioning it -- "a label we did not check" reading as nothing at all.
            #
            # The asymmetry was mine: the nested-claim loop below handles the identical
            # case with an identity-nested SKIP, and I wrote that one without noticing
            # this one was right above it. The existing identification FLAGs do not
            # cover it -- one fires only when the caller says identification never ran,
            # the other only when NOTHING matched anywhere.
            findings.append(Finding(
                "identity-unchecked", SKIP,
                f'The label "{b.claim_label.strip()}" was not checked',
                loc=f"{b.start}-{b.end}",
                detail=("Nothing in the reference set matches these bases well enough to "
                        "compare the label against, so whether the label is true is not "
                        "known. Not checked -- which is not the same as checked and "
                        "fine. A part not in the public set, a heavily diverged one, or "
                        "a region with an indel too large to bridge all look like this."),
                fix=("If this part is one of yours, pass your library with --library so "
                     "it can be compared. If it should be a public part, check it "
                     "against the Registry by hand.")))
            continue
        claim = b.claim_label.strip()
        claim_norm = claim.replace("BBa_", "").upper()
        ident = b.ident_id.upper()

        # The Registry holds the same bases under several part numbers. B0034's twelve
        # bases are also J34801, J70591, K1325011 and six more -- not similar sequences,
        # the SAME sequence. Which of them the identifier names is arbitrary, so a claim
        # naming any of them is correct, and must not be called a mislabel.
        #
        # Before this check consulted the alternatives, an honest B0034 label was told
        # 'Block labelled "B0034" is actually K1325011' -- an accusation of the one thing
        # this tool exists to detect, levelled at a construct that was right, with B0034
        # sitting in the block's own alternatives list at the time. A student who sees
        # that on every correct RBS learns that mislabel FLAGs are noise, and the one
        # real mislabel goes past them too.
        # ORIENTATION. The direction is part of the claim, and it is checked before the
        # name, because a part pointed the wrong way is wrong whatever it is called. An
        # RBS annotated on the reverse strand does not initiate translation of the CDS
        # after it; a promoter annotated backwards does not drive it.
        #
        # There was no orientation check here at all, so this went unreported for any
        # label. Accepting identical-sequence synonyms then made the worst case worse:
        # B0034 and K1045010 are reverse complements of each other, so forward B0034
        # bases annotated complement(1..12) were told 'Block labelled "B0034" is correct'.
        claimed_strand = getattr(b, "claim_strand", None)
        wrong_way = (claimed_strand is not None and b.strand is not None
                     and claimed_strand != b.strand)
        if wrong_way:
            said = "reverse" if claimed_strand == -1 else "forward"
            found = "reverse" if b.strand == -1 else "forward"
            findings.append(Finding(
                "identity-orientation", FLAG,
                f'Block labelled "{claim}" is annotated on the {said} strand, but its '
                f"bases read as {found}",
                loc=f"{b.start}-{b.end}",
                detail=f"The sequence matches {b.ident_id} on the {found} strand at "
                       f"{b.pident}% identity. The annotation claims {said}. One of the "
                       f"two is wrong, and which one matters: a part pointed the wrong "
                       f"way does not do its job -- an RBS on the reverse strand will not "
                       f"initiate translation of the CDS after it.",
                fix=f"Either correct the feature's orientation to {found}, or reverse-"
                    f"complement those bases if {said} is what you meant."))

        synonyms = {str(a).replace("BBa_", "").upper()
                    for a in (getattr(b, "alternatives", None) or [])}
        if claim_norm and claim_norm in synonyms and not wrong_way:
            others = sorted(s for s in synonyms if s != claim_norm)
            findings.append(Finding(
                "identity-synonym", NOTE,
                f'Block labelled "{claim}" is correct',
                loc=f"{b.start}-{b.end}",
                detail=f"These bases are registered under {len(others) + 1} part "
                       f"numbers, which are the same sequence, not similar ones: "
                       f"{', '.join([claim_norm] + others[:8])}"
                       f"{' and more' if len(others) > 8 else ''}. The identifier "
                       f"reported {b.ident_id}; naming any of them is right."))
            continue

        if claim_norm and claim_norm in synonyms:
            # Named correctly, pointed the wrong way. The orientation FLAG above already
            # says what is wrong; naming it a mislabel on top of that would be a second,
            # false accusation about the identity, which is the thing that IS right.
            continue

        # SUB-PART RE-DEPOSITS. The Registry holds fragments of parts as parts in their
        # own right, so the bases of one part are often also, exactly, some shorter
        # part's whole sequence. K1799015 is PyeaR[13:113]: 100 of PyeaR's 162 bases,
        # deposited separately.
        #
        # Both are hits, and _tile prefers completeness -- for a documented reason, so a
        # composite cannot swallow the real parts inside it -- so the 100%-complete
        # FRAGMENT won over the 96.9%-complete real part. The identifier then named the
        # fragment, and this check compared the claim against that name and accused it:
        #
        #   [FLAG] Block labelled "PyeaR" is actually K1799015
        #          fix: Re-label to K1799015, or swap in the real PyeaR sequence
        #
        # on a construct Katana had just built from its own sealed PyeaR. A false
        # accusation on the one check this software exists to perform, against the
        # engine's own output, telling a student to change a label that was right.
        #
        # So the claim is checked against the FACTS, not against the identifier's choice
        # of name: do these bases actually occur in the sequence the claim names? If they
        # do, the label is true whichever re-deposit got reported.
        #
        # Measured before settling on this: reporting every longer reference that
        # contains an identified one was tried and abandoned. The Registry is a dense web
        # of nested re-deposits and it fired dozens of times per construct -- the
        # cry-wolf failure this is meant to remove.
        # by_id() keys keep the Registry's own casing, and claim_norm is upper-cased, so
        # a direct lookup silently misses every part whose id is not already upper --
        # which is most of them, PyeaR included. Match case-insensitively, once.
        claim_ref_seq = ""
        if claim_norm:
            _byid = kg_refs.by_id()
            _cr = _byid.get(claim_norm)
            if _cr is None:
                for _k, _v in _byid.items():
                    if _k.replace("BBa_", "").upper() == claim_norm:
                        _cr = _v
                        break
            claim_ref_seq = kg_refs.normalise((_cr or {}).get("seq") or "")
        block_bases = kg_refs.normalise(record.seq[b.start - 1:b.end]) if b.start else ""
        if b.strand == -1 and block_bases:
            block_bases = revcomp(block_bases)
        if (claim_ref_seq and block_bases and claim_norm != ident
                and len(block_bases) >= 12 and block_bases in claim_ref_seq):
            at = claim_ref_seq.find(block_bases)
            # How much of the CLAIMED part is here. Accepting the containment and
            # stopping there traded a false accusation for a missed truncation: with only
            # PyeaR[13:113] present and labelled PyeaR, the report said 'Block labelled
            # "PyeaR" is correct', gave a verdict of PASS, and never mentioned that 62 of
            # PyeaR's 162 bases were absent. The label IS correct; the part is not all
            # there, and both of those have to be said.
            _claim_missing = len(claim_ref_seq) - len(block_bases)
            if _claim_missing > 0.05 * len(claim_ref_seq):
                findings.append(Finding(
                    "truncation", FLAG,
                    f'{claim}: {len(block_bases)} of its {len(claim_ref_seq)} bases are '
                    f'here, {_claim_missing} are not',
                    loc=f"{b.start}-{b.end}",
                    detail=(f"The label is right -- these bases are {claim}'s own, at its "
                            f"position {at + 1}-{at + len(block_bases)}. But {claim} is "
                            f"{len(claim_ref_seq)} bases long and only "
                            f"{len(block_bases)} of them are in this sequence. The "
                            f"identifier named {b.ident_id}, a shorter Registry entry "
                            f"that covers exactly the part that IS here, which is why "
                            f"the coverage beside it reads as complete."),
                    fix=(f"If you meant the whole of {claim}, the rest of it is missing. "
                         f"If you meant only this stretch, label it {b.ident_id} -- that "
                         f"is what it is.")))
            findings.append(Finding(
                "identity-subpart", NOTE,
                f'Block labelled "{claim}" names the right part — {b.ident_id} is a '
                f'shorter re-deposit of part of it',
                loc=f"{b.start}-{b.end}",
                detail=(f"These {len(block_bases)} bases are {claim}'s own, at its "
                        f"position {at + 1}-{at + len(block_bases)} of "
                        f"{len(claim_ref_seq)}. The Registry also holds them as "
                        f"{b.ident_id}, a separate shorter entry, and the identifier "
                        f"named that one because it matches end to end while {claim} "
                        f"here does not. Both readings are of the same bases."),
                fix=None))
            continue

        if claim_norm and claim_norm != ident and ident not in claim_norm:
            # Is it a strength-class change? (the B0032→B0034 case)
            claim_ref = kg_refs.by_id().get(claim_norm)
            extra = ""
            if claim_ref and claim_ref.get("variant") and b.ident_variant:
                cs = kg_refs.STRENGTH_ORDER.get(claim_ref["variant"].split("-")[-1])
                is_ = kg_refs.STRENGTH_ORDER.get(b.ident_variant.split("-")[-1])
                if cs and is_ and cs != is_:
                    direction = "stronger" if is_ > cs else "weaker"
                    extra = (f" Strength class differs "
                             f"({claim_ref['variant']}→{b.ident_variant}): expression "
                             f"will be markedly {direction} than the label implies.")
            findings.append(Finding(
                "identity-mislabel", FLAG,
                f'Block labelled "{claim}" is actually {b.ident_id} '
                f'({b.ident_name})',
                loc=f"{b.start}-{b.end}",
                detail=f"Sequence matches {b.ident_registry or b.ident_id} at "
                       f"{b.pident}% identity, {int(b.coverage*100)}% coverage." + extra,
                fix=f"Re-label to {b.ident_id}, or swap in the real {claim} sequence "
                    f"from the Registry via katana-parts-library."))

    # ---- NESTED CLAIMS: a claim inside an identified block is still a claim ----
    # These do not appear as blocks, because the blocks table is a decomposition and two
    # rows over the same bases makes it unreadable. They are checked here instead, with
    # the same rules the block claims get: a synonym is correct, a different sequence is
    # a mislabel, and an unidentifiable region is a SKIP rather than silence.
    for nc in ((identify_status or {}).get("nested_claims") or []):
        claim = str(nc.get("claim") or "").strip()
        if not claim:
            continue
        claim_norm = claim.replace("BBa_", "").upper()
        loc = f"{nc['start']}-{nc['end']}"
        ident = (nc.get("ident_id") or "").upper()
        syns = {str(a).replace("BBa_", "").upper()
                for a in (nc.get("alternatives") or [])}

        if not ident:
            findings.append(Finding(
                "identity-nested", SKIP,
                f'The region labelled "{claim}" inside a larger part was not identified',
                loc=loc,
                detail="Nothing in the reference set matches those bases well enough to "
                       "compare the label against. Not checked -- which is not the same "
                       "as checked and fine."))
            continue

        if claim_norm == ident or ident in claim_norm or claim_norm in syns:
            continue                      # named correctly, including by a synonym

        findings.append(Finding(
            "identity-mislabel", FLAG,
            f'Block labelled "{claim}" is actually {nc["ident_id"]}',
            loc=loc,
            detail=(f"This region sits inside a larger identified part, and its own "
                    f"label was checked separately. Its bases match {nc['ident_id']} at "
                    f"{nc['pident']}% identity over {int((nc['coverage'] or 0) * 100)}% "
                    f"of that reference."),
            fix=f"Re-label it {nc['ident_id']}, or put the real {claim} sequence there."))

    # ---- INDEL: a base inserted or deleted inside a part ----
    # The identifier joins the two pieces a shifted diagonal produces, so the part is
    # found and its label still gets checked. That must not make the indel INVISIBLE:
    # before the join the part vanished from the report entirely, and a silent repair is
    # the other half of the same mistake. One base out of a CDS shifts every codon after
    # it.
    for b in blocks:
        n_indel = getattr(b, "indel", 0) or 0
        if not b.ident_id or not n_indel:
            continue
        # The NET shift decides the frame, not the total bases involved. An insertion of
        # one base and a deletion of one base are two events and three bases of evidence
        # to go and look at, but the frame downstream is untouched. Two deletions are also
        # two events and the frame IS shifted. Using the total said "a multiple of three,
        # so a reading frame survives it" about a CDS carrying two deleted bases.
        n_net = abs(getattr(b, "indel_net", 0) or 0) or n_indel
        n_events = getattr(b, "indel_events", 0) or 1
        frame = n_net % 3 != 0
        role = (b.ident_role or b.claim_role or "").lower()
        coding = "cds" in role or "orf" in role
        findings.append(Finding(
            "indel", (FLAG if _accuses_the_design(b) else NOTE),
            (f"{b.ident_id} has {n_indel} base(s) inserted or deleted inside it"
             + (f", at {n_events} separate places" if n_events > 1 else "")),
            loc=f"{b.start}-{b.end}",
            detail=("The sequence matches this reference on either side of "
                    + (f"{n_events} shifts totalling {n_indel} base(s)"
                       if n_events > 1 else f"a shift of {n_indel} base(s)")
                    + ": the part is all there, but one stretch does not "
                    "line up with the rest. An insertion or deletion is the commonest "
                    "cloning and synthesis artifact."
                    + (f" The net shift is {n_net}, which is not a multiple of three, so "
                       "in a reading frame it SHIFTS every codon after it and the protein "
                       "downstream is wrong."
                       if frame else
                       f" The net shift is {n_net}, a multiple of three, so a reading "
                       "frame survives it -- but the residues there are not what the "
                       "reference encodes.")
                    + (" This block is coding, so that applies directly."
                       if coding else "")),
            fix=("Compare this region against the Registry entry base by base. If the "
                 "indel is deliberate, say so in the design; if it came from a primer or "
                 "a synthesis error, the part is not the one the label claims.")))

    # ---- FULL-LENGTH: truncated reference parts ----
    for b in blocks:
        # How many of the reference's bases are MISSING, which is the question the
        # threshold is about. Coverage alone counted a deleted base as covered, so a part
        # with 652 of 687 bases present -- a true 0.9491 -- came through at 0.9520 and
        # the FLAG never fired. A net deletion removes reference positions that have no
        # query base at all, so they are missing on top of whatever coverage did not
        # reach.
        _missing = None
        if b.ident_id and b.ref_len and b.matched_len is not None:
            _del = max(0, -(getattr(b, "indel_net", 0) or 0))
            _missing = max(0, b.ref_len - b.matched_len + _del)
        _short_frac = (_missing / float(b.ref_len)) if (_missing is not None
                                                        and b.ref_len) else None
        # A FLAG accuses the design. With no CLAIM on this block there is nothing to
        # accuse: nobody said this region was K1758101 -- the identifier offered it as
        # its best guess, and a guess that is partial is a fact about the guess, not a
        # truncation of anyone's part.
        #
        # Measured on a circular record, which is where it shows: the same molecule
        # rotated to five different origins gave five different sets of findings. The
        # parts that were actually there -- sfGFP and B0015 -- were identified completely
        # at every rotation, coverage 1.000, which is correct and rotation-invariant. What
        # moved was which incidental composite the tiler picked up in the filler
        # (K1758101 at one origin, K157005 at another, K112903 at a third), each then
        # FLAGGED for being incomplete. A circular plasmid has no canonical start, so a
        # FLAG that depends on where someone chose to cut the file is noise with a
        # severity attached.
        #
        # Unclaimed shortfalls are still reported, as a NOTE, which never changes a PASS.
        _claimed = _accuses_the_design(b)
        if b.ident_id and ((_short_frac is not None and _short_frac > 0.05)
                           or (_short_frac is None and b.coverage is not None
                               and b.coverage < 0.95)):
            _sev = FLAG if _claimed else NOTE
            # "Only N/M bp present" was a claim about the construct, and for a part whose
            # END has diverged it is a false one: those bases ARE present, they just
            # differ. Sequence alone cannot tell a substitution from a replacement -- they
            # are the same observation -- so this reports the shape of the disagreement
            # and names both readings instead of picking one.
            short = (b.ref_len - b.matched_len) if (b.ref_len and b.matched_len) else None
            present = getattr(b, "ref_in_query", None)
            if present and short:
                head = (f"{b.ident_id}: {b.matched_len} of its {b.ref_len} bases match "
                        f"here, {short} do not")
                detail = (f"The other {short} base(s) are present in your sequence but "
                          f"differ from the reference. Over the whole part that is "
                          f"{b.pident}% identity; over the matching stretch alone it is "
                          f"{getattr(b, 'core_pident', None)}%. Two things look like "
                          f"this: the part was truncated and something else sits in the "
                          f"gap, or the part is full length with a diverged end -- a "
                          f"cloning scar or a primer tail does exactly that. Which one "
                          f"it is cannot be told from the sequence.")
                fix = (f"Compare those {short} bases against the Registry entry for "
                       f"{b.ident_id}. If they are a scar you added deliberately, say so "
                       f"in the design; if not, the part is not the one you think it is.")
            else:
                head = (f"{b.ident_id} is incomplete here "
                        f"({int(b.coverage * 100)}% of the reference)")
                detail = (f"Only {b.matched_len} of {b.ref_len} bases of {b.ident_id} "
                          f"are in this sequence at all -- the rest runs past its end.")
                fix = ("Confirm the full-length part; a truncated "
                       "promoter/terminator/RBS often loses function.")
            if not _claimed:
                head = head + " (nothing here was labelled, so this is the "
                head = head + "identifier's own partial match)"
            findings.append(Finding("truncation", _sev, head,
                                    loc=f"{b.start}-{b.end}",
                                    detail=detail, fix=fix))

    # ---- ORF CLEAN: CDS-like blocks translate without internal stop ----
    for b in blocks:
        role = b.ident_role or b.claim_role
        if role not in ("cds", "reporter", "gene"):
            continue
        sub = record.sub(b.start, b.end)
        if b.strand == -1:
            sub = revcomp(sub)
        clean_len = len(sub) - (len(sub) % 3)
        codons = [sub[i:i+3] for i in range(0, clean_len, 3)]
        # A RUN of stops at the end is terminal, not internal. Two in a row (…AAA TGA TGA) is
        # deliberate practice, belt and braces against readthrough, and sfGFP ships that way.
        # Treating only the final codon as terminal reported the first of the pair as a premature
        # stop, and failed a construct that was simultaneously reported as 100% identical to the
        # reference - a contradiction that is the tell for this bug.
        last = len(codons)
        while last > 0 and codons[last - 1] in STOPS:
            last -= 1
        internal_stop = any(c in STOPS for c in codons[:last])
        starts_atg = sub[:3] in ("ATG", "GTG", "TTG")
        if internal_stop:
            findings.append(Finding(
                "orf", (FAIL if _accuses_the_design(b) else NOTE),
                f"CDS at {b.start}-{b.end} has an internal stop codon"
                + ("" if _accuses_the_design(b)
                   else " (nothing here was labelled -- this frame is the identifier's "
                        "own partial match, not a declared CDS)"),
                loc=f"{b.start}-{b.end}",
                detail="Premature stop in the reading frame — the protein is truncated.",
                fix="Fix the frame/sequence; re-check the upstream junction and any scar."))
        elif len(sub) % 3 != 0:
            findings.append(Finding(
                "orf", (FLAG if _accuses_the_design(b) else NOTE),
                f"CDS at {b.start}-{b.end} length is not a multiple of 3",
                loc=f"{b.start}-{b.end}",
                detail=f"{len(sub)} bp — frame is ambiguous; a fusion may be out of frame.",
                fix="Confirm the intended frame and junction spacing."))
        elif not starts_atg:
            findings.append(Finding(
                "orf", FLAG, f"CDS at {b.start}-{b.end} does not start with a start codon",
                loc=f"{b.start}-{b.end}",
                fix="Confirm the CDS boundary / start codon."))
        else:
            findings.append(Finding("orf", PASS,
                                    f"CDS at {b.start}-{b.end} ORF-clean",
                                    loc=f"{b.start}-{b.end}"))

    # ---- JUNCTIONS: RBS → ATG spacing 5–9 nt ----
    for i, b in enumerate(blocks):
        role = b.ident_role or b.claim_role
        if role != "rbs":
            continue
        # Look downstream of the RBS IN THE STRAND'S OWN DIRECTION. The search ran
        # along the forward sequence whatever the strand was, so for a reverse-strand
        # RBS it looked away from the CDS, which sits upstream in forward coordinates.
        # The same construct, reverse-complemented, changed its own verdict: forward
        # gave a spacing, reverse gave "No ATG found downstream of RBS" -- a false
        # finding about a junction that is correct, on a strand that is perfectly
        # ordinary in a real plasmid. The SD search a few lines below already
        # reverse-complements the block, so the two halves of the same check disagreed
        # about which way the gene runs.
        #
        # Spacing is a distance, so it is measured in the strand's own frame and the
        # result is the same number either way round.
        if b.strand == -1:
            # Walk the reverse complement of everything upstream of the block, which is
            # that strand's downstream. Offsets stay inside the reversed window; the
            # spacing arithmetic below uses them symmetrically.
            _up = revcomp(record.sub(max(1, b.start - 40), b.start - 1)) \
                if b.start > 1 else ""
            _rel = _next_atg(_up, 0, window=40)
            atg_rel = None if _rel < 0 else _rel
        else:
            _rel = _next_atg(seq, b.end, window=40)
            atg_rel = None if _rel < 0 else _rel
        if atg_rel is None:
            findings.append(Finding(
                "junction", FLAG, f"No ATG found downstream of RBS at {b.start}-{b.end}",
                loc=f"{b.start}-{b.end}",
                fix="Check the RBS is followed by a CDS start within ~30 nt."))
            continue
        # Measure from the SD core, not the block edge. A sealed RBS is a functional unit (SD plus
        # spacer), so its block ends flush against the ATG and edge-arithmetic always returns 0.
        block_seq = record.sub(b.start, b.end)
        if b.strand == -1:
            block_seq = revcomp(block_seq)
        sd = _find_sd(block_seq)
        if sd is None:
            findings.append(Finding(
                "junction", FLAG,
                f"No Shine-Dalgarno core found in the RBS at {b.start}-{b.end}",
                loc=f"{b.start}-{b.end}",
                detail="Spacing is measured from the SD, so it cannot be checked without one.",
                fix="Confirm this block really is a ribosome binding site."))
            continue
        # Both axes in the strand's own frame. On the forward strand the SD's last base
        # is at b.start + sd[1] - 1 and the ATG is at an absolute index; on the reverse
        # strand `block_seq` was already reverse-complemented for the SD search, so
        # sd[1] counts from the block's 3'-in-strand edge and atg_rel counts from the
        # block's upstream edge. The gap is the bases between them either way.
        if b.strand == -1:
            # distance from the SD's end to the block's strand-downstream edge, plus
            # however far into the upstream window the ATG sits
            spacing = (len(block_seq) - sd[1]) + atg_rel
        else:
            sd_abs_end = b.start + sd[1] - 1      # 1-based, inclusive, SD's last base
            spacing = (atg_rel + 1) - sd_abs_end - 1  # nt between the SD and the A of ATG
        if 5 <= spacing <= 9:
            findings.append(Finding("junction", PASS,
                                    f"RBS→ATG spacing {spacing} nt (in 5–9, measured from the SD)",
                                    loc=f"{b.start}-{b.end}"))
        else:
            findings.append(Finding(
                "junction", FLAG,
                f"RBS→ATG spacing {spacing} nt (outside 5–9, measured from the SD)",
                loc=f"{b.start}-{b.end}",
                detail="Spacing outside the 5–9 nt window shifts translation initiation "
                       "rate and can silence the downstream CDS.",
                fix="Adjust the spacer between the RBS and the start codon to 5–9 nt."))

    # ---- RESTRICTION SITES ----
    # Severity depends on the assembly context. A site is only an actionable FLAG when the chosen
    # method's enzyme would actually cut it; otherwise it is surfaced as a NOTE, so a SapI site does
    # not read as a problem for a construct that is simply being synthesised.
    biobrick_context = (assembly == "BioBrick")
    for name in BIOBRICK_FORBIDDEN:
        pos = _site_positions(seq, SIXCUT[name])
        if pos:
            loc = ", ".join(str(p) for p in pos[:6]) + ("…" if len(pos) > 6 else "")
            if biobrick_context:
                findings.append(Finding(
                    "restriction", FLAG,
                    f"{name} site present ({len(pos)}×) — breaks the BioBrick (RFC[10]) assembly "
                    f"you selected",
                    loc=loc,
                    detail=f"{name} = {SIXCUT[name]}. RFC[10] assembly and any digest using {name} "
                           f"will cut here.",
                    fix=f"Domesticate the {name} site (synonymous change in a CDS)."))
            else:
                findings.append(Finding(
                    "restriction", NOTE,
                    f"{name} site present ({len(pos)}×)",
                    loc=loc,
                    detail=f"{name} = {SIXCUT[name]}. Only matters if you submit this as a BioBrick "
                           f"(RFC[10]) part or digest with {name}; irrelevant for synthesis.",
                    fix=""))
    chosen_ts = assembly if assembly in TYPEIIS else None
    for name, site in TYPEIIS.items():
        n = len(_site_positions(seq, site))
        if not n:
            continue
        if name == chosen_ts:
            findings.append(Finding(
                "restriction", FLAG,
                f"{name} site present ({n}×) — cuts inside your Golden Gate / MoClo assembly",
                detail=f"{name} = {site}. You selected {name} assembly, so an internal {name} site "
                       f"will fragment the construct.",
                fix=f"Domesticate the internal {name} site(s) before assembly."))
        else:
            findings.append(Finding(
                "restriction", NOTE,
                f"Type IIS site present: {name}×{n}",
                detail=f"{name} = {site}. Only matters if you assemble with {name} "
                       f"(Golden Gate / MoClo); irrelevant for synthesis or a non-{name} method.",
                fix=""))

    # ---- GC / HOMOPOLYMER ----
    gc = _gc(seq)
    if gc < 25 or gc > 65:
        findings.append(Finding("composition", FLAG, f"Global GC {gc:.1f}% (outside 25–65%)",
                                fix="Extreme GC raises synthesis failure risk; consider "
                                    "codon/UTR adjustment."))
    else:
        findings.append(Finding("composition", PASS, f"Global GC {gc:.1f}%"))
    hp = _longest_homopolymer(seq)
    if hp >= 9:
        findings.append(Finding("composition", FLAG,
                                f"Homopolymer run of {hp} nt",
                                detail="Runs ≥9 nt are a synthesis + sequencing-slippage risk.",
                                fix="Break up the homopolymer (synonymous change if in a CDS)."))

    # ---- REPEATS / MISASSEMBLY ----
    dr = _longest_direct_repeat(seq)
    if dr >= 30:
        findings.append(Finding("repeats", FLAG,
                                f"Direct repeat ≥{dr} bp detected",
                                detail="Long exact repeats misassemble and cannot be PCR'd "
                                       "across reliably (cf. the bARGSer tandem gvpA region).",
                                fix="Excise on non-repeat flanks; avoid PCR across the repeat."))

    # ---- HOST OFF-TARGET (optional; >~40 bp exact = recombination substrate) ----
    if host_seq:
        longest = _longest_shared(seq, host_seq)
        # Past the cap the number is a floor, so it has to be said as one. "60 bp exact
        # match" for 3000 bases of verbatim chromosome was true of the measurement and
        # false about the construct.
        _howmuch = (f"at least {longest} bp" if longest >= HOST_MATCH_CAP
                    else f"{longest} bp")
        if longest >= HOST_GENE_SCALE:
            # recA governs whether a homology is a recombination substrate. At this
            # length that is no longer the interesting question, so the recA− strain
            # does not downgrade it: a cloning strain makes a short homology inert, not
            # a kilobase of verbatim chromosome correct.
            findings.append(Finding(
                "host-homology", FLAG,
                f"{_howmuch} of this sequence is verbatim host chromosome",
                detail="That is gene scale, not a chance homology, so the question is which "
                       "part is here rather than how it recombines. Either a host-derived "
                       "part is intended — in which case this is expected and worth "
                       "confirming — or the construct carries a stretch of the host genome "
                       "nobody asked for. This audit has no Design Spec to tell the two "
                       "apart, so it reports both readings rather than choosing one."
                       + ("" if host_reca is not False else
                          " The recA− strain you selected does not settle it: recA decides "
                          "whether a homology recombines, not whether the right part is here."),
                fix="Check this stretch against the Design Spec: confirm it is the "
                    "host-derived part you intended, at the coordinates you intended."))
        elif longest > HOST_SUBSTRATE_MIN:
            if host_reca is False:
                findings.append(Finding(
                    "host-homology", NOTE,
                    f"{_howmuch} exact match to the host chromosome "
                    f"(>{HOST_SUBSTRATE_MIN} bp)",
                    detail="A recombination substrate only in a recA+ background; inert in the "
                           "recA− cloning strain you selected.",
                    fix="Fine as-is in a recA− strain; recode the stretch before moving it into "
                        "a recA+ host."))
            else:
                findings.append(Finding(
                    "host-homology", FLAG,
                    f"{_howmuch} exact match to the host chromosome "
                    f"(>{HOST_SUBSTRATE_MIN} bp)",
                    detail="A recombination substrate in the recA+ host selected — the plasmid can "
                           "recombine into the chromosome across this stretch.",
                    fix="Recode/replace the host-identical stretch, or clone/propagate in a "
                        "recA− strain."))
        else:
            # Deliberately not "longest host match N bp". The search samples 12-base
            # seeds, so below 12 it cannot state the longest match and reports 0 — and
            # "longest host match 0 bp" was a precision the method does not have. What
            # the check establishes is the threshold, so that is what it says.
            findings.append(Finding(
                "host-homology", PASS,
                f"No exact match over {HOST_SUBSTRATE_MIN} bp to the host chromosome"))
    else:
        findings.append(Finding("host-homology", SKIP,
                                "Host off-target scan not run (no host selected)",
                                detail="Choose a host to run the >40 bp recombination-substrate "
                                       "check. This is a not-run check, not a problem with the "
                                       "sequence.",
                                fix=""))

    # ---- SIZE vs VENDOR ----
    if fragment_bp_max:
        if len(seq) > fragment_bp_max:
            findings.append(Finding("size", FLAG,
                                    f"{len(seq)} bp > {vendor or 'vendor'} cap "
                                    f"{fragment_bp_max} bp/fragment",
                                    fix="Split into fragments under the cap for synthesis."))
        else:
            findings.append(Finding("size", PASS,
                                    f"{len(seq)} bp ≤ {vendor or 'vendor'} cap {fragment_bp_max} bp"))

    return findings


def _longest_shared(seq, genome, cap=HOST_MATCH_CAP):
    """Longest exact contiguous match (either strand), measured up to `cap` bp.

    This used to binary-search the length while SAMPLING the start positions at
    `max(1, k // 2)` intervals. That predicate is not monotonic in k -- a longer k
    samples more coarsely, so a match present at one start could be stepped over -- and a
    binary search over a non-monotonic predicate returns an arbitrary answer.

    Measured: sfGFP against its own [1:42] slice as the host reported 28, with an exact
    41-base match sitting there. The host-homology check compares this against its
    threshold, so a recombination substrate over the actionable length came back as a
    PASS. A false all-clear on the one off-target check an audit can run without a
    genome.

    Seed-and-extend fixed that, but `cap` was 60 and the caller printed the returned
    number as the length. So the fix stopped at the point where the answer stopped being
    informative: measured against the bundled MG1655 chromosome, 3000 bp of verbatim host
    DNA reported 60, and pSense-Nit -- this branch's own ORACLE construct -- reported 60
    for a match that is 147. Nothing was missed; the number was simply wrong, and the
    advice attached to it ("recode the stretch") was advice for a 60 bp stretch.

    Two changes. `cap` is HOST_MATCH_CAP, far above every threshold that changes a
    verdict, and the caller says "at least" when it is reached. And the extension now
    walks the genome POSITION rather than asking `s[lo:hi] in genome` after each step:
    the old form re-searched 4.6 Mb per base added, which is why a low cap was load
    bearing. Extending at the locus the seed was found at is both exact and cheap, and
    every occurrence of every seed is tried, so nothing of length >= SEED can hide.
    """
    SEED = 12                 # every match of 12+ bases starts at some multiple of 12
    best = 0
    glen = len(genome)
    for s in (seq, revcomp(seq)):
        n = len(s)
        limit = min(n, cap)
        if limit <= 0:
            continue
        # Short sequences: just ask directly, which is exact and cheap.
        if n <= SEED * 4:
            for k in range(limit, best, -1):
                if any(s[i:i + k] in genome for i in range(0, n - k + 1)):
                    best = k
                    break
            continue
        for i in range(0, n - SEED + 1, SEED):
            seed = s[i:i + SEED]
            at = genome.find(seed)
            # Every occurrence, not just the first: the first one may sit in a short
            # island while a later one runs for a kilobase.
            while at != -1:
                lo, hi, glo, ghi = i, i + SEED, at, at + SEED
                while lo > 0 and glo > 0 and s[lo - 1] == genome[glo - 1]:
                    lo -= 1
                    glo -= 1
                while hi < n and ghi < glen and s[hi] == genome[ghi]:
                    hi += 1
                    ghi += 1
                if hi - lo > best:
                    best = min(hi - lo, limit)
                    if best >= limit:
                        return best
                at = genome.find(seed, at + 1)
    return best


def count(findings, status):
    return sum(1 for f in findings if f.status == status)


def verdict_kind(findings):
    """Canonical verdict token for logic, colour and exit codes: FAIL / REVIEW / PASS.

    NOTE and SKIP never make a construct anything other than a PASS here, and that is a
    DELIBERATE difference from core.result.BuildResult, where a SKIP does force REVIEW.
    The two are not drifting; they answer different questions.

    An audit's SKIP is usually unavoidable and usually the default: host-homology skips
    whenever no host is given, which is most runs. Escalating that would put REVIEW on
    nearly every audit, and a warning that fires every time is the cry-wolf failure, not
    a safeguard. A BUILD's SKIP means a gate the engine was supposed to run did not --
    rare, and worth holding back the verdict for.

    What makes this safe is that the audit's SKIPs stay VISIBLE: rendered as [----] and
    never a tick, sorted above PASS, and counted in the verdict line's tally as
    "N NOT CHECKED". Not escalating is not the same as not saying."""
    if any(f.status == FAIL for f in findings):
        return FAIL
    if any(f.status == FLAG for f in findings):
        return "REVIEW"
    # NOT every SKIP, but not none of them either. The blanket rule above was written for
    # host-homology, which skips whenever no host is given -- the default, on most runs --
    # and escalating that would put REVIEW on nearly every audit, which is the cry-wolf
    # failure rather than a safeguard.
    #
    # An unchecked IDENTITY is a different thing. It means the tool could not establish
    # what a labelled region actually is, which is the one question this software exists
    # to answer. Measured before this: an eight-base region labelled "B0015" produced
    # `identity-unchecked: SKIP` in the findings and a verdict of PASS, so the browser
    # showed a green pill over a label nobody had verified. Surfacing a row does not stop
    # a consumer treating the overall verdict as approval -- the verdict has to carry it.
    if any(f.status == SKIP and f.category.startswith("identity") for f in findings):
        return "REVIEW"
    return PASS


def verdict(findings):
    """Human verdict line. A construct with only contextual notes reads PASS — N notes, so a
    genuinely clean design is not dressed up as a problem, while the notes stay visible below."""
    kind = verdict_kind(findings)
    if kind == FAIL:
        return FAIL
    if kind == "REVIEW":
        n = count(findings, FLAG)
        if n == 0:
            # REVIEW with nothing to resolve means an identity could not be checked.
            # "REVIEW - 0 to resolve" would read like a bug; say what it is instead.
            _unchecked = [f for f in findings if f.status == SKIP
                          and f.category.startswith("identity")]
            return ("REVIEW — %d label(s) could not be checked" % len(_unchecked)
                    if _unchecked else "REVIEW")
        return f"REVIEW — {n} to resolve"
    n = count(findings, NOTE)
    if n:
        return f"PASS — {n} note" + ("s" if n != 1 else "")
    return PASS
