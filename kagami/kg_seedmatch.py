"""kg_seedmatch.py — pure-Python part identification.

Replaces the blastn shell-out in kg_identify. Measured against blastn on real
constructs: identical identity and coverage on full-length parts, on truncation at
60/80/90%, and on point mutations at 94.4/97.7/98.6% identity; roughly twice as fast,
because kg_identify rebuilt a BLAST database of 18,538 sequences on every single run.

Hybrid, because the reference length distribution forces it. Of 18,538 references, 894
are shorter than 25 bp and 123 are shorter than 12 bp, and short RBS parts are exactly
where the B0032/B0034 class of mislabel happens:

  references <= SHORT_MAX   direct exact match, C-speed via str.find, stricter than blastn
  references >  SHORT_MAX   k-mer seeding with stride filtering, then FULL-REFERENCE
                            alignment on the best diagonal

The full-reference alignment is load-bearing. Extending only to the last seed
under-reports coverage when a mutation falls near a part's end, which raised a FALSE
truncation FLAG on a full-length part (B0015 with 3 mutations read 89% coverage instead
of 100%). Aligning the whole reference against its query window fixes it exactly.

Output shape is identical to kg_identify._blast() so _tile() and everything downstream
are unchanged.
"""
import sys
from collections import defaultdict

try:
    import kg_refs
except ImportError:                      # running from another directory
    import os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import kg_refs

K = 12                  # seed length for the long path
SHORT_MAX = 25          # references this long or shorter take the exact path
SHORT_MIN = 12          # shorter than this carries no identifying information
MIN_HIT = 25            # shortest alignment the long path guarantees to find
MIN_IDENT = 80.0        # below this, not reported at all
MAX_LOCI = 4            # distinct diagonals extended per reference per strand

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def low_complexity(s):
    """True when a sequence carries no identifying information.

    Deliberately NOT "two or fewer distinct bases": B0034 is AAAGAGGAGAAA, which
    contains only A and G, and purine richness is the Shine-Dalgarno sequence's
    FUNCTION — it base-pairs with the pyrimidine-rich 16S rRNA anti-SD. A rule that
    rejects two-base sequences deletes the most important part class for the mislabel
    bug this project exists to catch.
    """
    if not s:
        return True
    if len(set(s)) == 1:
        return True                      # pure homopolymer: AAAAAAAAAA
    top = max(s.count(b) for b in set(s))
    return top / float(len(s)) > 0.8     # one base overwhelmingly dominant


def _usable(ref):
    """A reference that can identify something."""
    return len(ref) >= SHORT_MIN and not low_complexity(ref)


def _short_hits(sid, ref, query, qrc, n):
    """Exact occurrences of a short reference, on both strands.

    Exact rather than approximate on purpose. For a 12 bp RBS, 'these are the same
    bases' and 'these are nearly the same bases' are different claims, and the audit's
    job is the first one.
    """
    out = []
    L = len(ref)
    for strand, hay in ((1, query), (-1, qrc)):
        start = hay.find(ref)
        while start >= 0:
            if strand == 1:
                qs, qe = start + 1, start + L
            else:
                qs, qe = n - (start + L) + 1, n - start
            out.append(dict(sid=sid, pident=100.0, length=L, qstart=qs, qend=qe,
                            strand=strand, cov=1.0, bit=2.0 * L))
            start = hay.find(ref, start + 1)
    return out


def _query_index(q, k):
    """k-mer -> positions in q. Built once per query; the query is small."""
    idx = defaultdict(list)
    for i in range(len(q) - k + 1):
        idx[q[i:i + k]].append(i)
    return idx


def _long_hits(sid, ref, q, qidx, n, k, stride):
    """Find where a long reference occurs in the query. Returns a list of hits.

    Two stages. A stride prefilter rejects references with no seed at all, then every
    surviving reference is scanned fully and its seeds grouped by diagonal; each of the
    best MAX_LOCI diagonals is extended and scored.

    The stride prefilter's guarantee is narrower than it looks, and the limit is worth
    stating rather than glossing: an EXACT run of >= MIN_HIT bases contains
    (MIN_HIT - k + 1) consecutive k-mer starts, so sampling every (MIN_HIT - k + 1)//2
    positions must land inside one. That guarantee does NOT extend to approximate
    matches. A 40 bp reference at 92.5% identity with substitutions at positions 11, 24
    and 37 has a substitution inside every sampled window (0, 7, 14, 21, 28) while the
    exactly-matching windows at 12 and 25 are never sampled -- so it is rejected before
    the full scan can find it. Halving the stride improves sensitivity; it does not make
    high-identity detection certain. Use --deep for that.
    """
    out = []
    L = len(ref)
    get = qidx.get
    for strand, s in ((1, ref), (-1, revcomp(ref))):
        # cheap rejection: does any sampled seed hit at all?
        seeded = False
        for i in range(0, L - k + 1, stride):
            if get(s[i:i + k]):
                seeded = True
                break
        if not seeded:
            continue

        # full scan of this one reference, grouping seeds by diagonal
        diags = defaultdict(int)
        for i in range(L - k + 1):
            h = get(s[i:i + k])
            if h:
                for qp in h:
                    diags[qp - i] += 1
        if not diags:
            continue

        # Several diagonals, not just the best one. A part used TWICE in one construct
        # sits on two different diagonals, and architecture.order explicitly permits a
        # repeated id -- pAP-Logic's own version history is about a repeated RBS that
        # formed a 43 bp direct repeat. Taking only max(diags) reported the second
        # occurrence as an unidentified block, and a mislabel on it would then never be
        # checked, because the identity comparison skips blocks with no ident_id.
        #
        # Bounded at MAX_LOCI by seed count: unbounded would mean a full extension per
        # diagonal per reference across 18,538 references.
        #
        # The cap used to drop the rest IN SILENCE. Five tandem copies of B0015 returned
        # four hits and the fifth, at 517-645, was simply missing -- so a wrong label on
        # that copy received no identity comparison at all. Raising the cap would not fix
        # that; any cap drops the next one. What fixes it is saying so, which is this
        # project's own SKIP principle: "we did not check" must never read as "checked
        # and fine". loci_capped rides on the hits so the audit can report it.
        ranked = sorted(diags.items(), key=lambda kv: kv[1], reverse=True)
        capped = len(ranked) > MAX_LOCI
        _found = []
        for d, _seeds in ranked[:MAX_LOCI]:
            hit = _extend(sid, q, s, d, L, n, strand, allow_short=True)
            if hit is not None:
                if capped:
                    hit["loci_capped"] = len(ranked)
                _found.append(hit)

        # A provisional (sub-MIN_HIT) piece survives ONLY if it could actually merge:
        # there has to be a real piece of this same reference on a diagonal within
        # MAX_INDEL of it. Without this condition every diagonal of every reference
        # yielded a short piece, and the cost was measured rather than guessed -- 40
        # random 900 bp sequences went from 37 raw hits to 10,611, with one false
        # positive reaching coverage 0.800 on pure noise. That is the cry-wolf failure
        # this whole module is careful about, arriving through a fix for something else.
        _real = [h for h in _found if not h.get("provisional")]
        for h in _found:
            if not h.get("provisional"):
                out.append(h)
                continue
            _d = h.get("diag", 0)
            if any(0 < abs(_d - r.get("diag", 0)) <= MAX_INDEL for r in _real):
                out.append(h)
    return out


def _best_clean_segment(qseg, sseg, lo, hi):
    """The longest sub-segment of [lo,hi) whose identity clears MIN_IDENT.

    Kadane hands back the highest-SCORING segment, which for a part carrying two indels
    spans the shifted middle and so fails the identity floor -- taking the whole diagonal
    with it even though two clean stretches are plainly inside. This recovers the longer
    of them.

    Prefix sums over matches, then the longest window whose match ratio clears the floor.
    Returns (lo, hi, matches) or None.
    """
    m = [1 if a == b else 0 for a, b in zip(qseg[lo:hi], sseg[lo:hi])]
    n = len(m)
    if n == 0:
        return None
    pre = [0] * (n + 1)
    for i, v in enumerate(m):
        pre[i + 1] = pre[i] + v
    need = MIN_IDENT / 100.0
    best = None
    # Longest first: the first length that works anywhere is the answer, and lengths are
    # tried descending so this stops early on the common case.
    for ln in range(n, MIN_HIT - 1, -1):
        for a in range(0, n - ln + 1):
            mm = pre[a + ln] - pre[a]
            if mm >= need * ln:
                best = (lo + a, lo + a + ln, mm)
                break
        if best:
            break
    return best


def _extend(sid, q, s, d, L, n, strand, allow_short=False):
    """Score one diagonal and return a hit dict, or None. d = query pos - reference pos.

    The diagonal fixes the correspondence between reference and query positions, so the
    whole reference is laid against its query window -- and then the part of it that is
    actually SUPPORTED is found.

    Taking the whole window was wrong in a way that matters more than it looks. The span
    then measured how much CONSTRUCT was available rather than how much REFERENCE
    aligned, so a part present at 103 of its 129 bases followed by unrelated sequence
    reported cov=1.000 and never raised the truncation FLAG; at 77 of 129 the unrelated
    tail dragged identity under the 80% floor and the part vanished from the report
    altogether. A truncated part reported as full-length, or not reported at all, is
    this project's cardinal failure.

    Maximal-scoring segment (+1 match, -1 mismatch) -- ungapped local extension.
    Unrelated flanks score about -0.5 per base and are clipped, while a full-length part
    with scattered mutations scores positive throughout and keeps its whole span, so a
    mutated part is not mistaken for a truncated one. That would be the opposite error
    and just as wrong.
    """
    ws, we = d, d + L
    cs, ce = max(0, ws), min(n, we)
    if ce - cs < MIN_HIT:
        return None
    qseg = q[cs:ce]
    sseg = s[cs - ws:ce - ws]

    best_score = cur = 0
    best_lo = best_hi = cur_lo = 0
    for i, (a, b) in enumerate(zip(qseg, sseg)):
        cur += 1 if a == b else -1
        if cur <= 0:
            cur, cur_lo = 0, i + 1
        elif cur > best_score:
            best_score, best_lo, best_hi = cur, cur_lo, i + 1
    if best_score <= 0:
        return None

    span = best_hi - best_lo
    # A segment shorter than MIN_HIT is not an identification on its own -- 20 bases of
    # agreement is not evidence that a part is here. But it can be EVIDENCE FOR A MERGE,
    # and destroying it here meant any indel within 25 bases of either end of a part was
    # structurally unmergeable: B0015 with one deletion near each end lost both end
    # pieces, came back at coverage 0.752 with no indel at all, lost the tile to an
    # unrelated reference, and collected three confident false FLAGs -- mislabel, indel
    # and truncation -- on a correctly labelled part.
    #
    # So it survives marked `provisional`, usable only as a merge partner. _merge_indels
    # drops any that did not merge, so nothing short reaches a report by itself.
    _provisional = False
    if span < MIN_HIT:
        if not allow_short:
            return None
        _provisional = True
    matches = sum(1 for a, b in zip(qseg[best_lo:best_hi], sseg[best_lo:best_hi])
                  if a == b)
    core_pident = 100.0 * matches / span
    # Acceptance is judged on the clipped segment, so a genuinely truncated part still
    # clears the floor and gets reported with low coverage instead of disappearing --
    # which is what happened before Phase 1, at 77 of 129 bases.
    if core_pident < MIN_IDENT:
        # DO NOT throw the diagonal away. Kadane maximises SCORE, not identity, so when
        # two well-matching stretches sit either side of a shifted middle the maximal
        # segment spans the middle: for B0015 with an insertion and a deletion 40 bases
        # apart, a 68-point segment at ~76% identity beat either clean block's 40 or 48
        # points, failed this floor, and the whole diagonal was discarded -- so the part
        # fell to coverage 0.31, dropped out of the decomposition entirely, and the audit
        # filled the gap with three confident falsehoods about a correctly labelled part.
        #
        # The evidence a stricter reading would have kept is still here. Fall back to the
        # best segment that DOES clear the floor, found by sliding a window inward from
        # each end of the maximal one. That is the sub-segment the score hid, not a new
        # guess.
        _alt = _best_clean_segment(qseg, sseg, best_lo, best_hi)
        if _alt is None:
            return None
        best_lo, best_hi, matches = _alt
        span = best_hi - best_lo
        if span < MIN_HIT and not allow_short:
            return None
        _provisional = _provisional or span < MIN_HIT
        core_pident = 100.0 * matches / span

    # Identity over the reference's own span, as far as the query holds it. This is the
    # number a person reads, so it is the one that must not hide a diverged end.
    window = len(qseg)
    window_matches = sum(1 for a, b in zip(qseg, sseg) if a == b)
    pident = 100.0 * window_matches / window if window else 0.0
    # Did the WHOLE reference fit inside the query? It decides what a shortfall means.
    # Bases that are present but different is a diverged end; bases that run off the end
    # of the sequence are genuinely absent. Reporting "only 121/129 present" for the
    # first is a false statement about the construct.
    # rstart/rend are the slice of the REFERENCE this piece aligned, 0-based half-open.
    # Merging two pieces of one gapped alignment needs them: the pieces' query intervals
    # say where they sit, and only the reference intervals say how much of the part is
    # actually accounted for.
    return dict(sid=sid, pident=round(pident, 1), core_pident=round(core_pident, 1),
                length=span, qstart=cs + best_lo + 1, qend=cs + best_hi, strand=strand,
                cov=span / float(L), bit=2.0 * matches,
                ref_in_query=(ws >= 0 and we <= n),
                provisional=_provisional,
                diag=d, matches=matches, aligned=span,
                rstart=cs + best_lo - ws, rend=cs + best_hi - ws)


# The largest shift the merge will bridge. An indel from cloning or synthesis is one to a
# few bases; a few tens covers a sloppy primer or a small deletion. Two SEPARATE
# occurrences of a part differ in diagonal by the distance between them -- 129 for
# adjacent B0015 copies -- which is far above this, so they are never merged into one.
MAX_INDEL = 30


def _group_by_sid_strand(hits):
    """(sid, strand) -> its hits. Merging only ever joins pieces of the same reading."""
    groups = {}
    for h in hits:
        groups.setdefault((h["sid"], h["strand"]), []).append(h)
    return list(groups.items())


def _merge_indels(hits, reflen, qseq="", refseq=""):
    """Join hits of one reference that are pieces of a single GAPPED alignment.

    _extend scores one diagonal, and an insertion or deletion SHIFTS the diagonal, so the
    true alignment arrives here as two pieces -- each covering about half the reference.
    _tile then dropped both at its 0.6 coverage floor, and the part vanished from the
    report entirely.

    Measured before this existed, on real AmCyan bases labelled "sfGFP": intact, the
    mislabel FLAG fired at coverage 1.000; with ONE base deleted the block came back with
    no identity at all and the mislabel finding count went to zero. A block with no
    identity is never compared against its claim, so the one check this software exists to
    perform was silently skipped -- by the single most common cloning artifact there is.

    Two pieces are one alignment when they are the same reference and strand, their
    diagonals differ by no more than MAX_INDEL, their reference intervals barely overlap
    (they are different parts of the same part), and they sit next to each other in the
    query. The merged reading reports the union of the reference covered and carries the
    shift as `indel`, so the audit can say a base is missing rather than implying the part
    is short.
    """
    out = []
    by_group = {}
    for h in hits:
        by_group.setdefault((h["sid"], h["strand"]), []).append(h)

    for _key, group in by_group.items():
        if len(group) < 2:
            out.extend(group)
            continue
        group = sorted(group, key=lambda h: h["qstart"])
        used = [False] * len(group)
        for i, a in enumerate(group):
            if used[i]:
                continue
            merged = dict(a)
            used[i] = True
            for j in range(i + 1, len(group)):
                b = group[j]
                if used[j]:
                    continue
                # Against the PREVIOUS piece's diagonal, kept on the merged dict as
                # `diag`, and SIGNED. Measuring against the growing merged dict made n
                # single-base deletions report n(n+1)/2 -- two reported 3 -- and the
                # audit decides whether a reading frame survives with `indel % 3`, so a
                # CDS whose frame IS shifted was told "a reading frame survives it". A
                # false reassurance about a frameshift is worse than the miscount
                # underneath it.
                #
                # The sign matters for the same reason: d = qstart - refstart, so a
                # deletion lowers it and an insertion raises it. One of each cancels and
                # the frame really does survive; two deletions do not.
                delta = b.get("diag", 0) - merged.get("diag", 0)
                shift = abs(delta)
                if not 0 < shift <= MAX_INDEL:
                    continue
                # Reference intervals must be largely disjoint: two pieces of one part,
                # not the same stretch found twice on neighbouring diagonals.
                lo = max(merged["rstart"], b["rstart"])
                hi = min(merged["rend"], b["rend"])
                overlap = max(0, hi - lo)
                shorter = min(merged["rend"] - merged["rstart"],
                              b["rend"] - b["rstart"])
                if shorter <= 0 or overlap > 0.5 * shorter:
                    continue
                # And they must be adjacent in the query, allowing for the shift itself.
                gap = max(merged["qstart"], b["qstart"]) - min(merged["qend"], b["qend"])
                if gap > MAX_INDEL + 1:
                    continue

                rs = min(merged["rstart"], b["rstart"])
                re_ = max(merged["rend"], b["rend"])
                # TRIM AND RECOUNT at the junction, rather than estimate.
                #
                # Kadane extends each piece PAST the indel so both pieces cover the
                # junction, which means they overlap on BOTH axes and the two overlaps
                # DIFFER. For AmCyan with one base deleted: reference overlap 10, query
                # overlap 11, and the number of matches actually counted twice is 7.
                # Charging `overlap x density` -- the first attempt at this -- subtracted
                # 10 and landed 3 matches too low, so across 200 references with one base
                # deleted the error ran from -3.81 to +0.50 and 64 of them reported
                # EXACTLY 100.0% for a part with a base missing. The clamp removed the
                # ">" and left the "=", and "=" is what decides which reference wins a
                # tile.
                #
                # So the later piece is cut back to where the earlier one stops, on both
                # axes, and its contribution is counted against the real bases. One zip
                # over the trimmed piece; no heuristic, and the clamp below becomes an
                # assertion rather than a correction.
                # The cut has to satisfy BOTH axes, and they do not coincide: on this
                # piece's diagonal q0 = r0 + diag, and the two pieces' diagonals differ
                # by exactly the indel -- which is why the reference overlap (10) and the
                # query overlap (11) were different numbers in the measured case.
                # Trimming on the reference axis alone still let one query base be
                # counted by both pieces, so a part with a base deleted came back at
                # exactly 100.0%.
                #
                # Smallest r0 on this diagonal that is past where the earlier piece
                # stopped on each axis: r0 >= merged.rend, and r0 + diag >= merged.qend.
                _b_rs, _b_re = b["rstart"], b["rend"]
                _d_b = b.get("diag", 0)
                _cut_r = max(_b_rs, merged["rend"], merged["qend"] - _d_b)
                _cut_q = _cut_r + _d_b
                _tail_len = _b_re - _cut_r
                _tail_matches = 0
                if _tail_len > 0 and qseq and refseq:
                    _qs = qseq[_cut_q:_cut_q + _tail_len]
                    _rs_ = refseq[_cut_r:_cut_r + _tail_len]
                    _tail_matches = sum(1 for _x, _y in zip(_qs, _rs_) if _x == _y)
                    _tail_len = min(_tail_len, len(_qs), len(_rs_))
                elif _tail_len > 0:
                    # No sequences to recount against (a caller that passed none): fall
                    # back to the piece's own density, which is what this replaced.
                    _tail_matches = int(round(_tail_len * b["matches"]
                                              / float(max(1, b["aligned"]))))
                matches = merged["matches"] + _tail_matches
                aligned = merged["aligned"] + max(0, _tail_len)
                span = float(max(1, re_ - rs))
                # And clamp. The arithmetic above is an estimate over two local
                # alignments, not a global one, so it must not be able to express an
                # impossible answer even if a future change gets the estimate wrong.
                _pid = min(100.0, max(0.0, 100.0 * matches / span))
                _core = min(100.0, max(0.0, 100.0 * matches / float(max(1, aligned))))
                merged.update(
                    qstart=min(merged["qstart"], b["qstart"]),
                    qend=max(merged["qend"], b["qend"]),
                    rstart=rs, rend=re_,
                    matches=matches, aligned=aligned,
                    # Identity over the reference accounted for, with the shifted bases
                    # counted as the mismatches they are.
                    pident=round(_pid, 1),
                    core_pident=round(_core, 1),
                    # Coverage is how much of the reference was actually COMPARED, not
                    # the interval spanned. The span counted every base deleted inside it
                    # as covered, and the error always ran in the direction that
                    # suppresses the truncation FLAG: AmCyan with its last 34 bases
                    # trimmed AND one deleted has 652 of 687 present -- a true 0.9491,
                    # below the 0.95 threshold -- and reported 0.9520, so the student was
                    # not told 35 bases of their reporter were missing.
                    #
                    # `aligned` is the exact recount above: reference positions that had
                    # a query base against them. The remaining judgement -- how many
                    # bases the student is MISSING -- belongs where the threshold lives,
                    # in kg_audit, and it is made there from matched_len and indel_net.
                    cov=min(1.0, aligned / float(reflen)),
                    length=aligned, bit=2.0 * matches,
                    # indel: how many bases were inserted or deleted in total, which is
                    # what the headline says. indel_net: the signed sum, which is what
                    # decides whether a reading frame survives. indel_events: how many
                    # separate shifts, because "3 bases at one place" and "3 bases at
                    # three places" are different things to go and look at.
                    indel=merged.get("indel", 0) + shift,
                    indel_net=merged.get("indel_net", 0) + delta,
                    indel_events=merged.get("indel_events", 0) + 1,
                    # The merged piece's own diagonal becomes the LAST one absorbed, so
                    # the next comparison is against its neighbour and not against where
                    # this chain started.
                    diag=b.get("diag", merged.get("diag", 0)))
                used[j] = True
            # A provisional piece that found no partner is not a hit: MIN_HIT applies to
            # anything that reaches a report on its own. One that DID merge is part of a
            # full-length alignment now, so the mark comes off.
            if merged.get("provisional") and not merged.get("indel_events"):
                continue
            merged.pop("provisional", None)
            out.append(merged)
    return out


def identify_hits(query, refs, k=K, circular=False):
    """Identify which references occur in `query`. Returns kg_identify._blast()'s shape.

    circular: the query is a plasmid, so a part may sit ACROSS its origin.

    Where the origin falls is an arbitrary choice made by whoever exported the file --
    often the cloning site, which sits right next to the parts. Read linearly, a part
    spanning it is split in two: B0015 as a 129 bp circle rotated by 26 bases gave hits
    at coverage 0.798 and 0.202, and the audit raised a truncation FLAG saying 80% of the
    part was there when all of it was. A false warning on a correct plasmid is the
    cry-wolf failure -- once a student learns the warnings are noise, the real one goes
    past them too.

    So a circular query is searched with its own beginning appended, long enough for the
    longest reference to be found whole. Hits are then mapped back: one that lies wholly
    in the appended tail is a duplicate of one already found and is dropped, and one that
    crosses the seam is reported once, with an end coordinate past the sequence length --
    the same convention GenBank's join(104..129,1..103) expresses, and marked
    wraps_origin so nothing downstream has to infer it from arithmetic.
    """
    q = kg_refs.normalise(query)
    n = len(q)
    if n < SHORT_MIN:
        return []                        # nothing in the set can identify this

    # The overlap only needs to cover the longest reference that could straddle the
    # origin, and can never usefully exceed the plasmid itself.
    wrap = 0
    if circular:
        longest = max((len(kg_refs.normalise(r.get("seq") or "")) for r in refs),
                      default=0)
        wrap = min(n, max(0, longest - 1))
        if wrap:
            q = q + q[:wrap]
            n = len(q)
    qrc = revcomp(q)

    qidx = _query_index(q, k) if n >= k else {}
    stride = max(1, (MIN_HIT - k + 1) // 2)

    hits = []
    for r in refs:
        ref = kg_refs.normalise(r.get("seq") or "")
        if not _usable(ref):
            continue
        if len(ref) <= SHORT_MAX:
            hits.extend(_short_hits(r["id"], ref, q, qrc, n))
        elif qidx:
            hits.extend(_long_hits(r["id"], ref, q, qidx, n, k, stride))

    # One reference can seed on both strands at the same place (a palindrome, or a
    # self-complementary terminator). Keep the better reading rather than reporting
    # the same block twice.
    # Join pieces of one gapped alignment BEFORE anything judges coverage. Doing it after
    # the dedup, or after _tile's floor, is too late: the pieces are already gone.
    refseqs = {}
    for r in refs:
        refseqs[r["id"]] = kg_refs.normalise(r.get("seq") or "")
    merged = []
    for (sid, strand), grp in _group_by_sid_strand(hits):
        _ref = refseqs.get(sid) or ""
        # The reference as the hit's strand read it, and the query it was read against.
        # Both are needed to RECOUNT at a junction instead of estimating -- see
        # _merge_indels.
        merged.extend(_merge_indels(grp, len(_ref) or 1,
                                    q if strand == 1 else qrc,
                                    _ref if strand == 1 else revcomp(_ref)))
    hits = merged

    keep = {}
    for h in hits:
        key = (h["sid"], h["qstart"] // 10, h["qend"] // 10)
        cur = keep.get(key)
        if cur is None or (h["pident"], h["cov"]) > (cur["pident"], cur["cov"]):
            keep[key] = h
    out = list(keep.values())

    if circular and wrap:
        real = n - wrap                  # the plasmid's own length
        mapped = []
        for h in out:
            if h["qstart"] > real:
                continue                 # wholly inside the appended tail: a duplicate
            h = dict(h, wraps_origin=h["qend"] > real)
            mapped.append(h)
        # A part found whole across the seam is ALSO found as the two halves either side
        # of it, so it would be reported three times. Only those two halves are dropped.
        #
        # Keying this by sid alone -- which is what I wrote first -- discards every other
        # copy of that reference once any one copy reaches full coverage. Measured: a
        # plasmid with B0015 across the origin and a complete B0015 in the middle reported
        # one of them, and a plasmid with a whole copy across the origin and a TRUNCATED
        # copy in the middle reported only the whole one, so the truncated part vanished
        # from the report entirely. A truncated part not reported at all is this project's
        # cardinal failure, and the fix for one defect had introduced it.
        #
        # So the test is overlap, not identity of reference: a hit is superseded only if
        # it lies within the span of a whole wrapping hit of the same reference.
        def _pieces(h):
            """The hit's span in plasmid coordinates, as one or two intervals."""
            if h["qend"] <= real:
                return [(h["qstart"], h["qend"])]
            return [(h["qstart"], real), (1, h["qend"] - real)]

        def _covered_by(inner, outer):
            for a, b in _pieces(inner):
                if not any(c <= a and b <= d for c, d in _pieces(outer)):
                    return False
            return True

        wholes = [h for h in mapped
                  if h.get("wraps_origin") and h["cov"] > 0.99]
        out = [h for h in mapped
               if not any(w is not h and w["sid"] == h["sid"] and _covered_by(h, w)
                          for w in wholes)]
    return out
