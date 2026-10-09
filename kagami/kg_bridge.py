"""
kg_bridge.py — close the loop back to FORWARD Katana.

The question this answers: "can Kagami save the parts to the library, hash them,
and rebuild the construct into a working one?"

The tempting shortcut is: take each block's bytes straight from the pasted
construct, hash them, and seal them as library parts. Katana FORBIDS this. It is
the circular-provenance anti-pattern (the Katana spec's
law 3): "Every part has an independent, primary source that predates the construct.
Verifying a construct against the same bank it was built from is circular and
forbidden." Sealing a part from the very construct you are auditing would launder an
unverified sequence into a trusted one — exactly the failure Katana exists to kill.

So Kagami does NOT seal or hash parts into the library itself. Instead it emits two
handoff artifacts, and the REAL sealing/hashing/rebuild is done by the existing
forward skills, which already do it correctly:

  1. INTAKE REQUESTS  → katana-parts-library
       For each identified block, a request that names the PRIMARY SOURCE to fetch
       (Registry id / NCBI accession), never the construct bytes. The library's
       intake gate fetches fresh, verifies, seals, hashes (seq_sha256/file_sha256),
       and LOCKs — independently of this construct.

  2. A DRAFT DESIGN SPEC (reverse Stage 0) → katana-spec → katana-assemble
       INTENT reconstructed from the decomposition: ordered part ids + roles +
       architecture. Feed it to forward Katana and it regenerates a clean,
       deterministic .gb from the freshly-sealed parts — the "working one" — with
       the audit's fixes applied (correct label, domesticated sites).

  3. Then katana-diff (Stage 6) compares regenerate(draft spec) against the pasted
     input and shows exactly what was fixed.

That is the honest round-trip: Kagami finds and explains; forward Katana seals,
hashes, and rebuilds — no law bent.
"""


def _class_for(block):
    """reference if it resolved to a public reference part; else needs a primary
    source decision (designed/synthesised) — never sealed from the construct."""
    if block.ident_registry or block.ident_id:
        return "reference"
    return "unresolved"


def intake_requests(blocks):
    reqs = []
    for b in blocks:
        if b.ident_id:
            reqs.append(dict(
                id=b.ident_id,
                cls="reference",
                primary_source=b.ident_registry or b.ident_id,
                action="fetch from primary source → verify → seal → hash → LOCK "
                       "(katana-parts-library intake gate)",
                note=(b.provenance or "") +
                     "  |  DO NOT seal from this construct (circular-provenance).",
            ))
        else:
            role = b.ident_role or b.claim_role or "unknown"
            reqs.append(dict(
                id=f"UNRESOLVED_{b.start}_{b.end}",
                cls="unresolved",
                primary_source="NONE — no reference match",
                action=f"IDENTIFY the primary source for this {role} "
                       f"({b.length} bp) before it can enter the library. "
                       f"A block with no independent source CANNOT be sealed.",
                note=b.note or "",
            ))
    return reqs


def draft_spec(record, blocks, findings, vendor="Twist"):
    """Emit a reverse-engineered Katana Design Spec (YAML text). It holds INTENT
    (ordered part ids), never base pairs — same contract katana-spec produces."""
    lines = []
    lines.append("# DRAFT Katana Design Spec — reverse-engineered by Kagami.")
    lines.append("# INTENT only, NO base pairs. Every part must still earn its seal via")
    lines.append("# katana-parts-library intake (fetch primary source; never from this .gb).")
    lines.append("# Review, resolve any UNRESOLVED parts, then run forward Katana to rebuild.")
    lines.append(f"id:        {record.name}-recovered")
    lines.append("version:   1")
    lines.append("track:     acoustic          # community / WIST iGEM — public parts only")
    lines.append('purpose:   "recovered from an audited sequence; verify intent"')
    lines.append("host:      E_coli_MG1655")
    lines.append("assembly:  { method: null, decide_with: assembly-strategy-advisor }")
    lines.append(f"vendor:    {vendor}")
    lines.append("constraints:")
    lines.append("  forbid_sites: [EcoRI, XbaI, SpeI, PstI, BsaI, BsmBI, SapI]")
    lines.append("  host_context: E_coli_MG1655")
    lines.append("  output_gate:  acoustic-only")
    lines.append("")
    lines.append("parts:                        # ORDERED; pins added at intake (id@version@seq_sha12)")
    order = []
    for b in blocks:
        role = b.ident_role or b.claim_role or "misc"
        if b.ident_id:
            pid = b.ident_id
            src = f"{{ registry: iGEM, part: {b.ident_registry} }}" if b.ident_registry \
                else "{ db: NCBI, accession: TBD }"
            note = ""
            if b.claim_label and b.claim_label.replace("BBa_", "").upper() != b.ident_id.upper():
                note = f"   # NOTE: input labelled '{b.claim_label}' — corrected to {b.ident_id}"
            # One key per line. These used to be column-aligned onto a single line, which reads
            # nicely and is not YAML: `- id: X  role: Y  class: reference` is one scalar, not
            # three keys. Handing the emitted Spec to forward Katana died on a scanner error, so
            # the documented recovery path did not work. Alignment is not worth that.
            lines.append(f"  - id: {pid}")
            lines.append(f"    role: {role}")
            lines.append("    class: reference")
            lines.append("    pin: TBD@v1@TBD")
            lines.append(f"    source: {src}{note}")
            order.append(pid)
        else:
            pid = f"UNRESOLVED_{b.start}_{b.end}"
            lines.append(f"  - id: {pid}")
            lines.append(f"    role: {role}")
            lines.append("    class: unresolved")
            lines.append(f"    # {b.note or 'needs a primary source before it can be sealed'}")
            order.append(pid)
    lines.append("")
    lines.append("architecture:")
    lines.append(f"  order:    [{', '.join(order)}]")
    lines.append(f"  topology: {record.topology}")
    lines.append('  junction_rules: [ "RBS→ATG spacing 5-9 nt", "no spurious internal ATG/stop",')
    lines.append('                    "no forbidden RE site across any junction" ]')
    lines.append("")
    lines.append("records:   # PIPELINE-OWNED — stamped by forward Katana, not by Kagami")
    lines.append("  spec_hash: null")
    lines.append("  ref_parts: []")
    lines.append("  seal_certificate: null")
    return "\n".join(lines) + "\n"
