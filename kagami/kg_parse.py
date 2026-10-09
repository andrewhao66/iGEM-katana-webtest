"""
kg_parse.py — minimal, dependency-free FASTA / GenBank reader for Kagami.

Returns a Record: sequence (str, upper A/C/G/T/N), topology, and a list of
declared Features (the construct's *claims* about itself — label + role + span).
Kagami treats these claims as assertions to be checked, never as ground truth.
"""
import os
import re

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


class Feature:
    __slots__ = ("start", "end", "strand", "kind", "label")

    def __init__(self, start, end, strand, kind, label):
        self.start = start      # 1-based inclusive
        self.end = end          # 1-based inclusive
        self.strand = strand    # +1 / -1
        self.kind = kind        # feature key, e.g. CDS, promoter, RBS, terminator
        self.label = label      # human label / gene / product / note

    def span(self):
        return (self.start, self.end)


class Record:
    def __init__(self, name, seq, topology, features):
        self.name = name
        self.seq = seq
        self.topology = topology            # "circular" | "linear"
        self.features = features            # list[Feature]

    def sub(self, start, end):
        """1-based inclusive slice; wraps the origin when the record is circular.

        A circular molecule has no end, so a block that crosses the origin has
        `end > len(seq)` and a plain slice silently returns the TRUNCATED front half.
        Measured: sfGFP identified at 650-1369 on a 1049-base circular record came back
        400 bases long, so the ORF check reported "length is not a multiple of 3" for a
        complete, intact 720-base gene -- and the finding appeared or vanished depending
        on where the file happened to have been cut, which for a plasmid is an arbitrary
        choice. The same molecule rotated to five origins gave three different verdicts.

        A linear record keeps the old behaviour exactly: there is nothing past its end,
        and a request for bases beyond it is clamped rather than wrapped.
        """
        if self.topology == "circular" and self.seq and end > len(self.seq):
            n = len(self.seq)
            s = (start - 1) % n
            want = end - start + 1
            if want >= n:                 # a block longer than the molecule: one lap
                return (self.seq * (want // n + 2))[s:s + want]
            tail = self.seq[s:]
            return (tail + self.seq)[:want]
        return self.seq[start - 1:end]


def _clean(seq: str) -> str:
    return "".join(c for c in seq.upper() if c.isalpha())


# feature keys we map onto Kagami roles for downstream audit
_ROLE_KEYS = {
    "promoter": "promoter",
    "rbs": "rbs",
    "ribosome_binding_site": "rbs",
    "cds": "cds",
    "gene": "cds",
    "terminator": "terminator",
    "misc_feature": "misc",
    "protein_bind": "misc",
}

_LABEL_QUALS = ("label", "gene", "product", "standard_name", "note")



# ── loose input: csv / tsv / xlsx / txt / pasted ────────────────────────────

_DNA_OK = set("ACGTUNRYKMSWBDHV")


def _dna_fraction(tok):
    if not tok:
        return 0.0
    up = tok.upper()
    return sum(1 for c in up if c in _DNA_OK) / len(up)


def _dna_tokens(cells, min_len=8, min_frac=0.9):
    """Cells that are plausibly DNA, in document order.

    min_frac is deliberately below 1.0: exports carry stray spaces, line breaks and the odd digit
    from a position ruler. It is high enough that prose, part names and dates do not qualify.
    """
    out = []
    for c in cells:
        tok = "".join(ch for ch in str(c) if not ch.isspace() and not ch.isdigit())
        if len(tok) >= min_len and _dna_fraction(tok) >= min_frac:
            out.append("".join(ch for ch in tok.upper() if ch.isalpha()))
    return out


def _xlsx_cells(path):
    """Cell text from an .xlsx, in sheet order. stdlib only - an xlsx is a zip of XML."""
    import zipfile
    import xml.etree.ElementTree as ET

    NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    with zipfile.ZipFile(path) as z:
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            root = ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in root.findall(f"{NS}si"):
                shared.append("".join(node.text or "" for node in si.iter(f"{NS}t")))
        sheets = sorted(n for n in z.namelist()
                        if n.startswith("xl/worksheets/sheet") and n.endswith(".xml"))
        cells = []
        for name in sheets:
            root = ET.fromstring(z.read(name))
            for c in root.iter(f"{NS}c"):
                v = c.find(f"{NS}v")
                if c.get("t") == "s" and v is not None:
                    try:
                        cells.append(shared[int(v.text)])
                    except (ValueError, IndexError):
                        pass
                elif c.get("t") == "inlineStr":
                    cells.append("".join(n.text or "" for n in c.iter(f"{NS}t")))
                elif v is not None and v.text:
                    cells.append(v.text)
    return cells


class NotASequenceFile(Exception):
    """This file could not be read as DNA, and saying so is the only honest answer.

    Returning an empty Record instead is the fourth rule's failure in miniature: it
    resolves "I cannot read this" into "this is empty", silently, and destroys the
    information that there was a discrepancy. A report headed

        KAGAMI . sequence audit  .  photo  .  0 bp  . linear
          [FAIL] invariant    Empty sequence

    is a claim that this is a sequence, that it is 0 bp, and that it is linear. None of
    that was read from the file. And it sends a student to debug their sequence when what
    they need to debug is which file they picked -- which, for somebody dropping a photo
    or a Word document, is the whole of the problem.
    """


# Magic numbers of the files people actually drop by mistake. Only used to say WHICH kind
# of wrong file it was, which is the difference between a person fixing it and giving up.
_MAGIC = (
    (b"\x89PNG\r\n\x1a\n", "a PNG image"),
    (b"%PDF", "a PDF"),
    (b"\xff\xd8\xff", "a JPEG image"),
    (b"GIF8", "a GIF image"),
    (b"\x1f\x8b", "a gzip archive"),
    (b"BM", "a BMP image"),
    (b"\x00\x01\x00\x00", "a font"),
    (b"RIFF", "a media file"),
    (b"\xd0\xcf\x11\xe0", "an old Office document (.doc/.xls)"),
)

_ZIP_KINDS = {".docx": "a Word document", ".pptx": "a PowerPoint file",
              ".odt": "an OpenDocument file", ".zip": "a zip archive",
              ".jar": "a Java archive", ".epub": "an EPUB book"}


def _diagnose(path):
    """Why this file holds no DNA, in words that name the next thing to do."""
    try:
        size = os.path.getsize(path)
    except OSError:
        size = -1
    if size == 0:
        return ("%s is empty -- it has no bytes in it at all. Check you saved the file, "
                "and that you are pointing at the one you meant."
                % os.path.basename(str(path)))

    try:
        with open(path, "rb") as f:
            head = f.read(4096)
    except OSError:
        head = b""

    what = None
    for magic, name in _MAGIC:
        if head.startswith(magic):
            what = name
            break
    if what is None and head.startswith(b"PK\x03\x04"):
        # .xlsx and .xlsm are zips too, and they are genuinely supported -- so a zip is
        # only the wrong kind of file once the spreadsheet reader has found nothing.
        ext = os.path.splitext(str(path).lower())[1]
        what = _ZIP_KINDS.get(ext, "a zip-based document (.docx, .pptx or similar)")
    if what is None and b"\x00" in head:
        what = "a binary file of some kind"

    base = os.path.basename(str(path))
    if what:
        return ("%s is %s, not a sequence file. Katana reads FASTA (.fasta, .fa), "
                "GenBank (.gb, .gbk), a spreadsheet (.xlsx, .csv), or a plain text file "
                "with the bases in it. If your sequence is inside this document, export "
                "or copy it out first." % (base, what))
    return ("No DNA was found in %s. Katana reads FASTA (.fasta, .fa), GenBank "
            "(.gb, .gbk), a spreadsheet (.xlsx, .csv), or a plain text file with the "
            "bases in it -- A, C, G and T. If the bases are in there, check they are not "
            "split across a column the reader did not look in." % base)


def _parse_loose(path, text=None) -> Record:
    """Last resort: find the DNA in whatever this file is."""
    if str(path).lower().endswith((".xlsx", ".xlsm")):
        cells = _xlsx_cells(path)
    else:
        cells = re.split(r'[,;\t\r\n"]+', text or "")
    pieces = _dna_tokens(cells)
    if not pieces:
        # One last try: the whole file might be bare bases with no delimiters at all.
        whole = "".join(ch for ch in (text or "") if ch.isalpha())
        if len(whole) >= 8 and _dna_fraction(whole) >= 0.95:
            pieces = [whole.upper()]
    rec = Record(os.path.basename(str(path)).rsplit(".", 1)[0], _clean("".join(pieces)),
                 "linear", [])
    # Surfaced by the CLI. Joining in the wrong order yields a different construct that would still
    # audit cleanly, so the join must be visible, never assumed.
    rec.assembled_from = [len(p) for p in pieces] if len(pieces) > 1 else []
    return rec

def parse(path: str) -> Record:
    """FASTA, GenBank, or whatever the student actually saved.

    Dispatch is on CONTENT, not on the extension, because a file called .txt is as likely to hold
    GenBank as anything else - and a beginner's ".csv" is often a sequence pasted into one cell.
    """
    record = _parse_any(path)
    # A record with no bases is not a sequence; it is a file we could not read. Say which,
    # rather than handing an empty Record to an auditor that will dutifully report
    # "0 bp, linear" and "Empty sequence" as though those were findings about DNA.
    if not record.seq:
        raise NotASequenceFile(_diagnose(path))
    return record


def _parse_any(path):
    """Dispatch to a reader. May return a record with no bases; parse() judges that."""
    if str(path).lower().endswith((".xlsx", ".xlsm")):
        return _parse_loose(path)                      # binary; do not read it as text
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        text = f.read()
    stripped = text.lstrip()
    if stripped[:1] == ">":
        return _parse_fasta(path, text)
    if "LOCUS" in text[:2000] or "\nORIGIN" in text:
        return _parse_genbank(path, text)
    # Not a format with a spec. Go and find the DNA.
    return _parse_loose(path, text)


def _parse_fasta(path, text) -> Record:
    name = "sequence"
    chunks = []
    for line in text.splitlines():
        if line.startswith(">"):
            name = line[1:].strip().split()[0] if line[1:].strip() else name
        else:
            chunks.append(line)
    seq = _clean("".join(chunks))
    return Record(name, seq, "linear", [])


_LOC_RE = re.compile(r"(complement\()?\s*<?(\d+)\.\.>?(\d+)\)?")


def _parse_location(loc):
    """Return (start, end, strand) for simple/complement/join spans (1-based)."""
    strand = -1 if "complement" in loc else 1
    coords = re.findall(r"(\d+)\.\.(\d+)", loc)
    if not coords:
        one = re.findall(r"(\d+)", loc)
        if not one:
            return None
        p = int(one[0])
        return (p, p, strand)
    starts = [int(a) for a, _ in coords]
    ends = [int(b) for _, b in coords]
    return (min(starts), max(ends), strand)


def _parse_genbank(path, text) -> Record:
    lines = text.splitlines()

    # name + topology + DECLARED LENGTH from LOCUS
    name, topology, declared = "sequence", "linear", None
    for ln in lines:
        if ln.startswith("LOCUS"):
            parts = ln.split()
            if len(parts) >= 2:
                name = parts[1]
            if "circular" in ln.lower():
                topology = "circular"
            # The LOCUS line declares the length, and it was read for the name and the
            # topology and not for that. A file whose ORIGIN block had been cut off --
            # LOCUS saying 129 bp with 60 bases present -- parsed as a 60-base sequence
            # and was audited as though it were whole. Every per-base finding below is
            # then a finding about a fragment, reported as a finding about the construct,
            # with nothing saying which one it is.
            for _i, _tok in enumerate(parts):
                if _tok == "bp" and _i >= 1:
                    try:
                        declared = int(parts[_i - 1].replace(",", ""))
                    except ValueError:
                        declared = None
                    break
            break

    # sequence from ORIGIN..//
    seq_lines, in_origin = [], False
    for ln in lines:
        if ln.startswith("ORIGIN"):
            in_origin = True
            continue
        if in_origin:
            if ln.startswith("//"):
                break
            seq_lines.append(ln)
    seq = _clean("".join(seq_lines))

    # features
    features = []
    in_feat = False
    cur_kind = None
    cur_loc = ""
    cur_quals = {}

    def flush():
        if not cur_kind:
            return
        loc = _parse_location(cur_loc)
        if not loc:
            return
        role = _ROLE_KEYS.get(cur_kind.lower(), cur_kind.lower())
        label = ""
        for q in _LABEL_QUALS:
            if q in cur_quals and cur_quals[q]:
                label = cur_quals[q]
                break
        features.append(Feature(loc[0], loc[1], loc[2], role, label))

    for ln in lines:
        if ln.startswith("FEATURES"):
            in_feat = True
            continue
        if not in_feat:
            continue
        if ln.startswith("ORIGIN") or ln.startswith("//"):
            flush()
            break
        # feature key line: 5 spaces, key in cols ~6-20, then location
        if len(ln) > 5 and ln[5] != " " and not ln.lstrip().startswith("/"):
            flush()
            cur_kind = ln[5:21].strip()
            cur_loc = ln[21:].strip()
            cur_quals = {}
        elif ln.lstrip().startswith("/"):
            m = re.match(r'\s*/([^=]+)=?"?([^"]*)"?', ln)
            if m:
                cur_quals[m.group(1).strip().lower()] = m.group(2).strip()
        elif cur_kind and cur_loc and ln.strip() and not cur_quals:
            # location continuation (join spanning lines)
            cur_loc += ln.strip()

    # A declared length that does not match the bases present means the file is not the
    # record it says it is. Refusing is the only honest answer: the audit cannot report
    # on a construct from a fragment of it, and reporting on the fragment under the
    # construct's name is exactly the drift this project exists to prevent. Which two
    # things disagree, and what each says, per rule four.
    if declared is not None and seq and declared != len(seq):
        raise NotASequenceFile(
            "%s: the LOCUS line declares %d bp but the ORIGIN block holds %d base(s).\n"
            "  Two readings of the same file disagree about how long the sequence is, so "
            "neither can be\n"
            "  trusted. This is what a truncated download, an interrupted write or a "
            "partial copy looks\n"
            "  like. Nothing was audited. Fetch the file again and compare it against "
            "its source."
            % (os.path.basename(str(path)), declared, len(seq)))

    return Record(name, seq, topology, features)
