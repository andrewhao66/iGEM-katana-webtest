"""core.lock — the one implementation of the Parts Library manifest.

There were eight, in katana_lock.py, katana_build.py, add_part.py, find_part.py,
katana_init.py, kagami/kg_refs.py, kagami/kg_rebuild.py, kagami/kg_katana_tabs.py and
kagami/build_refs.py. katana_lock.py was meant to be the shared one and was imported by
exactly two files, both on the verifier path; everything else wrote its own. Nothing
forced convergence, and two of the copies disagreed about the two things that matter:
whether a row's hash is recomputed before the root is checked, and whether a Spec's pin
or the highest version number decides which part is loaded.

Standard library only. Loaded by Pyodide in the browser front end.
"""
import os

from . import hashing


class LockError(Exception):
    """The manifest cannot be trusted. Never swallowed -- it ends the caller's turn."""


def read(path):
    """Read a manifest. Returns (header, rows). Raises LockError if it is unusable.

    Blank lines are skipped. Every row is a dict keyed by the header's own names, so a
    reordered column is harmless; a MISSING one is not, and is refused here rather than
    surfacing as a KeyError three stages later.
    """
    with open(path, "r", encoding="utf-8") as f:
        lines = [ln.rstrip("\n") for ln in f if ln.strip() != ""]
    if not lines:
        raise LockError("%s is empty" % path)
    header = lines[0].split("\t")
    missing = [k for k in hashing.FIELDS if k not in header]
    if missing:
        raise LockError("%s is missing the column(s) %s, so its rows cannot be hashed"
                        % (path, ", ".join(missing)))
    rows = []
    for ln in lines[1:]:
        cells = ln.split("\t")
        if len(cells) < len(header):
            cells = cells + [""] * (len(header) - len(cells))
        rows.append(dict(zip(header, cells)))
    return header, rows


def verify_root(lock_path, root_path, pinned=None):
    """Verify the manifest is self-consistent, and optionally that it is THE library.

    Two independent checks, and the first is the one CODE-REPORT finding A was about:

    1. SELF-CONSISTENCY, always on. Every row's hash is RECOMPUTED from its fields and
       compared to the row_sha256 column, and the root is recomputed from those
       recomputed hashes. The engine used to hash the column as written, so editing a
       recorded accession without touching its row hash passed the gate -- the column
       still agreed with itself, and a root over unverified values is a root over
       nothing.
    2. EXTERNAL PIN, only when `pinned` is given. A library can be perfectly
       self-consistent and still be the wrong library: a stale sync, an old checkout, a
       second machine. Only a hash carried in from outside catches that.

    Returns (ok, message). The message names what disagrees, because "the manifest
    changed" is not a diagnosis.
    """
    if not os.path.exists(root_path):
        return False, "LOCK.root file missing"
    with open(root_path, "r", encoding="utf-8") as f:
        disk_root = f.read().strip()

    try:
        _header, rows = read(lock_path)
    except LockError as exc:
        return False, str(exc)

    # 1. Each row, from its fields.
    for row in rows:
        recomputed = hashing.row_sha256(row)
        written = row.get("row_sha256", "")
        if recomputed != written:
            _who = "%s v%s" % (row.get("id", "?"), row.get("version", "?"))
            # A TRUNCATED hash column is not an edited trust field, and saying so sent
            # the reader to audit source/version/class/outfile -- all of which were fine.
            # Measured: cutting the last 40 bytes off LOCK.tsv left this row's
            # row_sha256 at 26 of 64 characters, and the message said a trust field had
            # been edited.
            if written and len(written) != 64:
                return False, (
                    "%s: row_sha256 is %d characters, not 64 -- the manifest is "
                    "truncated or was written incompletely, so this row cannot be "
                    "checked at all.\n       Take a fresh copy of the library. Nothing "
                    "here can be trusted while a row is cut short."
                    % (_who, len(written)))
            if not written:
                return False, (
                    "%s: the row_sha256 column is empty, so there is nothing to check "
                    "this row against." % _who)
            # And show ENOUGH of each hash to see that they differ. Both were printed at
            # [:12], and when the difference fell past character 12 the message told the
            # reader two things disagreed and then showed them as identical.
            _n = 12
            while _n < 64 and recomputed[:_n] == written[:_n]:
                _n += 4
            return False, (
                "%s: row_sha256 does not match the row's own fields "
                "(recomputed %s..., column says %s...). A trust field -- source, "
                "version, class or outfile -- was edited."
                % (_who, recomputed[:_n], written[:_n]))

    # 2. The root, from the recomputed hashes.
    computed = hashing.lock_root(rows)
    if computed != disk_root:
        return False, ("LOCK is not self-consistent: recomputed root %s... != LOCK.root "
                       "file %s... (rows added or removed?)"
                       % (computed[:16], disk_root[:16]))

    if pinned:
        pinned = pinned.strip()
        if disk_root != pinned:
            return False, ("LOCK.root mismatch: library=%s... pinned=%s... "
                           "(wrong or stale library)" % (disk_root[:16], pinned[:16]))
        return True, "self-consistent + matches pin %s..." % pinned[:16]

    return True, ("self-consistent (%s...); NO EXTERNAL PIN -- pass --expect-root to bind"
                  % disk_root[:16])


def filename_sha12(outfile):
    """The 12-hex fingerprint embedded in a sealed part's filename.

    Sealed parts are named <id>__v<N>__<first 12 of seq_sha256>.<ext>. Returns "" when
    the name does not carry one, so a caller can tell "absent" from "wrong".
    """
    base = os.path.basename(str(outfile)).rsplit(".", 1)[0]
    tail = base.split("__")[-1]
    if len(tail) == 12 and all(c in "0123456789abcdefABCDEF" for c in tail):
        return tail.lower()
    return ""


def resolve(rows, part_id, pin=None, version=None, lib=None, notes=None):
    """Select exactly one manifest row for `part_id`. Raises LockError otherwise.

    Selection is by what the CALLER asked for -- the Spec's `pin`, an explicit
    `version`, or the seal's `lib` filename -- and NOT by the highest version number.
    That difference is CODE-REPORT finding B. resolve_parts() used to take max(version)
    and compare the pin against only that row, so a Spec pinned to an older sealed
    version was refused even though its row and its file were both still present, with
    a message advising the reader to update the Spec to match the library. Following
    that advice changes the construct, and it contradicts the append-only history the
    architecture promises: HrpS.Ec-opt has v1, v2 and v3 sealed, and pAP-Logic v5 pinned
    v2.

    A bare id with no pin resolves only when the library holds exactly one version of
    it. Guessing is what this function exists not to do.
    """
    candidates = [r for r in rows if r.get("id") == part_id]
    if not candidates:
        raise LockError("no sealed row for '%s'" % part_id)

    if lib is not None:
        want = os.path.basename(str(lib))
        narrowed = [r for r in candidates
                    if os.path.basename(r.get("outfile", "")) == want]
        # A lib filename that names nothing is a hint, not a command: the pin below is
        # the authority. Narrowing to nothing here would refuse a Spec whose seal block
        # carries a stale filename alongside a correct fingerprint.
        #
        # But not refusing is not the same as not SAYING. Preferring the pin in silence
        # resolves a discrepancy between two things the Spec asserts, and resolving a
        # discrepancy silently destroys the information that there was one -- which is
        # this project's fourth rule. So the disagreement is reported, and the build
        # continues: a stale filename is a documentation error, not a wrong part, and the
        # fingerprint is what binds the bases.
        if narrowed:
            candidates = narrowed
        elif notes is not None:
            notes.append(
                "%s: the Spec's seal names the file %s, which this library does not "
                "hold. The library has %s. The build used the fingerprint, which is the "
                "authority, so the part itself is right -- but the two disagree and one "
                "of them is out of date."
                % (part_id, want,
                   ", ".join(sorted(os.path.basename(r.get("outfile", ""))
                                    for r in candidates)) or "no file for this part"))

    if version is not None:
        want = str(version)
        candidates = [r for r in candidates if str(r.get("version", "")) == want]
        if not candidates:
            raise LockError("no sealed row for %s v%s" % (part_id, want))

    if pin:
        pin = pin.lower()
        matched = [r for r in candidates
                   if r.get("seq_sha256", "").lower().startswith(pin)]
        if not matched:
            held = ", ".join("v%s=%s" % (r.get("version", "?"),
                                         r.get("seq_sha256", "")[:12])
                             for r in rows if r.get("id") == part_id)
            raise LockError(
                "'%s': no sealed version matches the pin %s. The library holds %s. "
                "One of the two has moved on; this is the check working, not a bug."
                % (part_id, pin, held or "nothing"))
        candidates = matched

    if len(candidates) > 1:
        seen = sorted(str(r.get("version", "?")) for r in candidates)
        raise LockError(
            "'%s' is ambiguous: %d rows match (versions %s). Pin it in the Spec's seal "
            "block, or give a version." % (part_id, len(candidates), ", ".join(seen)))

    row = candidates[0]

    # Numeric versions are a convention the rest of the code relies on; an int() three
    # stages away would raise ValueError instead of saying what is wrong.
    if not str(row.get("version", "")).isdigit():
        raise LockError("'%s': version %r is not a number"
                        % (part_id, row.get("version")))

    # The filename carries the sequence hash. A disagreement here IS failure 5 from the
    # README: a part file named ...__47c4687cca62.gb whose sequence hashed to 3c840d2b,
    # same length and different bases, with a prepared manifest row claiming the name.
    fn12 = filename_sha12(row.get("outfile", ""))
    if fn12 and fn12 != row.get("seq_sha256", "")[:12].lower():
        raise LockError(
            "%s v%s: the filename says %s but the row's seq_sha256 is %s. Two sources "
            "disagree about which sequence this is -- report it, do not pick one."
            % (part_id, row.get("version"), fn12, row.get("seq_sha256", "")[:12]))

    return row
