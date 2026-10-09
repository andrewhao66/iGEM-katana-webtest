"""core.hashing — the single definition of this project's hashing conventions.

Reverse-engineered and confirmed against the live library on 2026-07-05, and frozen
since: every sealed part and every manifest row in existence was produced with these
exact rules, so changing one invalidates the library rather than improving it.

  seq_sha256   sha256 of the sequence, uppercased, letters only
  file_sha256  sha256 of the raw file bytes
  row_sha256   sha256 over the nine trust-bearing fields, "k=v" joined by newlines
  lock_root    sha256 over the row hashes, joined by newlines

These were previously defined in five places -- katana_lock.py, add_part.py,
katana_init.py, kagami/build_refs.py and test_determinism.py -- each carrying a comment
saying the duplication was deliberate and load-bearing, because a drift would make every
part that tool admitted unbuildable. The comments were right about the danger and wrong
about the remedy: one definition cannot drift from itself.

Standard library only. This module is loaded by Pyodide in the browser front end.
"""
import hashlib

# The order is part of the convention: row_sha256 is computed over these keys in this
# sequence. It must not contain row_sha256 itself, while the manifest HEADER does --
# a distinction that cost a build once.
FIELDS = ["id", "version", "seq_sha256", "file_sha256", "length",
          "source", "date", "class", "outfile"]

HEADER = FIELDS + ["row_sha256"]


def sha256_hex(data):
    """sha256 of bytes, as lowercase hex."""
    return hashlib.sha256(data).hexdigest()


def seq_sha256(seq):
    """The canonical sequence hash: uppercase, ASCII, nothing else.

    KATANA_SPEC v2 section 3.4 specifies UPPER+"|"+topology, but the existing LOCK and
    every sealed part and construct use plain UPPER. The engine follows the established
    convention, because the alternative is re-sealing the library.
    """
    return sha256_hex(seq.upper().encode("ascii"))


def file_sha256(path):
    """sha256 of a file's raw bytes. Never decoded -- line endings are part of the seal."""
    with open(path, "rb") as f:
        return sha256_hex(f.read())


def row_manifest(row):
    """The canonical string a row's hash is taken over."""
    return "\n".join("%s=%s" % (k, row.get(k, "")) for k in FIELDS)


def row_sha256(row):
    """Hash of a manifest row's trust-bearing fields.

    This is what seals the ROW rather than the sequence: editing a recorded accession,
    bumping a version, or repointing an outfile changes this value.
    """
    return sha256_hex(row_manifest(row).encode())


def lock_root(rows):
    """Merkle-style root over the rows.

    RECOMPUTES each row hash from its fields rather than reading the row_sha256 column.
    That difference is CODE-REPORT finding A: hashing the column as written let an
    edited `source` field pass the build engine's gate, because the column still agreed
    with itself. A root over values you did not verify is a root over nothing.
    """
    return sha256_hex("\n".join(row_sha256(r) for r in rows).encode())
