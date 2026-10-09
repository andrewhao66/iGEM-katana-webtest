"""core — the shared, dependency-free heart of Katana.

Everything in this package is standard library only (plus the vendored YAML parser),
because it is the module set the browser front end loads through Pyodide. An import of
subprocess or shutil.which here would silently cost the web version, so there is a test
that asserts it does not happen.

  core.hashing   the hashing conventions, defined once
  core.lock      LOCK.tsv: read, resolve, verify
  core.parts     sealed part files: read a sequence, render a record, audit a library
"""
