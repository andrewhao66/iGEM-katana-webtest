"""core.result — the structured outcome of a build.

Why this exists. The GUI used to run the engine as a subprocess and grep its English
prose: kg_verdict.classify looked for "BLOCK", "SEALED:" and "NOT enforced",
kg_rebuild looked for ".gb:" and ".csv:", and test_determinism looked for "seq_sha256:"
and "not self-consistent". The engine's human-readable output was therefore a machine
interface, and rewording any message broke the window SILENTLY -- classify failing to
find "SEALED:" fell through to "exited cleanly but printed no SEALED line" and reported
REVIEW, presenting a successful build as a problem.

kg_katana_tabs chose the subprocess for a good reason, stated in its own docstring: one
implementation of every gate, so the window cannot disagree with the command line.
Returning a structured object achieves that goal more strongly -- one literal code path
whose result is data rather than prose -- and it is what the browser front end needs,
since there is no stdout to grep in Pyodide.

The five status tokens are the ones kg_audit already uses, so the forward engine and the
reverse auditor speak one verdict language instead of two that must be kept in step.

Standard library only. Loaded by Pyodide in the browser front end.
"""

PASS = "PASS"
FLAG = "FLAG"
FAIL = "FAIL"
NOTE = "NOTE"
SKIP = "SKIP"


class Finding(object):
    """One thing the build observed, with why it matters and what to do."""

    __slots__ = ("category", "status", "summary", "loc", "detail", "fix")

    def __init__(self, category, status, summary, loc="", detail="", fix=""):
        self.category = category
        self.status = status
        self.summary = summary
        self.loc = loc
        self.detail = detail
        self.fix = fix

    def to_dict(self):
        return {"category": self.category, "status": self.status,
                "summary": self.summary, "loc": self.loc,
                "detail": self.detail, "fix": self.fix}

    def __repr__(self):
        return "<Finding %s %s: %s>" % (self.status, self.category, self.summary)


class StageResult(object):
    """What one pipeline stage did. `ok=False` means it refused to continue."""

    __slots__ = ("name", "ok", "findings", "data")

    def __init__(self, name, ok, findings=None, data=None):
        self.name = name
        self.ok = bool(ok)
        self.findings = list(findings or [])
        self.data = dict(data or {})

    def to_dict(self):
        return {"name": self.name, "ok": self.ok, "data": self.data,
                "findings": [f.to_dict() for f in self.findings]}


class BuildResult(object):
    """Everything a build produced: stages, findings, the seal, and the files.

    A renderer turns this into text, widgets, HTML or JSON. It is the only thing a
    caller needs, and nothing downstream has to parse a sentence.
    """

    def __init__(self):
        self.construct_id = ""
        self.version = None
        self.insert_len = 0
        self.seq_sha256 = ""
        self.library = ""
        self.library_root = ""
        self.pinned = False
        self.dry_run = False
        self.outputs = {}
        self.features = []
        self.consumed = {}
        self.stages = []
        # The text the pipeline printed. A renderer that captured it can
        # replay it; a wrapper that wants only data can ignore it.
        self.log = ""

    # ---- building it up -------------------------------------------------------
    def add(self, stage):
        """Record a stage. Stages accumulate in the order they ran."""
        self.stages.append(stage)
        return stage

    def stage(self, name):
        """The LAST stage recorded under `name`, or None.

        Last rather than first: a stage can report more than once -- validate and the
        dry-lab gate both do -- and a caller asking for it wants the latest word.
        """
        for s in reversed(self.stages):
            if s.name == name:
                return s
        return None

    # ---- reading it off ------------------------------------------------------
    @property
    def findings(self):
        out = []
        for s in self.stages:
            out.extend(s.findings)
        return out

    def count(self, status):
        return sum(1 for f in self.findings if f.status == status)

    @property
    def verdict(self):
        """PASS / REVIEW / FAIL, by the same rules the reverse auditor uses."""
        # A result with no stages has not been computed. It used to read PASS with exit
        # code 0 -- "we ran nothing" presented as "we ran everything and it was fine",
        # which is the exact confusion the SKIP tier exists to prevent, one level up.
        #
        # Nothing reaches it today: build() records a stage for a refusal and for a
        # SystemExit before anything else can happen. But a default that is safe only
        # because no caller has hit it yet is a defect waiting for one, and the cost of
        # closing it is this comment and two lines.
        if not self.stages:
            return FAIL
        for f in self.findings:
            if f.status == FAIL:
                return FAIL
        for s in self.stages:
            if s.ok is False:
                return FAIL
        for f in self.findings:
            if f.status == FLAG:
                return "REVIEW"
        # A SKIP means a check did not run, and that is not a pass. The two readings of
        # the same result used to disagree: verdict / exit_code / --json said PASS while
        # kg_verdict.from_result escalated the same SKIP to REVIEW, so a CI step gated on
        # the exit code was told PASS on a build where a gate never ran, with the window
        # beside it saying one had not. One of the two had to be wrong, and this project's
        # rule says which: "we did not check" must never read as "checked and fine".
        #
        # This is only safe to enforce because the off-target gate now actually finds the
        # genome the bundle ships. While it did not, every single build carried a SKIP and
        # making that REVIEW would have been noise on every run -- which is the cry-wolf
        # failure, not a fix for it.
        for f in self.findings:
            if f.status == SKIP:
                return "REVIEW"
        return PASS

    @property
    def exit_code(self):
        """0 PASS, 5 REVIEW, 1 FAIL -- the codes the reverse auditor already returns."""
        return {PASS: 0, "REVIEW": 5, FAIL: 1}[self.verdict]

    def blocked_stage(self):
        """The name of the stage that refused, or None."""
        for s in self.stages:
            if s.ok is False:
                return s.name
        for s in self.stages:
            for f in s.findings:
                if f.status == FAIL:
                    return s.name
        return None

    def not_run(self):
        """The categories of check that did not run.

        A caller can say so out loud, which is the whole reason SKIP is a separate tier:
        a gate that vanished silently is worse than no gate.
        """
        return [f.category for f in self.findings if f.status == SKIP]

    def to_dict(self):
        return {
            "construct_id": self.construct_id,
            "version": self.version,
            "insert_len": self.insert_len,
            "seq_sha256": self.seq_sha256,
            "library": self.library,
            "library_root": self.library_root,
            "pinned": self.pinned,
            "dry_run": self.dry_run,
            "verdict": self.verdict,
            "exit_code": self.exit_code,
            "blocked_stage": self.blocked_stage(),
            "not_run": self.not_run(),
            "outputs": dict(self.outputs),
            "features": list(self.features),
            "consumed": dict(self.consumed),
            "stages": [s.to_dict() for s in self.stages],
            "findings": [f.to_dict() for f in self.findings],
            "log": self.log,
        }
