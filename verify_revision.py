#!/usr/bin/env python3
"""Check revision evidence, printed equations, and generated numerical assets."""

import os
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys
import tempfile

PAPER = Path(__file__).resolve().parent
ROOT = PAPER.parents[1]


def main():
    for filename, schema in [
        ("crn-review-experiments.json", "synkit.crn-review-experiments/1"),
        ("crn-supplement-audit.json", "synkit.crn-supplement-audit/1"),
        ("crn-si-donor-control.json", "synkit.crn-si-donor-control/1"),
    ]:
        report = json.loads((PAPER / "analysis" / filename).read_text())
        assert report["schema"] == schema, filename
        assert report["status"] == "PASS" and all(report["checks"].values()), filename
    review = json.loads((PAPER / "analysis/crn-review-experiments.json").read_text())
    for filename, digest in review["source_sha256"].items():
        assert (
            hashlib.sha256((ROOT / filename).read_bytes()).hexdigest() == digest
        ), f"Stale experiment source: {filename}"
    supplement = json.loads((PAPER / "analysis/crn-supplement-audit.json").read_text())
    donors = json.loads((PAPER / "analysis/crn-si-donor-control.json").read_text())
    for filename, digest in donors["source_sha256"].items():
        assert (
            hashlib.sha256((ROOT / filename).read_bytes()).hexdigest() == digest
        ), filename
    assert donors["equations"] == supplement["equations"]
    assert (
        hashlib.sha256(
            (ROOT / "Experiment/CRN/supplement_audit.py").read_bytes()
        ).hexdigest()
        == supplement["script_sha256"]
    )
    printed = dict(
        (f"r{n}", eq.replace(r"\(\rightarrow\)", ">>"))
        for n, eq in re.findall(
            r"^r(\d+) & (.*?) \\\\", (PAPER / "SI/si.tex").read_text(), re.M
        )
    )
    assert (
        printed == supplement["equations"]
    ), "SI reaction table differs from audited equations"
    with tempfile.TemporaryDirectory(prefix="synkit-crn-assets-") as temp:
        env = {**os.environ, "PYTHONPATH": str(PAPER / "analysis/scripts")}
        subprocess.run(
            [
                sys.executable,
                str(PAPER / "analysis/scripts/review_figures.py"),
                "--evidence-dir",
                str(PAPER / "analysis"),
                "--output-dir",
                temp,
                "--tables-dir",
                temp,
            ],
            env=env,
            check=True,
        )
        main_text = (PAPER / "main.tex").read_text()
        si_text = (PAPER / "SI/si.tex").read_text()
        for line in (Path(temp) / "crn-review-macros.tex").read_text().splitlines():
            if line.startswith(r"\newcommand"):
                assert line in main_text, f"Stale or missing macro in main.tex: {line}"
        for line in (Path(temp) / "crn-formose-growth.tex").read_text().splitlines():
            if line.strip() and not line.startswith("%"):
                assert line.strip() in si_text, f"Stale or missing formose row in SI/si.tex: {line}"
        for line in (Path(temp) / "crn-numerical-audit.tex").read_text().splitlines():
            if line.strip() and not line.startswith("%"):
                assert line.strip() in si_text, f"Stale or missing audit row in SI/si.tex: {line}"
    text = (PAPER / "main.tex").read_text()
    abstract = text.split("\\abstract{", 1)[1].split("\\keywords", 1)[0]
    assert len(abstract.split()) <= 350, "Abstract exceeds 350 words"
    assert "\\cite" not in abstract, "Abstract contains citations"
    print(
        "Revision evidence, source hashes, SI equations, generated tables, and abstract: PASS"
    )


if __name__ == "__main__":
    main()
