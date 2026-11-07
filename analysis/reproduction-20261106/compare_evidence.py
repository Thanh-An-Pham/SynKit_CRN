#!/usr/bin/env python3
"""Compare reproduced scientific payloads with the paper's saved evidence."""

import argparse
import json
from pathlib import Path


RENAMES = {
    "crn-review-experiments.json": "crn-representation.json",
    "crn-supplement-audit.json": "crn-supplementary.json",
    "crn-si-donor-control.json": "crn-glycolysis-donors.json",
}
EXCLUDED_KEYS = {
    "environment", "source_sha256", "script_sha256", "source",
    "external_execution", "claim_boundary", "native_backend",
    "verification_source", "independent_search_source", "speedup",
    "samples", "median", "q1", "q3", "ratio_of_medians",
    "median_speedup", "closure_median", "exhaustive_median",
    "closure_iqr", "exhaustive_iqr",
}


def differences(before, after, path=""):
    if isinstance(before, dict) and isinstance(after, dict):
        for key in sorted(before.keys() | after.keys()):
            if key in EXCLUDED_KEYS or "seconds" in key:
                continue
            location = path + "/" + key
            if key not in before or key not in after:
                yield {"path": location, "difference": "missing key"}
            else:
                yield from differences(before[key], after[key], location)
    elif isinstance(before, list) and isinstance(after, list):
        if len(before) != len(after):
            yield {"path": path, "before_length": len(before), "after_length": len(after)}
        else:
            for index, (a, b) in enumerate(zip(before, after)):
                yield from differences(a, b, path + "/" + str(index))
    elif before != after:
        yield {"path": path, "paper": before, "reproduced": after}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("reproduced", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for old in sorted(args.baseline.glob("crn-*.json")):
        new = args.reproduced / RENAMES.get(old.name, old.name)
        changes = list(differences(json.loads(old.read_text()), json.loads(new.read_text())))
        rows.append({"paper_report": old.name, "reproduced_report": new.name,
                     "differences": changes})
    for name in ("structural-verification.json", "glycolysis_analysis.json"):
        before = json.loads((args.baseline / name).read_text())
        after = json.loads((args.reproduced / name).read_text())
        for value in (before, after):
            for key in ("environment", "native_backend", "verification_source"):
                value.pop(key, None)
        rows.append({"paper_report": name, "reproduced_report": name,
                     "scientific_payload_identical": before == after})
    status = "PASS" if len(rows) == 12 and all(
        not row.get("differences") and row.get("scientific_payload_identical", True)
        for row in rows
    ) else "FAIL"
    result = {"status": status, "reports": rows,
              "excluded_keys": sorted(EXCLUDED_KEYS), "excluded_patterns": ["*seconds*"],
              "note": "Environment, provenance and runtime values are excluded from study comparisons. Structural sets and the full glycolysis trajectory are also compared as exact payloads after removing environment and verification-script text."}
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(status, "-", len(rows), "evidence payloads compared")
    raise SystemExit(status != "PASS")


if __name__ == "__main__":
    main()
