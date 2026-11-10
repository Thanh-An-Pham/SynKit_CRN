#!/usr/bin/env python3
"""Build manuscript figures from validated scientific results."""

from pathlib import Path
import argparse
import json
from visual_abstract import render_visual_abstract
from results_panels import representation, scaling, formose_figure


def load(path, schema):
    data = json.loads(path.read_text())
    if data.get("schema") != schema or data.get("status") != "PASS":
        raise ValueError(f"Invalid or failed evidence: {path}")
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence-dir", type=Path, default=Path("analysis"))
    parser.add_argument("--output-dir", type=Path, default=Path("Figure"))
    args = parser.parse_args()
    e = args.evidence_dir
    review = load(e / "crn-review-experiments.json", "synkit.crn-representation/1")
    archived = load(e / "crn-scaling.json", "synkit.crn-scaling/1")
    formose = load(e / "crn-formose.json", "synkit.crn-formose/1")
    supplement = load(e / "crn-supplement-audit.json", "synkit.crn-supplementary/1")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    render_visual_abstract(args.output_dir)
    representation(args.output_dir, review)
    formose_figure(args.output_dir, formose, supplement)
    scaling(args.output_dir, review, archived)


if __name__ == "__main__":
    main()
