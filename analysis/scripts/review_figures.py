#!/usr/bin/env python3
"""Build review figures and numerical macros from validated JSON evidence."""

from pathlib import Path
import argparse
import json
from visual_abstract import render_visual_abstract
from results_panels import representation, scaling, formose_figure
from biochemical_schemes import glycolysis_donors


def load(path, schema):
    data = json.loads(path.read_text())
    if data.get("schema") != schema or data.get("status") != "PASS":
        raise ValueError(f"Invalid or failed evidence: {path}")
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence-dir", type=Path, default=Path("analysis"))
    parser.add_argument("--output-dir", type=Path, default=Path("Figure"))
    parser.add_argument("--tables-dir", type=Path, default=None)
    args = parser.parse_args()
    e = args.evidence_dir
    review = load(e / "crn-review-experiments.json", "synkit.crn-review-experiments/1")
    archived = load(e / "crn-scaling.json", "synkit.crn-scaling/1")
    formose = load(e / "crn-formose.json", "synkit.crn-formose/1")
    supplement = load(e / "crn-supplement-audit.json", "synkit.crn-supplement-audit/1")
    donors = load(e / "crn-si-donor-control.json", "synkit.crn-si-donor-control/1")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.tables_dir is not None:
        args.tables_dir.mkdir(parents=True, exist_ok=True)
    render_visual_abstract(args.output_dir)
    representation(args.output_dir, review)
    formose_figure(args.output_dir, formose, supplement)
    glycolysis_donors(args.output_dir, donors)
    scaling(args.output_dir, review, archived)
    last = {
        r["algorithm"]: r
        for r in review["paired_siphon_timings"]
        if r["n_species"] == 21
    }
    macros = "% Generated from crn-review-experiments.json; do not edit.\n"
    macros += (
        "\\newcommand{\\ReviewExhaustiveSeconds}{%.2f}\n"
        % last["exhaustive"]["median_seconds"]
    )
    macros += "\\newcommand{\\ReviewClosureMilliseconds}{%.2f}\n" % (
        1000 * last["closure"]["median_seconds"]
    )
    failures = sum(
        r["rounded_invalid_vectors"] > 0 for r in supplement["numerical_comparison"]
    )
    if args.tables_dir is not None:
        (args.tables_dir / "crn-review-macros.tex").write_text(macros)
        rows = [
            f"{r['iteration']} & {r['n_species']} & {r['n_reactions']} & {r['seconds']:.3f}"
            for r in formose["growth"]
        ]
        (args.tables_dir / "crn-formose-growth.tex").write_text(
            "% Generated from crn-formose.json.\n"
            + (" " + chr(92) * 2 + "\n").join(rows)
            + "\n"
        )
        aliases = {
            "reversible_isomerization": "Reversible isomerization",
            "irreversible_isomerization": "Irreversible isomerization",
            "isomerization_cycle": "Isomerization cycle",
            "disjoint_reversible_blocks": "Disjoint pairs",
            "disjoint_blocks_with_stoichiometry": "Non-unit disjoint pairs",
            "michaelis_menten": "Michaelis-Menten",
            "reversible_michaelis_menten": "Reversible Michaelis-Menten",
            "edelstein": "Edelstein",
            "futile_cycle_1_site": "One-site futile cycle",
            "futile_cycle_2_site": "Two-site futile cycle",
            "shinar_feinberg_acr": "Shinar-Feinberg",
            "open_inflow_outflow": "Open reactor",
        }
        rows = [
            f"{aliases[r['network']]} & {r['svd_max_abs_residual']:.1e} & {r['rounded_invalid_vectors']} & 0"
            for r in supplement["numerical_comparison"]
        ]
        (args.tables_dir / "crn-numerical-audit.tex").write_text(
            "% Generated from crn-supplement-audit.json.\n"
            + (" " + chr(92) * 2 + "\n").join(rows)
            + "\n"
        )


if __name__ == "__main__":
    main()
