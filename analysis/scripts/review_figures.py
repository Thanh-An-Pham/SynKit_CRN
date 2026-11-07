#!/usr/bin/env python3
"""Build review figures and render manuscript values from validated JSON evidence."""

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


def render_tables(e):
    """Return numerical macro and table text for embedding in the manuscripts."""
    review = load(e / "crn-review-experiments.json", "synkit.crn-representation/1")
    archived = load(e / "crn-scaling.json", "synkit.crn-scaling/1")
    supplement = load(e / "crn-supplement-audit.json", "synkit.crn-supplementary/1")
    tables = {}
    last = {
        r["algorithm"]: r
        for r in review["paired_siphon_timings"]
        if r["n_species"] == 21
    }
    macros = "% Generated from crn-review-experiments.json; do not edit.\n"
    macros += "\\newcommand{\\ReviewExhaustiveMilliseconds}{%.2f}\n" % (
        1000 * last["exhaustive"]["median_seconds"]
    )
    macros += "\\newcommand{\\ReviewClosureMilliseconds}{%.2f}\n" % (
        1000 * last["closure"]["median_seconds"]
    )
    failures = sum(
        r["rounded_invalid_vectors"] > 0 for r in supplement["numerical_comparison"]
    )
    macros += "\\newcommand{\\ReviewRoundedFailures}{%d}\n" % failures
    first = {
        r["algorithm"]: r
        for r in review["paired_siphon_timings"]
        if r["n_species"] == 11
    }
    for name, rows in (("ReviewRatioSmall", first), ("ReviewRatioLarge", last)):
        ratio = rows["exhaustive"]["median_seconds"] / rows["closure"]["median_seconds"]
        macros += "\\newcommand{\\%s}{%.1f}\n" % (name, ratio)
    # Keep every timing quoted in prose tied to the same reports as the plots.
    timings = {
        (row["family"], row["size"], row["task"]): row["seconds"]
        for row in archived["measurements"]
    }
    for name, task in (
        ("ReviewReversibleSiphonSeconds", "minimal_siphons"),
        ("ReviewReversibleCanonicalSeconds", "canonical_form"),
        ("ReviewReversibleMoietySeconds", "conserved_moieties"),
        ("ReviewReversibleRankSeconds", "rank"),
    ):
        macros += "\\newcommand{\\%s}{%.3g}\n" % (
            name, timings["reversible_chain", 400, task]
        )
    for name, task in (
        ("ReviewRandomSiphonSeconds", "minimal_siphons"),
        ("ReviewRandomMoietySeconds", "conserved_moieties"),
    ):
        macros += "\\newcommand{\\%s}{%.3g}\n" % (
            name, timings["random_sparse", 400, task]
        )
    biomodels = load(e / "crn-biomodels.json", "synkit.crn-biomodels/1")
    largest = next(
        row for row in biomodels["models"] if row["model_id"] == "BIOMD0000000255"
    )
    for name, field in (
        ("ReviewBioModelSiphonSeconds", "siphon_seconds"),
        ("ReviewBioModelMoietySeconds", "moiety_seconds"),
    ):
        macros += "\\newcommand{\\%s}{%.3g}\n" % (name, largest[field])
    tables["crn-review-macros.tex"] = (macros)
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
    tables["crn-numerical-audit.tex"] = (
        "% Generated from crn-supplement-audit.json.\n"
        + "".join(row + " " + chr(92) * 2 + "\n" for row in rows)
    )
    random_rows = sorted(
        (
            r
            for r in archived["measurements"]
            if r["family"] == "random_sparse" and r["task"] == "minimal_siphons"
        ),
        key=lambda r: r["size"],
    )
    if len(random_rows) != 4 or not all(r["complete"] for r in random_rows):
        raise ValueError(
            "The paper requires all four complete random-network searches"
        )
    rows = [
        f"{r['n_species']} & {r['n_reactions']} & {r['result']} & {r['seconds']:.3g}"
        for r in random_rows
    ]
    tables["crn-random-siphons.tex"] = (
        "% Generated from crn-scaling.json.\n"
        + "".join(row + " " + chr(92) * 2 + "\n" for row in rows)
    )
    return tables


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
