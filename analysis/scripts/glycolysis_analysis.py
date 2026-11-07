#!/usr/bin/env python
"""Generate the numerical inputs used by the glycolysis figure."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

def build_analysis() -> dict[str, object]:
    """Run the structural analysis and replay the returned certificate."""
    from Experiment.CRN.benchmarks import (
        GLYCOLYSIS_FLUX,
        analyze_kegg_module,
        build_kegg_crn,
    )
    from Experiment.CRN.benchmarks.kegg import _species_labels
    from Experiment.CRN.tools.common import environment
    from synkit.CRN.Pathway.realizability import PathwayRealizability
    from synkit.CRN.Props.stoich import conserved_moieties

    crn = build_kegg_crn("M00001", drop_currency=True)
    summary = analyze_kegg_module("M00001")

    checker = PathwayRealizability().load_syncrn_and_flow(
        crn,
        flow=GLYCOLYSIS_FLUX,
        initial_marking={"alpha-D-Glucose": 1},
        species="label",
        reaction="label",
    )
    checker.build_petri_net_from_flow()
    outcome = checker.realizability_result()
    certificate = list(outcome.certificate or [])

    labels = _species_labels(crn)
    laws = [list(map(int, vector)) for vector in conserved_moieties(crn)]
    backbone = max(laws, key=lambda vector: sum(value > 0 for value in vector))
    weights = dict(zip(labels, backbone))
    backbone_species = {name for name, value in weights.items() if value > 0}

    net = checker.petri
    marking = dict(checker.initial_marking)
    trajectory: list[dict[str, object]] = []

    def record(step: int, reaction: str | None) -> None:
        species_marking = {
            name: int(marking.get(name, 0))
            for name in labels
        }
        trajectory.append(
            {
                "step": step,
                "reaction": reaction,
                "molecule_tokens": sum(
                    species_marking[name] for name in backbone_species
                ),
                "occupied_species": sum(
                    species_marking[name] > 0 for name in backbone_species
                ),
                "backbone_invariant": sum(
                    weights[name] * species_marking[name]
                    for name in backbone_species
                ),
                "nonnegative": all(value >= 0 for value in species_marking.values()),
                "marking": {
                    name: value for name, value in species_marking.items() if value
                },
            }
        )

    record(0, None)
    for step, reaction in enumerate(certificate, start=1):
        if not net.enabled(marking, reaction):
            raise RuntimeError(f"certificate fails before step {step}: {reaction}")
        marking = net.fire(marking, reaction)
        record(step, reaction)

    fired: dict[str, int] = {}
    for reaction in certificate:
        fired[reaction] = fired.get(reaction, 0) + 1

    invariant_values = {row["backbone_invariant"] for row in trajectory}
    if fired != dict(GLYCOLYSIS_FLUX):
        raise RuntimeError("certificate multiplicities differ from requested flux")
    if invariant_values != {2}:
        raise RuntimeError(f"backbone invariant changed during replay: {invariant_values}")
    if not all(bool(row["nonnegative"]) for row in trajectory):
        raise RuntimeError("certificate produced a negative species marking")

    return {
        "environment": environment(),
        "module": {
            "species": summary.n_species,
            "reactions": summary.n_reactions,
            "complexes": summary.n_complexes,
            "linkage_classes": summary.n_linkage_classes,
            "rank": summary.rank,
            "deficiency": summary.deficiency,
            "weakly_reversible": summary.weakly_reversible,
        },
        "structure": {
            "moiety_support_sizes": sorted(
                (sum(value > 0 for value in law) for law in laws), reverse=True
            ),
            "minimal_siphons": summary.siphons,
            "uncovered_siphons": summary.uncovered_siphons,
            "persistence_condition": summary.persistent,
            "backbone_weights": {
                name: weights[name] for name in labels if weights[name] > 0
            },
        },
        "realizability": {
            "status": outcome.status.value,
            "certificate": certificate,
            "certificate_length": len(certificate),
            "explored_states": outcome.explored_states,
            "max_depth": outcome.max_depth_reached,
            "termination_reason": outcome.termination_reason,
            "upper_firings": sum(count for count in fired.values() if count == 1),
            "lower_firings": sum(count for count in fired.values() if count == 2),
            "trajectory": trajectory,
        },
    }


def render_tex(data: dict[str, object]) -> str:
    """Render compact TeX macros consumed by ``glycolysis.tex``."""
    module = data["module"]
    structure = data["structure"]
    replay = data["realizability"]
    points = ",".join(
        f"{row['step']}/{row['molecule_tokens']}/{row['backbone_invariant']}"
        for row in replay["trajectory"]
    )
    split_step = next(
        row["step"]
        for row in replay["trajectory"]
        if row["molecule_tokens"] > replay["trajectory"][0]["molecule_tokens"]
    )
    support_sizes = structure["moiety_support_sizes"]
    return "\n".join(
        [
            "% Generated by glycolysis_analysis.py; do not edit by hand.",
            f"\\def\\GlySpecies{{{module['species']}}}",
            f"\\def\\GlyReactions{{{module['reactions']}}}",
            f"\\def\\GlyComplexes{{{module['complexes']}}}",
            f"\\def\\GlyLinkage{{{module['linkage_classes']}}}",
            f"\\def\\GlyRank{{{module['rank']}}}",
            f"\\def\\GlyDeficiency{{{module['deficiency']}}}",
            f"\\def\\GlyMoieties{{{len(support_sizes)}}}",
            f"\\def\\GlyBackboneSupport{{{support_sizes[0]}}}",
            f"\\def\\GlySiphons{{{len(structure['minimal_siphons'])}}}",
            f"\\def\\GlyUncovered{{{len(structure['uncovered_siphons'])}}}",
            f"\\def\\GlyCertificateLength{{{replay['certificate_length']}}}",
            f"\\def\\GlyExplored{{{replay['explored_states']}}}",
            f"\\def\\GlyDepth{{{replay['max_depth']}}}",
            f"\\def\\GlyUpperFirings{{{replay['upper_firings']}}}",
            f"\\def\\GlyLowerFirings{{{replay['lower_firings']}}}",
            f"\\def\\GlySplitStep{{{split_step}}}",
            f"\\def\\GlyTrajectory{{{points}}}",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--tex", type=Path)
    args = parser.parse_args()

    data = build_analysis()
    args.json.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    if args.tex is not None:
        args.tex.write_text(render_tex(data), encoding="utf-8")
    return 0


if __name__ == "__main__":
    from Experiment.CRN.tools.common import native_experiment

    raise SystemExit(native_experiment(main)())
