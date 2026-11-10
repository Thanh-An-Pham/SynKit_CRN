#!/usr/bin/env python3
"""Keep numerical values in main.tex and SI/si.tex tied to the saved evidence."""

import json
import os
from pathlib import Path
import re
import sys

PAPER = Path(__file__).resolve().parents[2]
ROOT = Path(os.environ.get("SYNKIT_ROOT", PAPER.parents[1] / "SynEco/SynKit"))
sys.path.insert(0, str(ROOT))

from analysis_overview import reference_table
from glycolysis_analysis import render_tex
from review_figures import render_tables
from formose_dpo_assets import table as dpo_table, macros as dpo_macros


def biomodels_table(report):
    rows = ["% Generated from crn-biomodels.json."]
    for model in report["models"]:
        values = [rf"\texttt{{{model['model_id']}}}", model["n_species"], model["n_reactions"]]
        values += [model.get(key, "--") for key in (
            "rank", "deficiency", "n_conserved_moieties", "n_minimal_siphons"
        )]
        for key in ("roundtrip_ok", "persistent"):
            value = model.get(key)
            values.append("--" if value is None else "yes" if value else "no")
        rows.append(" & ".join(map(str, values)) + r" \\")
    return "\n".join(rows) + "\n"


def biomodels_modifier_table(report):
    rows = ["% Generated from crn-biomodels-modifiers.json."]
    for model in report["models"]:
        if model.get("reaction_free"):
            continue
        values = [rf"\texttt{{{model['model_id']}}}", model["n_modifier_refs"]]
        values += [model[key] for key in (
            "deficiency_kept", "deficiency_cancelled", "deficiency_modifiers"
        )]
        rows.append(" & ".join(map(str, values)) + r" \\")
    return "\n".join(rows) + "\n"


def kegg_table(report):
    rows = ["% Generated from crn-kegg.json."]
    for model in report["modules"]:
        values = [rf"\texttt{{{model['module_id']}}}"]
        values += [model[key] for key in (
            "n_species", "n_reactions", "n_complexes", "n_linkage_classes", "rank", "deficiency"
        )]
        values += [len(model["moieties"]), len(model["siphons"])]
        rows.append(" & ".join(map(str, values)) + r" \\")
    return "\n".join(rows) + "\n"


def glycolysis_flux_table(report):
    flux = report["glycolysis_flux"]
    order = list(dict.fromkeys(flux["certificate"]))
    assert len(order) == 10 and sum(flux["flux"].values()) == 15
    rows = ["% Generated from crn-kegg.json."]
    for first, second in zip(order[:5], order[5:]):
        rows.append(
            rf"\texttt{{{first}}} & {flux['flux'][first]} & "
            rf"\texttt{{{second}}} & {flux['flux'][second]} \\"
        )
    return "\n".join(rows) + "\n"


def glycolysis_currency_table(report):
    rows = ["% Generated from crn-review-experiments.json."]
    for arm in sorted(report["representation"]["currency_arms"], key=lambda row: not row["drop_currency"]):
        values = ["Removed" if arm["drop_currency"] else "Retained", arm["n_species"], arm["n_reactions"]]
        outcome = arm["outcome"]
        values += [outcome["status"], len(outcome["certificate"]) if outcome["realizable"] else "--"]
        rows.append(" & ".join(map(str, values)) + r" \\")
    return "\n".join(rows) + "\n"


def formose_route_table(report):
    rows = ["% Generated from crn-supplement-audit.json."]
    qualitative = sorted(report["qualitative_routes"], key=lambda row: row["target"].count("C"))
    for name, entry in zip(("T4", "T6"), qualitative):
        route = next(row for row in report["formose_routes"] if row["target"] == name and row["pool"] == "C2-only")
        flow = route["outcome"]["requested_flux"]
        terms = [((str(count) if count > 1 else "") + "r_{" + rid.split("_")[1] + "}")
                 for rid, count in sorted(flow.items())]
        rows.append(
            "$" + name[0] + "_{" + name[1:] + "}$ & "
            + str(entry["outcome"]["max_depth_reached"]) + " & "
            + str(sum(flow.values())) + " & $" + "+".join(terms) + r"$ \\")
    return "\n".join(rows) + "\n"


def formose_inventory_table(report):
    rows = ["% Generated from crn-supplement-audit.json."]
    for index, pool in enumerate(("M1", "M2", "M3", "M4", "C2-only"), 1):
        routes = {row["target"]: row["outcome"] for row in report["formose_routes"] if row["pool"] == pool}
        marking = routes["T4"]["initial_marking"]
        c1, c2 = marking["C=O"], marking["O=CCO"]
        values = [rf"$M_{index}$", c1, c2, c1 + 2 * c2]
        values += ["yes" if routes[target]["realizable"] else "no" for target in ("T4", "T6")]
        rows.append(" & ".join(map(str, values)) + r" \\")
    return "\n".join(rows) + "\n"


def generated_blocks(paper=PAPER):
    """Return the expected inline blocks without writing source or figures."""
    evidence = paper / "analysis"
    tables = render_tables(evidence)
    def read(name):
        return json.loads((evidence / name).read_text())
    macros = tables["crn-review-macros.tex"]
    si_macros = "\n".join(
        line for line in macros.splitlines()
        if line.startswith(r"\newcommand{\ReviewRandomMoietySeconds}")
    ) + "\n"
    return {
        "main.tex": {
            "dpo-generation-macros": dpo_macros(),
            "review-macros": macros,
            "glycolysis-macros": render_tex(read("glycolysis_analysis.json")),
            "reference-networks": reference_table(
                read("crn-validation.json"), read("crn-external-crnt4sbml.json")
            ),
        },
        "SI/si.tex": {
            "formose-dpo-comparison": dpo_table(),
            "review-macros": si_macros,
            "biomodels-models": biomodels_table(read("crn-biomodels.json")),
            "biomodels-modifiers": biomodels_modifier_table(read("crn-biomodels-modifiers.json")),
            "kegg-modules": kegg_table(read("crn-kegg.json")),
            "glycolysis-flux": glycolysis_flux_table(read("crn-kegg.json")),
            "glycolysis-currency-control": glycolysis_currency_table(read("crn-review-experiments.json")),
            "formose-routes": formose_route_table(read("crn-supplement-audit.json")),
            "formose-inventories": formose_inventory_table(read("crn-supplement-audit.json")),
        },
    }


def block_pattern(name):
    return re.compile(
        r"^% BEGIN GENERATED: " + re.escape(name) + r"\n.*?"
        r"^% END GENERATED: " + re.escape(name) + r"$",
        re.M | re.S,
    )


def marked_block(name, body):
    return (
        f"% BEGIN GENERATED: {name}\n" + body.rstrip() + "\n"
        f"% END GENERATED: {name}"
    )


def main():
    for filename, blocks in generated_blocks().items():
        path = PAPER / filename
        previous = path.read_text()
        text = previous
        for name, body in blocks.items():
            pattern = block_pattern(name)
            if len(pattern.findall(text)) != 1:
                raise ValueError(f"Expected one {name} block in {filename}")
            text = pattern.sub(lambda match: marked_block(name, body), text)
        if text != previous:
            path.write_text(text)
            print(f"Updated numerical values in {filename}")


if __name__ == "__main__":
    main()
