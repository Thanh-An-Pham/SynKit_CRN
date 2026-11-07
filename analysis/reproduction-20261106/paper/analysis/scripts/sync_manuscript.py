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


def _count(value):
    text = f"{int(value):d}"
    return text if len(text) < 5 else text[:-3] + r"\," + text[-3:]


def hydrogen_macros(report):
    """Macros for the hydrogen-model comparison in the formose section."""
    deepest = max(report["formose"], key=lambda row: row["depth"])
    cases = {case["name"]: case for case in report["balance_cases"]}
    values = {
        "HydrogenDepth": str(deepest["depth"]),
        "HydrogenSpecies": _count(deepest["lewis"]["n_species"]),
        "HydrogenReactions": _count(deepest["lewis"]["n_reactions"]),
        "HydrogenExplicitSeconds": f"{deepest['explicit']['seconds']:.1f}",
        "HydrogenLewisSeconds": f"{deepest['lewis']['seconds']:.1f}",
        "HydrogenLewisReactions": str(sum(c["lewis"]["n_reactions"] for c in cases.values())),
        "HydrogenInsertionExplicit": str(cases["incomplete_insertion"]["explicit"]["n_reactions"]),
    }
    lines = ["% Generated from crn-hydrogen.json; do not edit."]
    lines += [rf"\newcommand{{\{name}}}{{{value}}}" for name, value in values.items()]
    return "\n".join(lines) + "\n"


def hydrogen_table(report):
    """SI table rows comparing explicit hydrogen vertices with hydrogen counts."""
    lines = ["% Generated from crn-hydrogen.json."]
    for row in report["formose"]:
        e, l = row["explicit"], row["lewis"]
        lines.append(
            f"Formose, {row['depth']} iterations & {_count(l['n_species'])} & "
            f"{_count(e['n_reactions'])} & {e['seconds']:.3f} & {_count(l['n_reactions'])} & "
            f"{l['seconds']:.3f} & {'yes' if row['identical'] else 'no'} \\\\"
        )
    lines.append(r"\midrule")
    lines.append(
        r"Control & Seeds & \multicolumn{2}{c}{Balanced / reactions} & "
        r"\multicolumn{2}{c}{Balanced / reactions} & \\"
    )
    lines.append(r"\midrule")
    for case in report["balance_cases"]:
        e, l = case["explicit"], case["lewis"]
        seeds = ", ".join(r"\texttt{" + s.replace("#", r"\#") + "}" for s in case["seeds"])
        lines.append(
            f"{case['label']} & {seeds} & \\multicolumn{{2}}{{c}}{{{e['n_balanced']} / {e['n_reactions']}}} & "
            f"\\multicolumn{{2}}{{c}}{{{l['n_balanced']} / {l['n_reactions']}}} & \\\\"
        )
    return "\n".join(lines) + "\n"


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
            "review-macros": macros,
            "glycolysis-macros": render_tex(read("glycolysis_analysis.json")),
            "hydrogen-macros": hydrogen_macros(read("crn-hydrogen.json")),
            "reference-networks": reference_table(
                read("crn-validation.json"), read("crn-external-crnt4sbml.json")
            ),
        },
        "SI/si.tex": {
            "review-macros": si_macros,
            "formose-growth": tables["crn-formose-growth.tex"],
            "hydrogen-models": hydrogen_table(read("crn-hydrogen.json")),
            "numerical-audit": tables["crn-numerical-audit.tex"],
            "random-siphons": tables["crn-random-siphons.tex"],
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
