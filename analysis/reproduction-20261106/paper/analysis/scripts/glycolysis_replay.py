#!/usr/bin/env python3
"""Render every recorded marking of the glycolysis certificate as a token map."""

import argparse
import json
from pathlib import Path
import matplotlib.pyplot as plt
from result_style import configure, BLUE, GRAY, HAIR, save

SPECIES = [
    ("alpha-D-Glucose", "Glc"),
    ("alpha-D-Glucose 6-phosphate", "G6P"),
    ("D-Fructose 6-phosphate", "F6P"),
    ("D-Fructose 1,6-bisphosphate", "FBP"),
    ("Glycerone phosphate", "DHAP"),
    ("D-Glyceraldehyde 3-phosphate", "GAP"),
    ("3-Phospho-D-glyceroyl phosphate", "BPG"),
    ("3-Phospho-D-glycerate", "3PG"),
    ("2-Phospho-D-glycerate", "2PG"),
    ("Phosphoenolpyruvate", "PEP"),
    ("Pyruvate", "Pyr"),
]
# Enzyme of the step that consumes each intermediate, as labelled in panel a.
ENZYME = dict(zip((key for key, _ in SPECIES[:-1]),
                  ("HK", "PGI", "PFK", "ALD", "TPI", "GAPDH", "PGK", "PGM",
                   "ENO", "PK")))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "glycolysis_analysis.json",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("Figure"))
    args = parser.parse_args()
    report = json.loads(args.input.read_text())
    run = report["realizability"]
    assert run["status"] == "realizable"
    rows = run["trajectory"]
    assert len(rows) == len(run["certificate"]) + 1
    configure()  # serif labels, matching Figure/glycolysis.tex
    fig, ax = plt.subplots(figsize=(7.2, 2.6))
    fig.subplots_adjust(left=0.025, right=0.98, top=0.95, bottom=0.04)
    ax.set(xlim=(-3.0, 15.6), ylim=(13.3, -2.6))
    ax.axis("off")
    ax.text(-2.9, -2.02, "(c)", fontsize=9.5, color="black")
    weights = report["structure"]["backbone_weights"]
    index = {key: i for i, (key, _) in enumerate(SPECIES)}
    for i in range(len(SPECIES)):
        ax.plot([-0.4, len(rows) - 0.6], [i, i], color="black", lw=0.3,
                ls=(0, (1, 2)), zorder=1)
    for j, row in enumerate(rows):
        assert row["step"] == j and row["nonnegative"]
        assert row["reaction"] == (run["certificate"][j - 1] if j else None)
        marking = row["marking"]
        assert (
            sum(marking.get(k, 0) * v for k, v in weights.items())
            == row["backbone_invariant"]
        )
        for key, value in marking.items():
            if key not in index:
                continue
            offsets = [0.0] if value == 1 else [-0.15, 0.15][:value]
            for dx in offsets:
                ax.scatter(j + dx, index[key], s=16, color="black", lw=0, zorder=4)
        for yy, key in ((11.45, "molecule_tokens"), (12.55, "backbone_invariant")):
            ax.text(j, yy, str(row[key]), ha="center", va="center",
                    fontsize=8, color="black")
        if not j:
            continue
        # Link each consumed species to what the fired reaction produced, and
        # carry untouched tokens straight across.
        before = rows[j - 1]["marking"]
        names = set(before) | set(marking)
        consumed = [k for k in names if marking.get(k, 0) < before.get(k, 0)]
        produced = [k for k in names if marking.get(k, 0) > before.get(k, 0)]
        assert len(consumed) == 1 and consumed[0] in index
        for k in names & set(index):
            if min(before.get(k, 0), marking.get(k, 0)):
                ax.plot([j - 1, j], [index[k]] * 2, color=GRAY, lw=0.6,
                        zorder=2)
        for k in produced:
            ax.plot([j - 1, j], [index[consumed[0]], index[k]], color="black",
                    lw=0.7, zorder=3)
        ax.text(j - 0.5, -0.95, ENZYME[consumed[0]], ha="center", va="bottom",
                fontsize=6.5, color="black")
    for i, (_, label) in enumerate(SPECIES):
        ax.text(-0.68, i, label, ha="right", va="center", fontsize=8)
    ax.plot([-2.9, len(rows) - 0.4], [10.75, 10.75], color="black", lw=0.5)
    ax.text(-0.68, 11.45, "Tokens", ha="right", va="center", fontsize=8)
    ax.text(-0.68, 12.55, "Weighted pool", ha="right", va="center", fontsize=8)
    # Aldolase is the observed first increase in molecule count.
    split = next(
        r["step"] for r in rows if r["molecule_tokens"] > rows[0]["molecule_tokens"]
    )
    ax.plot([split - 0.5, split - 0.5], [-0.5, 10.5], color="black", lw=0.5, ls="--")
    save(fig, args.output_dir, "glycolysis_replay")


if __name__ == "__main__":
    main()
