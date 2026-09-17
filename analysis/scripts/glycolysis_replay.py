#!/usr/bin/env python3
"""Render every recorded marking of the glycolysis certificate as a token map."""

import argparse
import json
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle
from result_style import configure, INK, BLUE, GREEN, ORANGE, GRAY, HAIR, PALE, save

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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(__file__).with_name("glycolysis_analysis.json"),
    )
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent)
    args = parser.parse_args()
    report = json.loads(args.input.read_text())
    run = report["realizability"]
    assert run["status"] == "realizable"
    rows = run["trajectory"]
    assert len(rows) == len(run["certificate"]) + 1
    configure()
    # This panel is embedded in Figure/glycolysis.tex; match its Arial labels.
    font_manager.findfont("Arial", fallback_to_default=False)
    plt.rcParams.update({"font.family": "Arial", "text.color": "black"})
    fig, ax = plt.subplots(figsize=(7.2, 2.8))
    fig.subplots_adjust(left=0.025, right=0.98, top=0.95, bottom=0.04)
    ax.set(xlim=(-3.0, 15.6), ylim=(13.6, -2.6))
    ax.axis("off")
    ax.text(-2.9, -2.02, "c", fontsize=9, color="black", weight="bold")
    ax.text(
        -1.95,
        -2.02,
        "Certificate replay: every marking and firing",
        fontsize=9,
        weight="bold",
    )
    ax.text(
        15.4,
        -2.02,
        f"{len(run['certificate'])} firings · {run['explored_states']} states searched",
        fontsize=8,
        ha="right",
        color="black",
    )
    for j, row in enumerate(rows):
        assert row["step"] == j and row["nonnegative"]
        assert row["reaction"] == (run["certificate"][j - 1] if j else None)
        marking = row["marking"]
        weights = report["structure"]["backbone_weights"]
        assert (
            sum(marking.get(k, 0) * v for k, v in weights.items())
            == row["backbone_invariant"]
        )
        ax.text(j, -0.97, str(j), ha="center", fontsize=8, color="black")
        for i, (key, label) in enumerate(SPECIES):
            value = marking.get(key, 0)
            ax.add_patch(
                Rectangle(
                    (j - 0.47, i - 0.45),
                    0.94,
                    0.90,
                    facecolor=PALE if value == 0 else "#D6E8F2" if value == 1 else "#80B8D4",
                    edgecolor="none",
                )
            )
            if value:
                ax.text(
                    j,
                    i,
                    str(value),
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="black",
                    weight="bold",
                )
        for yy, key, color in (
            (11.50, "molecule_tokens", ORANGE),
            (12.68, "backbone_invariant", GREEN),
        ):
            ax.text(
                j,
                yy,
                str(row[key]),
                ha="center",
                va="center",
                fontsize=8,
                color=color,
                weight="bold",
            )
    for i, (_, label) in enumerate(SPECIES):
        ax.text(-0.68, i, label, ha="right", va="center", fontsize=8)
    ax.text(-0.68, 11.5, "Tokens", ha="right", va="center", fontsize=8, color=ORANGE)
    ax.text(
        -0.68, 12.68, "Weighted pool", ha="right", va="center", fontsize=8, color=GREEN
    )
    ax.text(-0.68, -0.97, "Step", ha="right", fontsize=8, color="black")
    # Aldolase is the observed first increase in molecule count.
    split = next(
        r["step"] for r in rows if r["molecule_tokens"] > rows[0]["molecule_tokens"]
    )
    ax.plot([split - 0.5, split - 0.5], [-0.5, 10.5], color=ORANGE, lw=1, ls="--")
    save(fig, args.output_dir, "glycolysis_replay")


if __name__ == "__main__":
    main()
