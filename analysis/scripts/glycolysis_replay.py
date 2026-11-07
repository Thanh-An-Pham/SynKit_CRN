#!/usr/bin/env python3
"""Render every recorded marking of the glycolysis certificate as a token map."""

import argparse
import json
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from result_style import configure, GRAY, HAIR, PALE, save

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
    weights = report["structure"]["backbone_weights"]
    index = {key: i for i, (key, _) in enumerate(SPECIES)}
    n_steps = len(rows) - 1
    n_hex = sum(1 for key, _ in SPECIES if weights[key] == 2)

    fig = plt.figure(figsize=(7.2, 2.85))
    # Boxed token map; the counter table sits under it on the same x scale.
    ax = fig.add_axes([0.105, 0.255, 0.83, 0.60])
    ax.set(xlim=(-0.6, n_steps + 0.6), ylim=(len(SPECIES) - 0.5, -0.5))
    ax.set_yticks(range(len(SPECIES)), [label for _, label in SPECIES])
    ax.set_xticks(range(len(rows)), [str(j) for j in range(len(rows))])
    ax.tick_params(axis="both", labelsize=7.5, length=2.5)
    ax.tick_params(axis="x", top=False, pad=2)
    ax.tick_params(axis="y", right=False)
    # Flush left, so the letter lines up with "(a)" in Figure/glycolysis.tex.
    fig.text(0.0, 0.985, "(b)", fontsize=9.5, ha="left", va="top")
    # Hexose rows carry pool weight two: shade them as one band.
    ax.axhspan(-0.5, n_hex - 0.5, color=PALE, lw=0, zorder=0)
    for i in range(len(SPECIES)):
        ax.axhline(i, color=HAIR, lw=0.35, zorder=1)
    # Pool weight of each species, read off the conservation vector.
    right = ax.secondary_yaxis("right")
    right.set_yticks(range(len(SPECIES)),
                     [str(weights[key]) for key, _ in SPECIES])
    right.tick_params(labelsize=7.5, length=2.5, direction="in")
    ax.text(1.012, 1.015, "$w$", transform=ax.transAxes, fontsize=8,
            ha="left", va="bottom", style="italic")

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
            assert value in (1, 2)
            ax.plot(j, index[key], "o", ms=3.6, color="black", mew=0, zorder=5)
            if value == 2:  # two tokens: dot inside a ring
                ax.plot(j, index[key], "o", ms=7.0, mfc="white", mec="black",
                        mew=0.7, zorder=4)
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
                ax.plot([j - 1, j], [index[k]] * 2, color=GRAY, lw=0.9,
                        zorder=2, solid_capstyle="butt")
        for k in produced:
            ax.plot([j - 1, j], [index[consumed[0]], index[k]], color="black",
                    lw=0.8, zorder=3)
        # Enzyme fired between markings j-1 and j, set above the frame.
        ax.text(j - 0.5, -0.72, ENZYME[consumed[0]], ha="center", va="bottom",
                fontsize=6.5, rotation=90)
    # Aldolase is the observed first increase in molecule count.
    split = next(
        r["step"] for r in rows if r["molecule_tokens"] > rows[0]["molecule_tokens"]
    )
    ax.axvline(split - 0.5, color="black", lw=0.5, ls=(0, (4, 2)), zorder=1)
    ax.text(split - 0.38, len(SPECIES) - 0.75, "cleavage", fontsize=6.5,
            ha="left", va="bottom", style="italic")
    for yy, va, label in ((n_hex - 0.62, "bottom", "hexoses"),
                          (n_hex - 0.38, "top", "trioses")):
        ax.text(n_steps + 0.45, yy, label, fontsize=6.5, ha="right", va=va,
                style="italic")

    # Legend in the empty lower-left corner of the map.
    handles = [
        Line2D([], [], ls="", marker="o", ms=3.6, color="black", mew=0),
        Line2D([], [], ls="", marker="o", ms=7.0, mfc="white", mec="black",
               mew=0.7),
        Line2D([], [], color="black", lw=0.8),
        Line2D([], [], color=GRAY, lw=0.9),
    ]
    handles[1] = (handles[1], handles[0])
    two = ax.legend(handles, ["one token", "two tokens", "reaction fired",
                              "token carried"],
                    loc="lower left", bbox_to_anchor=(0.005, 0.01), fontsize=6.5,
                    handlelength=1.6, borderpad=0.45, labelspacing=0.3,
                    handletextpad=0.5, framealpha=1, borderaxespad=0.3,
                    handler_map={tuple: HandlerTuple(ndivide=1)})
    two.get_frame().set_linewidth(0.5)
    two.set_zorder(6)

    # Counter table: total tokens and the conserved weighted pool per marking.
    table = fig.add_axes([0.105, 0.02, 0.83, 0.12], sharex=ax)
    table.set(ylim=(1.6, -0.6))
    table.axis("off")
    for yy, key, label in ((0, "molecule_tokens", "Tokens"),
                           (1, "backbone_invariant", "Weighted pool")):
        table.text(-0.9, yy, label, ha="right", va="center", fontsize=7.5)
        for j, row in enumerate(rows):
            table.text(j, yy, str(row[key]), ha="center", va="center",
                       fontsize=7.5)
    table.plot([-0.6, n_steps + 0.6], [-0.6, -0.6], color="black", lw=0.5,
               clip_on=False)
    ax.set_xlabel("Firing step", fontsize=7.5, labelpad=1.5)
    save(fig, args.output_dir, "glycolysis_replay")


if __name__ == "__main__":
    main()
