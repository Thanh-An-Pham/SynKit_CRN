#!/usr/bin/env python3
"""Show validation coverage, BioModels cost, and KEGG conservation structure."""

import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.path import Path as MarkerPath
from matplotlib.patches import Rectangle
from result_style import (
    BLUE,
    GREEN,
    ORANGE,
    RED,
    HAIR,
    PALE,
    configure,
    plain_heading as heading,
    save,
)


def load(path):
    data = json.loads(path.read_text())
    if data.get("status") != "PASS" or not data.get("schema", "").startswith(
        "synkit.crn-"
    ):
        raise ValueError(f"Invalid evidence: {path}")
    return data


ALIASES = {
    "reversible_isomerization": "Reversible isomerization",
    "irreversible_isomerization": "Irreversible isomerization",
    "isomerization_cycle": "Isomerization cycle",
    "disjoint_reversible_blocks": "Disconnected reversible pairs",
    "disjoint_blocks_with_stoichiometry": "Pairs with stoichiometry",
    "michaelis_menten": "Michaelis–Menten",
    "reversible_michaelis_menten": "Reversible Michaelis–Menten",
    "edelstein": "Edelstein switch",
    "futile_cycle_1_site": "Phosphorylation: one site",
    "futile_cycle_2_site": "Phosphorylation: two sites",
    "shinar_feinberg_acr": "Shinar–Feinberg motif",
    "open_inflow_outflow": "Open reactor",
}
SHARED = [
    "n_species",
    "n_reactions",
    "n_complexes",
    "n_linkage_classes",
    "rank",
    "deficiency",
    "weakly_reversible",
    "linkage_class_deficiencies",
]


def validation_panel(ax, validation, external):
    networks = validation["networks"]
    n = len(networks)
    numeric = [
        "n_species",
        "n_reactions",
        "n_complexes",
        "n_linkage_classes",
        "rank",
        "deficiency",
        "n_conservation_laws",
    ]
    raw = np.array([[r["computed"][k] for k in numeric] for r in networks])
    cmap = LinearSegmentedColormap.from_list("counts", ["#FFFFFF", "#BBDDED"])
    scaled = raw / np.maximum(raw.max(axis=0), 1)
    ax.set(xlim=(-4.8, 10.7), ylim=(n - 0.4, -1.8))
    ax.axis("off")
    ax.imshow(
        scaled,
        cmap=cmap,
        aspect="auto",
        extent=(-0.5, 6.5, n - 0.5, -0.5),
        vmin=0,
        vmax=1,
    )
    headers = [
        "Species",
        "Rxns",
        "$n$",
        "$\ell$",
        "$s$",
        "$\delta$",
        "Laws",
        "WR",
        "Tests",
        "External",
    ]
    xs = list(range(8)) + [8.4, 9.9]
    for x, label in zip(xs, headers):
        ax.text(x, -1.03, label, ha="center", va="center", fontsize=8, weight="bold")
    external_by_name = {r["name"]: r for r in external["networks"]}
    for i, row in enumerate(networks):
        ax.text(-4.7, i, ALIASES[row["name"]], va="center", fontsize=8)
        for j, v in enumerate(raw[i]):
            ax.text(j, i, str(v), va="center", ha="center", fontsize=8)
        # A filled/open dot avoids confusing weak reversibility with a test pass.
        wr = row["computed"]["weakly_reversible"]
        ax.scatter(
            7, i, s=26, facecolor=BLUE if wr else "white", edgecolor=BLUE, lw=0.8
        )
        ok = row["passed"] and not row["mismatches"] and not row["crosscheck_failures"]
        # Draw the check as a vector: Arial does not contain this glyph.
        if ok:
            check = MarkerPath([(-0.45, 0), (-0.12, -0.3), (0.5, 0.4)])
            ax.scatter(8.4, i, marker=check, s=70, facecolor="none",
                       edgecolor=GREEN, linewidths=1.2)
        else:
            ax.scatter(8.4, i, marker="x", s=42, color=RED, linewidths=1.2)
        ext = external_by_name[row["name"]]
        matches = sum(ext["synkit"][k] == ext["crnt4sbml"][k] for k in SHARED)
        ax.text(
            9.9,
            i,
            f"{matches}/{len(SHARED)}",
            ha="center",
            va="center",
            fontsize=8,
            color="black",
            weight="bold",
        )
    heading(ax, "a", "Reference checks and independent SBML comparison", y=1.02)


def biomodels_panel(ax, evidence):
    models = sorted(
        (r for r in evidence["models"] if r["n_reactions"] > 0),
        key=lambda r: r["n_species"],
    )
    for i, row in enumerate(models):
        siph = row["siphon_seconds"]
        moiety = row["moiety_seconds"]
        if moiety is not None:
            ax.plot([siph, moiety], [i, i], color=HAIR, lw=1.2, ls="--", zorder=1)
        ax.scatter(siph, i, s=23, marker="s", color=GREEN, zorder=3)
        if moiety is not None:
            ax.scatter(moiety, i, s=25, marker="o", color=BLUE, zorder=3)
        else:
            ax.scatter(
                120, i, s=35, marker=">", facecolor="white", edgecolor=RED, zorder=4
            )
            ax.annotate(
                ">120 s",
                (120, i),
                xytext=(-2, 8),
                textcoords="offset points",
                ha="right",
                fontsize=8,
                color="black",
            )
    ax.set(
        xscale="log",
        xlim=(5e-5, 400),
        ylim=(len(models) - 0.5, -0.6),
        yticks=range(len(models)),
        yticklabels=[],
        xlabel="Wall time (s)",
    )
    ax.set_xticks([1e-4, 1e-2, 1, 100])
    ax.minorticks_off()
    ax.tick_params(axis="y", length=0, pad=4)
    # Separate columns keep identifiers and differently sized counts aligned.
    id_x, species_x = -0.29, -0.10
    for i, row in enumerate(models):
        ax.text(id_x, i, row["model_id"][-4:],
                transform=ax.get_yaxis_transform(), ha="left", va="center",
                fontsize=8, color="black", clip_on=False)
        ax.text(species_x, i, str(row["n_species"]),
                transform=ax.get_yaxis_transform(), ha="center", va="center",
                fontsize=8, color="black", clip_on=False)
    ax.grid(False, which="both", axis="both")
    heading(ax, "b", "BioModels analysis costs", y=1.20)
    for x, title, alignment in [(id_x, "ID\nsuffix", "left"),
                                 (species_x, "Species", "center")]:
        ax.text(x, 1.04, title, transform=ax.transAxes, ha=alignment,
                fontsize=7.5, color="black", clip_on=False)
    ax.scatter([], [], color=GREEN, marker="s", s=23, label="Siphons")
    ax.scatter([], [], color=BLUE, marker="o", s=25, label="Moieties")
    ax.scatter([], [], marker=">", s=35, facecolor="white", edgecolor=RED,
               label="Timeout (120 s)")
    ax.legend(
        loc="lower left",
        bbox_to_anchor=(0, 1.005),
        fontsize=7.5,
        ncol=3,
        frameon=False,
        handlelength=1.2,
        handletextpad=0.4,
        columnspacing=0.8,
    )
    ax.text(
        0.5,
        -0.33,
        f"{sum(r['roundtrip_ok'] for r in models)}/{len(models)} supported structural round trips",
        transform=ax.transAxes,
        ha="center",
        fontsize=8,
        color="black",
        weight="bold",
    )


def kegg_panel(ax, evidence):
    names = {
        "M00001": "Glycolysis",
        "M00307": "Pyruvate oxidation",
        "M00009": "Citrate cycle",
        "M00004": "Pentose phosphate",
    }
    modules = evidence["modules"]
    ax.set(xlim=(0, 11.5), ylim=(-2.75, 9.0))
    ax.axis("off")
    heading(ax, "c", "Conserved pools in KEGG", y=1.1111)
    for i, row in enumerate(modules):
        y = 8.9 - i * 3.0
        ax.text(0, y + 0.20, names[row["module_id"]], fontsize=9, weight="bold")
        ax.text(11.5, y + 0.20, row["module_id"], ha="right", fontsize=7.5, color="black")
        # Pool counts alone hide how much network each pool has to cover, so
        # name the reduced model size that the supports are drawn from.
        ax.text(
            0,
            y - 0.45,
            f"{row['n_species']} species · {row['n_reactions']} reactions · "
            f"δ = {row['deficiency']}",
            fontsize=7.5,
            color="black",
        )
        for j, pool in enumerate(sorted(row["moieties"], key=len, reverse=True)):
            x = 0.50 + j * 0.96
            ax.scatter(
                x, y - 1.15, s=125, color=GREEN, alpha=0.20, edgecolor=GREEN, lw=0.7
            )
            ax.text(
                x,
                y - 1.15,
                str(len(pool)),
                ha="center",
                va="center",
                fontsize=7.5,
                color="black",
            )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence-dir", type=Path, default=Path("analysis"))
    parser.add_argument(
        "--output", type=Path, default=Path("Figure/analysis_overview.pdf")
    )
    args = parser.parse_args()
    configure()
    font_manager.findfont("Arial", fallback_to_default=False)
    plt.rcParams.update({
        "font.family": "Arial",
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "mathtext.sf": "Arial",
        "text.color": "black",
        "axes.labelcolor": "black",
        "axes.edgecolor": "black",
        "xtick.color": "black",
        "ytick.color": "black",
        "axes.grid": False,
    })
    e = args.evidence_dir
    validation = load(e / "crn-validation.json")
    external = load(e / "crn-external-crnt4sbml.json")
    biomodels = load(e / "crn-biomodels.json")
    kegg = load(e / "crn-kegg.json")
    fig = plt.figure(figsize=(7.2, 6.9))
    a = fig.add_axes([0.045, 0.525, 0.925, 0.40])
    b = fig.add_axes([0.145, 0.115, 0.355, 0.275])
    c = fig.add_axes([0.585, 0.045, 0.385, 0.360])
    validation_panel(a, validation, external)
    biomodels_panel(b, biomodels)
    kegg_panel(c, kegg)
    save(fig, args.output.parent, args.output.stem)


if __name__ == "__main__":
    main()
