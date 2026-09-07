#!/usr/bin/env python3
"""Show validation coverage, BioModels cost, and KEGG conservation structure."""

import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from result_style import (
    INK,
    BLUE,
    GREEN,
    ORANGE,
    RED,
    GRAY,
    HAIR,
    PALE,
    configure,
    heading,
    clean,
    verdict,
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
        verdict(ax, 8.4, i, ok, 11)
        ext = external_by_name[row["name"]]
        matches = sum(ext["synkit"][k] == ext["crnt4sbml"][k] for k in SHARED)
        ax.text(
            9.9,
            i,
            f"{matches}/{len(SHARED)}",
            ha="center",
            va="center",
            fontsize=8,
            color=GREEN if matches == len(SHARED) else RED,
            weight="bold",
        )
        ax.plot([-4.7, 10.55], [i + 0.5, i + 0.5], color=HAIR, lw=0.45)
    heading(ax, "A", "Reference checks and independent SBML comparison", y=1.02)
    ax.text(
        0,
        -0.09,
        "Colour scaled within each count column; WR: ● yes / ○ no; external: shared quantities agreeing.",
        transform=ax.transAxes,
        fontsize=7.5,
        color=GRAY,
    )


def biomodels_panel(ax, evidence):
    models = sorted(
        (r for r in evidence["models"] if r["n_reactions"] > 0),
        key=lambda r: r["n_species"],
    )
    for i, row in enumerate(models):
        siph = row["siphon_seconds"]
        moiety = row["moiety_seconds"]
        if moiety is not None:
            ax.plot([siph, moiety], [i, i], color=HAIR, lw=1.2, zorder=1)
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
                color=RED,
            )
    ax.set(
        xscale="log",
        xlim=(5e-5, 400),
        ylim=(len(models) - 0.5, -0.6),
        yticks=range(len(models)),
        yticklabels=[f"{r['model_id'][-4:]} · {r['n_species']}" for r in models],
        xlabel="Wall time (s)",
    )
    ax.set_xticks([1e-4, 1e-2, 1, 100])
    ax.tick_params(axis="y", length=0, pad=4)
    clean(ax, "x")
    heading(ax, "B", "BioModels analysis costs", y=1.20)
    ax.text(
        -0.28,
        1.04,
        "ID suffix · species",
        transform=ax.transAxes,
        fontsize=7.5,
        color=GRAY,
    )
    ax.scatter([], [], color=GREEN, marker="s", s=23, label="Siphons")
    ax.scatter([], [], color=BLUE, marker="o", s=25, label="Moieties")
    ax.legend(
        loc="lower left",
        bbox_to_anchor=(0.07, 1.005),
        fontsize=7.5,
        ncol=2,
        frameon=False,
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
        color=GREEN,
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
    ax.set(xlim=(0, 11.5), ylim=(-2.15, 9.0))
    ax.axis("off")
    heading(ax, "C", "Conserved pools in KEGG", y=1.15)
    ax.text(0, 8.8, "One circle = one minimal conserved pool", fontsize=7.5, color=GRAY)
    for i, row in enumerate(modules):
        y = 7.4 - i * 2.6
        ax.text(0, y + 0.20, names[row["module_id"]], fontsize=9, weight="bold")
        ax.text(11.5, y + 0.20, row["module_id"], ha="right", fontsize=7.5, color=GRAY)
        # Pool counts alone hide how much network each pool has to cover, so
        # name the reduced model size that the supports are drawn from.
        ax.text(
            0,
            y - 0.45,
            f"{row['n_species']} species · {row['n_reactions']} reactions · "
            f"δ = {row['deficiency']}",
            fontsize=7.5,
            color=GRAY,
        )
        for j, pool in enumerate(sorted(row["moieties"], key=len, reverse=True)):
            x = 0.50 + j * 0.88
            ax.scatter(
                x, y - 1.15, s=105, color=GREEN, alpha=0.20, edgecolor=GREEN, lw=0.7
            )
            ax.text(
                x,
                y - 1.15,
                str(len(pool)),
                ha="center",
                va="center",
                fontsize=7.5,
                color=INK,
            )
        ax.plot([0, 11.5], [y - 1.55, y - 1.55], color=HAIR, lw=0.6)
    # Aligned with panel B's closing note now that this panel's box is taller.
    ax.text(
        0.5,
        -0.19,
        "Numbers = species in each pool",
        transform=ax.transAxes,
        ha="center",
        fontsize=8,
        color=GRAY,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence-dir", type=Path, default=Path("analysis"))
    parser.add_argument(
        "--output", type=Path, default=Path("Figure/analysis_overview.pdf")
    )
    args = parser.parse_args()
    configure()
    e = args.evidence_dir
    validation = load(e / "crn-validation.json")
    external = load(e / "crn-external-crnt4sbml.json")
    biomodels = load(e / "crn-biomodels.json")
    kegg = load(e / "crn-kegg.json")
    fig = plt.figure(figsize=(7.2, 6.9))
    a = fig.add_axes([0.045, 0.525, 0.925, 0.40])
    b = fig.add_axes([0.145, 0.115, 0.355, 0.275])
    c = fig.add_axes([0.585, 0.085, 0.385, 0.320])
    validation_panel(a, validation, external)
    biomodels_panel(b, biomodels)
    kegg_panel(c, kegg)
    save(fig, args.output.parent, args.output.stem)


if __name__ == "__main__":
    main()
