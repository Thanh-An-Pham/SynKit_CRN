#!/usr/bin/env python3
"""Render a cross-cohort structural landscape from committed evidence."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import result_style

BLUE = "#0072B2"
SKY = "#56B4E9"
GREEN = "#009E73"
ORANGE = "#E69F00"
VERMILLION = "#D55E00"
PURPLE = "#CC79A7"
INK = "#1A1A1A"
MUTED = "#666666"
HAIR = "#D9D9D9"

# Black and white: marker shape and fill separate the network sets.
COHORTS = {
    "validation": ("validation (12)", "black", "o"),
    "KEGG": ("KEGG (4)", "white", "s"),
    "BioModels": ("BioModels (10)", "white", "D"),
    "formose": ("formose (1)", "black", "^"),
}


def load(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        data = json.load(handle)
    if data.get("status") != "PASS" or not data.get("schema", "").startswith(
        "synkit.crn-"
    ):
        raise ValueError(f"Invalid evidence: {path}")
    return data


def configure() -> None:
    result_style.configure()
    mpl.rcParams.update({"font.size": 8.5, "axes.labelsize": 8.5,
                         "xtick.labelsize": 8, "ytick.labelsize": 8,
                         "legend.fontsize": 8})


def panel_title(ax: plt.Axes, letter: str, title: str) -> None:
    """Plain "(a)" above the axes; the caption names the panel."""
    result_style.plain_heading(ax, letter, y=1.03)


def collect(evidence_dir: Path) -> dict[str, list[dict]]:
    validation = load(evidence_dir / "crn-validation.json")
    kegg = load(evidence_dir / "crn-kegg.json")
    biomodels = load(evidence_dir / "crn-biomodels.json")
    formose = load(evidence_dir / "crn-formose.json")

    records: dict[str, list[dict]] = {key: [] for key in COHORTS}
    for item in validation["networks"]:
        row = dict(item["computed"])
        row["id"] = item["name"]
        records["validation"].append(row)
    for item in kegg["modules"]:
        records["KEGG"].append(
            {
                "id": item["module_id"],
                "n_species": item["n_species"],
                "n_reactions": item["n_reactions"],
                "n_complexes": item["n_complexes"],
                "n_linkage_classes": item["n_linkage_classes"],
                "rank": item["rank"],
                "deficiency": item["deficiency"],
                "weakly_reversible": item["weakly_reversible"],
                "persistent": item["persistent"],
            }
        )
    for item in biomodels["models"]:
        if item["n_species"] <= 0:
            continue
        records["BioModels"].append(
            {
                "id": item["model_id"],
                "n_species": item["n_species"],
                "n_reactions": item["n_reactions"],
                "n_complexes": item["n_complexes"],
                "n_linkage_classes": item["n_linkage_classes"],
                "rank": item["rank"],
                "deficiency": item["deficiency"],
                "weakly_reversible": item["weakly_reversible"],
                "persistent": item["persistent"],
            }
        )
    structure = formose["structure"]
    records["formose"].append(
        {
            "id": "iteration 4",
            "n_species": structure["n_species"],
            "n_reactions": structure["n_reactions"],
            "n_complexes": structure["n_complexes"],
            "n_linkage_classes": structure["n_linkage_classes"],
            "rank": structure["rank"],
            "deficiency": structure["deficiency"],
            "weakly_reversible": structure["weakly_reversible"],
            "persistent": structure["persistent"],
        }
    )
    return records


def scale_panel(ax: plt.Axes, records: dict[str, list[dict]]) -> None:
    for cohort, (label, color, marker) in COHORTS.items():
        x = [item["n_species"] for item in records[cohort]]
        y = [item["n_reactions"] for item in records[cohort]]
        ax.scatter(
            x,
            y,
            s=18,
            facecolor=color,
            marker=marker,
            label=label,
            edgecolors="black",
            linewidths=0.6,
            zorder=3,
        )
    limits = np.asarray([1.5, 1200.0])
    ax.plot(limits, limits, color="black", linewidth=0.5, linestyle="--", zorder=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*limits)
    ax.set_ylim(*limits)
    ax.set_xlabel("Species $|X|$")
    ax.set_ylabel("Directed reactions $|R|$")
    ax.legend(
        loc="lower right",
        handletextpad=0.35,
        borderaxespad=0.3,
    )
    panel_title(ax, "a", "Network size")


def decomposition_panel(ax: plt.Axes, records: dict[str, list[dict]]) -> None:
    cohorts = list(COHORTS)
    linkage, rank, deficiency = [], [], []
    for cohort in cohorts:
        total = sum(item["n_complexes"] for item in records[cohort])
        linkage.append(
            100 * sum(item["n_linkage_classes"] for item in records[cohort]) / total
        )
        rank.append(100 * sum(item["rank"] for item in records[cohort]) / total)
        deficiency.append(
            100 * sum(item["deficiency"] for item in records[cohort]) / total
        )
    y = np.arange(len(cohorts))
    bar = dict(height=0.58, edgecolor="black", linewidth=0.5)
    ax.barh(y, linkage, color="white", label=r"linkage $\ell$", **bar)
    ax.barh(y, rank, left=linkage, color="#BFBFBF", label="rank $s$", **bar)
    ax.barh(
        y,
        deficiency,
        left=np.asarray(linkage) + np.asarray(rank),
        color="#404040",
        label=r"deficiency $\delta$",
        **bar,
    )
    ax.set_yticks(y, labels=[COHORTS[key][0] for key in cohorts])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_ylim(len(cohorts) - 0.5, -0.5)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Share of aggregate complex count (%)")
    ax.tick_params(axis="y", which="both", left=False, right=False)
    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, -0.30),
        ncol=3,
        handlelength=1.2,
        columnspacing=0.8,
    )
    panel_title(ax, "b", "Complex-count composition")
    ax.text(
        1.0, -0.13, r"", transform=ax.transAxes, ha="right", fontsize=8.3, color=MUTED
    )


def deficiency_panel(ax: plt.Axes, records: dict[str, list[dict]]) -> None:
    cohorts = list(COHORTS)
    for row, cohort in enumerate(cohorts):
        values = np.asarray(
            [item["deficiency"] / item["n_complexes"] for item in records[cohort]]
        )
        offsets = (
            np.linspace(-0.13, 0.13, len(values))
            if len(values) > 1
            else np.asarray([0.0])
        )
        _, color, marker = COHORTS[cohort]
        ax.scatter(
            values,
            row + offsets,
            facecolor=color,
            marker=marker,
            s=16,
            edgecolors="black",
            linewidths=0.6,
            zorder=3,
        )
        ax.plot(
            [np.median(values), np.median(values)],
            [row - 0.23, row + 0.23],
            color="black",
            linewidth=1.0,
            zorder=4,
        )
    ax.set_yticks(np.arange(len(cohorts)), labels=[COHORTS[key][0] for key in cohorts])
    ax.invert_yaxis()
    ax.set_xlim(-0.015, 0.43)
    ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4], labels=["0", "0.1", "0.2", "0.3", "0.4"])
    ax.set_xlabel(r"Normalized deficiency $\delta/n$")
    ax.tick_params(axis="y", which="both", left=False, right=False)
    panel_title(ax, "c", "Normalized deficiency")
    ax.text(
        1.0, -0.20, "", transform=ax.transAxes, ha="right", fontsize=8.3, color=MUTED
    )


def verdict_panel(ax: plt.Axes, records: dict[str, list[dict]]) -> None:
    labels, yes, no, unknown = [], [], [], []
    for cohort in COHORTS:
        for key, short in (("weakly_reversible", "WR"), ("persistent", "pers.")):
            values = [item.get(key) for item in records[cohort]]
            total = len(values)
            labels.append(f"{cohort} · {short}")
            yes.append(100 * sum(value is True for value in values) / total)
            no.append(100 * sum(value is False for value in values) / total)
            unknown.append(100 * sum(value is None for value in values) / total)
    y = np.arange(len(labels))
    bar = dict(height=0.58, edgecolor="black", linewidth=0.5)
    ax.barh(y, yes, color="#404040", label="yes", **bar)
    ax.barh(y, no, left=yes, color="white", label="no", **bar)
    if any(unknown):
        ax.barh(
            y,
            unknown,
            left=np.asarray(yes) + np.asarray(no),
            color="#BFBFBF",
            label="not assessed",
            **bar,
        )
    ax.set_yticks(y, labels=labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_ylim(len(labels) - 0.5, -0.5)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Networks in set (%)")
    ax.tick_params(axis="y", which="both", left=False, right=False)
    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, -0.30),
        ncol=3,
        handlelength=1.2,
        columnspacing=0.9,
    )
    panel_title(ax, "d", "Structural verdicts")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence-dir", type=Path, default=Path("analysis"))
    parser.add_argument(
        "--output", type=Path, default=Path("Figure/structural_landscape.pdf")
    )
    args = parser.parse_args()

    configure()
    records = collect(args.evidence_dir)
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 7.0), constrained_layout=False)
    fig.subplots_adjust(
        left=0.145, right=0.975, top=0.95, bottom=0.12, wspace=0.38, hspace=0.48
    )
    scale_panel(axes[0, 0], records)
    decomposition_panel(axes[0, 1], records)
    deficiency_panel(axes[1, 0], records)
    verdict_panel(axes[1, 1], records)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
