#!/usr/bin/env python3
"""Render topology-sensitive performance diagnostics from committed evidence."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

BLUE = "#0072B2"
GREEN = "#009E73"
ORANGE = "#E69F00"
VERMILLION = "#D55E00"
INK = "#1A1A1A"
MUTED = "#666666"
HAIR = "#D9D9D9"

FAMILIES = {
    "chain": ("directed chain", BLUE, "o"),
    "reversible_chain": ("reversible chain", GREEN, "s"),
    "random_sparse": ("random sparse", ORANGE, "D"),
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
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.3,
            "axes.labelsize": 8.3,
            "xtick.labelsize": 8.3,
            "ytick.labelsize": 8.3,
            "text.color": INK,
            "axes.labelcolor": INK,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "pdf.fonttype": 42,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.025,
        }
    )


def panel_title(ax: plt.Axes, letter: str, title: str) -> None:
    ax.text(
        0.0,
        1.12,
        letter.upper(),
        transform=ax.transAxes,
        ha="left",
        va="center",
        color="white",
        fontsize=8.3,
        fontweight="bold",
        bbox={"boxstyle": "circle,pad=0.19", "facecolor": INK, "edgecolor": "none"},
        clip_on=False,
    )
    ax.text(
        0.085,
        1.12,
        title,
        transform=ax.transAxes,
        ha="left",
        va="center",
        fontsize=9.5,
        fontweight="bold",
        clip_on=False,
    )


def task_panel(
    ax: plt.Axes,
    evidence: dict,
    task: str,
    letter: str,
    title: str,
    ylim: tuple[float, float],
) -> None:
    measurements = evidence["measurements"]
    for family, (label, color, marker) in FAMILIES.items():
        rows = sorted(
            (
                item
                for item in measurements
                if item["family"] == family
                and item["task"] == task
                and item["seconds"] is not None
            ),
            key=lambda item: item["n_species"],
        )
        if not rows:
            continue
        x = np.asarray([item["n_species"] for item in rows], dtype=float)
        y = np.asarray([item["seconds"] for item in rows], dtype=float)
        ax.plot(
            x,
            y,
            color=color,
            marker=marker,
            markersize=3.1,
            linewidth=0.9,
            markeredgewidth=0,
            label=label,
            zorder=3,
        )

    if task == "conserved_moieties":
        timed_out = [
            item
            for item in measurements
            if item["task"] == task and item.get("timed_out")
        ]
        for item in timed_out:
            ax.scatter(
                item["n_species"],
                evidence["parameters"]["time_budget"],
                marker="^",
                s=28,
                color=VERMILLION,
                zorder=5,
            )
            ax.annotate(
                ">300 s",
                (item["n_species"], evidence["parameters"]["time_budget"]),
                xytext=(4, 3),
                textcoords="offset points",
                fontsize=8.3,
                color=VERMILLION,
                fontweight="bold",
            )
        ax.text(
            0.98,
            0.04,
            "random sparse: 395-species run omitted\nafter the 197-species timeout",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8.3,
            color=MUTED,
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(42, 560)
    ax.set_ylim(*ylim)
    ax.set_xticks([50, 100, 200, 400], labels=["50", "100", "200", "400"])
    ax.get_xaxis().set_minor_formatter(mpl.ticker.NullFormatter())
    ax.set_xlabel("species $|X|$")
    ax.set_ylabel("wall clock (s)")
    ax.yaxis.grid(True, which="major", color=HAIR, linewidth=0.45)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(HAIR)
    panel_title(ax, letter, title)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--evidence", type=Path, default=Path("analysis/crn-scaling.json")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("Figure/performance_diagnostics.pdf")
    )
    args = parser.parse_args()

    configure()
    evidence = load(args.evidence)
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.4), constrained_layout=False)
    fig.subplots_adjust(
        left=0.11, right=0.985, top=0.91, bottom=0.17, wspace=0.31, hspace=0.50
    )

    task_panel(
        axes[0, 0], evidence, "minimal_siphons", "a", "Siphon enumeration", (5e-5, 80)
    )
    task_panel(axes[0, 1], evidence, "crnt_summary", "b", "CRNT summary", (7e-2, 100))
    task_panel(
        axes[1, 0], evidence, "canonical_form", "c", "Canonicalization", (7e-3, 40)
    )
    task_panel(
        axes[1, 1],
        evidence,
        "conserved_moieties",
        "d",
        "Moiety enumeration",
        (1e-1, 700),
    )

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=3,
        frameon=False,
        fontsize=8.3,
        handlelength=1.8,
        columnspacing=1.8,
        bbox_to_anchor=(0.5, 0.015),
    )
    fig.text(
        0.5,
        0.06,
        "One archived measurement per size; lines connect observations",
        ha="center",
        color=MUTED,
        fontsize=8.3,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
