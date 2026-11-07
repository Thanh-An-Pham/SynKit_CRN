#!/usr/bin/env python3
"""Render topology-sensitive performance diagnostics from committed evidence."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import result_style
from schematics import arc, place, transition

BLUE = "#0072B2"
GREEN = "#009E73"
ORANGE = "#E69F00"
VERMILLION = "#D55E00"
INK = "#1A1A1A"
MUTED = "#666666"
HAIR = "#D9D9D9"

# Black and white: marker fill and dash pattern separate the families.
FAMILIES = {
    "chain": ("directed chain", ("black", "-"), "o"),
    "reversible_chain": ("reversible chain", ("white", "--"), "s"),
    "random_sparse": ("random sparse", ("white", ":"), "^"),
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
            color="black",
            linestyle=color[1],
            marker=marker,
            markerfacecolor=color[0],
            markeredgecolor="black",
            markersize=3.6,
            linewidth=0.8,
            markeredgewidth=0.6,
            label=label,
            zorder=3,
        )

    limited = [item for item in measurements if item["task"] == task
               and not item.get("complete", item["seconds"] is not None)
               and not item.get("error")]
    for item in limited:
        elapsed = item.get("elapsed_seconds")
        if elapsed is None:
            continue
        reason = item.get("termination_reason", "time_limit")
        ax.scatter(item["n_species"], elapsed, marker="x", s=35,
                   color="black", zorder=5)
        ax.annotate(reason.replace("_", " "), (item["n_species"], elapsed),
                    xytext=(5, 5), textcoords="offset points", fontsize=7.5,
                    color="black")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(42, 560)
    observed = [item["seconds"] for item in measurements
                if item["task"] == task and item["seconds"] is not None]
    observed += [item["elapsed_seconds"] for item in limited if item.get("elapsed_seconds")]
    if observed:
        ax.set_ylim(max(min(observed) * 0.45, 1e-7), max(observed) * 3)
    else:
        ax.set_ylim(*ylim)
    ax.set_xticks([50, 100, 200, 400], labels=["50", "100", "200", "400"])
    ax.get_xaxis().set_minor_formatter(mpl.ticker.NullFormatter())
    ax.set_xlabel("Species $|X|$")
    ax.set_ylabel("Wall time (s)")
    panel_title(ax, letter, title)


def family_key(ax: plt.Axes) -> None:
    """Draw each benchmark family beside the line style that plots it."""
    ax.set(xlim=(0, 30), ylim=(-1.2, 1.25), aspect="equal")
    ax.axis("off")
    r = 0.22
    for x0, family in ((0.4, "chain"), (10.6, "reversible_chain"), (20.8, "random_sparse")):
        label, (face, dash), marker = FAMILIES[family]
        ax.plot([x0, x0 + 1.6], [0, 0], color="black", ls=dash, lw=0.8)
        ax.plot([x0 + 0.8], [0], marker=marker, mfc=face, mec="black", mew=0.6,
                ms=3.6, color="black")
        nx = x0 + 2.6
        if family == "random_sparse":
            # bimolecular reactions S_a + S_b -> S_c
            pts = {"a": (nx, 0.45), "b": (nx, -0.45), "c": (nx + 2.0, 0.0),
                   "d": (nx + 2.0, 0.9), "e": (nx + 4.0, 0.45)}
            for x, y in pts.values():
                place(ax, x, y, r=r)
            for ins, t, out in ((("a", "b"), (nx + 1.0, 0.0), "c"),
                                (("c", "d"), (nx + 3.0, 0.45), "e")):
                transition(ax, *t, h=0.3)
                for k in ins:
                    arc(ax, pts[k], t, shrink_a=4, shrink_b=1.5)
                arc(ax, t, pts[out], shrink_a=1.5, shrink_b=4)
            ax.text(nx + 4.5, 0.45, "...", va="center", fontsize=8)
        else:
            for i in range(4):
                place(ax, nx + 1.35 * i, 0, r=r)
            for i in range(3):
                a, b = (nx + 1.35 * i, 0), (nx + 1.35 * (i + 1), 0)
                if family == "chain":
                    arc(ax, a, b, shrink_a=4, shrink_b=4)
                else:
                    arc(ax, (a[0], 0.07), (b[0], 0.07), shrink_a=4, shrink_b=4)
                    arc(ax, (b[0], -0.07), (a[0], -0.07), shrink_a=4, shrink_b=4)
            ax.text(nx + 4.4, 0.0, "...", va="center", fontsize=8)
        ax.text(x0 + 4.6, -0.85, label, ha="center", va="top", fontsize=8,
                style="italic")


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
        left=0.11, right=0.985, top=0.94, bottom=0.13, wspace=0.28, hspace=0.38
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

    family_key(fig.add_axes([0.06, 0.0, 0.90, 0.085]))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
