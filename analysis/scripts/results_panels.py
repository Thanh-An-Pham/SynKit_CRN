"""Evidence-driven result plots and controlled-experiment illustrations."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Circle, Rectangle, FancyBboxPatch
from matplotlib.ticker import NullFormatter
from biochemical_schemes import carbon_chain, formose_route
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
    plain_heading,
    clean,
    arrow,
    verdict,
    save,
)


def configure_arial():
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


def representation(out, report):
    configure_arial()
    rep = report["representation"]
    assert (
        rep["net_matrix_unchanged"]
        and rep["species_order_unchanged"]
        and rep["reaction_order_unchanged"]
    )
    arms = {(r["arm"], r["catalyst_tokens"]): r for r in rep["catalyst_arms"]}
    fig = plt.figure(figsize=(7.2, 4.9))
    a = fig.add_axes([0.045, 0.55, 0.91, 0.37])
    a.set(xlim=(0, 10), ylim=(0, 3))
    a.axis("off")
    plain_heading(a, "a", "Glucose phosphorylation with a shared PolyP species")
    # Display the actual reaction equations; orange participant is retained twice.
    for xx, key, title in ((0, "intact", "Retained"), (5.4, "cancelled", "Cancelled")):
        a.text(xx + 2.18, 2.6, title, ha="center", weight="bold", fontsize=10)
        carbon_chain(a, xx + 0.63, 1.94, 6, spacing=0.18)
        carbon_chain(a, xx + 3.85, 1.94, 6, spacing=0.18, phosphates=1)
        a.text(xx + 0.63, 1.52, "Glc", fontsize=9, ha="center")
        a.text(xx + 3.85, 1.52, "G6P", fontsize=9, ha="center")
        arrow(a, (xx + 1.35, 1.86), (xx + 3.0, 1.86), BLUE, 1.4)
        if key == "intact":
            for px in (xx + 0.62, xx + 3.82):
                a.text(px, 1.23, "+ PolyP", color=ORANGE, fontsize=10, ha="center")
            a.text(xx + 2.18, 2.07, "R02189", color="black", fontsize=8, ha="center")
        else:
            a.text(
                xx + 2.18,
                1.20,
                "PolyP incidence removed",
                color="black",
                fontsize=8,
                ha="center",
            )
        row = arms[(key, 0)]
        a.text(
            xx + 2.18,
            0.53,
            f"{row['n_complexes']} complexes   ·   {row['n_linkage_classes']} classes   ·   δ = {row['deficiency']}",
            ha="center",
            fontsize=9,
        )
    # Anchor the separator halfway between the two diagram column centres.
    midpoint = 2.18 + 5.4 / 2
    a.text(midpoint - 0.24, -0.02,
           r"$S_{\mathrm{retained}}=S_{\mathrm{cancelled}}$",
           ha="right", va="baseline", fontsize=9)
    a.text(midpoint, -0.02, r"$\bullet$",
           ha="center", va="baseline", fontsize=9)
    a.text(midpoint + 0.24, -0.02, "14 species × 15 reactions",
           ha="left", va="baseline", fontsize=9, weight="bold")

    b = fig.add_axes([0.045, 0.13, 0.43, 0.28])
    b.set(xlim=(0, 4), ylim=(0, 2.2))
    b.axis("off")
    plain_heading(b, "b", "Single-reaction enabling", y=1.16)
    b.text(2.8, 2.28, "Initial PolyP tokens", ha="center", fontsize=8)
    for x, tokens in ((2.25, 0), (3.35, 1)):
        b.text(x, 1.94, str(tokens), ha="center", fontsize=9, weight="bold")
    for y, key, title in (
        (1.23, "intact", "Retained"),
        (0.42, "cancelled", "Cancelled"),
    ):
        b.text(0.03, y, title, va="center", fontsize=9)
        for x, tokens in ((2.25, 0), (3.35, 1)):
            ok = arms[(key, tokens)]["outcome"]["status"] == "realizable"
            b.add_patch(
                Rectangle(
                    (x - 0.40, y - 0.30),
                    0.8,
                    0.6,
                    facecolor="#E5F2ED" if ok else "#FAEAE5",
                    edgecolor="none",
                )
            )
            verdict(b, x, y, ok, 16, vector=True)

    c = fig.add_axes([0.565, 0.13, 0.39, 0.28])
    c.set(xlim=(0, 16), ylim=(0, 2.2))
    c.axis("off")
    plain_heading(c, "c", "Full glycolytic flux", y=1.16)
    c.text(0, 2.23, "Initial inventory: one glucose", fontsize=8, color="black")
    for y, row in zip((1.25, 0.35), rep["currency_arms"]):
        reduced = row["drop_currency"]
        ok = row["outcome"]["status"] == "realizable"
        c.text(
            0,
            y + 0.34,
            f"Currency metabolites {'removed' if reduced else 'retained'} · {row['n_species']}",
            fontsize=8,
        )
        # Indent each consequence beneath its condition and connect the two.
        c.plot([0.25, 0.25], [y + 0.22, y], color=GRAY, lw=0.8)
        arrow(c, (0.25, y), (0.95, y), GRAY, 0.8)
        n = len(row["outcome"]["certificate"] or [])
        if ok:
            for step in range(n):
                c.add_patch(
                    Rectangle(
                        (1.25 + step * 0.8, y - 0.11),
                        0.58,
                        0.22,
                        facecolor=GREEN,
                        edgecolor="none",
                    )
                )
            c.text(13.6, y, str(n), color="black", va="center", weight="bold")
        else:
            verdict(c, 1.54, y, False, 15, vector=True)
            c.text(2.5, y, "No first firing", color="black", va="center", fontsize=8)
    save(fig, out, "representation_sensitivity")


def scaling(out, review, archived):
    configure_arial()
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.0))
    fig.subplots_adjust(left=0.10, right=0.97, bottom=0.20, top=0.76, wspace=0.40)
    ax = axes[0]
    medians = {}
    for alg, label, color, marker in (
        ("exhaustive", "Exhaustive", RED, "o"),
        ("closure", "Closure", GREEN, "s"),
    ):
        rows = sorted(
            (r for r in review["paired_siphon_timings"] if r["algorithm"] == alg),
            key=lambda r: r["n_species"],
        )
        x = np.array([r["n_species"] for r in rows])
        y = np.array([r["median_seconds"] for r in rows])
        medians[alg] = dict(zip(x.tolist(), y.tolist()))
        errors = np.array(
            [
                [r["median_seconds"] - r["q1_seconds"] for r in rows],
                [r["q3_seconds"] - r["median_seconds"] for r in rows],
            ]
        )
        ax.errorbar(
            x,
            y,
            yerr=errors,
            color=color,
            marker=marker,
            ms=5,
            lw=1.4,
            capsize=4,
            label=label,
            zorder=3,
        )
        for row in rows:
            ax.scatter(
                row["n_species"]
                + np.linspace(-0.20, 0.20, len(row["samples_seconds"])),
                row["samples_seconds"],
                s=9,
                alpha=0.55,
                color=color,
                zorder=4,
            )
        last = rows[-1]
        labelval = (
            f"{last['median_seconds']:.2f} s"
            if alg == "exhaustive"
            else f"{last['median_seconds']*1000:.2f} ms"
        )
        ax.annotate(
            labelval,
            (x[-1], y[-1]),
            xytext=(-5, 10 if alg == "exhaustive" else -17),
            textcoords="offset points",
            ha="right",
            color="black",
            fontsize=9,
            weight="bold",
        )
    # The gap between the two curves is the quantity of interest, so name it
    # once per size instead of leaving readers to subtract log decades.
    for size in sorted(medians["exhaustive"]):
        slow = medians["exhaustive"][size]
        fast = medians["closure"][size]
        ratio = slow / fast
        digits = 0 if ratio < 100 else -(len(f"{ratio:.0f}") - 2)
        ax.text(
            size,
            np.sqrt(slow * fast),
            f"{round(ratio, digits):.0f}×",
            ha="center",
            va="center",
            fontsize=7.5,
            color="black",
            bbox=dict(facecolor="white", edgecolor="none", pad=0.8),
            zorder=5,
        )
    ax.set(
        xlim=(10, 22),
        xticks=[11, 15, 19, 21],
        yscale="log",
        xlabel="Species",
        ylabel="Wall time (s)",
    )
    plain_heading(ax, "a", "Paired siphon search", y=1.27)
    ax.legend(
        frameon=False,
        fontsize=8,
        ncol=2,
        loc="lower left",
        bbox_to_anchor=(-0.02, 1.01),
        handlelength=1.3,
        columnspacing=1,
    )
    ax = axes[1]
    for task, label, color, marker in (
        ("conserved_moieties", "Moieties", RED, "o"),
        ("canonical_form", "Canonical", ORANGE, "D"),
        ("crnt_summary", "CRNT", "#9063A5", "^"),
        ("minimal_siphons", "Siphons", GREEN, "s"),
        ("rank", "Num. rank", BLUE, "P"),
    ):
        rows = sorted(
            (
                r
                for r in archived["measurements"]
                if r["family"] == "reversible_chain"
                and r["task"] == task
                and r["seconds"] is not None
            ),
            key=lambda r: r["n_species"],
        )
        ax.plot(
            [r["n_species"] for r in rows],
            [r["seconds"] for r in rows],
            label=label,
            color=color,
            marker=marker,
            ms=4,
            lw=1.2,
        )
    ax.set(xscale="log", yscale="log", xlabel="Species", ylabel="Wall time (s)")
    ax.set_xticks([51, 101, 201, 401], labels=["51", "101", "201", "401"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    plain_heading(ax, "b", "Larger reversible chains", y=1.27)
    ax.legend(
        frameon=False,
        fontsize=7.6,
        ncol=3,
        loc="lower left",
        bbox_to_anchor=(-0.02, 1.0),
        handlelength=1,
        columnspacing=0.7,
    )
    for ax in axes:
        ax.grid(False, which="both", axis="both")
        ax.minorticks_off()
    save(fig, out, "scaling")


def formose_figure(out, evidence, supplement):
    configure_arial()
    fig = plt.figure(figsize=(7.2, 6.8))
    route = fig.add_axes([0.095, 0.75, 0.855, 0.19])
    formose_route(route, supplement)
    a = fig.add_axes([0.095, 0.405, 0.36, 0.28])
    b = fig.add_axes([0.59, 0.405, 0.36, 0.28])
    growth = evidence["growth"]
    for key, label, color, marker in [
        ("n_species", "Species", BLUE, "o"),
        ("n_reactions", "Reactions", ORANGE, "s"),
    ]:
        x = [r["iteration"] for r in growth]
        y = [r[key] for r in growth]
        a.plot(x, y, color=color, marker=marker, lw=1.4, label=label)
        a.annotate(
            str(y[-1]),
            (x[-1], y[-1]),
            xytext=(6, 2 if key == "n_reactions" else -9),
            textcoords="offset points",
            color="black",
            weight="bold",
            fontsize=10,
        )
    a.set(
        xlabel="Rule-expansion iteration",
        ylabel="Count",
        xticks=[1, 2, 3, 4],
        xlim=(0.8, 4.55),
        ylim=(0, 51),
    )
    a.legend(frameon=False, fontsize=8, loc="upper left")
    a.grid(False)
    a.minorticks_off()
    plain_heading(a, "b", "Network growth")
    entries = evidence["carbon_conservation"]["entries"]
    # A separate small point for each species avoids hiding multiplicities.
    # Horizontal offsets separate coincident observations; vertical values are exact.
    b.plot([0, 11], [0, 11], ls="--", lw=0.9, color=GRAY, zorder=1)
    for carbon in sorted(set(row["carbons"] for row in entries)):
        rows = [row for row in entries if row["carbons"] == carbon]
        offsets = np.linspace(-0.24, 0.24, len(rows)) if len(rows) > 1 else [0]
        b.scatter(
            carbon + np.array(offsets),
            [row["coefficient"] for row in rows],
            s=15,
            color=GREEN,
            edgecolor="white",
            lw=0.35,
            zorder=3,
        )
    b.set(
        xlabel="Carbon atoms per molecule",
        ylabel="Integer conservation coefficient",
        xlim=(0, 11),
        ylim=(0, 11),
        xticks=[2, 4, 6, 8, 10],
        yticks=[2, 4, 6, 8, 10],
    )
    b.text(
        0.04,
        0.95,
        f"{sum(r['coefficient'] == r['carbons'] for r in entries)}/{len(entries)} exact matches",
        transform=b.transAxes,
        va="top",
        fontsize=9,
        color=GREEN,
        weight="bold",
    )
    b.grid(False)
    b.minorticks_off()
    plain_heading(b, "c", "Carbon conservation")
    c = fig.add_axes([0.095, 0.065, 0.855, 0.23])
    c.set(xlim=(0, 10), ylim=(0, 4.2))
    c.axis("off")
    plain_heading(c, "d", "Equal carbon, different precursor supply", y=1.04)
    c.text(1.8, 3.63, "Initial inventory", fontsize=9, weight="bold", ha="center")
    c.text(5.0, 3.63, "Total C", fontsize=9, weight="bold", ha="center")
    c.text(7.0, 3.63, "Tetrose route", fontsize=9, weight="bold", ha="center")
    c.text(9.1, 3.63, "Hexose route", fontsize=9, weight="bold", ha="center")
    rows = supplement["formose_routes"]
    pools = list(dict.fromkeys(row["pool"] for row in rows))
    for i, pool in enumerate(pools):
        y = 3.03 - i * 0.59
        by_target = {r["target"]: r for r in rows if r["pool"] == pool}
        marking = by_target["T4"]["outcome"]["initial_marking"]
        c1 = marking["C=O"]
        c2 = marking["O=CCO"]
        if i % 2 == 0:
            c.add_patch(
                Rectangle(
                    (-0.2, y - 0.26), 10.3, 0.52, facecolor=PALE, edgecolor="none"
                )
            )
        # One dot = formaldehyde; a joined pair = glycolaldehyde. Counts are molecules.
        for j in range(c1):
            c.scatter(0.18 + j * 0.29, y, s=24, color=ORANGE, zorder=3)
        for j in range(c2):
            x = 1.80 + j * 0.54
            c.plot([x, x + 0.16], [y, y], color=BLUE, lw=1.5)
            c.scatter([x, x + 0.16], [y, y], s=24, color=BLUE, zorder=3)
        c.text(3.8, y, f"{c1} + {c2}", fontsize=8, ha="center", va="center", color="black")
        c.text(5.0, y, str(c1 + 2 * c2), ha="center", va="center", fontsize=9)
        for x, target in ((7.0, "T4"), (9.1, "T6")):
            outcome = by_target[target]["outcome"]
            ok = outcome["status"] == "realizable"
            assert outcome["independent_replay_ok"] if ok else outcome["exhaustive"]
            verdict(c, x, y, ok, 15, vector=True)
    save(fig, out, "formose")
