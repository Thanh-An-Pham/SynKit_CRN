"""Evidence-driven result plots and controlled-experiment illustrations."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import NullFormatter
from biochemical_schemes import carbon_chain, formose_route
from schematics import arc, minimal_siphons, place
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
    plain_heading,
    inset_label,
    clean,
    arrow,
    verdict,
    save,
)


def configure_arial():
    """Historical name: the result plots now share the serif print style."""
    configure()


def representation(out, report):
    configure_arial()
    rep = report["representation"]
    assert (
        rep["net_matrix_unchanged"]
        and rep["species_order_unchanged"]
        and rep["reaction_order_unchanged"]
    )
    # Only the currency-boundary panel remains. The former panels on
    # polyphosphate retention were dropped because KEGG's generic C00404
    # makes the polymer look like a retained catalyst when it is a donor.
    # Drawn at the full text width (7.2 in) with equal aspect, so places stay
    # circular and text prints at the size used in the other figures.
    W, H = 10.0, 2.8                      # data units; 0.72 in per unit
    fig = plt.figure(figsize=(7.2, 7.2 * H / W))
    c = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    c.set(xlim=(0, W), ylim=(0.05, 0.05 + H))
    c.set_aspect("equal")
    c.axis("off")
    tx, dx, dy, r = 3.0, 0.75, 0.33, 0.2  # transition x, place offsets, radius

    def past(xy, bar, extra=0.15):
        """Point beyond the bar on the line from xy, so that the 13 pt
        shrink leaves a gap at the circle but stops close to the bar."""
        vx, vy = bar[0] - xy[0], bar[1] - xy[1]
        n = (vx * vx + vy * vy) ** 0.5
        return (bar[0] + extra * vx / n, bar[1] + extra * vy / n)
    x0, x1 = 5.3, 9.2                     # search path
    c.text(tx, 2.6, "Hexokinase at step 0", ha="center", fontsize=10)
    c.text((x0 + x1) / 2, 2.6, "Search from one glucose", ha="center", fontsize=10)
    c.plot([1.6, 4.4], [2.47, 2.47], color=INK, lw=0.5)
    c.plot([4.9, 9.95], [2.47, 2.47], color=INK, lw=0.5)
    for cy, row in zip((1.75, 0.6), rep["currency_arms"]):
        kept = not row["drop_currency"]
        outcome = row["outcome"]
        ok = outcome["status"] == "realizable"
        assert outcome["initial_marking"] == {"alpha-D-Glucose": 1}
        assert outcome["requested_flux"]["R01786"] == 1
        assert ok != kept  # kept: unrealizable; removed: certified
        c.text(0.05, cy, "Cofactors\nkept" if kept else "Cofactors\nremoved",
               va="center", fontsize=9.5, linespacing=1.1)
        tr = (tx, cy)
        ins = [((tx - dx, cy + dy), 1, "Glc"), ((tx - dx, cy - dy), 0, "ATP")]
        outs = [((tx + dx + 0.2, cy + dy), "G6P"), ((tx + dx + 0.2, cy - dy), "ADP")]
        if not kept:
            ins, outs = [((tx - dx, cy), 1, "Glc")], [((tx + dx + 0.2, cy), "G6P")]
        for xy, tokens, name in ins:
            place(c, *xy, r=r, tokens=tokens)
            arrow(c, xy, past(xy, tr), GRAY, 0.8, shrink=13)
            c.text(xy[0] - r - 0.08, xy[1], name, ha="right", va="center",
                   fontsize=8.5)
        for xy, name in outs:
            place(c, *xy, r=r)
            arrow(c, past(xy, tr), xy, GRAY, 0.8, shrink=13)
            c.text(xy[0] + r + 0.08, xy[1], name, ha="left", va="center",
                   fontsize=8.5)
        c.plot(*tr, marker="|", ms=19, mew=3.2, color=INK, zorder=3)

        # State space actually reached: one node per marking on the path.
        start = dict(ms=8, mfc="white", mec=INK, mew=0.9, zorder=4)
        if not ok:
            # Caption: "exhaustive search ends after one state".
            assert outcome["exhaustive"] and outcome["max_depth_reached"] == 0
            assert outcome["explored_states"] == 1
            c.plot(x0, cy, "o", **start)
            c.text(x0, cy + 0.16, "Glc", ha="center", va="bottom", fontsize=8.5)
            # Blocked first firing: an arrow stub ending in a bar.
            c.plot([x0 + 0.12, x0 + 0.9], [cy, cy], color=GRAY, lw=0.8)
            c.plot([x0 + 0.9] * 2, [cy - 0.13, cy + 0.13], color=INK, lw=1.2)
            c.text(x0 + 0.55, cy + 0.05, "HK", ha="center", va="bottom",
                   fontsize=8.5)
            verdict(c, x0 + 1.3, cy, False, 10, vector=True)
        else:
            assert outcome["independent_replay_ok"]
            depth = len(outcome["certificate"])
            assert depth == outcome["max_depth_reached"]
            assert outcome["explored_states"] == 30  # quoted in the caption
            xs = np.linspace(x0, x1, depth + 1)
            c.plot(xs, [cy] * len(xs), color=INK, lw=0.8, zorder=2)
            c.plot(xs[1:-1], [cy] * (depth - 1), "o", ms=3.0, color=INK,
                   mew=0, zorder=3)
            c.plot(xs[0], cy, "o", **start)
            # Final marking: two pyruvate tokens (ringed dot, as in Fig. 5b).
            c.plot(xs[-1], cy, "o", ms=8, mfc="white", mec=INK, mew=0.9,
                   zorder=4)
            c.plot(xs[-1], cy, "o", ms=3.6, color=INK, mew=0, zorder=5)
            c.text(xs[-1], cy + 0.16, "2 Pyr", ha="center", va="bottom",
                   fontsize=8.5)
            c.text(xs[0], cy + 0.16, "Glc", ha="center", va="bottom",
                   fontsize=8.5)
            c.text((x0 + x1) / 2, cy - 0.28, f"{depth}-step certificate",
                   ha="center", va="top", fontsize=8.5, color=GRAY)
            verdict(c, x1 + 0.4, cy, True, 10, vector=True)
    save(fig, out, "representation_sensitivity")


def search_panel(ax):
    """Exhaustive subset testing against implication propagation (panel a).

    In the reversible chain S0 <=> S1 <=> ... every producer of a selected
    species consumes a neighbour, so one selection forces the whole chain,
    which is the only minimal siphon.
    """
    chain = [f"S{i}" for i in range(4)]
    reactions = [({a: 1}, {b: 1}) for a, b in zip(chain, chain[1:])]
    reactions += [({b: 1}, {a: 1}) for a, b in zip(chain, chain[1:])]
    assert minimal_siphons(chain, reactions) == [set(chain)]
    ax.set(xlim=(-0.3, 13.7), ylim=(-0.95, 1.38), aspect="equal")
    ax.axis("off")
    plain_heading(ax, "a", y=0.97)
    # Left: every non-empty subset of three species is a candidate.
    levels = [["0", "1", "2"], ["01", "02", "12"], ["012"]]
    xy = {}
    for row, names in enumerate(levels):
        for k, name in enumerate(names):
            xy[name] = (1.1 + (k - (len(names) - 1) / 2) * 1.15 + 1.2, -0.75 + row * 0.85)
    for name, (x, y) in xy.items():
        for other, (u, v) in xy.items():
            if len(other) == len(name) + 1 and set(name) <= set(other):
                ax.plot([x, u], [y, v], color="black", lw=0.5, zorder=1)
    for name, (x, y) in xy.items():
        place(ax, x, y, "", r=0.13, fill="black")
        ax.text(x + 0.2, y, "{" + ",".join(name) + "}", fontsize=6.5,
                va="center", ha="left")
    # Right: selecting S0 forces S1, which forces S2, and so on.
    x0 = 7.0
    for i, name in enumerate(chain):
        place(ax, x0 + 1.6 * i, -0.2, f"$S_{i}$", r=0.24, fill=PALE if i == 0 else "white",
              fontsize=8)
    for i in range(len(chain) - 1):
        a, b = (x0 + 1.6 * i, -0.2), (x0 + 1.6 * (i + 1), -0.2)
        arc(ax, (a[0], a[1] + 0.07), (b[0], b[1] + 0.07), shrink_a=8, shrink_b=8)
        arc(ax, (b[0], b[1] - 0.07), (a[0], a[1] - 0.07), shrink_a=8, shrink_b=8)
        ax.annotate("", (b[0] - 0.05, 0.12), (a[0] + 0.05, 0.12),
                    arrowprops=dict(arrowstyle="-|>", lw=0.5, ls="--", color="black",
                                    connectionstyle="arc3,rad=-0.55",
                                    shrinkA=4, shrinkB=4, mutation_scale=6))
    ax.text(x0 + 2.4, 1.15, "selecting $S_0$ forces each neighbour", ha="center",
            fontsize=7.5)


def scaling(out, review, archived):
    configure_arial()
    fig = plt.figure(figsize=(7.2, 4.5))
    search_panel(fig.add_axes([0.02, 0.70, 0.96, 0.28]))
    axes = [fig.add_axes([0.10, 0.10, 0.37, 0.55]),
            fig.add_axes([0.61, 0.10, 0.37, 0.55])]
    ax = axes[0]
    medians = {}
    for alg, label, color, marker, face, dash in (
        ("exhaustive", "Exhaustive", "black", "o", "white", "--"),
        ("closure", "Pruned", "black", "s", "black", "-"),
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
            mfc=face,
            mec="black",
            ms=4.5,
            lw=0.8,
            ls=dash,
            elinewidth=0.6,
            capsize=2.5,
            capthick=0.6,
            label=label,
            zorder=3,
        )
        for row in rows:
            ax.scatter(
                row["n_species"]
                + np.linspace(-0.20, 0.20, len(row["samples_seconds"])),
                row["samples_seconds"],
                s=5,
                color=GRAY,
                lw=0,
                zorder=2,
            )
        ax.annotate(
            label,
            (x[1], y[1]),
            xytext=(-9, 2) if alg == "exhaustive" else (0, -13),
            textcoords="offset points",
            ha="right" if alg == "exhaustive" else "center",
            fontsize=8.5,
        )
    # Shade the gap between the medians; one bracket names the final ratio.
    sizes = sorted(medians["exhaustive"])
    slow = np.array([medians["exhaustive"][s] for s in sizes])
    fast = np.array([medians["closure"][s] for s in sizes])
    ax.fill_between(sizes, fast, slow, color=PALE, lw=0, zorder=1)
    ratio = slow[-1] / fast[-1]
    ax.annotate(
        "",
        (21.55, slow[-1]),
        (21.55, fast[-1]),
        arrowprops=dict(arrowstyle="<->", color="black", lw=0.5, shrinkA=0, shrinkB=0),
    )
    ax.text(
        21.75, np.sqrt(slow[-1] * fast[-1]), f"{ratio:.2g}×", va="center", fontsize=8.5
    )
    ax.set(
        xlim=(10, 23.2),
        xticks=[11, 15, 19, 21],
        yscale="log",
        xlabel="Species",
        ylabel="Wall time (s)",
    )
    inset_label(ax, "b")
    ax = axes[1]
    ends = []
    for task, label, color, marker in (
        # Highlight conserved pools; direct labels identify every measured task.
        ("conserved_moieties", "Conserved pools", "black", "o"),
        ("canonical_form", "Canonical", "white", "D"),
        ("crnt_summary", "CRNT", "white", "^"),
        ("minimal_siphons", "Siphons", "white", "s"),
        ("rank", "Exact rank", "white", "v"),
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
            color="black",
            marker=marker,
            mfc=color,
            mec="black",
            mew=0.6,
            ms=4,
            lw=0.7,
            ls="-" if color == "black" else (0, (3, 1.5)),
        )
        ends.append([np.log10(rows[-1]["seconds"]), label, "black"])
    ax.set(xscale="log", yscale="log", xlabel="Species", ylabel="Wall time (s)")
    ax.set_xticks([51, 101, 201, 401], labels=["51", "101", "201", "401"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    inset_label(ax, "c")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.90), handlelength=2.6,
              borderpad=0.4, labelspacing=0.3, fontsize=8)
    ax.set_xlim(44, 470)
    for ax in axes:
        ax.grid(False, which="both", axis="both")
        ax.minorticks_off()
    save(fig, out, "scaling")


def formose_figure(out, evidence, supplement):
    configure_arial()
    fig = plt.figure(figsize=(7.2, 3.65))
    route = fig.add_axes([0.06, 0.62, 0.90, 0.34])
    formose_route(route, supplement)
    b = fig.add_axes([0.12, 0.14, 0.34, 0.42])
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
            s=10,
            color="black",
            lw=0,
            zorder=3,
        )
    b.set(
        xlabel="Carbon atoms in species",
        ylabel="Conservation law weight",
        xlim=(0, 11),
        ylim=(0, 11),
        xticks=[2, 4, 6, 8, 10],
        yticks=[2, 4, 6, 8, 10],
    )
    b.text(
        0.05,
        0.95,
        f"weight = carbons\nfor {sum(r['coefficient'] == r['carbons'] for r in entries)}/{len(entries)} species",
        transform=b.transAxes,
        va="top",
        fontsize=8,
    )
    b.minorticks_off()
    plain_heading(b, "b", y=1.02)

    # (c) Initial inventories in precursor space.  Dashed lines join inventories
    # of equal carbon; each inventory carries a T4 (left) and T6 (right) mark.
    c = fig.add_axes([0.64, 0.14, 0.32, 0.42])
    for carbon in (3, 4, 6):
        xs = np.array([0, carbon])
        c.plot(xs, (carbon - xs) / 2, ls="--", lw=0.5, color="black", zorder=1)
        lx = {3: 0.25, 4: 0.25, 6: 1.0}[carbon]
        c.text(
            lx,
            (carbon - lx) / 2 + 0.08,
            f"{carbon} C",
            fontsize=7.5,
            ha="left",
            va="bottom",
        )
    rows = supplement["formose_routes"]
    for pool in dict.fromkeys(row["pool"] for row in rows):
        by_target = {r["target"]: r for r in rows if r["pool"] == pool}
        marking = by_target["T4"]["outcome"]["initial_marking"]
        x, y = marking["C=O"], marking["O=CCO"]
        for dx, target in ((-0.17, "T4"), (0.17, "T6")):
            outcome = by_target[target]["outcome"]
            ok = outcome["status"] == "realizable"
            assert outcome["independent_replay_ok"] if ok else outcome["exhaustive"]
            c.scatter(
                x + dx,
                y,
                s=26,
                marker="s",
                zorder=3,
                facecolor="black" if ok else "white",
                edgecolor="black",
                lw=0.6,
            )
    c.set(
        xlim=(-0.6, 4.8),
        ylim=(-0.2, 3.75),
        xticks=[0, 1, 2, 3, 4],
        yticks=[0, 1, 2, 3],
        xlabel="Formaldehyde (C$_1$)",
        ylabel="Glycolaldehyde (C$_2$)",
    )
    # Name the two marks once, on the inventory that runs both routes.
    for dx, name, ha in ((-0.1, "T$_4$", "right"), (0.1, "T$_6$", "left")):
        c.text(dx, 3.24, name, ha=ha, va="bottom", fontsize=7.5)
    c.minorticks_off()
    plain_heading(c, "c", y=1.02)
    save(fig, out, "formose")
