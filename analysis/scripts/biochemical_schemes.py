"""Carbon-count and phosphoryl-donor illustrations, with evidence-derived outcomes.

The glyphs show molecular inventories, not atom maps or stereochemistry.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
from result_style import (
    BLUE,
    GREEN,
    ORANGE,
    RED,
    GRAY,
    INK,
    PALE,
    HAIR,
    arrow,
    configure,
    heading,
    plain_heading,
    save,
    verdict,
)


def carbon_chain(ax, x, y, n, spacing=0.18, phosphates=0, color=BLUE):
    """Draw n carbon units and zero, one or two attached phosphate symbols."""
    xx = x + (np.arange(n) - (n - 1) / 2) * spacing
    yy = y + np.array([(-1) ** i * spacing * 0.22 for i in range(n)])
    ax.plot(xx, yy, color=color, lw=1.3, zorder=3)
    ax.scatter(xx, yy, s=19, color=color, edgecolor="white", lw=0.3, zorder=4)
    for i in ([n - 1] if phosphates == 1 else [0, n - 1] if phosphates == 2 else []):
        px, py = xx[i], yy[i] + 0.36
        ax.plot([px, px], [yy[i], py], color=ORANGE, lw=1)
        ax.scatter(
            px, py, s=85, facecolor="#FFF0CF", edgecolor=ORANGE, lw=0.7, zorder=4
        )
        ax.text(
            px, py, "P", fontsize=6.5, ha="center", va="center", color=ORANGE, zorder=5
        )


def nucleotide(ax, x, y, n, size=7.5):
    """Adenosine plus its phosphate count; n=3 ATP, n=2 ADP, n=1 AMP."""
    ax.text(x, y, "Ado", fontsize=size, color=INK, ha="right", va="center")
    for i in range(n):
        px = x + 0.17 + 0.21 * i
        ax.plot([px - 0.19, px], [y, y], color=ORANGE, lw=1, zorder=1)
        ax.scatter(px, y, s=58, facecolor="#FFF0CF", edgecolor=ORANGE, lw=0.7, zorder=2)
        ax.text(px, y, "P", fontsize=5.8, ha="center", va="center", color=ORANGE)


def formose_route(ax, supplement):
    """Illustrate the actual four-firing hexose route in the audit."""
    outcome = next(
        r["outcome"]
        for r in supplement["formose_routes"]
        if r["pool"] == "C2-only" and r["target"] == "T6"
    )
    assert outcome["independent_replay_ok"]
    assert outcome["requested_flux"] == {"r_1": 2, "r_3": 1, "r_7": 1}
    ax.set(xlim=(0, 12), ylim=(0, 2.65))
    ax.axis("off")
    plain_heading(ax, "a", r"A selected route from $\mathbf{C_2}$ to $\mathbf{C_6}$", y=1.01)
    # Subscripts extend below the title; align its baseline with the panel letter.
    for text in ax.texts[-2:]:
        text.set_verticalalignment("baseline")
    # Condensed structures retain hydroxyl and carbonyl identities without
    # assigning a stereoisomer to the constitution-only generated targets.
    ax.text(0.0, 1.97, r"2 × $\mathbf{HOCH_2}$–CHO", fontsize=9, weight="bold")
    arrow(ax, (3.0, 2.07), (4.0, 2.07), ORANGE, 1.2)
    ax.annotate(
        r"$\mathrm{r}_1$ × 2", (3.5, 2.07), xytext=(0, 3),
        textcoords="offset points", fontsize=8, ha="center", va="bottom", color=ORANGE,
    )
    ax.text(4.15, 1.97, "2 × HO–CH=CH–OH", fontsize=9, weight="bold")
    ax.text(8.25, 1.97, "two enediol donors", fontsize=8, color="black")
    # Each aldol addition consumes one of the two enediols above.
    positions = [
        (1.1, 2, "Glycolaldehyde"),
        (5.6, 4, r"$\mathrm{T}_4$ aldotetrose"),
        (10.3, 6, r"$\mathrm{T}_6$ aldohexose"),
    ]
    for x, n, label in positions:
        carbon_chain(ax, x, 0.85, n, spacing=0.26)
        ax.text(x, 0.31, label, fontsize=8.5, ha="center")
    for xa, xb, rid in [(2.0, 4.7, r"$\mathrm{r}_3$"), (6.7, 9.3, r"$\mathrm{r}_7$")]:
        arrow(ax, (xa, 0.85), (xb, 0.85), GREEN, 1.4)
        ax.annotate(
            "+ one enediol", ((xa + xb) / 2, 0.85), xytext=(0, 3),
            textcoords="offset points", fontsize=8, ha="center", va="bottom", color=ORANGE,
        )
        ax.text(
            (xa + xb) / 2,
            0.38,
            f"Aldol addition · {rid}",
            fontsize=7.5,
            ha="center",
            color="black",
        )


def glycolysis_donors(out, supplement):
    configure()
    fig = plt.figure(figsize=(7.2, 5.6))
    ax = fig.add_axes([0.045, 0.36, 0.91, 0.56])
    ax.set(xlim=(0, 13.8), ylim=(0, 6.2))
    ax.axis("off")
    heading(ax, "A", "Two phosphoryl donors, two substrate requirements")
    # Shared upper pathway. Nucleotides retain their actual phosphate counts.
    for x, n, p, label in [
        (0.9, 6, 0, "Glucose"),
        (3.7, 6, 1, "G6P"),
        (6.3, 6, 1, "F6P"),
    ]:
        carbon_chain(ax, x, 3.25, n, phosphates=p)
        ax.text(x, 2.65, label, ha="center", fontsize=9, weight="bold")
    arrow(ax, (1.55, 3.25), (3.0, 3.25), BLUE, 1.5)
    ax.text(2.28, 2.82, "r10", fontsize=8, ha="center", color=GRAY)
    nucleotide(ax, 1.22, 4.85, 3)
    nucleotide(ax, 3.25, 4.85, 2)
    arrow(ax, (2.03, 4.85), (2.68, 4.85), ORANGE, 1)
    ax.text(1.45, 5.38, "ATP", fontsize=8, ha="center", color=ORANGE)
    ax.text(3.4, 5.38, "ADP", fontsize=8, ha="center", color=ORANGE)
    arrow(ax, (2.26, 4.50), (2.26, 3.58), ORANGE, 1)
    arrow(ax, (4.38, 3.25), (5.63, 3.25), BLUE, 1.5)
    ax.text(4.99, 2.82, "r15", fontsize=8, ha="center", color=GRAY)
    # Colours refer strictly to the one-ATP initial inventory.
    for yy, donor, product, rid, colour, ok in [
        (4.9, "ATP", "ADP", "r3", RED, False),
        (1.50, "ADP", "AMP", "r12", GREEN, True),
    ]:
        arrow(ax, (6.55, 3.6 if ok is False else 3.0), (7.55, yy), colour, 1.4)
        arrow(ax, (7.58, yy), (9.08, yy), colour, 1.4)
        ax.text(
            8.35,
            yy + 0.46,
            f"{donor} → {product}",
            fontsize=9,
            ha="center",
            color=colour,
            weight="bold",
        )
        ax.text(8.35, yy - 0.46, rid, fontsize=8, ha="center", color=GRAY)
        verdict(ax, 8.35, yy, ok, 16)
        carbon_chain(ax, 9.78, yy, 6, phosphates=2)
        ax.text(9.78, yy - 0.65, "FBP", ha="center", fontsize=9)
        arrow(ax, (10.47, yy), (11.37, yy), colour, 1.3)
        ax.text(10.91, yy - 0.46, "r7", fontsize=8, ha="center", color=GRAY)
        carbon_chain(ax, 12.15, yy + 0.3, 3, phosphates=1)
        carbon_chain(ax, 12.15, yy - 0.65, 3, phosphates=1)
        ax.text(12.73, yy + 0.3, "GAP", fontsize=7.5, va="center")
        ax.text(12.73, yy - 0.65, "DHAP", fontsize=7.5, va="center")
    ax.text(0.15, 1.4, "Starting with", fontsize=8, color=GRAY)
    ax.text(0.15, 0.91, "1 glucose + 1 ATP", fontsize=10, weight="bold")
    ax.text(0.15, 0.36, "ADP from r10 supplies r12", fontsize=8, color=GREEN)

    b = fig.add_axes([0.045, 0.035, 0.91, 0.24])
    b.set(xlim=(0, 13.8), ylim=(0, 3.2))
    b.axis("off")
    heading(b, "B", "Checked reaction sequences", y=1.12)
    b.text(0.15, 3.0, "Initial inventory", fontsize=9, weight="bold")
    for x, title in [(7.5, "ATP-dependent"), (11.4, "ADP-dependent")]:
        b.text(x, 3.0, title, ha="center", fontsize=9, weight="bold")
    for yy, pool, label in [
        (2.35, "one_atp", "1 glucose + 1 ATP"),
        (1.43, "two_atp", "1 glucose + 2 ATP"),
        (0.51, "added_f6p", "1 glucose + 1 ATP + 1 F6P"),
    ]:
        b.add_patch(
            Rectangle((0, yy - 0.36), 13.8, 0.72, facecolor=PALE, edgecolor="none")
        )
        b.text(0.15, yy, label, va="center", fontsize=9)
        for x, name in [(7.5, "ATP-dependent"), (11.4, "ADP-dependent")]:
            outcome = next(
                r["outcome"]
                for r in supplement["routes"]
                if r["name"] == name and r["pool"] == pool
            )
            ok = outcome["status"] == "realizable"
            assert outcome["independent_replay_ok"] if ok else outcome["exhaustive"]
            verdict(b, x - 0.7, yy, ok, 15)
            b.text(
                x - 0.3,
                yy,
                "4 steps" if ok else "Blocked",
                fontsize=8,
                va="center",
                color=GREEN if ok else RED,
            )
    save(fig, out, "glycolysis_donors")
