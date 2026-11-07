"""Line-art Petri-net glyphs for the explanatory panels of the result figures.

Species are open circles, reactions black bars, and arcs thin black arrows,
as in the plain black-and-white print style of result_style.py.
"""

from itertools import combinations, product

import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle

from result_style import PALE


def place(ax, x, y, label="", r=0.17, fill="white", lw=0.7, fontsize=8, tokens=0):
    """Species place; a gray fill marks a highlighted set."""
    ax.add_patch(Circle((x, y), r, facecolor=fill, edgecolor="black", lw=lw, zorder=3))
    if label:
        ax.text(x, y, label, ha="center", va="center", fontsize=fontsize, zorder=4)
    for dx in (np.linspace(-0.05, 0.05, tokens) if tokens > 1 else [0] * tokens):
        ax.add_patch(Circle((x + dx, y), 0.025, color="black", zorder=4))


def transition(ax, x, y, h=0.34, w=0.055, angle=0):
    """Reaction bar."""
    if angle:
        w, h = h, w
    ax.add_patch(Rectangle((x - w / 2, y - h / 2), w, h, color="black", zorder=3))


def arc(ax, a, b, rad=0.0, shrink_a=6, shrink_b=6, lw=0.6, style="-|>", ls="-"):
    """Arc between two points, shortened at both ends in points."""
    ax.add_patch(
        FancyArrowPatch(
            a, b, arrowstyle=style, mutation_scale=6, lw=lw, color="black",
            linestyle=ls, connectionstyle=f"arc3,rad={rad}",
            shrinkA=shrink_a, shrinkB=shrink_b, zorder=2,
        )
    )


def carbon_glyph(ax, x, y, n, phosphates=0, spacing=0.12, size=7):
    """Zigzag of n carbon dots, with circled phosphates on the chain ends."""
    xx = x + (np.arange(n) - (n - 1) / 2) * spacing
    yy = y + np.array([(-1) ** i * spacing * 0.25 for i in range(n)])
    ax.plot(xx, yy, color="black", lw=0.7, zorder=3)
    ax.scatter(xx, yy, s=size, color="black", lw=0, zorder=4)
    ends = [n - 1] if phosphates == 1 else [0, n - 1] if phosphates == 2 else []
    for i in ends:
        ax.plot([xx[i], xx[i]], [yy[i], yy[i] + 0.16], color="black", lw=0.5)
        ax.scatter(xx[i], yy[i] + 0.2, s=30, facecolor="white", edgecolor="black",
                   lw=0.5, zorder=4)
        ax.text(xx[i], yy[i] + 0.2, "P", fontsize=4.5, ha="center", va="center",
                zorder=5)


# ---------------------------------------------------------------------------
# Brute-force structure of a small net, used to check what a panel draws.
# ---------------------------------------------------------------------------

def minimal_siphons(species, reactions):
    """All inclusion-minimal siphons by enumeration (small nets only)."""
    def is_siphon(subset):
        return all(
            any(reactants.get(x, 0) for x in subset)
            for reactants, products in reactions
            if any(products.get(x, 0) for x in subset)
        )

    found = []
    for k in range(1, len(species) + 1):
        for subset in combinations(species, k):
            if is_siphon(subset) and not any(set(f) <= set(subset) for f in found):
                found.append(subset)
    return [set(s) for s in found]


def minimal_pools(species, reactions, bound=2):
    """Supports of minimal non-negative conservation vectors (small nets only)."""
    def conserved(y):
        return all(
            sum(w * (products.get(x, 0) - reactants.get(x, 0))
                for w, x in zip(y, species)) == 0
            for reactants, products in reactions
        )

    supports = []
    for y in product(range(bound + 1), repeat=len(species)):
        if any(y) and conserved(y):
            support = {x for w, x in zip(y, species) if w}
            if support not in supports:
                supports.append(support)
    return [s for s in supports if not any(t < s for t in supports)]


ENZYME_SPECIES = ("E", "S", "ES", "P")
ENZYME_REACTIONS = (
    ({"E": 1, "S": 1}, {"ES": 1}),
    ({"ES": 1}, {"E": 1, "S": 1}),
    ({"ES": 1}, {"E": 1, "P": 1}),
)


def enzyme_net(ax, x0, y0, highlight=(), scale=1.0):
    """E + S <=> ES -> E + P as a Petri net, highlighted places shaded."""
    s = scale
    pos = {
        "E": (x0 + 0.0 * s, y0 + 0.62 * s),
        "S": (x0 + 0.0 * s, y0 - 0.62 * s),
        "ES": (x0 + 2.0 * s, y0),
        "P": (x0 + 3.9 * s, y0),
    }
    t1, t2, t3 = (x0 + 1.0 * s, y0 + 0.32 * s), (x0 + 1.0 * s, y0 - 0.32 * s), \
        (x0 + 3.0 * s, y0)
    for name, (x, y) in pos.items():
        place(ax, x, y, name, r=0.2 * s, fill=PALE if name in highlight else "white",
              lw=1.0 if name in highlight else 0.7, fontsize=7.5)
    for t in (t1, t2, t3):
        transition(ax, *t, h=0.3 * s)
    arc(ax, pos["E"], t1, shrink_b=2)
    arc(ax, pos["S"], t1, shrink_b=2)
    arc(ax, t1, pos["ES"], shrink_a=2)
    arc(ax, pos["ES"], t2, shrink_b=2)
    arc(ax, t2, pos["E"], shrink_a=2)
    arc(ax, t2, pos["S"], shrink_a=2)
    arc(ax, pos["ES"], t3, shrink_b=2)
    arc(ax, t3, pos["P"], shrink_a=2)
    arc(ax, t3, pos["E"], rad=0.32, shrink_a=2)
    return pos
