#!/usr/bin/env python3
"""Illustrate the SynKit CRN-submodule workflow with editable vector primitives."""

from pathlib import Path
import argparse
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from result_style import use_sans
from matplotlib.transforms import Affine2D
from matplotlib.patches import (
    Circle,
    Rectangle,
    FancyBboxPatch,
    FancyArrowPatch,
    Ellipse,
    Polygon,
    Wedge,
    RegularPolygon,
)

INK = "#193342"
BLUE = "#0072B2"
GREEN = "#008564"
ORANGE = "#D89915"
MUTED = "#81939D"
PALE = "#F1F6F8"
# Color roles shared by all figures: ink and gray for structure and arrows,
# blue for species, green only for a verified result, orange for the enzyme
# (here) and for phosphorus (elsewhere).
LIGHT_BLUE = "#9CC5E3"


def render_visual_abstract(out):
    use_sans()
    plt.rcParams.update(
        {
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    # Lower the central diagram to the middle connector row, preserving its
    # internal layout and leaving room for the derived views and key below it.
    connector_rows = (6.95, 4.49, 2.2)
    connector_y = connector_rows[1]
    central_offset_y = connector_y - 5.65
    fig, ax = plt.subplots(figsize=(7.2, 4.65))
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    ax.set(xlim=(0, 16.5), ylim=(central_offset_y + 0.95, 8.55), aspect="equal")
    ax.axis("off")

    def label(x, y, value, size=10, color=INK, weight="normal", ha="center"):
        ax.text(
            x, y, value, ha=ha, va="center", fontsize=size, color=color, weight=weight
        )

    def line(a, b, color=BLUE, lw=1.4, **kw):
        ax.plot(
            [a[0], b[0]], [a[1], b[1]], color=color, lw=lw, solid_capstyle="round", **kw
        )

    def arrow(a, b, color=MUTED, lw=1.2, rad=0, shrink=0, shrink_end=None):
        ax.add_patch(
            FancyArrowPatch(
                a,
                b,
                arrowstyle="-|>",
                mutation_scale=9,
                linewidth=lw,
                color=color,
                connectionstyle=f"arc3,rad={rad}",
                shrinkA=shrink,
                shrinkB=shrink if shrink_end is None else shrink_end,
            )
        )

    def species(x, y, color=BLUE, r=0.24, token=False):
        ax.add_patch(
            Circle((x, y), r, facecolor="white", edgecolor=color, lw=1.4, zorder=5)
        )
        if token:
            ax.add_patch(
                Circle((x, y), r * 0.25, facecolor=INK, edgecolor="none", zorder=6)
            )

    def reaction(x, y, width=0.22, height=0.35):
        ax.add_patch(
            Rectangle(
                (x - width / 2, y - height / 2),
                width,
                height,
                facecolor=INK,
                edgecolor="none",
                zorder=5,
            )
        )

    def molecule(points, edges, color=BLUE, r=0.09):
        for i, j in edges:
            line(points[i], points[j], color, 1.8)
        for p in points:
            ax.add_patch(
                Circle(p, r, facecolor="white", edgecolor=color, lw=1.2, zorder=5)
            )

    def tick(x, y, color=GREEN):
        line((x - 0.09, y), (x - 0.015, y - 0.085), color, 1.7)
        line((x - 0.015, y - 0.085), (x + 0.14, y + 0.13), color, 1.7)

    # Three columns read left to right under plain bold headers.
    for x0, x1, head in (
        (0.25, 3.05, "Build"),
        (4.55, 11.05, "Represent"),
        (12.3, 16.35, "Analyze"),
    ):
        label((x0 + x1) / 2, 8.10, head, 11, INK, "bold")

    # Curated pathway database: a cylinder containing a branching pathway.
    x, y = 1.75, 6.95
    ax.add_patch(
        Rectangle(
            (x - 0.75, y - 0.45), 1.5, 0.9, facecolor=PALE, edgecolor=BLUE, lw=1.1
        )
    )
    ax.add_patch(
        Ellipse((x, y - 0.45), 1.5, 0.35, facecolor=PALE, edgecolor=BLUE, lw=1.1)
    )
    ax.add_patch(
        Rectangle((x - 0.73, y - 0.39), 1.46, 0.2, facecolor=PALE, edgecolor="none")
    )
    ax.add_patch(
        Ellipse((x, y + 0.45), 1.5, 0.35, facecolor="white", edgecolor=BLUE, lw=1.1)
    )
    for a, b in [
        ((1.28, 6.91), (1.77, 6.91)),
        ((1.77, 6.91), (2.19, 7.12)),
        ((1.77, 6.91), (2.19, 6.7)),
    ]:
        arrow(a, b, BLUE, lw=0.9, shrink=3)
    for px, py in [(1.28, 6.91), (1.77, 6.91), (2.19, 7.12), (2.19, 6.7)]:
        species(px, py, r=0.075)
    label(x, 5.96, "KEGG", 10, INK, "bold")

    # Structured model file: folded page and a small reaction inside it.
    x, y = 1.75, 4.49
    ax.add_patch(
        Polygon(
            [
                (x - 0.59, y - 0.56),
                (x + 0.59, y - 0.56),
                (x + 0.59, y + 0.31),
                (x + 0.32, y + 0.6),
                (x - 0.59, y + 0.6),
            ],
            closed=True,
            facecolor=PALE,
            edgecolor=BLUE,
            lw=1.1,
        )
    )
    line((x + 0.32, y + 0.6), (x + 0.32, y + 0.31), BLUE, 1)
    line((x + 0.32, y + 0.31), (x + 0.59, y + 0.31), BLUE, 1)
    for yy, xx in [(y + 0.25, x + 0.10), (y + 0.05, x + 0.29)]:
        line((x - 0.34, yy), (xx, yy), MUTED, 0.8)
    species(x - 0.3, y - 0.27, r=0.1)
    species(x + 0.32, y - 0.27, r=0.1)
    arrow((x - 0.16, y - 0.27), (x + 0.17, y - 0.27), BLUE, 0.9)
    label(x, 3.58, "SBML", 10, INK, "bold")

    # Graph rewrite: a bond is added to a molecular fragment.
    molecule([(0.63, 2.1), (0.94, 2.4), (1.25, 2.1)], [(0, 1), (1, 2)])
    arrow((1.43, 2.2), (1.91, 2.2), MUTED, 1)
    molecule(
        [(2.07, 2.1), (2.38, 2.4), (2.69, 2.1), (2.4, 1.83)], [(0, 1), (1, 2), (1, 3)]
    )
    ax.add_patch(
        Circle(
            (2.4, 1.83), 0.1, facecolor=BLUE, edgecolor=BLUE, lw=1.2, zorder=6
        )
    )
    label(x, 1.31, "Rules", 10, INK, "bold")

    # Three grey input branches join a vertical stem; only the combined
    # input is a blue arrow, aligned with the middle branch.
    for yy in connector_rows:
        line((3.12, yy), (4.03, yy), MUTED, 1)
    line((4.03, min(connector_rows)), (4.03, max(connector_rows)), MUTED, 1)
    arrow((4.03, connector_y), (4.66, connector_y), BLUE, 1.6)

    # Enzyme binding and product release: E+S <-> ES -> E+P.
    # Three transition nodes preserve the actual reactant/product incidence.
    existing_artists = set(ax.get_children())

    def enzyme(x, y, bound=False, r=0.39):
        ax.add_patch(
            Wedge(
                (x, y),
                r,
                35,
                325,
                facecolor="#FCE3A6",
                edgecolor=ORANGE,
                lw=1.4,
                zorder=5,
            )
        )
        if bound:
            ax.add_patch(
                RegularPolygon(
                    (x + r * 0.58, y),
                    6,
                    radius=r * 0.40,
                    orientation=np.pi / 6,
                    facecolor=BLUE,
                    edgecolor="white",
                    lw=0.7,
                    zorder=6,
                )
            )

    # Binding: E and S must both be present. Release regenerates E.
    positions = {
        "S": (5.5, 5.65),
        "E": (8.0, 7.06),
        "ES": (8.3, 5.65),
        "P": (10.42, 5.65),
    }
    bind, unbind, release = (6.65, 5.65), (6.65, 4.65), (9.5, 5.65)
    for a, b in [
        (positions["S"], bind),
        (positions["E"], bind),
        (bind, positions["ES"]),
        (positions["ES"], release),
        (release, positions["P"]),
        (release, positions["E"]),
        (positions["ES"], unbind),
        (unbind, positions["S"]),
    ]:
        arrow(
            a,
            b,
            MUTED,
            1.3,
            shrink=9,
            # Keep the incoming arrowhead outside the complex's circular edge.
            shrink_end=14 if b == positions["ES"] else None,
        )
    arrow(unbind, positions["E"], MUTED, 1.1, rad=-0.75, shrink=9)
    enzyme(*positions["E"])
    enzyme(*positions["ES"], bound=True)
    ax.add_patch(
        RegularPolygon(
            positions["S"],
            6,
            radius=0.23,
            orientation=np.pi / 6,
            facecolor=BLUE,
            edgecolor="white",
            zorder=6,
        )
    )
    ax.add_patch(
        RegularPolygon(
            positions["P"], 5, radius=0.23, facecolor=BLUE, edgecolor="white", zorder=6
        )
    )
    for pos in (bind, unbind, release):
        reaction(*pos)
    label(8.0, 7.60, "Enzyme", 8.5, ORANGE)
    label(5.50, 5.03, "Substrate", 8)
    label(positions["ES"][0] + 0.15, 5.04, "Complex", 8)
    label(10.40, 5.03, "Product", 8)
    label(8.0, 4.20, r"$E+S\rightleftharpoons ES\rightarrow E+P$", 10)

    # Three derived views are rendered as diagrams, not explanatory text cards.
    arrow((7.95, 3.84), (7.95, 3.3), BLUE, 1.1)
    line((5.8, 3.3), (10.05, 3.3), MUTED, 0.85)
    for xx in (5.8, 7.95, 10.05):
        arrow((xx, 3.3), (xx, 2.99), MUTED, 0.85)
    net = np.array([[-1, 1, 0], [-1, 1, 1], [1, -1, -1], [0, 0, 1]])
    for i in range(4):
        for j in range(3):
            ax.add_patch(
                Rectangle(
                    (5.46 + j * 0.22, 2.69 - i * 0.23),
                    0.18,
                    0.18,
                    facecolor=(
                        BLUE if net[i, j] < 0 else LIGHT_BLUE if net[i, j] > 0 else "#DFE7EB"
                    ),
                    edgecolor="none",
                )
            )
    label(5.79, 1.31, "Stoichiometry", 9)
    # One directed hyperedge: the tail set {E, S} is consumed jointly and the
    # head ES is produced, so the two tails meet before the single arrowhead.
    tails, junction, head = [(7.38, 2.82), (7.38, 1.98)], (7.98, 2.4), (8.62, 2.4)
    ax.add_patch(Ellipse((7.38, 2.4), 0.62, 1.36, facecolor="none",
                         edgecolor=MUTED, lw=0.8, ls=(0, (2, 1.5))))
    for a, rad in zip(tails, (0.3, -0.3)):
        ax.add_patch(FancyArrowPatch(a, junction, arrowstyle="-", lw=1.4,
                                     color=MUTED, shrinkA=4, shrinkB=0,
                                     connectionstyle=f"arc3,rad={rad}"))
    ax.add_patch(Circle(junction, 0.045, facecolor=MUTED, edgecolor="none", zorder=6))
    arrow(junction, head, MUTED, 1.4, shrink_end=5)
    species(*tails[0], color=ORANGE, r=0.13)
    species(*tails[1], color=BLUE, r=0.13)
    species(*head, color=BLUE, r=0.13)
    label(7.94, 1.31, "Hypergraph", 9)
    species(9.46, 2.4, token=True, r=0.18)
    reaction(10.04, 2.4, 0.12, 0.47)
    species(10.62, 2.4, r=0.18)
    arrow((9.68, 2.4), (9.94, 2.4), BLUE, 0.95)
    arrow((10.16, 2.4), (10.39, 2.4), BLUE, 0.95)
    label(10.05, 1.31, "Petri net", 9)

    # Translate the whole central group in data coordinates. The substrate,
    # complex, product, and forward reactions now share connector_y exactly.
    central_transform = Affine2D().translate(0, central_offset_y) + ax.transData
    for artist in ax.get_children():
        if artist not in existing_artists:
            artist.set_transform(central_transform)

    # A common output stem fans out to three distinct mathematical questions.
    line((11.14, connector_y), (11.62, connector_y), MUTED, 1)
    line((11.62, min(connector_rows)), (11.62, max(connector_rows)), MUTED, 1)
    for yy in connector_rows:
        arrow((11.62, yy), (12.27, yy), MUTED, 1.3)

    # A conserved enzyme pool redistributes between free and bound forms.
    enzyme(13.02, 7.02, r=0.31)
    label(13.72, 7.02, "+", 14, INK)
    enzyme(14.51, 7.02, bound=True, r=0.31)
    label(14.0, 6.34, r"$E+ES=E_{\mathrm{total}}$", 10, INK)
    label(14.0, 5.96, "Conserved pools", 10, INK, "bold")

    # CRNT: the complex graph of this example, E+S <-> ES -> E+P, one linkage
    # class, rank two, so the deficiency evaluates to zero.
    for xx, name in ((12.82, "E+S"), (14.0, "ES"), (15.18, "E+P")):
        label(xx, 4.68, name, 8.5, INK)
    arrow((13.2, 4.76), (13.72, 4.76), MUTED, 1.1)
    arrow((13.72, 4.60), (13.2, 4.60), MUTED, 1.1)
    arrow((14.28, 4.68), (14.8, 4.68), MUTED, 1.1)
    label(14.0, 3.99, r"$\delta=3-1-2=0$", 10.5, INK)
    label(14.0, 3.58, "CRNT tests", 10, INK, "bold")

    # A firing sequence from one E and one S: bind, then release.  Tokens keep
    # the species colours of the central scheme, and one check marks replay.
    markings = [(ORANGE, BLUE), ("ES",), (ORANGE, BLUE)]
    for i, (xx, tokens) in enumerate(zip((12.82, 14.0, 15.18), markings)):
        ax.add_patch(
            FancyBboxPatch(
                (xx - 0.35, 1.94),
                0.7,
                0.53,
                boxstyle="round,pad=.025,rounding_size=.12",
                facecolor="white",
                edgecolor=MUTED,
                lw=1.1,
            )
        )
        if tokens == ("ES",):
            ax.add_patch(Circle((xx, 2.205), 0.085, facecolor=ORANGE,
                                edgecolor="none", zorder=6))
            ax.add_patch(Circle((xx, 2.205), 0.04, facecolor=BLUE,
                                edgecolor="none", zorder=7))
        else:
            for dx, color in zip((-0.1, 0.1), tokens):
                ax.add_patch(Circle((xx + dx, 2.205), 0.07, facecolor=color,
                                    edgecolor="none", zorder=6))
        if i < 2:
            arrow((xx + 0.40, 2.21), (xx + 0.78, 2.21), MUTED, 1.1)
    tick(15.85, 2.21)
    label(14.0, 1.31, "Reaction sequences", 10, INK, "bold")

    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        out / "visual_abstract.pdf",
        dpi=300,
        facecolor="white",
        bbox_inches="tight",
        pad_inches=0.025,
    )
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=Path(__file__).resolve().parent
    )
    render_visual_abstract(parser.parse_args().output_dir)


if __name__ == "__main__":
    main()
