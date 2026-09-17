#!/usr/bin/env python3
"""Illustrate the SynKit CRN-submodule workflow with editable vector primitives."""

from pathlib import Path
import argparse
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
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


def render_visual_abstract(out):
    font_manager.findfont("Arial", fallback_to_default=False)
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "mathtext.sf": "Arial",
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
    ax.set(xlim=(0, 16.5), ylim=(central_offset_y, 9.25), aspect="equal")
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
                facecolor=ORANGE,
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

    # Three illustrated acquisition routes, with only the source names in text.
    label(1.7, 8.65, "Build", 12, BLUE, "bold")
    label(8.0, 8.75, "SynKit", 21, INK, "bold")
    label(8.0, 8.11, "CRN submodule", 9, MUTED)
    label(14.15, 8.65, "Analyse", 12, GREEN, "bold")

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
            (2.4, 1.83), 0.1, facecolor="#FFF0CD", edgecolor=ORANGE, lw=1.2, zorder=6
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
    ax.add_patch(
        FancyBboxPatch(
            (4.85, 3.95),
            6.15,
            3.78,
            boxstyle="round,pad=.08,rounding_size=.24",
            facecolor=PALE,
            edgecolor="none",
        )
    )

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
            ORANGE if b == positions["E"] or a == positions["E"] else BLUE,
            1.3,
            shrink=9,
            # Keep the incoming arrowhead outside the complex's circular edge.
            shrink_end=14 if b == positions["ES"] else None,
        )
    arrow(unbind, positions["E"], ORANGE, 1.1, rad=-0.75, shrink=9)
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
            positions["P"], 5, radius=0.23, facecolor=GREEN, edgecolor="white", zorder=6
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
                        BLUE if net[i, j] < 0 else GREEN if net[i, j] > 0 else "#DFE7EB"
                    ),
                    edgecolor="none",
                )
            )
    label(5.79, 1.31, "Stoichiometry", 9)
    # The shaded branch indicates a hyperedge joining two required substrates.
    ax.add_patch(
        Polygon(
            [
                (7.4, 2.74),
                (7.55, 2.95),
                (8.58, 2.39),
                (7.55, 1.89),
                (7.4, 2.07),
                (8.14, 2.4),
            ],
            closed=True,
            facecolor="#DCEFE9",
            edgecolor="none",
        )
    )
    for a in [(7.48, 2.79), (7.48, 2.0)]:
        line(a, (8.36, 2.4), GREEN, 1.4)
    for a in [(7.48, 2.79), (7.48, 2.0), (8.36, 2.4)]:
        species(*a, color=GREEN, r=0.135)
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
        arrow((11.62, yy), (12.27, yy), GREEN, 1.3)

    # A conserved enzyme pool redistributes between free and bound forms.
    enzyme(13.02, 7.02, r=0.31)
    label(13.72, 7.02, "+", 14, GREEN)
    enzyme(14.51, 7.02, bound=True, r=0.31)
    label(14.0, 6.34, r"$E+ES=E_{\mathrm{total}}$", 10, GREEN)
    label(14.0, 5.96, "Conserved pools", 10, INK, "bold")

    # CRNT: two complex-graph linkage classes and an exact deficiency expression.
    for xx in (13.28, 14.72):
        ax.add_patch(
            Ellipse((xx, 4.68), 1.13, 0.98, facecolor="#E6F2ED", edgecolor="none")
        )
        for aa, bb, rad in [
            ((xx - 0.29, 4.68), (xx + 0.29, 4.68), 0.55),
            ((xx + 0.29, 4.68), (xx - 0.29, 4.68), 0.55),
        ]:
            arrow(aa, bb, GREEN, 1, rad=rad, shrink=4)
        for dx in (-0.29, 0.29):
            species(xx + dx, 4.68, color=GREEN, r=0.12)
    label(14.0, 3.99, r"$\delta=n-\ell-s$", 11, GREEN)
    label(14.0, 3.58, "CRNT tests", 10, INK, "bold")

    # Executability: a token moves through enabled transitions, with replay checks.
    for i, xx in enumerate((12.88, 14.0, 15.12)):
        ax.add_patch(
            FancyBboxPatch(
                (xx - 0.35, 1.94),
                0.7,
                0.53,
                boxstyle="round,pad=.025,rounding_size=.12",
                facecolor="white",
                edgecolor=GREEN,
                lw=1.1,
            )
        )
        for j in range(2 if i == 0 else 1):
            ax.add_patch(
                Circle(
                    (xx + (j - 0.5) * 0.19 if i == 0 else xx, 2.21),
                    0.07,
                    facecolor=BLUE,
                    edgecolor="none",
                )
            )
        tick(xx, 2.78)
        if i < 2:
            arrow((xx + 0.42, 2.21), (xx + 0.70, 2.21), GREEN, 1.1)
    label(14.0, 1.31, "Reaction sequences", 10, INK, "bold")

    # A compact key defines the central graph symbols.
    key_y = 0.46 + central_offset_y
    species(6.17, key_y, r=0.1)
    label(6.42, key_y, "Species", 8, ha="left")
    reaction(8.60, key_y, 0.12, 0.2)
    label(8.85, key_y, "Reaction", 8, ha="left")
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
