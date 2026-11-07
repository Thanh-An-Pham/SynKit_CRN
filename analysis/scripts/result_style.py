"""Shared print typography and vector exports for manuscript results."""

from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyArrowPatch, Circle
from matplotlib.path import Path as MarkerPath
from matplotlib.transforms import offset_copy

# Result figures are black and white, after the plain style of the keriso
# paper: one ink, two grays, series told apart by marker shape and dashing.
# Only the visual abstract (visual_abstract.py) keeps its own colours.
INK = "black"
BLUE = "black"      # method / measured series
GREEN = "black"     # verified result (shape, not colour, carries it)
ORANGE = "#4D4D4D"  # phosphorus and highlighted participants
RED = "black"
GRAY = "#7F7F7F"    # baseline, structure, failed cases
HAIR = "#BFBFBF"
PALE = "#E6E6E6"
SERIF = "CMU Serif"


def configure():
    """Serif text as in the manuscript, boxed axes, inward ticks on all sides."""
    plt.rcParams.update(
        {
            "font.family": SERIF,
            "mathtext.fontset": "cm",
            "font.size": 9,
            "axes.labelsize": 9,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "legend.fontsize": 8.5,
            # CMU Serif is a CFF font: Type 3 embeds it cleanly, Type 42 does not.
            "pdf.fonttype": 3,
            "svg.fonttype": "none",
            "text.color": "black",
            "axes.labelcolor": "black",
            "axes.edgecolor": "black",
            "axes.linewidth": 0.5,
            "axes.spines.top": True,
            "axes.spines.right": True,
            "axes.grid": False,
            "axes.unicode_minus": True,
            "xtick.color": "black",
            "ytick.color": "black",
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 3,
            "ytick.major.size": 3,
            "xtick.minor.size": 1.5,
            "ytick.minor.size": 1.5,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "xtick.minor.width": 0.4,
            "ytick.minor.width": 0.4,
            "lines.linewidth": 0.8,
            "lines.markersize": 4,
            "legend.frameon": True,
            "legend.edgecolor": "black",
            "legend.fancybox": False,
            "legend.framealpha": 1,
            "patch.linewidth": 0.5,
            "hatch.linewidth": 0.5,
            "savefig.bbox": None,
        }
    )


def heading(ax, letter, title, y=1.06):
    """Panel header: an ink badge carrying the letter, then the title.

    The badge is the shared panel marker of the manuscript (paper_layout.tex):
    panel identity is encoded by the letter, not by colour alone, so the header
    survives greyscale printing.  Both pieces are offset from the axes origin in
    points, so the badge clears the figure edge whatever the axes width is.
    """
    figure = ax.get_figure()
    badge_at = offset_copy(ax.transAxes, fig=figure, x=7, y=0, units="points")
    title_at = offset_copy(ax.transAxes, fig=figure, x=20, y=0, units="points")
    ax.text(
        0.0,
        y,
        letter.upper(),
        transform=badge_at,
        ha="center",
        va="bottom",
        color="white",
        fontsize=8.3,
        fontweight="bold",
        bbox={"boxstyle": "circle,pad=0.19", "facecolor": INK, "edgecolor": "none"},
        clip_on=False,
    )
    ax.text(
        0.0,
        y,
        title,
        transform=title_at,
        ha="left",
        va="bottom",
        fontsize=10,
        weight="bold",
        clip_on=False,
    )


def sans_family():
    """Arial where installed, otherwise its metric-compatible Liberation Sans."""
    for name in ("Arial", "Liberation Sans"):
        try:
            font_manager.findfont(name, fallback_to_default=False)
            return name
        except ValueError:
            pass
    raise ValueError("Install Arial or Liberation Sans")


def use_sans():
    """Arial labels for the visual abstract, the one colour figure."""
    name = sans_family()
    plt.rcParams.update({
        "font.family": name,
        "mathtext.fontset": "custom",
        "mathtext.rm": name,
        "mathtext.it": f"{name}:italic",
        "mathtext.bf": f"{name}:bold",
        "mathtext.sf": name,
    })


def plain_heading(ax, letter, *args, y=1.06, x=0, **kwargs):
    """Plain parenthesized panel letter, "(a)". The caption names the panel."""
    if args and isinstance(args[0], (int, float)):
        y = args[0]
    if "y" in kwargs:
        y = kwargs["y"]
    if "x" in kwargs:
        x = kwargs["x"]
    ax.text(
        x, y, f"({letter.lower()})", transform=ax.transAxes,
        ha="left", va="bottom", fontsize=9.5,
        color="black", clip_on=False,
    )


def inset_label(ax, letter, x=0.03, y=0.96):
    """Panel letter inside the upper-left corner of a boxed plot."""
    ax.text(
        x, y, f"({letter.lower()})", transform=ax.transAxes,
        ha="left", va="top", fontsize=9.5, color="black",
    )


def clean(ax, grid="y"):
    """Boxed axes without grid lines."""
    ax.grid(False)
    ax.spines[:].set_visible(True)


def arrow(ax, a, b, color=GRAY, lw=1.1, rad=0, shrink=0):
    # Schematic arrows print as thin black lines; other colours are kept.
    if color == GRAY:
        color = "black"
    lw = min(lw, 0.75)
    ax.add_patch(
        FancyArrowPatch(
            a,
            b,
            arrowstyle="-|>",
            mutation_scale=8,
            linewidth=lw,
            color=color,
            connectionstyle=f"arc3,rad={rad}",
            shrinkA=shrink,
            shrinkB=shrink,
        )
    )


def verdict(ax, x, y, ok, size=11, vector=False):
    # Shape encodes status. The serif face lacks a check glyph, so the marks
    # are always drawn as paths.
    vector = True
    if vector:
        # Arial lacks a checkmark glyph, so draw the symbol directly.
        if ok:
            marker = MarkerPath([(-0.45, 0), (-0.12, -0.3), (0.5, 0.4)])
            ax.scatter(x, y, marker=marker, s=size**2 * 0.65,
                       facecolor="none", edgecolor=GREEN, linewidths=1.3)
        else:
            ax.scatter(x, y, marker="x", s=size**2 * 0.45,
                       color=GRAY, linewidths=1.3)
        return
    ax.text(
        x,
        y,
        "✓" if ok else "×",
        ha="center",
        va="center",
        color=GREEN if ok else GRAY,
        fontsize=size,
        weight="bold",
    )


def save(fig, out, name):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        out / f"{name}.pdf",
        dpi=300,
        facecolor="white",
        bbox_inches="tight",
        pad_inches=0.025,
    )
    plt.close(fig)
