"""Shared print typography and vector exports for manuscript results."""

from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle
from matplotlib.transforms import offset_copy

INK = "#193342"
BLUE = "#0072B2"
GREEN = "#008564"
ORANGE = "#B57A00"
RED = "#C34724"
GRAY = "#667680"
HAIR = "#DCE5E9"
PALE = "#F0F5F7"


def configure():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "text.color": INK,
            "axes.labelcolor": INK,
            "xtick.color": GRAY,
            "ytick.color": GRAY,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": HAIR,
            "savefig.bbox": None,
        }
    )


def heading(ax, letter, title, y=1.06):
    """Panel header: an ink badge carrying the letter, then the title.

    The badge is the shared panel marker of the manuscript (Figure/synstyle.tex):
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


def clean(ax, grid="y"):
    ax.grid(axis=grid, color=HAIR, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)


def arrow(ax, a, b, color=GRAY, lw=1.1, rad=0, shrink=0):
    ax.add_patch(
        FancyArrowPatch(
            a,
            b,
            arrowstyle="-|>",
            mutation_scale=10,
            linewidth=lw,
            color=color,
            connectionstyle=f"arc3,rad={rad}",
            shrinkA=shrink,
            shrinkB=shrink,
        )
    )


def verdict(ax, x, y, ok, size=11):
    # Shape and text encode status as well as colour.
    ax.text(
        x,
        y,
        "✓" if ok else "×",
        ha="center",
        va="center",
        color=GREEN if ok else RED,
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
