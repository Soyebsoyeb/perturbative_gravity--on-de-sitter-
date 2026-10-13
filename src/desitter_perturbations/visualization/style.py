"""Common plotting style.

Colours follow a fixed categorical order (never cycled); at most four series appear in
any one panel.  Text uses text tokens, never series colours; grids and axes are
recessive; series are identified by legends or direct labels, not by colour alone
(line styles differ where several series overlap).
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e4e3df"
# categorical slots 1..4 (blue, orange, aqua, yellow) + magenta, violet
C = {"blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a", "yellow": "#eda100",
     "magenta": "#e87ba4", "violet": "#6250d6"}
SERIES = [C["blue"], C["orange"], C["aqua"], C["yellow"]]
LINESTYLES = ["-", "--", "-.", ":"]


def apply_style():
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": INK2, "ytick.color": INK2,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.7,
        "axes.axisbelow": True,
        "lines.linewidth": 1.8, "lines.markersize": 5.5,
        "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "regular",
        "axes.labelsize": 10, "legend.frameon": False, "legend.fontsize": 9,
        "figure.dpi": 110, "savefig.dpi": 200, "mathtext.fontset": "dejavusans",
    })


apply_style()


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def figure_dir(sub: str) -> Path:
    base = Path(os.environ.get("DESITTER_FIGDIR", repo_root() / "figures"))
    d = base / sub
    d.mkdir(parents=True, exist_ok=True)
    return d


def savefig(fig, sub: str, name: str) -> str:
    path = figure_dir(sub) / f"{name}.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return str(path)


def direct_label(ax, x, y, text, color=INK2, dx=4, dy=0, ha="left", **kw):
    ax.annotate(text, (x, y), xytext=(dx, dy), textcoords="offset points",
                color=color, ha=ha, va="center", fontsize=8.5, **kw)
