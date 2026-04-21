"""Shared matplotlib style for figures in the report."""

from __future__ import annotations

import matplotlib as mpl

STYLE: dict = {
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "figure.autolayout": False,
    "font.family": "serif",
    "font.size": 10,
    "axes.labelsize": 10,
    "axes.titlesize": 11,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "legend.fontsize": 9,
    "lines.linewidth": 1.1,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
}


def apply() -> None:
    mpl.rcParams.update(STYLE)


# Colour palette shared across figures.
COLORS = {
    "cruithne": "#C0392B",
    "earth": "#2E86AB",
    "venus": "#E9B44C",
    "mars": "#9B2226",
    "sun": "#F1C40F",
    "neutral": "#2C3E50",
    "accent": "#1B9AAA",
    "muted": "#95A5A6",
}
