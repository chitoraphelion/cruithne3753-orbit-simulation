"""Recreation of Wiegert et al. (1998) Fig. 3: 3753 Cruithne in the Sun-Earth
corotating frame over roughly half a horseshoe cycle (~385 yr).

Usage: uv run python scripts/plot_fig3_horseshoe.py [--sim data/sim_main.npz]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle

from cruithne.constants import DAY
from cruithne.frames import corotating_xy
from cruithne.io import load_sim
from cruithne.plot_style import COLORS, apply as apply_style


def plot_horseshoe(
    res,
    year_from: float,
    year_to: float,
    out: Path,
    points_every_yr: float = 0.01,
) -> None:
    apply_style()
    years = res.years
    mask = (years >= year_from) & (years <= year_to)
    if not mask.any():
        raise ValueError(f"No samples in window [{year_from}, {year_to}]; have {years.min()}..{years.max()}")

    xa, ya, _xe, _ye = corotating_xy(res.Rh, res.REh)
    dt = float(np.mean(np.diff(res.t_sec)))
    step_pts = max(1, int(round((points_every_yr * 365.25 * DAY) / dt)))

    idx_all = np.where(mask)[0]
    idx_pts = idx_all[::step_pts]
    xs = xa[idx_pts]
    ys = ya[idx_pts]

    # Outer edge in angular bins (following the sense of Wiegert Fig. 3).
    theta = np.arctan2(ya[mask], xa[mask])
    rho = np.hypot(xa[mask], ya[mask])
    nbins = 720
    edges = np.linspace(-np.pi, np.pi, nbins + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    rmax = np.full(nbins, np.nan)
    bid = np.digitize(theta, edges) - 1
    bid[bid == nbins] = nbins - 1
    for b in range(nbins):
        m = bid == b
        if m.any():
            rmax[b] = np.nanmax(rho[m])
    valid = ~np.isnan(rmax)
    x_edge = rmax[valid] * np.cos(centres[valid])
    y_edge = rmax[valid] * np.sin(centres[valid])
    # Close the outline smoothly
    order = np.argsort(np.arctan2(y_edge, x_edge))
    x_edge, y_edge = x_edge[order], y_edge[order]
    x_edge = np.append(x_edge, x_edge[0])
    y_edge = np.append(y_edge, y_edge[0])

    fig, ax = plt.subplots(figsize=(6.6, 6.6))

    ax.scatter(xs, ys, s=5, alpha=0.35, linewidths=0, color=COLORS["cruithne"],
               label=f"Cruithne (каждые {points_every_yr:.2f} года)")
    ax.plot(x_edge, y_edge, lw=2.2, color=COLORS["neutral"], alpha=0.9,
            label="внешняя кромка")

    # Sun at origin (⊙ symbol)
    ax.scatter([0], [0], s=70, marker="o", facecolor=COLORS["sun"],
               edgecolor="k", lw=0.8, zorder=5)
    ax.scatter([0], [0], s=8, marker="o", color="k", zorder=6)
    ax.annotate("Солнце", (0, 0), xytext=(0.04, -0.02), fontsize=9,
                color=COLORS["neutral"])

    # Earth at (1,0) — circled plus
    earth_xy = (1.0, 0.0)
    ax.add_patch(Circle(earth_xy, 0.035, fill=True, facecolor=COLORS["earth"],
                        edgecolor="k", lw=0.8, zorder=5))
    ax.plot([earth_xy[0] - 0.03, earth_xy[0] + 0.03], [0, 0], lw=0.8, color="k", zorder=6)
    ax.plot([earth_xy[0], earth_xy[0]], [-0.03, 0.03], lw=0.8, color="k", zorder=6)
    ax.annotate("Земля (EMB)", earth_xy, xytext=(1.05, 0.06), fontsize=9,
                color=COLORS["neutral"])

    # Earth's orbit reference
    th = np.linspace(0, 2 * np.pi, 721)
    ax.plot(np.cos(th), np.sin(th), lw=0.6, ls="--", color=COLORS["muted"], alpha=0.6)

    ax.set_aspect("equal", "box")
    ax.set_xlim(-1.55, 1.55)
    ax.set_ylim(-1.55, 1.55)
    ax.set_xlabel(r"$x'$ (а.е.)")
    ax.set_ylabel(r"$y'$ (а.е.)")
    ax.set_title(f"3753 Cruithne в коротирующем кадре ({year_from:.0f}–{year_to:.0f})")
    ax.legend(loc="lower left")

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    print(f"[fig3] wrote {out}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--sim", default="data/sim_main.npz")
    p.add_argument("--from-year", type=float, default=1900.0)
    p.add_argument("--to-year", type=float, default=2285.0)
    p.add_argument("--out", default="plots/generated/fig3_horseshoe.png")
    a = p.parse_args()
    res = load_sim(a.sim)
    plot_horseshoe(res, a.from_year, a.to_year, Path(a.out))


if __name__ == "__main__":
    main()
