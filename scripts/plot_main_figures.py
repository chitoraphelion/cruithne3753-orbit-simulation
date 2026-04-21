"""Render the main-run figures other than Fig. 3 (that lives in its own script).

Produces:
  fig1_nonrotating.png  — 3-panel Sun-centered view (xy, xz, yz)
  fig2_a_t.png          — semimajor axis evolution (Wiegert Fig. 2)
  fig5_dlam.png         — resonant phase Δλ(t) over 2000 yr
  node_distances.png    — r±(t)
  close_earth.png       — yearly minimum |r_C − r_EMB|
  a_P_dP.png            — a(t), P(t), ΔP(t) combined
  a_vs_dlam.png         — phase portrait
  yearly_mean_horseshoe.png — averaged horseshoe (Wiegert Fig. 4 analogue)
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle

from cruithne.analysis import derive
from cruithne.constants import AU, DAY
from cruithne.frames import corotating_xy, yearly_mean_xy
from cruithne.io import load_sim
from cruithne.plot_style import COLORS, apply as apply_style


OUT = Path("plots/generated")


def fig1_nonrotating(res) -> None:
    x, y, z = res.Rh[:, 0] / AU, res.Rh[:, 1] / AU, res.Rh[:, 2] / AU
    xe, ye, ze = res.REh[:, 0] / AU, res.REh[:, 1] / AU, res.REh[:, 2] / AU
    mask = res.years <= res.years.min() + 1.0  # first year for Earth/planet outlines
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 4.0))
    for ax, (a, b, labels) in zip(
        axes,
        ((x, y, ("x", "y")), (x, z, ("x", "z")), (y, z, ("y", "z"))),
        strict=False,
    ):
        ax.plot(a, b, lw=0.4, color=COLORS["cruithne"], alpha=0.85, label="3753")
        if labels == ("x", "y"):
            ax.plot(xe[mask], ye[mask], lw=0.8, color=COLORS["earth"], label="Земля")
        elif labels == ("x", "z"):
            ax.plot(xe[mask], ze[mask], lw=0.8, color=COLORS["earth"])
        else:
            ax.plot(ye[mask], ze[mask], lw=0.8, color=COLORS["earth"])
        ax.scatter([0], [0], s=35, c=COLORS["sun"], edgecolor="k", lw=0.5, zorder=5)
        ax.set_aspect("equal", "box")
        ax.set_xlabel(f"{labels[0]} (а.е.)")
        ax.set_ylabel(f"{labels[1]} (а.е.)")
        ax.set_xlim(-2.2, 2.2)
        ax.set_ylim(-2.2, 2.2)
    axes[0].legend(loc="upper right")
    fig.suptitle("Орбита 3753 Cruithne в нерасположенной гелиоцентрической СК")
    fig.tight_layout()
    fig.savefig(OUT / "fig1_nonrotating.png")
    plt.close(fig)


def fig2_a_t(res, d) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 3.3))
    ax.plot(res.years, d.a_ast / AU, lw=0.7, color=COLORS["cruithne"])
    ax.axhline(1.0, color=COLORS["muted"], lw=0.6, ls="--")
    ax.set_xlabel("Год")
    ax.set_ylabel(r"$a$ (а.е.)")
    ax.set_title("Большая полуось 3753 Cruithne")
    fig.tight_layout()
    fig.savefig(OUT / "fig2_a_t.png")
    plt.close(fig)


def fig5_dlam(res, d, years_span: float = 2000.0) -> None:
    x0 = res.years.min()
    x1 = min(x0 + years_span, res.years.max())
    span = x1 - x0
    dt = float(np.mean(np.diff(res.t_sec)))
    step = max(1, int(round((0.1 * 365.25 * DAY) / dt)))
    m = (res.years >= x0) & (res.years <= x1)
    idx = np.where(m)[0][::step]
    fig, ax = plt.subplots(figsize=(8.5, 3.2))
    ax.scatter(res.years[idx], d.dlam[idx], s=5, alpha=0.6, linewidths=0,
               color=COLORS["cruithne"])
    ax.axhline(0, color=COLORS["muted"], lw=0.6)
    ax.set_xlabel("Год")
    ax.set_ylabel(r"$\Delta\lambda$ (°)")
    ax.set_title(rf"Резонансная фаза $\Delta\lambda(t)$ ({span:.0f} лет, выборка каждые 0.1 года)")
    fig.tight_layout()
    fig.savefig(OUT / "fig5_dlam.png")
    plt.close(fig)


def node_distances(res, d) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 3.2))
    ax.plot(res.years, d.r_plus / AU, lw=0.8, color=COLORS["accent"],
            label=r"$r_+$ (восходящий узел)")
    ax.plot(res.years, d.r_minus / AU, lw=0.8, color=COLORS["mars"],
            label=r"$r_-$ (нисходящий узел)")
    ax.axhline(1.0, color=COLORS["muted"], lw=0.6, ls="--", label="орбита Земли")
    ax.set_xlabel("Год")
    ax.set_ylabel(r"$r_\pm$ (а.е.)")
    ax.set_title("Нодальные расстояния Cruithne")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT / "node_distances.png")
    plt.close(fig)


def close_earth(res, d) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 3.2))
    ax.plot(res.years, d.dist_GE / AU, lw=0.6, color=COLORS["muted"], label="|r_C − r_⊕|")
    if len(d.peak_indices):
        ax.scatter(res.years[d.peak_indices], d.dist_GE[d.peak_indices] / AU,
                   s=10, c=COLORS["cruithne"], label="локальные минимумы")
    ax.set_xlabel("Год")
    ax.set_ylabel("Расстояние (а.е.)")
    ax.set_title("Сближения 3753 Cruithne с системой Земля–Луна")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT / "close_earth.png")
    plt.close(fig)


def a_P_dP(res, d) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(8.5, 5.0), sharex=True)
    axes[0].plot(res.years, d.a_ast / AU, lw=0.7, color=COLORS["cruithne"])
    axes[0].axhline(1.0, color=COLORS["muted"], lw=0.6, ls="--")
    axes[0].set_ylabel(r"$a$ (а.е.)")
    axes[1].plot(res.years, d.P_ast / DAY / 365.25, lw=0.7, color=COLORS["cruithne"],
                 label=r"$P_{\rm Cru}$")
    axes[1].plot(res.years, d.P_emb / DAY / 365.25, lw=0.7, color=COLORS["earth"],
                 alpha=0.6, label=r"$P_\oplus$")
    axes[1].plot(res.years, d.dP / DAY / 365.25, lw=0.7, color=COLORS["accent"],
                 label=r"$\Delta P$")
    axes[1].set_xlabel("Год")
    axes[1].set_ylabel("Период (годы)")
    axes[1].legend(loc="best")
    fig.suptitle(r"Эволюция $a(t)$, $P(t)$ и $\Delta P(t)$")
    fig.tight_layout()
    fig.savefig(OUT / "a_P_dP.png")
    plt.close(fig)


def a_vs_dlam(res, d) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.scatter(((d.lam_ast - d.lam_emb) % 360.0), d.a_ast / AU,
               s=1.2, alpha=0.3, linewidths=0, color=COLORS["cruithne"])
    ax.set_xlabel(r"$\Delta\lambda$ (°)")
    ax.set_ylabel(r"$a$ (а.е.)")
    ax.set_title(r"Фазовый портрет $a$ vs $\Delta\lambda$")
    fig.tight_layout()
    fig.savefig(OUT / "a_vs_dlam.png")
    plt.close(fig)


def yearly_mean_horseshoe(res) -> None:
    xa, ya, _xe, _ye = corotating_xy(res.Rh, res.REh)
    dt = float(np.mean(np.diff(res.t_sec)))
    xa_s, ya_s, _ = yearly_mean_xy(xa, ya, dt)
    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    ax.plot(xa, ya, lw=0.25, color=COLORS["muted"], alpha=0.6,
            label="траектория")
    ax.plot(xa_s, ya_s, lw=1.1, color=COLORS["cruithne"],
            label="среднее за 1 год")
    ax.scatter([0], [0], s=70, marker="o", facecolor=COLORS["sun"],
               edgecolor="k", lw=0.8, zorder=5)
    ax.add_patch(Circle((1.0, 0.0), 0.035, fill=True, facecolor=COLORS["earth"],
                        edgecolor="k", lw=0.8, zorder=5))
    th = np.linspace(0, 2 * np.pi, 721)
    ax.plot(np.cos(th), np.sin(th), lw=0.6, ls="--", color=COLORS["muted"], alpha=0.6)
    ax.set_aspect("equal", "box")
    ax.set_xlim(-1.55, 1.55)
    ax.set_ylim(-1.55, 1.55)
    ax.set_xlabel(r"$x'$ (а.е.)")
    ax.set_ylabel(r"$y'$ (а.е.)")
    ax.set_title("Годовое усреднение horseshoe в коротирующем кадре")
    ax.legend(loc="lower left")
    fig.tight_layout()
    fig.savefig(OUT / "yearly_mean_horseshoe.png")
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--sim", default="data/sim_main.npz")
    args = p.parse_args()

    apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    res = load_sim(args.sim)
    d = derive(res)

    fig1_nonrotating(res)
    fig2_a_t(res, d)
    fig5_dlam(res, d)
    node_distances(res, d)
    close_earth(res, d)
    a_P_dP(res, d)
    a_vs_dlam(res, d)
    yearly_mean_horseshoe(res)
    print(f"[plots] wrote all figures to {OUT}")


if __name__ == "__main__":
    main()
