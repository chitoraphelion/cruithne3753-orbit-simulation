"""Validate integrated a(t) against JPL Horizons for the modern era."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time

from cruithne.analysis import derive
from cruithne.constants import AU, GM_SUN
from cruithne.horizons import fetch_series
from cruithne.io import load_sim
from cruithne.kepler import rv_to_kepler
from cruithne.plot_style import COLORS, apply as apply_style


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--sim", default="data/sim_main.npz")
    p.add_argument("--start", default="1901-01-01")
    p.add_argument("--end", default="2299-01-01")
    p.add_argument("--step-days", type=int, default=7)
    p.add_argument("--out", default="plots/generated/a_validation.png")
    args = p.parse_args()

    apply_style()
    res = load_sim(args.sim)
    d = derive(res)

    jd_h, r_b, v_b = fetch_series("3753", args.start, args.end, args.step_days, "smallbody")
    jd_s, rs, vs = fetch_series("10", args.start, args.end, args.step_days, "majorbody")
    assert np.allclose(jd_h, jd_s)
    Rh_H = r_b - rs
    Vh_H = v_b - vs
    a_h, *_ = rv_to_kepler(Rh_H, Vh_H, GM_SUN)
    years_h = Time(jd_h, format="jd", scale="tdb").decimalyear

    fig, ax = plt.subplots(figsize=(9.0, 3.3))
    ax.plot(res.years, d.a_ast / AU, lw=0.9, color=COLORS["cruithne"], label="REBOUND/WHFast")
    ax.plot(years_h, a_h / AU, lw=0.9, color=COLORS["neutral"], alpha=0.9,
            label=f"JPL Horizons ({args.step_days} д.)")
    ax.set_xlabel("Год")
    ax.set_ylabel(r"$a$ (а.е.)")
    ax.set_title("Валидация: большая полуось Cruithne против JPL Horizons")
    ax.legend()
    fig.tight_layout()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out)
    plt.close(fig)
    print(f"[validation] wrote {args.out}")


if __name__ == "__main__":
    main()
