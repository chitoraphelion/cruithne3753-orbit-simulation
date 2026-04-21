"""Command-line entry points: ``cruithne-main`` and ``cruithne-clones``."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from .clones import CloneConfig, run_ensemble
from .io import save_sim
from .simulation import SimConfig, integrate


def run_main() -> None:
    p = argparse.ArgumentParser(description="Integrate 3753 Cruithne + planets.")
    p.add_argument("--fast", action="store_true", help="short smoke run (1900-2100, 7-day output)")
    p.add_argument("--start", default="1900-01-01")
    p.add_argument("--end", default="3900-01-01")
    p.add_argument("--out", default="data/sim_main.npz")
    args = p.parse_args()

    cfg = SimConfig(
        start_date=args.start,
        end_date="2100-01-01" if args.fast else args.end,
        output_step_days=7.0 if args.fast else 2.0,
        integrator_step_days=0.5,
    )
    print(f"[main] integrating {cfg.start_date} .. {cfg.end_date}")
    res = integrate(cfg)
    save_sim(args.out, res)
    print(f"[main] wrote {args.out}  (N={len(res.years)})")


def _write_csv(path: Path, result) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["clone", "planet", "past_SR", "past_RI", "future_SR", "future_RI"])
        for planet in result.R_I:
            for c in result.clone_names:
                w.writerow(
                    [
                        c,
                        planet,
                        result.fp_past_SR[planet][c],
                        result.fp_past_RI[planet][c],
                        result.fp_future_SR[planet][c],
                        result.fp_future_RI[planet][c],
                    ]
                )


def run_clones() -> None:
    p = argparse.ArgumentParser(description="Run Cruithne clone ensembles.")
    p.add_argument("--fast", action="store_true", help="short span, coarser sampling")
    p.add_argument("--start", default="1900-01-01")
    p.add_argument("--outdir", default="data/clones")
    p.add_argument("--eps", nargs="*", type=float, default=None)
    args = p.parse_args()

    cfg = CloneConfig(
        start_date=args.start,
        years_past=2000.0 if args.fast else 13000.0,
        years_future=2000.0 if args.fast else 13000.0,
        sample_years=0.5 if args.fast else 0.05,
        epsilon_list=tuple(args.eps) if args.eps else (1e-5, 1e-6),
    )
    outdir = Path(args.outdir)
    for eps in cfg.epsilon_list:
        print(f"[clones] eps={eps:.0e}")
        res = run_ensemble(cfg, eps)
        _write_csv(outdir / f"clones_eps{eps:.0e}_first_passages.csv", res)
        np.savez_compressed(
            outdir / f"clones_eps{eps:.0e}_meta.npz",
            R_I=np.array([res.R_I[k] for k in ("Venus", "Earth", "Mars")]),
            SR=np.array([res.SR[k] for k in ("Venus", "Earth", "Mars")]),
            planets=np.array(["Venus", "Earth", "Mars"]),
        )
    print(f"[clones] done -> {outdir}")
