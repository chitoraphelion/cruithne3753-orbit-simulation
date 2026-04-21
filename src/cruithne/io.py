"""Save/load ``SimResult`` and clone results to compressed .npz."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .simulation import SimResult


def save_sim(path: str | Path, res: SimResult) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        p,
        years=res.years,
        t_sec=res.t_sec,
        Rh=res.Rh,
        Vh=res.Vh,
        REh=res.REh,
        VEh=res.VEh,
    )


def load_sim(path: str | Path) -> SimResult:
    d = np.load(path)
    return SimResult(
        years=d["years"],
        t_sec=d["t_sec"],
        Rh=d["Rh"],
        Vh=d["Vh"],
        REh=d["REh"],
        VEh=d["VEh"],
        initial_states={},
    )
