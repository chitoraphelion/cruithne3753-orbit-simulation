"""JPL Horizons fetchers for barycentric (ICRF/equatorial) state vectors.

Results are cached to ``data/horizons_cache/`` so repeated runs are instant
and don't hammer the server.
"""

from __future__ import annotations

import hashlib
import os
import warnings
from pathlib import Path

import numpy as np
from astropy.time import Time
from astroquery.jplhorizons import Horizons

from .constants import AU, DAY

warnings.filterwarnings("ignore", message=".*id_type.*")
warnings.filterwarnings("ignore", message=".*dubious year.*")

CACHE_DIR = Path(os.environ.get("CRUITHNE_CACHE", "data/horizons_cache"))


def _col(tab, name: str) -> np.ndarray:
    c = tab[name]
    arr = getattr(c, "filled", lambda v=None: c)(np.nan)
    return np.asarray(arr, dtype=float)


def _cache_path(key: str) -> Path:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    h = hashlib.sha1(key.encode()).hexdigest()[:16]
    return CACHE_DIR / f"{h}.npz"


def fetch_vector(
    obj_id: str,
    epoch: str,
    id_type: str = "smallbody",
) -> tuple[np.ndarray, np.ndarray]:
    """Single state vector (r, v) in SI, barycentric, equatorial (ICRF)."""
    key = f"vec|{obj_id}|{id_type}|{epoch}"
    path = _cache_path(key)
    if path.exists():
        d = np.load(path)
        return d["r"], d["v"]

    jd = Time(epoch, format="iso").jd
    q = Horizons(id=obj_id, id_type=id_type, location="@ssb", epochs=[jd])
    tab = q.vectors(refplane="earth")
    if len(tab) == 0:
        tab = Horizons(id=obj_id, id_type=id_type, location="@ssb", epochs=[jd + 1]).vectors(
            refplane="earth"
        )
        if len(tab) == 0:
            raise RuntimeError(f"Horizons returned no data for id={obj_id} at {epoch}")

    r = np.array([_col(tab, "x")[0], _col(tab, "y")[0], _col(tab, "z")[0]]) * AU
    v = np.array([_col(tab, "vx")[0], _col(tab, "vy")[0], _col(tab, "vz")[0]]) * AU / DAY
    np.savez(path, r=r, v=v)
    return r, v


def fetch_series(
    obj_id: str,
    start: str,
    end: str,
    step_days: int = 7,
    id_type: str = "smallbody",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """State-vector time series (jd, r, v) in SI."""
    key = f"ser|{obj_id}|{id_type}|{start}|{end}|{step_days}"
    path = _cache_path(key)
    if path.exists():
        d = np.load(path)
        return d["jd"], d["r"], d["v"]

    tab = Horizons(
        id=obj_id,
        id_type=id_type,
        location="@ssb",
        epochs={"start": start, "stop": end, "step": f"{step_days}d"},
    ).vectors(refplane="earth")
    if len(tab) == 0:
        raise RuntimeError("Horizons series empty")
    jd = np.asarray(tab["datetime_jd"], dtype=float)
    r = np.column_stack([_col(tab, "x"), _col(tab, "y"), _col(tab, "z")]) * AU
    v = np.column_stack([_col(tab, "vx"), _col(tab, "vy"), _col(tab, "vz")]) * AU / DAY
    np.savez(path, jd=jd, r=r, v=v)
    return jd, r, v
