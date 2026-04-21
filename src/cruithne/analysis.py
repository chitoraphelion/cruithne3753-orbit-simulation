"""Derived quantities (Kepler elements, Δλ, periods, close approaches)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.signal import find_peaks

from .constants import DAY, GM_SUN
from .kepler import mean_longitude, nodal_distances, rv_to_kepler, wrap_pm180
from .simulation import SimResult


@dataclass
class Derived:
    a_ast: np.ndarray
    e_ast: np.ndarray
    i_ast: np.ndarray
    Om_ast: np.ndarray
    om_ast: np.ndarray
    M_ast: np.ndarray
    a_emb: np.ndarray
    lam_ast: np.ndarray
    lam_emb: np.ndarray
    dlam: np.ndarray
    P_ast: np.ndarray
    P_emb: np.ndarray
    dP: np.ndarray
    r_plus: np.ndarray
    r_minus: np.ndarray
    dist_GE: np.ndarray
    peak_indices: np.ndarray


def derive(res: SimResult, min_approach_days: float = 300.0) -> Derived:
    a_ast, e_ast, i_ast, Om_ast, om_ast, M_ast = rv_to_kepler(res.Rh, res.Vh, GM_SUN)
    a_emb, e_emb, i_emb, Om_emb, om_emb, M_emb = rv_to_kepler(res.REh, res.VEh, GM_SUN)

    lam_ast = mean_longitude(Om_ast, om_ast, M_ast)
    lam_emb = mean_longitude(Om_emb, om_emb, M_emb)
    dlam = wrap_pm180(lam_ast - lam_emb)

    P_ast = 2 * np.pi * np.sqrt((a_ast**3) / GM_SUN)
    P_emb = 2 * np.pi * np.sqrt((a_emb**3) / GM_SUN)
    dP = P_ast - P_emb

    r_plus, r_minus = nodal_distances(a_ast, e_ast, om_ast)
    dist_GE = np.linalg.norm(res.Rh - res.REh, axis=1)

    dt = float(np.mean(np.diff(res.t_sec)))
    peaks, _ = find_peaks(-dist_GE, distance=int((min_approach_days * DAY) / dt))

    return Derived(
        a_ast=a_ast,
        e_ast=e_ast,
        i_ast=i_ast,
        Om_ast=Om_ast,
        om_ast=om_ast,
        M_ast=M_ast,
        a_emb=a_emb,
        lam_ast=lam_ast,
        lam_emb=lam_emb,
        dlam=dlam,
        P_ast=P_ast,
        P_emb=P_emb,
        dP=dP,
        r_plus=r_plus,
        r_minus=r_minus,
        dist_GE=dist_GE,
        peak_indices=peaks,
    )
