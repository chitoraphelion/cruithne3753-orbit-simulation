"""Osculating Kepler elements from position/velocity vectors.

All angles are returned in degrees, wrapped to [0, 360).
"""

from __future__ import annotations

import numpy as np


def rv_to_kepler(
    R: np.ndarray, V: np.ndarray, mu: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Vectorized r,v → (a, e, i, Ω, ω, M). Inputs shape (N, 3). Units SI."""
    R = np.asarray(R, float)
    V = np.asarray(V, float)

    rn = np.linalg.norm(R, axis=1)
    vn = np.linalg.norm(V, axis=1)

    h = np.cross(R, V)
    hn = np.linalg.norm(h, axis=1)

    k = np.array([0.0, 0.0, 1.0])
    n = np.cross(np.broadcast_to(k, h.shape), h)
    nn = np.linalg.norm(n, axis=1)

    e_vec = (np.cross(V, h) / mu) - (R / rn[:, None])
    e = np.linalg.norm(e_vec, axis=1)

    eps = 0.5 * vn * vn - mu / rn
    a = np.where(eps < 0, -mu / (2 * eps), np.inf)

    ci = np.clip(h[:, 2] / np.maximum(hn, 1e-30), -1, 1)
    i = np.degrees(np.arccos(ci)) % 360.0

    Omega = np.zeros_like(nn)
    mask_n = nn >= 1e-16
    Omega[mask_n] = np.degrees(np.arctan2(n[mask_n, 1], n[mask_n, 0])) % 360.0

    omega = np.zeros_like(e)
    mask_e = e >= 1e-12
    mask_ew = mask_e & mask_n
    if np.any(mask_ew):
        cosw = np.clip(
            np.sum(n[mask_ew] * e_vec[mask_ew], axis=1) / (nn[mask_ew] * e[mask_ew]), -1, 1
        )
        sinw = np.sum(np.cross(n[mask_ew], e_vec[mask_ew]) * h[mask_ew], axis=1) / (
            nn[mask_ew] * e[mask_ew] * hn[mask_ew]
        )
        omega[mask_ew] = np.degrees(np.arctan2(sinw, cosw)) % 360.0

    nu = np.zeros_like(e)
    if np.any(mask_e):
        cosnu = np.clip(
            np.sum(e_vec[mask_e] * R[mask_e], axis=1) / (e[mask_e] * rn[mask_e]), -1, 1
        )
        sinnu = np.sum(np.cross(e_vec[mask_e], R[mask_e]) * h[mask_e], axis=1) / (
            e[mask_e] * rn[mask_e] * hn[mask_e]
        )
        nu[mask_e] = np.degrees(np.arctan2(sinnu, cosnu)) % 360.0
    mask_circ = ~mask_e
    if np.any(mask_circ):
        nu[mask_circ] = np.degrees(np.arctan2(R[mask_circ, 1], R[mask_circ, 0])) % 360.0

    M = np.zeros_like(e)
    mask_ell = e < (1 - 1e-12)
    if np.any(mask_ell):
        nu_r = np.radians(nu[mask_ell])
        em = e[mask_ell]
        denom = 1 + em * np.cos(nu_r)
        cosE = (em + np.cos(nu_r)) / np.maximum(denom, 1e-30)
        sinE = (np.sqrt(np.maximum(1 - em * em, 0.0)) * np.sin(nu_r)) / np.maximum(denom, 1e-30)
        E = np.arctan2(sinE, cosE)
        M[mask_ell] = np.degrees(E - em * np.sin(E)) % 360.0
    mask_hyp = e > (1 + 1e-12)
    if np.any(mask_hyp):
        nu_r = np.radians(nu[mask_hyp])
        eh = e[mask_hyp]
        coshH = np.maximum((eh + np.cos(nu_r)) / (1 + eh * np.cos(nu_r)), 1.0)
        H = np.arccosh(coshH)
        H = np.where(np.sin(nu_r) < 0, -H, H)
        M[mask_hyp] = np.degrees(eh * np.sinh(H) - H) % 360.0
    mask_par = ~(mask_ell | mask_hyp)
    if np.any(mask_par):
        D = np.tan(np.radians(nu[mask_par]) / 2)
        M[mask_par] = np.degrees(D + (D**3) / 3) % 360.0

    return a, e, i, Omega, omega, M


def mean_longitude(Omega_deg: np.ndarray, omega_deg: np.ndarray, M_deg: np.ndarray) -> np.ndarray:
    return (Omega_deg + omega_deg + M_deg) % 360.0


def wrap_pm180(x_deg: np.ndarray) -> np.ndarray:
    """Wrap angle difference into (-180, 180]."""
    return ((np.asarray(x_deg) + 540.0) % 360.0) - 180.0


def nodal_distances(a: np.ndarray, e: np.ndarray, omega_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """r_+ (ascending node) and r_− (descending node) heliocentric distances."""
    w = np.radians(omega_deg)
    r_plus = a * (1.0 - e**2) / (1.0 + e * np.cos(w))
    r_minus = a * (1.0 - e**2) / (1.0 - e * np.cos(w))
    return r_plus, r_minus
