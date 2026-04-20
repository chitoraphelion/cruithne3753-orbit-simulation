"""Reference-frame transforms: corotating Sun-Earth frame, yearly averaging."""

from __future__ import annotations

import numpy as np

from .constants import AU, DAY


def corotating_xy(Rh: np.ndarray, REh: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Rotate (Rh, REh) so Earth sits on +x. Returns Cruithne (xa, ya) and EMB (xe, ye) in AU."""
    phi = np.arctan2(REh[:, 1], REh[:, 0])
    c, s = np.cos(phi), np.sin(phi)
    xa = (c * Rh[:, 0] + s * Rh[:, 1]) / AU
    ya = (-s * Rh[:, 0] + c * Rh[:, 1]) / AU
    xe = (c * REh[:, 0] + s * REh[:, 1]) / AU
    ye = (-s * REh[:, 0] + c * REh[:, 1]) / AU
    return xa, ya, xe, ye


def rolling_mean(x: np.ndarray, window_steps: int) -> np.ndarray:
    """Centered rolling mean with a flat kernel (edge effects via same-mode)."""
    if window_steps <= 1:
        return x
    kernel = np.ones(window_steps) / window_steps
    return np.convolve(x, kernel, mode="same")


def yearly_mean_xy(
    xa: np.ndarray, ya: np.ndarray, dt_sec: float
) -> tuple[np.ndarray, np.ndarray, int]:
    """Apply a ~1 yr centered mean to (xa, ya); returns smoothed arrays and window size."""
    window = max(1, int(round((365.25 * DAY) / dt_sec)))
    return rolling_mean(xa, window), rolling_mean(ya, window), window
