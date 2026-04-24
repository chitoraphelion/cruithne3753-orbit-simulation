"""Ensemble-of-clones first-passage statistics.

Each ensemble consists of 2^6 = 64 clones obtained by perturbing the nominal
Cruithne state by ±ε in each of (x, y, z, vx, vy, vz), normalized by the
magnitudes of r0 and v0. We scan forwards and backwards in time, recording
the first passage of each clone into the Laplace sphere of influence R_I
and the "large sphere" S·R_I of each terrestrial planet.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import product

import numpy as np
import rebound
from astropy.time import Time
from tqdm import tqdm

from .constants import (
    CRUITHNE_ID,
    DAY,
    G,
    GM_SUN,
    M_EARTH,
    M_MARS,
    M_MOON,
    M_SUN,
    M_VENUS,
    PLANETS,
    YEAR,
)
from .horizons import fetch_vector

ORIGIN_YEAR = 1997.0
TERRESTRIAL = (
    ("Venus", "Venus", M_VENUS),
    ("Earth", "EMB", M_EARTH + M_MOON),
    ("Mars", "Mars", M_MARS),
)


def laplace_sphere_of_influence(a_m: float, m_p: float) -> float:
    return a_m * (m_p / M_SUN) ** (2.0 / 5.0)


def instantaneous_a(r: np.ndarray, v: np.ndarray, r_sun: np.ndarray, v_sun: np.ndarray) -> float:
    R = np.linalg.norm(r - r_sun)
    V = np.linalg.norm(v - v_sun)
    eps = 0.5 * V * V - GM_SUN / R
    return -GM_SUN / (2.0 * eps)


@dataclass
class CloneConfig:
    start_date: str = "1900-01-01"
    integrator_step_days: float = 0.5
    integrator_name: str = "whfast"
    years_past: float = 13000.0
    years_future: float = 13000.0
    sample_years: float = 0.05
    epsilon_list: tuple[float, ...] = (1e-5, 1e-6)
    s_factor: float = 5.0


@dataclass
class CloneResult:
    eps: float
    clone_names: list[str]
    R_I: dict[str, float]
    SR: dict[str, float]
    fp_future_SR: dict[str, dict[str, float]] = field(default_factory=dict)
    fp_future_RI: dict[str, dict[str, float]] = field(default_factory=dict)
    fp_past_SR: dict[str, dict[str, float]] = field(default_factory=dict)
    fp_past_RI: dict[str, dict[str, float]] = field(default_factory=dict)


def build_base_sim(cfg: CloneConfig) -> tuple[rebound.Simulation, dict[str, tuple[np.ndarray, np.ndarray]]]:
    sim = rebound.Simulation()
    sim.G = G
    sim.integrator = cfg.integrator_name
    sim.dt = cfg.integrator_step_days * DAY
    sim.units = ("m", "s", "kg")
    if cfg.integrator_name == "mercurius":
        # Switch to IAS15 inside r_crit = 3 Hill radii of any planet
        sim.ri_mercurius.r_crit_hill = 3.0
    states: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for name, obj_id, id_type, mass in PLANETS:
        r, v = fetch_vector(obj_id, cfg.start_date, id_type=id_type)
        states[name] = (r, v)
        sim.add(m=mass, x=r[0], y=r[1], z=r[2], vx=v[0], vy=v[1], vz=v[2], hash=name)
    return sim, states


def add_clones(sim: rebound.Simulation, eps: float, r0: np.ndarray, v0: np.ndarray) -> list[str]:
    rn = float(np.linalg.norm(r0))
    vn = float(np.linalg.norm(v0))
    names: list[str] = []
    for sig in product((+1, -1), repeat=6):
        dr = eps * rn * np.array(sig[:3], float) / np.sqrt(3.0)
        dv = eps * vn * np.array(sig[3:], float) / np.sqrt(3.0)
        name = "c" + "".join("+" if s > 0 else "-" for s in sig)
        sim.add(
            m=0.0,
            x=r0[0] + dr[0],
            y=r0[1] + dr[1],
            z=r0[2] + dr[2],
            vx=v0[0] + dv[0],
            vy=v0[1] + dv[1],
            vz=v0[2] + dv[2],
            hash=name,
        )
        names.append(name)
    return names


def spheres_of_influence(
    states: dict[str, tuple[np.ndarray, np.ndarray]], s_factor: float
) -> tuple[dict[str, float], dict[str, float]]:
    r_sun, v_sun = states["Sun"]
    R_I: dict[str, float] = {}
    for key, state_name, mass in TERRESTRIAL:
        r, v = states[state_name]
        a = instantaneous_a(r, v, r_sun, v_sun)
        R_I[key] = laplace_sphere_of_influence(a, mass)
    SR = {k: s_factor * v for k, v in R_I.items()}
    return R_I, SR


def _scan(
    sim: rebound.Simulation,
    clone_names: list[str],
    years_span: float,
    sample_years: float,
    direction: str,
    t_origin_offset: float,
    R_I: dict[str, float],
    SR: dict[str, float],
) -> tuple[dict[str, dict[str, float]], dict[str, dict[str, float]]]:
    """Walk the simulation in steps of ``sample_years`` and record first passages.

    ``t_origin_offset`` is the decimal-year value of ``sim.t = 0`` relative to
    ``ORIGIN_YEAR`` (so that the recorded time for a clone is years-from-1997).
    """
    assert direction in ("future", "past")
    sign = +1.0 if direction == "future" else -1.0
    t_start = sim.t
    nsamples = int(np.ceil(years_span / sample_years))

    fp_SR = {p: {c: np.nan for c in clone_names} for p in R_I}
    fp_RI = {p: {c: np.nan for c in clone_names} for p in R_I}
    planet_hashes = {p: sn for p, sn, _ in TERRESTRIAL}

    for i in tqdm(range(1, nsamples + 1), desc=f"scan {direction}", leave=False):
        t_sim = t_start + sign * i * sample_years * YEAR
        sim.integrate(t_sim)
        t_rel_1997 = t_origin_offset + sign * i * sample_years
        planets = {p: sim.particles[planet_hashes[p]] for p in R_I}
        for cname in clone_names:
            pc = sim.particles[cname]
            for pname, pp in planets.items():
                dx = pc.x - pp.x
                dy = pc.y - pp.y
                dz = pc.z - pp.z
                d2 = dx * dx + dy * dy + dz * dz
                if np.isnan(fp_SR[pname][cname]) and d2 < SR[pname] ** 2:
                    fp_SR[pname][cname] = t_rel_1997
                if np.isnan(fp_RI[pname][cname]) and d2 < R_I[pname] ** 2:
                    fp_RI[pname][cname] = t_rel_1997
    return fp_SR, fp_RI


def run_ensemble(cfg: CloneConfig, eps: float) -> CloneResult:
    sim_base, states = build_base_sim(cfg)
    r_ast, v_ast = fetch_vector(CRUITHNE_ID[0], cfg.start_date, id_type=CRUITHNE_ID[1])
    names = add_clones(sim_base, eps, r_ast, v_ast)
    sim_base.move_to_com()
    R_I, SR = spheres_of_influence(states, cfg.s_factor)

    year0 = Time(cfg.start_date, format="iso", scale="tdb").decimalyear
    offset = year0 - ORIGIN_YEAR  # at sim.t=0 we are at START_DATE

    sim_future = sim_base.copy()
    sim_past = sim_base.copy()

    fp_fSR, fp_fRI = _scan(
        sim_future, names, cfg.years_future, cfg.sample_years, "future", offset, R_I, SR
    )
    fp_pSR, fp_pRI = _scan(
        sim_past, names, cfg.years_past, cfg.sample_years, "past", offset, R_I, SR
    )
    return CloneResult(
        eps=eps,
        clone_names=names,
        R_I=R_I,
        SR=SR,
        fp_future_SR=fp_fSR,
        fp_future_RI=fp_fRI,
        fp_past_SR=fp_pSR,
        fp_past_RI=fp_pRI,
    )
