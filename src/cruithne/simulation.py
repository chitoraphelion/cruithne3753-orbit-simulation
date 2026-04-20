"""REBOUND/WHFast N-body setup and integration for the Cruithne problem.

The Earth-Moon system is collapsed into its barycentre (EMB). Cruithne is a
massless test particle. Output is the heliocentric state history of Cruithne
and EMB sampled on a uniform time grid.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import rebound
from astropy.time import Time
from tqdm import tqdm

from .constants import CRUITHNE_ID, DAY, G, PLANETS
from .horizons import fetch_vector


@dataclass
class SimConfig:
    start_date: str = "1900-01-01"
    end_date: str = "3900-01-01"
    output_step_days: float = 2.0
    integrator_step_days: float = 0.5


@dataclass
class SimResult:
    years: np.ndarray  # decimal years
    t_sec: np.ndarray  # seconds since START_DATE
    Rh: np.ndarray  # (N,3) heliocentric Cruithne position, m
    Vh: np.ndarray  # (N,3) heliocentric Cruithne velocity, m/s
    REh: np.ndarray  # (N,3) heliocentric EMB position, m
    VEh: np.ndarray  # (N,3) heliocentric EMB velocity, m/s
    initial_states: dict[str, tuple[np.ndarray, np.ndarray]]


def build_sim(cfg: SimConfig) -> tuple[rebound.Simulation, dict[str, tuple[np.ndarray, np.ndarray]]]:
    states: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for name, obj_id, id_type, _ in PLANETS:
        r, v = fetch_vector(obj_id, cfg.start_date, id_type=id_type)
        states[name] = (r, v)
    r_ast, v_ast = fetch_vector(CRUITHNE_ID[0], cfg.start_date, id_type=CRUITHNE_ID[1])
    states["Cruithne"] = (r_ast, v_ast)

    sim = rebound.Simulation()
    sim.G = G
    sim.integrator = "whfast"
    sim.dt = cfg.integrator_step_days * DAY
    sim.units = ("m", "s", "kg")

    for name, _, _, mass in PLANETS:
        r, v = states[name]
        sim.add(m=mass, x=r[0], y=r[1], z=r[2], vx=v[0], vy=v[1], vz=v[2], hash=name)
    sim.add(
        m=0.0,
        x=r_ast[0],
        y=r_ast[1],
        z=r_ast[2],
        vx=v_ast[0],
        vy=v_ast[1],
        vz=v_ast[2],
        hash="Cruithne",
    )
    sim.move_to_com()
    return sim, states


def integrate(cfg: SimConfig, progress: bool = True) -> SimResult:
    sim, states = build_sim(cfg)

    jd0 = Time(cfg.start_date, format="iso", scale="tdb").jd
    jd1 = Time(cfg.end_date, format="iso", scale="tdb").jd
    jd = jd0 + np.arange(0, int(np.floor(jd1 - jd0)) + 1, cfg.output_step_days, dtype=float)
    t_sec = (jd - jd0) * DAY
    years = t_sec / (365.25 * DAY) + Time(cfg.start_date, format="iso", scale="tdb").decimalyear

    N = len(t_sec)
    Rh = np.empty((N, 3))
    Vh = np.empty((N, 3))
    REh = np.empty((N, 3))
    VEh = np.empty((N, 3))

    it = tqdm(enumerate(t_sec), total=N, disable=not progress, desc="integrate")
    for k, t in it:
        sim.integrate(t)
        sun = sim.particles["Sun"]
        emb = sim.particles["EMB"]
        ast = sim.particles["Cruithne"]
        Rh[k] = [ast.x - sun.x, ast.y - sun.y, ast.z - sun.z]
        Vh[k] = [ast.vx - sun.vx, ast.vy - sun.vy, ast.vz - sun.vz]
        REh[k] = [emb.x - sun.x, emb.y - sun.y, emb.z - sun.z]
        VEh[k] = [emb.vx - sun.vx, emb.vy - sun.vy, emb.vz - sun.vz]

    return SimResult(
        years=years, t_sec=t_sec, Rh=Rh, Vh=Vh, REh=REh, VEh=VEh, initial_states=states
    )
