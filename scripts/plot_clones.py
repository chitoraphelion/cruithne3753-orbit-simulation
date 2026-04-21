"""Histograms of first-passage times for the clone ensembles.

Reads CSVs produced by ``cruithne-clones`` (one per epsilon) and emits:
  clones_epsX_future_SR.png  / _future_RI.png / _past_SR.png / _past_RI.png
  hist_first_entries_{Venus,Earth,Mars}_{SR,RI}.png  — pooled across epsilons
"""

from __future__ import annotations

import argparse
import glob
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cruithne.plot_style import COLORS, apply as apply_style

OUT = Path("plots/generated")


PLANET_COLOR = {"Venus": COLORS["venus"], "Earth": COLORS["earth"], "Mars": COLORS["mars"]}


def _canon(df: pd.DataFrame) -> pd.DataFrame:
    mapping = {}
    for c in df.columns:
        lc = c.lower().replace(" ", "")
        if re.search(r"^past_*sr", lc):
            mapping[c] = "past_SR"
        elif re.search(r"^past_*ri", lc):
            mapping[c] = "past_RI"
        elif re.search(r"^future_*sr", lc):
            mapping[c] = "future_SR"
        elif re.search(r"^future_*ri", lc):
            mapping[c] = "future_RI"
    return df.rename(columns=mapping)


def plot_per_clone(df: pd.DataFrame, eps_label: str) -> None:
    for phase, col in (("future", "future_SR"), ("future", "future_RI"),
                        ("past", "past_SR"), ("past", "past_RI")):
        sub = df.pivot(index="clone", columns="planet", values=col)
        fig, ax = plt.subplots(figsize=(9.0, 3.3))
        x = np.arange(len(sub.index))
        for i, planet in enumerate(("Venus", "Earth", "Mars")):
            if planet in sub.columns:
                ax.scatter(x + (i - 1) * 0.2, sub[planet], s=18,
                           color=PLANET_COLOR[planet], label=planet, alpha=0.9)
        ax.axhline(0, lw=0.6, color=COLORS["muted"])
        ax.set_xlabel("индекс клона")
        ax.set_ylabel("время 1-го входа (лет от 1997)")
        ax.set_title(f"Первые входы: {eps_label}, {col}")
        ax.legend()
        fig.tight_layout()
        tag = col
        fig.savefig(OUT / f"clones_eps{eps_label}_{tag}.png")
        plt.close(fig)


def plot_pooled_histograms(all_df: pd.DataFrame) -> None:
    for planet in ("Venus", "Earth", "Mars"):
        sub = all_df[all_df["planet"] == planet]
        for col in ("past_SR", "future_SR", "past_RI", "future_RI"):
            pass  # handled in combined histogram below

        for metric in ("SR", "RI"):
            past = sub[f"past_{metric}"].dropna().to_numpy()
            future = sub[f"future_{metric}"].dropna().to_numpy()
            data = np.concatenate([past, future])
            fig, ax = plt.subplots(figsize=(8.0, 3.3))
            if len(data) > 0:
                ax.hist(data, bins=40, color=PLANET_COLOR[planet], alpha=0.85)
            ax.axvline(0, lw=0.6, color=COLORS["muted"])
            ax.set_xlabel("время 1-го входа (лет от 1997)")
            ax.set_ylabel("число клонов")
            label = "5·R_I" if metric == "SR" else "R_I"
            ax.set_title(f"{planet}: время 1-го входа в сферу {label}")
            fig.tight_layout()
            fig.savefig(OUT / f"hist_first_entries_{planet}_{metric}.png")
            plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--dir", default="data/clones")
    a = p.parse_args()

    apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    files = sorted(glob.glob(str(Path(a.dir) / "clones_eps*_first_passages.csv")))
    if not files:
        raise FileNotFoundError(f"no clone CSVs in {a.dir}; run `cruithne-clones` first")
    frames = []
    for f in files:
        df = _canon(pd.read_csv(f))
        m = re.search(r"eps([0-9eE+\-.]+)", Path(f).stem)
        eps_label = m.group(1) if m else "x"
        plot_per_clone(df, eps_label)
        frames.append(df)
    plot_pooled_histograms(pd.concat(frames, ignore_index=True))
    print(f"[clones-plot] wrote figures to {OUT}")


if __name__ == "__main__":
    main()
