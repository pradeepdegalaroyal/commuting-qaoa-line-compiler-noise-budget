"""Black-and-grey figures (markers and line styles carry the distinctions)."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .analysis import BASELINE, LADDER, PROPOSED

DPI = 200
LABEL = {
    "qiskit_L1": "Qiskit level 1",
    "qiskit_L3": "Qiskit level 3",
    "qiskit_L3_best16": "Qiskit level 3, best of 16 seeds",
    "swap_network": "swap network",
    "ordered_route": "fixed-order routing",
    "commute_order": "+ free term order",
    "commute_fuse": "+ ZZ-SWAP fusion",
    "commute_fuse_multistart": "+ 16 random restarts",
}
FULL_LABEL = dict(LABEL, commute_fuse_multistart="commutation-aware router (full)")
STYLE = {  # (marker, grey, linestyle)
    "qiskit_L3_best16": ("o", "0.45", "-"),
    "swap_network": ("^", "0.7", "-."),
    "commute_fuse_multistart": ("s", "0.0", "-"),
    "qiskit_L3": ("v", "0.3", ":"),
}


def _save(fig: plt.Figure, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def _waterfall(ax, labels: list[str], values: list[float], ylabel: str, title: str, ref: float | None,
               ref_label: str, fmt: str) -> None:
    start, deltas = values[0], np.diff(values)
    xs = np.arange(len(values))
    ax.bar(0, start, color="0.55", edgecolor="k")
    cum = start
    for i, d in enumerate(deltas, start=1):
        ax.bar(i, d, bottom=cum, color="white", edgecolor="k", hatch="//")
        ax.text(i, cum + d + (0.01 * max(values) if d > 0 else -0.01 * max(values)), f"{d:+{fmt}}",
                ha="center", va="bottom" if d > 0 else "top", fontsize=8)
        cum += d
    ax.bar(len(values), cum, color="0.15", edgecolor="k")
    ax.text(0, start, f"{start:{fmt}}", ha="center", va="bottom", fontsize=8)
    ax.text(len(values), cum, f"{cum:{fmt}}", ha="center", va="bottom", fontsize=8)
    if ref is not None:
        ax.axhline(ref, color="k", ls="--", lw=1, label=ref_label)
        ax.legend(fontsize=8, loc="best")
    ax.set_xticks(list(xs) + [len(values)])
    ax.set_xticklabels([labels[0]] + labels[1:] + ["full compiler"], rotation=25, ha="right", fontsize=8)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)


def fig_waterfall(main: pd.DataFrame, lam: float, path: Path) -> None:
    sub = main[(main.p == 2) & (main.lam2 == lam)]
    cx = [sub[sub.compiler == c].cx.mean() for c in LADDER]
    gain = [sub[sub.compiler == c].gain.mean() for c in LADDER]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    labels = [LABEL[c] for c in LADDER]
    ref_cx = sub[sub.compiler == BASELINE].cx.mean()
    ref_g = sub[sub.compiler == BASELINE].gain.mean()
    _waterfall(axes[0], labels, cx, "mean number of cx gates (p = 2)",
               "Compiled two-qubit gates", ref_cx, LABEL[BASELINE], ".0f")
    _waterfall(axes[1], labels, gain, f"mean normalised cut gain (p = 2, noise {lam:g})",
               "Noisy solution quality", ref_g, LABEL[BASELINE], ".3f")
    n = int(main.n.iloc[0])
    fig.suptitle(f"Ablation of the router components, n = {n}, {sub[sub.compiler == LADDER[0]].shape[0]} instances",
                 fontsize=11)
    _save(fig, path)


def fig_dumbbell(main: pd.DataFrame, lam: float, path: Path) -> None:
    sub = main[(main.p == 2) & (main.lam2 == lam)]
    a = sub[sub.compiler == BASELINE].sort_values(["density", "seed"]).reset_index(drop=True)
    b = sub[sub.compiler == PROPOSED].sort_values(["density", "seed"]).reset_index(drop=True)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6), sharey=True)
    y = np.arange(len(a))
    for ax, col, xlabel, title in [
        (axes[0], "cx", "compiled cx gates", "Two-qubit gate count"),
        (axes[1], "gain", "normalised cut gain", f"Noisy solution quality (error {lam:g})"),
    ]:
        for i in y:
            ax.plot([a[col][i], b[col][i]], [i, i], color="0.6", lw=1)
        ax.scatter(a[col], y, marker="o", facecolor="white", edgecolor="k", label=LABEL[BASELINE], zorder=3)
        ax.scatter(b[col], y, marker="s", color="k", label=FULL_LABEL[PROPOSED], zorder=3)
        ax.set_xlabel(xlabel)
        ax.set_title(title, fontsize=10)
        ax.grid(axis="x", color="0.9")
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([f"density {d:g}, graph {s}" for d, s in zip(a.density, a.seed)], fontsize=7)
    axes[0].set_ylabel("instance (p = 2)")
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def fig_violin(main: pd.DataFrame, lams: list[float], path: Path) -> None:
    sub = main[main.p == 2]
    data, labels = [], []
    for lam in lams:
        s = sub[sub.lam2 == lam]
        a = s[s.compiler == BASELINE].sort_values(["density", "seed"]).gain.to_numpy()
        b = s[s.compiler == PROPOSED].sort_values(["density", "seed"]).gain.to_numpy()
        data.append(b - a)
        labels.append(f"{lam:g}")
    fig, ax = plt.subplots(figsize=(7, 4.2))
    parts = ax.violinplot(data, showmeans=False, showextrema=False)
    for body in parts["bodies"]:
        body.set_facecolor("0.85")
        body.set_edgecolor("k")
        body.set_alpha(1)
    ax.boxplot(data, widths=0.15, showfliers=False, medianprops={"color": "k"})
    rng = np.random.default_rng(0)
    for i, d in enumerate(data, start=1):
        ax.scatter(i + rng.uniform(-0.12, 0.12, len(d)), d, s=14, marker="o", facecolor="white",
                   edgecolor="k", zorder=3)
    ax.axhline(0, color="k", ls="--", lw=1)
    ax.set_xticks(range(1, len(lams) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlabel("two-qubit depolarising error per cx")
    ax.set_ylabel("gain(commutation-aware) − gain(Qiskit), paired")
    ax.set_title("Paired gain difference per instance (p = 2, dots = instances)", fontsize=10)
    _save(fig, path)


def fig_gain_vs_error(main: pd.DataFrame, path: Path) -> None:
    dens = sorted(main.density.unique())
    fig, axes = plt.subplots(1, len(dens), figsize=(4.2 * len(dens), 4), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, d in zip(axes, dens):
        sd = main[main.density == d]
        for comp, (mk, grey, _) in STYLE.items():
            if comp == "qiskit_L3":
                continue
            for p, ls in [(1, "--"), (2, "-")]:
                g = sd[(sd.compiler == comp) & (sd.p == p)].groupby("lam2").gain
                m, s = g.mean(), g.std(ddof=1)
                ax.errorbar(np.arange(len(m)), m.values, yerr=s.values, marker=mk, ls=ls, color=grey, ms=4, lw=1.2,
                            capsize=2, label=f"{FULL_LABEL[comp]}, p = {p}" if d == dens[0] else None,
                            mfc="white" if p == 1 else grey)
        ax.axhline(sd.classical_local_search_gain.mean(), color="k", ls=":", lw=1.5,
                   label="classical 1-flip local search" if d == dens[0] else None)
        lams = sorted(sd.lam2.unique())
        ax.set_xticks(range(len(lams)))
        ax.set_xticklabels([f"{l:g}" for l in lams], rotation=60, fontsize=7)
        ax.set_xlabel("two-qubit depolarising error per cx (grid, not to scale)")
        ax.set_title(f"edge density {d:g}, n = {int(sd.n.iloc[0])}", fontsize=10)
        ax.grid(color="0.92")
    axes[0].set_ylabel("normalised cut gain (0 = random, 1 = optimum)", fontsize=8)
    fig.legend(loc="lower center", ncol=3, fontsize=7.5, frameon=False)
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def fig_scaling(scal: pd.DataFrame, lam: float, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8))
    for comp in ["qiskit_L3_best16", "swap_network", "commute_fuse_multistart"]:
        mk, grey, ls = STYLE[comp]
        s = scal[(scal.compiler == comp) & (scal.lam2 == lam)]
        for ax, col in zip(axes, ["cx", "gain"]):
            g = s.groupby("n")[col]
            ax.errorbar(g.mean().index, g.mean().values, yerr=g.std(ddof=1).fillna(0).values, marker=mk, ls=ls,
                        color=grey, capsize=2, label=FULL_LABEL[comp])
    axes[0].set_ylabel("compiled cx gates (p = 2)")
    axes[1].set_ylabel(f"normalised cut gain (error {lam:g})")
    for ax in axes:
        ax.set_xlabel("number of qubits n (edge density 0.5)")
        ax.set_xticks(sorted(scal.n.unique()))
        ax.grid(color="0.92")
    axes[0].legend(fontsize=8)
    fig.suptitle("Scaling with qubit count (mean ± std over graphs)", fontsize=11)
    _save(fig, path)
