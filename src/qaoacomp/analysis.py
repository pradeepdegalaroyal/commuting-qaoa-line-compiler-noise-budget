"""Summary statistics for the compiler comparison (means, spreads, paired tests, break-even error rates)."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

BASELINE = "qiskit_L3_best16"
PROPOSED = "commute_fuse_multistart"
LADDER = ["ordered_route", "commute_order", "commute_fuse", "commute_fuse_multistart"]


def crossing(lams: np.ndarray, g1: np.ndarray, g2: np.ndarray) -> float:
    """Smallest two-qubit error at which the p = 2 gain falls to (or below) the p = 1 gain (linear interpolation).
    Returns ``nan`` if p = 2 stays better over the whole grid (right-censored)."""
    d = g2 - g1
    for i in range(1, len(lams)):
        if d[i] <= 0:
            if d[i - 1] <= 0:
                return float(lams[i - 1])
            return float(lams[i - 1] + (lams[i] - lams[i - 1]) * d[i - 1] / (d[i - 1] - d[i]))
    return float("nan")


def half_gain_error(lams: np.ndarray, g: np.ndarray) -> float:
    """Two-qubit error at which the gain falls to half of its noiseless value (linear interpolation, nan if never)."""
    target = 0.5 * g[0]
    for i in range(1, len(lams)):
        if g[i] <= target:
            return float(lams[i - 1] + (lams[i] - lams[i - 1]) * (g[i - 1] - target) / (g[i - 1] - g[i]))
    return float("nan")


def per_instance_budgets(df: pd.DataFrame) -> pd.DataFrame:
    """Break-even (p = 2 vs p = 1) and half-gain error rates for each compiler and instance."""
    rows = []
    for (comp, dens, seed), g in df.groupby(["compiler", "density", "seed"]):
        piv = g.pivot(index="lam2", columns="p", values="gain").sort_index()
        lams = piv.index.to_numpy()
        be = crossing(lams, piv[1].to_numpy(), piv[2].to_numpy())
        hg1 = half_gain_error(lams, piv[1].to_numpy())
        hg2 = half_gain_error(lams, piv[2].to_numpy())
        top = float(lams.max())  # right-censoring value: the event did not happen on the grid
        rows.append({"compiler": comp, "density": dens, "seed": seed,
                     "breakeven_p2_vs_p1": be, "half_gain_p1": hg1, "half_gain_p2": hg2,
                     "breakeven_censored": bool(np.isnan(be)), "half_gain_p2_censored": bool(np.isnan(hg2)),
                     "breakeven_c": top if np.isnan(be) else be,
                     "half_gain_p2_c": top if np.isnan(hg2) else hg2})
    return pd.DataFrame(rows)


def mean_std(x: pd.Series) -> str:
    return f"{x.mean():.3f} ± {x.std(ddof=1):.3f}"


def paired(df: pd.DataFrame, col: str, a: str, b: str, **filters) -> dict:
    """Paired comparison of ``col`` between compilers ``a`` and ``b`` over instances (b - a)."""
    sub = df
    for k, v in filters.items():
        sub = sub[sub[k] == v]
    x = sub[sub.compiler == a].sort_values(["density", "seed"])[col].to_numpy()
    y = sub[sub.compiler == b].sort_values(["density", "seed"])[col].to_numpy()
    d = y - x
    out = {"n_pairs": int(len(d)), "mean_a": float(x.mean()), "mean_b": float(y.mean()),
           "mean_diff": float(d.mean()), "std_diff": float(d.std(ddof=1)) if len(d) > 1 else float("nan"),
           "wins_b": int((d > 0).sum()), "ties": int((d == 0).sum())}
    if len(d) > 5 and np.any(d != 0):
        out["wilcoxon_p"] = float(stats.wilcoxon(d).pvalue)
    return out
