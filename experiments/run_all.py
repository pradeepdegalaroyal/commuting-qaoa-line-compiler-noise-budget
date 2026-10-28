"""Reproduce every result and figure:  python experiments/run_all.py   (add --quick for a < 1 minute smoke run)."""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
import time
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from qaoacomp import analysis as an  # noqa: E402
from qaoacomp import plots  # noqa: E402
from qaoacomp.experiments import FULL, QUICK, run_main, run_scaling  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="tiny smoke run (n = 6, 2 graphs per density)")
    ap.add_argument("--reuse", action="store_true", help="re-analyse and re-plot the saved raw rows in results/")
    ap.add_argument("--out", default=None, help="output root (default: repository root; quick runs use a temp dir)")
    args = ap.parse_args()
    cfg = QUICK if args.quick else FULL
    out = Path(args.out) if args.out else (Path(tempfile.mkdtemp(prefix="qaoacomp_quick_")) if args.quick else ROOT)
    res, fig = out / "results", out / "figures"
    res.mkdir(parents=True, exist_ok=True)
    fig.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    prev_runtime = None
    if args.reuse:
        main_df, scal_df = pd.read_csv(res / "main_runs.csv"), pd.read_csv(res / "scaling_runs.csv")
        prev_runtime = json.loads((res / "key_results.json").read_text()).get("runtime_seconds")
    else:
        main_df = run_main(cfg)
        main_df.to_csv(res / "main_runs.csv", index=False)
        scal_df = run_scaling(cfg)
        scal_df.to_csv(res / "scaling_runs.csv", index=False)

    lam_ref = 0.01
    budgets = an.per_instance_budgets(main_df)
    budgets.to_csv(res / "error_budgets_per_instance.csv", index=False)

    # ---- summary tables -------------------------------------------------------------------------------
    rows = []
    for (comp, dens, p), g in main_df[main_df.lam2 == lam_ref].groupby(["compiler", "density", "p"]):
        rows.append({"compiler": comp, "density": dens, "p": p, "n_graphs": len(g), "cx_mean": g.cx.mean(),
                     "cx_std": g.cx.std(ddof=1), "depth_mean": g.depth.mean(), "depth_std": g.depth.std(ddof=1),
                     "gain_mean": g.gain.mean(), "gain_std": g.gain.std(ddof=1),
                     "ideal_gain_mean": g.ideal_gain.mean()})
    pd.DataFrame(rows).to_csv(res / "summary_by_density.csv", index=False)

    brows = []
    for (comp, dens), g in budgets.groupby(["compiler", "density"]):
        brows.append({"compiler": comp, "density": dens, "n_graphs": len(g),
                      "breakeven_mean_censored_at_grid_max": g.breakeven_c.mean(),
                      "breakeven_std": g.breakeven_c.std(ddof=1),
                      "breakeven_n_censored": int(g.breakeven_censored.sum()),
                      "half_gain_p2_mean_censored_at_grid_max": g.half_gain_p2_c.mean(),
                      "half_gain_p2_std": g.half_gain_p2_c.std(ddof=1),
                      "half_gain_p2_n_censored": int(g.half_gain_p2_censored.sum())})
    pd.DataFrame(brows).to_csv(res / "summary_error_budgets.csv", index=False)

    lam_grid = cfg["lams"]
    key = {"lam_ref": lam_ref, "baseline": an.BASELINE, "proposed": an.PROPOSED,
           "classical_local_search_gain_mean": float(main_df.drop_duplicates(["density", "seed"])
                                                     .classical_local_search_gain.mean())}
    ref = main_df[(main_df.lam2 == lam_ref) & (main_df.p == 2)]
    key["paired_cx"] = {str(d): an.paired(ref, "cx", an.BASELINE, an.PROPOSED, density=d) for d in cfg["densities"]}
    key["paired_cx"]["all"] = an.paired(ref, "cx", an.BASELINE, an.PROPOSED)
    key["paired_gain"] = {str(d): an.paired(ref, "gain", an.BASELINE, an.PROPOSED, density=d) for d in cfg["densities"]}
    key["paired_gain"]["all"] = an.paired(ref, "gain", an.BASELINE, an.PROPOSED)
    key["paired_gain_vs_swap_network"] = an.paired(ref, "gain", "swap_network", an.PROPOSED)
    key["ladder_p2"] = {c: {"cx_mean": float(ref[ref.compiler == c].cx.mean()),
                            "gain_mean": float(ref[ref.compiler == c].gain.mean())} for c in an.LADDER}
    key["breakeven_overall_censored_at_grid_max"] = {
        c: {"mean": float(g.breakeven_c.mean()), "std": float(g.breakeven_c.std(ddof=1)),
            "n_censored": int(g.breakeven_censored.sum()), "n": int(len(g)),
            "half_gain_p2_mean": float(g.half_gain_p2_c.mean()), "half_gain_p2_std": float(g.half_gain_p2_c.std(ddof=1)),
            "half_gain_p2_n_censored": int(g.half_gain_p2_censored.sum())}
        for c, g in budgets.groupby("compiler")}
    key["paired_breakeven"] = {str(d): an.paired(budgets, "breakeven_c", an.BASELINE, an.PROPOSED, density=d)
                               for d in cfg["densities"]}
    key["paired_breakeven"]["all"] = an.paired(budgets, "breakeven_c", an.BASELINE, an.PROPOSED)
    key["paired_half_gain_error"] = an.paired(budgets, "half_gain_p2_c", an.BASELINE, an.PROPOSED)
    a_ = ref[ref.compiler == an.BASELINE].sort_values(["density", "seed"])
    b_ = ref[ref.compiler == an.PROPOSED].sort_values(["density", "seed"])
    rho = stats.spearmanr(b_.cx.to_numpy() - a_.cx.to_numpy(), b_.gain.to_numpy() - a_.gain.to_numpy())
    key["spearman_dcx_vs_dgain"] = {"rho": float(rho.statistic), "p": float(rho.pvalue)}
    cx_red = 1 - b_.cx.to_numpy() / a_.cx.to_numpy()
    key["cx_reduction_fraction"] = {"mean": float(cx_red.mean()), "std": float(cx_red.std(ddof=1)),
                                    "by_density": {str(d): float(np.mean(1 - b_[b_.density == d].cx.to_numpy()
                                                                         / a_[a_.density == d].cx.to_numpy()))
                                                   for d in cfg["densities"]}}
    key["noiseless_check_max_abs_diff"] = float(
        (main_df[main_df.lam2 == 0].gain - main_df[main_df.lam2 == 0].ideal_gain).abs().max())
    key["lam_grid"] = lam_grid
    key["runtime_seconds"] = None
    key["quick"] = bool(args.quick)

    # ---- figures --------------------------------------------------------------------------------------
    plots.fig_waterfall(main_df, lam_ref, fig / "fig1_waterfall_ablation.png")
    plots.fig_dumbbell(main_df, lam_ref, fig / "fig2_dumbbell_instances.png")
    plots.fig_violin(main_df, [l for l in lam_grid if l in (0.002, 0.005, 0.01, 0.02) or args.quick and l > 0],
                     fig / "fig3_violin_gain_difference.png")
    plots.fig_gain_vs_error(main_df, fig / "fig4_gain_vs_error.png")
    plots.fig_scaling(scal_df, cfg["scaling_lams"][-1], fig / "fig5_scaling.png")

    srows = []
    for (comp, n, lam), g in scal_df.groupby(["compiler", "n", "lam2"]):
        srows.append({"compiler": comp, "n": n, "lam2": lam, "n_graphs": len(g), "cx_mean": g.cx.mean(),
                      "cx_std": g.cx.std(ddof=1), "depth_mean": g.depth.mean(), "gain_mean": g.gain.mean(),
                      "gain_std": g.gain.std(ddof=1)})
    pd.DataFrame(srows).to_csv(res / "summary_scaling.csv", index=False)

    key["runtime_seconds"] = prev_runtime if args.reuse else round(time.time() - t0, 1)
    (res / "key_results.json").write_text(json.dumps(key, indent=2))
    print(f"done in {key['runtime_seconds']} s -> {res}")


if __name__ == "__main__":
    main()
