#!/usr/bin/env python3
"""
Aggregate and plot outputs from simulate.py.
Generates:
- Time series panels across regimes
- Summary curves vs kappa
- Variance peak (Var(phi) vs kappa) with CIs
- Lead–lag cross-correlation C_{r,phi}(tau)
"""
import argparse, glob, json, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

def find_runs(pattern="runs/*_summary.json"):
    return sorted(glob.glob(pattern))

def load_summary(path):
    with open(path, "r") as f:
        return json.load(f)

def load_timeseries_from_summary_path(summary_path):
    stem = summary_path.replace("_summary.json", "")
    csv_path = f"{stem}_timeseries.csv"
    return pd.read_csv(csv_path)

def plot_timeseries(df, outpath):
    # one plot per variable to keep it clean
    for col in ["r", "phi", "K", "H_phase"]:
        plt.figure()
        plt.plot(df["t"].values, df[col].values, linewidth=1.3)
        plt.xlabel("time (s)")
        plt.ylabel(col)
        plt.tight_layout()
        plt.savefig(outpath.format(col=col))
        plt.close()

def aggregate_by_kappa(summaries):
    rows = []
    for s in summaries:
        rows.append({
            "kappa": s["kappa"], "seed": s["seed"], "N": s["N"],
            "mean_r": s["mean_r"], "mean_phi": s["mean_phi"], "var_phi": s["var_phi"]
        })
    return pd.DataFrame(rows).sort_values(["kappa","seed"])

def ci_mean(x):
    mu = np.mean(x)
    se = np.std(x, ddof=1) / np.sqrt(len(x))
    return mu, 1.96*se

def plot_summary_vs_kappa(df, out_prefix):
    grp = df.groupby("kappa")
    kappas = sorted(df["kappa"].unique())
    mu_r, ci_r = [], []
    mu_phi, ci_phi = [], []
    mu_v, ci_v = [], []
    for k in kappas:
        g = df[df["kappa"]==k]
        m,c = ci_mean(g["mean_r"].values); mu_r.append(m); ci_r.append(c)
        m,c = ci_mean(g["mean_phi"].values); mu_phi.append(m); ci_phi.append(c)
        m,c = ci_mean(g["var_phi"].values); mu_v.append(m); ci_v.append(c)

    # r
    plt.figure()
    plt.errorbar(kappas, mu_r, yerr=ci_r, fmt='o-')
    plt.xlabel("kappa")
    plt.ylabel("⟨r⟩")
    plt.tight_layout()
    plt.savefig(out_prefix + "_mean_r_vs_kappa.png")
    plt.close()

    # phi
    plt.figure()
    plt.errorbar(kappas, mu_phi, yerr=ci_phi, fmt='o-')
    plt.xlabel("kappa")
    plt.ylabel("⟨Φ_c⟩")
    plt.tight_layout()
    plt.savefig(out_prefix + "_mean_phi_vs_kappa.png")
    plt.close()

    # var phi (variance peak)
    plt.figure()
    plt.errorbar(kappas, mu_v, yerr=ci_v, fmt='o-')
    plt.xlabel("kappa")
    plt.ylabel("Var(Φ_c)")
    plt.tight_layout()
    plt.savefig(out_prefix + "_var_phi_vs_kappa.png")
    plt.close()

def xcorr(a, b, max_lag):
    # normalized cross-correlation for lags in [-max_lag, max_lag]
    a = (a - a.mean()) / (a.std() + 1e-12)
    b = (b - b.mean()) / (b.std() + 1e-12)
    lags = np.arange(-max_lag, max_lag+1)
    corr = []
    for L in lags:
        if L < 0:
            corr.append(np.mean(a[-L:] * b[:L or None]))
        elif L > 0:
            corr.append(np.mean(a[:-L] * b[L:]))
        else:
            corr.append(np.mean(a * b))
    return lags, np.array(corr)

def plot_lead_lag(summary_path, dt, outpath, max_lag_steps=2000):
    stem = summary_path.replace("_summary.json", "")
    r_full = np.load(stem + "_r_full.npy")
    phi_full = np.load(stem + "_phi_full.npy")
    lags, corr = xcorr(r_full, phi_full, max_lag_steps)
    taus = lags * dt
    plt.figure()
    plt.plot(taus, corr, linewidth=1.3)
    plt.xlabel("lag τ (s), positive = r leads Φ_c")
    plt.ylabel("C_{r,Φ_c}(τ)")
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()
    # return argmax lag in seconds
    return taus[np.argmax(corr)], float(np.max(corr))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs_glob", default="runs/*_summary.json")
    ap.add_argument("--outdir", default="figs")
    ap.add_argument("--plot_example_timeseries", action="store_true")
    ap.add_argument("--leadlag_from", default=None, help="path to one summary.json to compute lead–lag")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    paths = find_runs(args.runs_glob)
    if not paths:
        print("No runs found.", file=sys.stderr)
        return

    # optional time series for first path
    if args.plot_example_timeseries:
        df_ts = load_timeseries_from_summary_path(paths[0])
        plot_timeseries(df_ts, os.path.join(args.outdir, "timeseries_{col}.png"))

    # aggregate summaries
    summaries = [load_summary(p) for p in paths]
    df = aggregate_by_kappa(summaries)
    df.to_csv(os.path.join(args.outdir, "aggregate.csv"), index=False)
    # summary curves
    out_prefix = os.path.join(args.outdir, "summary")
    plot_summary_vs_kappa(df, out_prefix)

    # lead–lag on a chosen run
    if args.leadlag_from is not None:
        s = load_summary(args.leadlag_from)
        dt = s["dt"]
        tau_star, c_star = plot_lead_lag(args.leadlag_from, dt, os.path.join(args.outdir, "leadlag.png"))
        print(json.dumps({"tau_star_seconds": tau_star, "corr_peak": c_star}, indent=2))

if __name__ == "__main__":
    main()
