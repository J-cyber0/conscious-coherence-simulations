#!/usr/bin/env python3
"""
Hybrid Kuramoto + coherence-state simulation.
- Phases: Euler–Maruyama
- Phi_c:  Heun (RK2)
Saves time series (thin-sampled) and summary CSVs.
"""
import argparse, json, math, os, sys
import numpy as np
import pandas as pd

def wrap_pi(x):
    # wrap to (-pi, pi]
    return (x + np.pi) % (2*np.pi) - np.pi

def phase_entropy_thetas(thetas, bins=36):
    hist, _ = np.histogram(thetas, bins=bins, range=(-np.pi, np.pi), density=True)
    p = hist + 1e-12
    p = p / p.sum()
    return -np.sum(p * np.log(p))

def simulate(
    N=200, T=100.0, dt=5e-3, burn=10.0,
    K0=0.5, kappa=0.6, alpha=1.0, beta=0.5, gamma=0.0,
    D=0.02, sigma_omega=0.5, phi0=0.2,
    freq_dist="normal", seed=42,
    thin=10, outdir="runs", tag=None,
    do_ablate_no_feedback=False, do_shuffle_surrogate=False
):
    rng = np.random.default_rng(seed)
    steps = int(T/dt)
    burn_steps = int(burn/dt)

    # intrinsic frequencies
    if freq_dist == "normal":
        omegas = rng.normal(loc=0.0, scale=sigma_omega, size=N)
    elif freq_dist == "lorentz":
        # Cauchy with scale sigma_omega
        omegas = rng.standard_cauchy(size=N) * sigma_omega
    else:
        raise ValueError("freq_dist must be 'normal' or 'lorentz'")

    thetas = rng.uniform(-np.pi, np.pi, size=N)
    phi = float(np.clip(phi0, 0.0, 1.0))

    # pre-alloc thin-sampled collectors
    K_list, r_list, phi_list, H_list, t_list = [], [], [], [], []
    # for lead–lag
    r_store, phi_store = [], []

    # integration
    sqrt_2Ddt = math.sqrt(2.0 * D * dt)
    for n in range(steps):
        # mean field
        z = np.exp(1j * thetas).mean()
        r = np.abs(z)
        psi = np.angle(z)

        # adaptive coupling
        Kt = K0 if do_ablate_no_feedback else (K0 + kappa * phi)

        # SDE: Euler–Maruyama
        drift = omegas + Kt * r * np.sin(psi - thetas)
        noise = rng.normal(0.0, 1.0, size=N) * sqrt_2Ddt
        thetas = wrap_pi(thetas + drift * dt + noise)

        # ODE: Heun (explicit RK2) for phi
        S = 1.0 - r
        Deff = D + gamma * (1.0 - r)
        f1 = alpha * (1.0 - phi) * S - beta * phi * Deff
        phi_pred = np.clip(phi + dt * f1, 0.0, 1.0)
        # evaluate slope at predictor
        # re-use S, Deff at same step (quasi-static within dt)
        f2 = alpha * (1.0 - phi_pred) * S - beta * phi_pred * Deff
        phi = np.clip(phi + 0.5 * dt * (f1 + f2), 0.0, 1.0)

        # optional surrogate: break coupling in the measurement channel
        r_meas = r
        if do_shuffle_surrogate:
            # shuffle phases before measuring order; not affecting dynamics
            shuf = np.array(thetas)
            rng.shuffle(shuf)
            z_surr = np.exp(1j * shuf).mean()
            r_meas = np.abs(z_surr)

        # diagnostics
        if (n % thin) == 0:
            t = n * dt
            H = phase_entropy_thetas(thetas)
            K_list.append(Kt)
            r_list.append(r_meas)
            phi_list.append(phi)
            H_list.append(H)
            t_list.append(t)
        # store for lead–lag
        r_store.append(r)
        phi_store.append(phi)

    # discard burn-in for summaries
    keep = np.array(t_list) >= burn
    r_arr = np.array(r_list)[keep]
    phi_arr = np.array(phi_list)[keep]
    H_arr = np.array(H_list)[keep]
    K_arr = np.array(K_list)[keep]
    t_arr = np.array(t_list)[keep]

    # summary metrics
    summary = {
        "N": N, "T": T, "dt": dt, "burn": burn,
        "K0": K0, "kappa": kappa, "alpha": alpha, "beta": beta, "gamma": gamma,
        "D": D, "sigma_omega": sigma_omega, "phi0": phi0,
        "freq_dist": freq_dist, "seed": seed,
        "mean_r": float(r_arr.mean()),
        "mean_phi": float(phi_arr.mean()),
        "var_phi": float(phi_arr.var(ddof=1)),
        "mean_H": float(H_arr.mean()),
        "do_ablate_no_feedback": bool(do_ablate_no_feedback),
        "do_shuffle_surrogate": bool(do_shuffle_surrogate),
    }

    os.makedirs(outdir, exist_ok=True)
    run_id = tag or f"N{N}_kap{kappa:.3f}_a{alpha:.3f}_b{beta:.3f}_K0{K0:.3f}_D{D:.3f}_seed{seed}"
    # time series
    df_ts = pd.DataFrame({
        "t": t_arr, "r": r_arr, "phi": phi_arr, "H_phase": H_arr, "K": K_arr
    })
    df_ts.to_csv(os.path.join(outdir, f"{run_id}_timeseries.csv"), index=False)
    # summary
    with open(os.path.join(outdir, f"{run_id}_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    # lead–lag arrays (full resolution, no burn removal here)
    np.save(os.path.join(outdir, f"{run_id}_r_full.npy"), np.array(r_store, dtype=np.float32))
    np.save(os.path.join(outdir, f"{run_id}_phi_full.npy"), np.array(phi_store, dtype=np.float32))

    return summary

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, default=200)
    p.add_argument("--T", type=float, default=100.0)
    p.add_argument("--dt", type=float, default=5e-3)
    p.add_argument("--burn", type=float, default=10.0)
    p.add_argument("--K0", type=float, default=0.5)
    p.add_argument("--kappa", type=float, default=0.6)
    p.add_argument("--alpha", type=float, default=1.0)
    p.add_argument("--beta", type=float, default=0.5)
    p.add_argument("--gamma", type=float, default=0.0)
    p.add_argument("--D", type=float, default=0.02)
    p.add_argument("--sigma_omega", type=float, default=0.5)
    p.add_argument("--phi0", type=float, default=0.2)
    p.add_argument("--freq_dist", choices=["normal","lorentz"], default="normal")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--thin", type=int, default=10)
    p.add_argument("--outdir", type=str, default="runs")
    p.add_argument("--tag", type=str, default=None)
    p.add_argument("--ablate_no_feedback", action="store_true")
    p.add_argument("--shuffle_surrogate", action="store_true")
    args = p.parse_args()

    summary = simulate(
        N=args.N, T=args.T, dt=args.dt, burn=args.burn,
        K0=args.K0, kappa=args.kappa, alpha=args.alpha, beta=args.beta, gamma=args.gamma,
        D=args.D, sigma_omega=args.sigma_omega, phi0=args.phi0,
        freq_dist=args.freq_dist, seed=args.seed, thin=args.thin,
        outdir=args.outdir, tag=args.tag,
        do_ablate_no_feedback=args.ablate_no_feedback,
        do_shuffle_surrogate=args.shuffle_surrogate
    )
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
