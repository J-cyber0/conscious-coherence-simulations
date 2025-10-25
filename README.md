# Conscious Coherence Simulation Framework

Implements the hybrid ODE–Kuramoto ensemble described in  
*The Conscious Coherence Principle: Toward a Third Constant in Relativistic Physics* (Martinez, 2025).

---

## 🚀 Quick Start

### 1️⃣ Install dependencies
```bash
pip install -r requirements.txt
```

---

### 2️⃣ Run parameter sweeps

The following commands sweep the coupling gain `κ` from 0.40 → 1.00  
and run 5 random seeds for each value.  
This creates a `runs/` directory containing CSV and JSON outputs.

#### 🐧 macOS / Linux (bash)
```bash
for k in 0.40 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.90 1.00; do
  for s in 1 2 3 4 5; do
    python simulate.py --kappa $k --seed $s --tag "kap${k}_seed${s}"
  done
done
```

#### 🪟 Windows (PowerShell)
```powershell
foreach ($k in 0.40,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.90,1.00) {
  foreach ($s in 1,2,3,4,5) {
    python simulate.py --kappa $k --seed $s --tag "kap${k}_seed${s}"
  }
}
```

### 3️⃣ Generate plots
Once the sweep completes, analyze the results:

```bash
python analyze.py --runs_glob "runs/*_summary.json" --outdir figs --plot_example_timeseries
```

This will create the following figures inside the `figs/` directory:

- `summary_mean_r_vs_kappa.png`  
- `summary_mean_phi_vs_kappa.png`  
- `summary_var_phi_vs_kappa.png`  
- `timeseries_r.png`, `timeseries_phi.png`, `timeseries_K.png`, `timeseries_H_phase.png`

---

### 4️⃣ (Optional) Lead–lag analysis
Compute the cross-correlation between \( r(t) \) and \( \Phi_c(t) \):

```bash
python analyze.py --leadlag_from runs/kap0.70_seed3_summary.json --outdir figs
```

This produces:

- `figs/leadlag.png`
- Console output showing the lag (in seconds) where \( r(t) \) best predicts \( \Phi_c(t) \).

---

## 🧩 Directory structure
```
├── simulate.py        # Core simulation script
├── analyze.py         # Analysis and plotting utilities
├── requirements.txt   # Dependencies
├── runs/              # Generated CSV + JSON outputs
└── figs/              # Generated figures
```

---

## 🔧 Customization

| Parameter | Description | Default |
|------------|-------------|----------|
| `--N` | Number of oscillators | 200 |
| `--T` | Total simulation time (s) | 100 |
| `--dt` | Integration step size (s) | 0.005 |
| `--alpha`, `--beta`, `--gamma` | Coherence dynamics coefficients | 1.0, 0.5, 0.0 |
| `--K0`, `--kappa` | Baseline and feedback coupling gains | 0.5, variable |
| `--D` | Diffusion constant | 0.02 |
| `--sigma_omega` | Frequency dispersion | 0.5 |

Run `python simulate.py --help` for the full list of options.

---

## 📊 Output Summary

| File | Description |
|------|--------------|
| `*_summary.json` | Run-level summary statistics (means, variances) |
| `*_timeseries.csv` | Thin-sampled time series of \(r(t)\), \(\Phi_c(t)\), \(K(t)\), and \(H_{\mathrm{phase}}(t)\) |
| `aggregate.csv` | Combined summary data across all runs |
| `leadlag.png` | Cross-correlation showing temporal lag between coherence and synchronization |

---

## 🧠 Reproducibility

All figures in the manuscript can be regenerated using the commands above.  
Set `--seed 42` for deterministic initialization.

---

## 📄 License
MIT License © 2025 Juan-Carlos Martinez
