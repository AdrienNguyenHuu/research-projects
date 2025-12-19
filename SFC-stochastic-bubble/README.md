# SFC-Stochastic-Bubble

This repository contains Julia code to simulate a **stock–flow consistent (SFC) macro-financial model**
with endogenous asset-price bubbles, stochastic dynamics, and crisis events.

The model combines:
- SFC macroeconomic accounting,
- nonlinear real–financial feedbacks,
- stochastic volatility and jump processes,
- Monte Carlo estimation of crisis probabilities.

The code accompanies ongoing research on endogenous financial instability and bubble dynamics.

---

## Model overview

State variables include:
- wage share,
- employment rate,
- firm money balances,
- firm leverage,
- asset price,
- market trend indicator.

The asset price follows a **jump-diffusion process** with state-dependent intensities,
capturing speculative booms and crashes.

Crises are identified through:
- employment collapses,
- excessive leverage imbalances,
- numerical blow-ups.

---

## Repository structure

SFC-stochastic-bubble/

├─ Project.toml

├─ Manifest.toml

├─ src/

│ ├─ SFCStochasticBubble.jl

│ ├─ model.jl

│ ├─ simulate.jl

│ ├─ montecarlo.jl

│ ├─ reporting.jl

│ └─ io.jl

├─ scripts/

│ ├─ run_single_trajectory.jl

│ ├─ run_scenarios.jl

│ ├─ run_sweeps_1d.jl

│ └─ run_heatmaps_2d.jl

├─ output/ (generated, not tracked)

└─ README.md


- `src/` contains the package code.
- `scripts/` reproduces figures and Monte Carlo experiments.
- `output/` is created automatically.

---

## Installation

Requires **Julia ≥ 1.9**.

From the repository root:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()

Usage
Load the package
using SFCStochasticBubble

Run a single trajectory
julia --project=. scripts/run_single_trajectory.jl

Run benchmark scenarios
julia --project=. scripts/run_scenarios.jl

Run Monte Carlo sweeps
julia --project=. scripts/run_sweeps_1d.jl

Generate crisis probability heatmaps
julia --project=. scripts/run_heatmaps_2d.jl

```

All figures and CSV outputs are written to output/.

## Reproducibility

- Random seeds are controlled at the script level.

- The package itself never sets the RNG seed.

- Manifest.toml is provided to ensure exact reproducibility.

## Dependencies

Main dependencies include:

- DifferentialEquations.jl

- JumpProcesses.jl

- Plots.jl

- DataFrames.jl

- CSV.jl

- See Project.toml for the full list.

- Citation

If you use this code, please cite the associated paper or working paper
(links to be added).

## Contact

Adrien Nguyen-Huu
University of Montpellier


---

## 3) `.gitignore`

```gitignore
/output/
/data/
.vscode/
*.log
```
