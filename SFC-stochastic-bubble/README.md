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


## Installation

Requires **Julia ≥ 1.9**.


## Contact

Adrien Nguyen-Huu
University of Montpellier


---
