# SFC Stochastic Bubble Model (Jump-SDE)

Julia code to simulate a stock-flow consistent macro model with an endogenous speculative/bubble component
and stochastic asset market dynamics (diffusion + state-dependent jumps).

## Model overview
State vector (baseline):
- ω : wage share
- e : employment rate
- m : money/firm deposits (definition in code)
- ℓ : loans/debt (definition in code)
- S : asset price index (jump-diffusion)
- μ : latent factor driving the interest/risk-premium channel (OU-type with jump terms)

Core blocks:
- Real sector: Phillips curve + investment function + growth and nominal growth
- Financial sector: endogenous risk-premium channel + speculative demand Ψ(·)
- Asset market: diffusion + asymmetric jump mechanism with intensities depending on Ψ

Crisis event definition:
- e(t) < crisis_threshold  OR  (ℓ(t) - m(t)) < leverage_gap_threshold

## Installation
Requirements:
- Julia >= 1.10 (recommended)

Clone and instantiate:
```bash
git clone https://github.com/AdrienNguyenHuu/research-projects.git
cd research-projects/SFC-stochastic-bubble
julia --project=. -e 'using Pkg; Pkg.instantiate()'
