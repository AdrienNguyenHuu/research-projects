Stock–Flow Consistent Model for Asset Price Bubbles

This repository accompanies an ongoing research project on macro–financial instability, speculative dynamics, and asset price bubbles within a stock–flow consistent (SFC) framework.

The work studies how credit dynamics, investment behavior, wages, employment, and speculative financial flows interact with a stochastic asset market. The model includes diffusion, jumps with state-dependent intensity, and feedback mechanisms linking real activity to speculative pressure.

The objective is to understand:

when and why asset price bubbles emerge endogenously;

how credit-financed speculation amplifies financial instability;

how crashes propagate to the real economy;

which parameters or institutional features increase (or reduce) systemic fragility.

The model integrates a Keen-type SFC core with a stochastic price process. Asset prices evolve through a diffusion term and two types of jumps (upward and downward), whose intensity depends on the macroeconomic environment and on speculative pressure. This creates a feedback loop: speculation fuels price acceleration, higher prices relax financing conditions, and adverse shifts trigger abrupt corrections.

To study these dynamics, the project relies on numerical simulation, including:
– forward simulations of the coupled system;
– Monte-Carlo estimation of crisis probabilities;
– detection of critical transitions through thresholds on employment, leverage, and speculative flows;
– sensitivity analysis on key financial and behavioral parameters.

What this repository contains:
– A description of the problem, modeling choices, and motivation.
– Julia scripts implementing the full coupled SFC–jump-diffusion system.
– Tools for trajectory simulation, crisis detection, and Monte-Carlo analysis.
– Numerical illustrations of bubble growth, crash events, and macro-financial feedback.
– Plots and stored simulation outputs for selected parameter sets.
– A link to the Overleaf document or PDF containing the full research draft (narrative description, model structure, figures, references).
– Notes, drafts, intermediate results, and work-in-progress materials.

This is not a software package but a research inquiry. The aim is to refine the modeling assumptions, explore stability and instability regimes, and better understand the conditions under which financial speculation can destabilize the real economy.

If you wish to propose ideas, suggest improvements, or discuss theoretical or numerical aspects, feel free to open an issue or contact the authors.

Citation

Grasselli, M. R., & Nguyen-Huu, A. (2025).
"A Stock–Flow Consistent Model for Asset Price Bubbles."
