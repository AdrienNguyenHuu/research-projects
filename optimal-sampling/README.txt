Optimal Sample Size for Multiple Discoveries
===========================================

This repository accompanies an ongoing research project on optimal sample
size and optimal allocation of sampling effort in sequential multiple testing.

The work studies how to guarantee:
1) strict family-wise error control (FWER), at any time;
2) high probability of detecting all true alternatives (FWPD);
3) minimal total sampling budget T.

The problem is formulated in continuous time. Evidence in each arm is
represented by an e-process, which provides anytime-valid testing without
inflating the error rate. The goal is to determine the smallest sampling
horizon for which it is possible to detect all non-null arms with probability
at least 1-beta, while ensuring FWER <= alpha.

To achieve this, we cast the question as a backward stochastic target
problem: the target consists of “all true alternatives declared discoveries”.
The control is the allocation of sampling effort across arms. The horizon T
is the variable to be minimized. This leads naturally to a dynamic
programming principle and to a characterization of optimal allocations. In
particular, optimal allocations tend to be “bang-bang”: at each moment,
sampling should focus on a single arm.

What this repository contains:
- A short description of the problem and motivation.
- Notes, drafts, and intermediate results.
- Numerical illustrations of detection behavior under various allocations.
- A link to the Overleaf document detailing the full mathematical
  formulation (state dynamics, definitions, proofs, and the corresponding
  dynamic programming equation).

This is not a software package but a research inquiry. Comments, suggestions,
and remarks are welcome. The goal is to refine the formulation, clarify the
dynamic structure of the problem, and explore implications for optimal
multiple testing under strict error control.

If you wish to contribute ideas, point out connections to existing work, or
discuss extensions, feel free to open an issue or contact the author.


Citation
--------

Nguyen-Huu, A. (2025).
"Optimal Sample Size for Multiple Discoveries".
