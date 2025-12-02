# Optimal Sample Size for Multiple Discoveries  
**Sequential Multiple Testing • E-Processes • Stochastic Target Problems**

This repository accompanies an ongoing **math/stat research project** on
optimal sample size and optimal allocation of sampling effort in a
multiple-testing setting with strict error control.

All detailed proofs and derivations are available in the Overleaf document
(link in the repo). This README provides an accessible mathematical overview.


---

# 1. Model Overview

We observe **N independent arms**, each generating a controlled continuous-time
Gaussian signal:

<img src="https://latex.codecogs.com/png.latex?dX_t^i%20=%20\theta^i%20a_t^i%20dt%20+%20\sigma%20\sqrt{a_t^i}\,dW_t^i" />

Here:

- \(a_t^i\) is the **sampling allocation** at time \(t\),
- \(\sum_i a_t^i = 1\),
- each arm satisfies a **simple hypothesis test**:

<img src="https://latex.codecogs.com/png.latex?H_0^i:\%20\theta^i%20=%20\theta_0%20\quad\text{vs.}\quad%20H_1^i:\%20\theta^i%20=%20\theta_1%20=\theta_0+\Delta" />

For exposition, we adopt the homogeneous case:

<img src="https://latex.codecogs.com/png.latex?\sigma^i=\sigma,\quad\theta_0^i=\theta_0,\quad\theta_1^i=\theta_1" />


---

# 2. Likelihood Ratio and E-Processes

The sequential log-likelihood ratio for arm \(i\) is:

<img src="https://latex.codecogs.com/png.latex?L_t^i%20=%20\int_0^t%20\frac{\Delta\,a_s^i}{\sigma^2}\,dX_s^i%20-%20\frac12%20\int_0^t%20\frac{\Delta^2(a_s^i)^2}{\sigma^2}\,ds" />

The associated **e-process** is:

<img src="https://latex.codecogs.com/png.latex?E_t^i%20=%20e^{L_t^i}" />

Under the null:

<img src="https://latex.codecogs.com/png.latex?dE_t^i%20=%20E_t^i%20\frac{\Delta}{\sigma}\,a_t^i\,dW_t^i" />

so \(E_t^i\) is a **nonnegative martingale**, giving **anytime-valid** Type I error control.


---

# 3. Multiple Testing: FWER and FWPD

## 3.1 FWER-\(\alpha\)

We reject \(H_0^i\) as soon as:

<img src="https://latex.codecogs.com/png.latex?E_t^i%20\ge%20\frac1{\alpha}" />

This ensures:

<img src="https://latex.codecogs.com/png.latex?\mathbb{P}_0(\text{any false discovery})%20\le%20\alpha" />

at **all times** (Ville’s inequality).


## 3.2 FWPD-(1−β)

We require that, by time \(T\):

<img src="https://latex.codecogs.com/png.latex?\mathbb{P}(\min_{i\in\mathcal{H}_1}E_T^i%20\ge%201/\alpha)%20\ge%201-\beta" />

This is the **Family-Wise Probability of Detection**:  
all true alternatives must be detected with probability ≥ \(1-\beta\).


---

# 4. Stochastic Target Formulation

Let \(\mathbf{E}_t = (E_t^1,\dots,E_t^N)\).  
Define the target event:

<img src="https://latex.codecogs.com/png.latex?\mathcal{E}^\alpha(\mathbf{E}_T)%20=%20\{\min_{i\in\mathcal{H}_1}E_T^i%20\ge%201/\alpha\}" />

The value function is:

<img src="https://latex.codecogs.com/png.latex?V(t,e)%20=%20\inf\Big\{T%20\ge%20t:%20\exists\,a_s\in\Delta_N,\;\mathbb{P}_1(\mathcal{E}^\alpha(\mathbf{E}_T))%20\ge%201-\beta\Big\}" />

Thus:

- **control**: the allocation process \(a_t\),
- **state**: the e-values,
- **target**: all alternatives detected,
- **objective**: minimize the horizon \(T\).

This is a **time-optimal stochastic target problem** with a probabilistic constraint.


---

# 5. Dynamic Programming and HJB Equation

Each e-process under \(\mathbb{P}_1\) follows:

<img src="https://latex.codecogs.com/png.latex?dE_t^i%20=%20E_t^i\,\frac{\Delta}{\sigma}\,a_t^i\,dW_t^i%20+%20\text{drift}" />

The generator acting on a smooth test function \(\phi\) includes terms:

<img src="https://latex.codecogs.com/png.latex?\frac12%20\sum_{i=1}^N%20\left(\frac{\Delta}{\sigma}\right)^2(a^i)^2e_i^2\,\partial_{e_ie_i}^2\phi" />

The dynamic programming yields the HJB equation:

<img src="https://latex.codecogs.com/png.latex?\partial_tV%20+%20\inf_{a\in\Delta_N}\{\mathcal{L}^aV\}%20=%200" />

A key qualitative result:

> **Optimal allocation is bang–bang:**  
> at each time, sample one arm at full rate.


---

# 6. Literature

This project sits at the intersection of:

- **Continuous-time bandits**  
  Mandelbaum (1987), El Karoui & Karatzas (1994)
- **Stochastic target problems**  
  Soner, Touzi & Zhang; Bouchard & Vu
- **Anytime-valid tests and e-values**  
  Ramdas, Wang, Xu et al.

It addresses the question:

> *How much sampling is needed to achieve multiple discoveries under strict
> FWER control, using e-processes, in a continuous-time Gaussian setting?*


---

# 7. Contents of the Repository

- Conceptual notes about the model  
- Derivation sketches  
- HJB structure and intuition  
- Simple numerical experiment scripts (optional)  
- Link to the full Overleaf document


---

# 8. Citation

Nguyen-Huu, A. (2025).  
**Optimal Sample Size for Multiple Discoveries**.  
Working paper.
