<p align="center">
<img src="assets/BABI_logo.svg" alt="BABI logo" width="450"/>
</p>

# `BABI`: Barrier Certificates via Bayesian Inference
[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.17036521-blue)](https://doi.org/10.5281/zenodo.17036521)

`BABI` is a software package to identify and verify forward-invariant sets of unknown dynamical systems with latent states and polynomial dynamics, using only output data. To enable data-efficient inference with rigorous uncertainty quantification, a Bayesian approach is employed: a prior in state-space representation, which allows incorporating structural insights or domain knowledge, is updated with a targeted marginal Metropolis–Hastings (MMH) sampler. Posterior samples are then used to synthesize a barrier certificate and corresponding invariant set via a sum-of-squares (SOS) optimization problem. An additional test set of posterior samples provides finite-sample probabilistic guarantees that the certificate holds for the true, unknown system.

The approach is explained in the paper "Barrier Certificates for Unknown Systems with Latent States and Polynomial Dynamics using Bayesian Inference", available on [IEEExplore](https://doi.org/10.1109/CDC57313.2025.11312207) and as a preprint on [arXiv](https://doi.org/10.48550/arXiv.2504.01807).

## Installation and Requirements
To run the code, the following software is required:
* [MATLAB](https://mathworks.com/products/matlab.html) (tested with R2024b), including:
    - _Parallel Computing Toolbox_ (for Monte Carlo simulations)
    -  _Symbolic Math Toolbox_ (for plotting forward-invariant sets)
* [SOSTOOLS](https://github.com/oxfordcontrol/SOSTOOLS) (tested with v4.03): download and add to the MATLAB path, e.g., `addpath(genpath('Path\to\SOSTOOLS'))`.
* [MOSEK](https://www.mosek.com/) (tested with v11.0.13): default SDP solver. A free academic license is available. Install according to the [documentation](https://docs.mosek.com/latest/install/installation.html) and add MOSEK to the MATLAB path.

MOSEK is the recommended SDP solver. Other solvers (e.g., [SeDuMi](https://github.com/sqlp/sedumi)) may work with SOSTOOLS, but this has not been tested and results may differ.

All experiments reported in the paper were run under Windows 11 with the software versions listed above. The code is expected to work on other platforms with MATLAB support.

> **⚠️ Warning**
>
> With older MATLAB releases (e.g., R2021b), we observed compatibility issues between SOSTOOLS and MOSEK that may cause the validation of barrier certificates to fail.

## Usage and Experiments
The `/experiments` folder contains scripts that reproduce the results presented in Section V of the paper. These serve as examples of how to use the code. Users are encouraged to start by exploring these examples to understand how the code works. Since the implementation heavily relies on [SOSTOOLS](https://github.com/oxfordcontrol/SOSTOOLS), we recommend consulting its documentation for details on the SOS formulation.
1. `autocorrelation.m`: reproduces the autocorrelation analysis of the MMH sampler (Figure 2 in the paper).
2. `exemplary_run.m`: generates the forward-invariant set and barrier function for a single run (Figure 3 in the paper).
3. `monte_carlo_posterior.m`: runs the Monte Carlo simulation of the proposed method (Table 1; column _Posterior (proposed method)_).
4. `monte_carlo_prior.m`: runs the Monte Carlo simulation of the baseline method (Table 1; column _Prior (baseline)_).

Below is an example output from `experiments/exemplary_run.m`, showing the synthesized forward-invariant set $\mathcal{C}$ (red) enclosing the system’s limit cycle (yellow), with circular regions $\mathcal{X}\_{\text{min}}$ (green) and $\mathcal{X}\_{\text{max}}$ (grey).

<img src="assets/forward-invariant_set.svg" alt="forward-invariant set" width="400"/>

## Erratum
Algorithm 1 in the paper (up to and including v0.1.2) contains a bookkeeping error in the marginal Metropolis–Hastings sampler: the chain index was only advanced on acceptance, so rejected proposals were not recorded as repeated states. As a result, the stored samples do not follow the posterior, as states from which proposals are frequently accepted are over-represented. Accordingly, thinning in the paper is described as discarding samples between accepted samples, whereas it should refer to MH iterations. In the corrected version below, every iteration produces a sample (line 8), and burn-in and thinning refer to MH iterations.

```text
Algorithm 1 (corrected): Marginal Metropolis–Hastings sampler

Input:  observations D, model f_θ(·), g_θ(·), p_W(· | θ),
        priors p(θ), p(x_0), proposal distribution q(θ', x_0' | θ[k], x_0[k]),
        initial values θ[1], x_0[1], number of samples K,
        burn-in K_b, thinning parameter k_d
Output: K samples from p(θ, x_{0:T} | D)

 1: set N ← K_b + 1 + (K - 1)(k_d + 1)
 2: get x_{0:T}[1] and compute likelihood
 3: for k = 1, ..., N - 1 do
 4:     propose θ', x_0' ~ q(θ', x_0' | θ[k], x_0[k])
 5:     perform forward simulation to obtain x_{0:T}'
 6:     compute likelihood p(D | θ', x_0')
 7:     compute acceptance probability α
 8:     set (θ, x_{0:T})[k+1] ← (θ', x_{0:T}') with probability α,
        otherwise (θ, x_{0:T})[k+1] ← (θ, x_{0:T})[k]
 9: end for
10: discard burn-in and perform thinning
```

The fix is included from v0.2.0. Since burn-in and thinning now refer to MH iterations rather than accepted proposals, the autocorrelation analysis in Fig. 2 was repeated with the corrected sampler and the sampler settings were tuned accordingly (burn-in `K_b = 1000`, samples per stage `K_stage = 5000`, thinning `k_d = 600`). All experiments were then rerun with these settings.

The qualitative results and conclusions of the paper are unchanged. The SOS program is still feasible and the resulting certificate still valid for the true system in 100/100 runs, while the 99.9 % confidence lower bound 1 - δ changes only slightly, from 98.75 % ± 0.42 % (range 97.41 % to 99.08 %) to 98.80 % ± 0.42 % (range 97.01 % to 99.08 %). The barrier function of the exemplary run in Fig. 3 has different coefficients, but the resulting forward-invariant set is qualitatively the same. The prior baseline is unaffected by the error, and its SOS program remains infeasible in all runs. The code used for the published results is available under tag [v0.1.2](https://github.com/TUM-ITR/BABI/tree/v0.1.2).

## Reference
If you use this software in your research, please cite the following paper:
```
@article{BABI2025,
 title={Barrier Certificates for Unknown Systems with Latent States and Polynomial Dynamics using Bayesian Inference},
 author={Lefringhausen, Robert and Noel Aziz Hanna, Sami Leon and August, Elias and Hirche, Sandra},
 journal={arXiv preprint arXiv:2504.01807},
 year={2025}
}
```
