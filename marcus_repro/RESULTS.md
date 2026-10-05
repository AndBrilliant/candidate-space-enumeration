# Marcus reference audit results

These are **independent replication diagnostics**, not evidential significances.

## Deterministic checks

Independent reconstruction reproduces the manuscript's heavy-sector headline
residuals:

- charm: **+0.255914%**
- bottom: **-0.081531%**
- (c b/G^2): **+0.174175%**

## Six-class grammar

The independent grammar rebuild gives exactly:

- P = 324
- T = 252
- A = 728
- S = 426
- D = 426
- K = 1584
- total = **3740**

For the manuscript spectrum the tolerance ladder is **[0, 1, 3, 6]** at
0.1%, 0.25%, 0.5%, and 1.0%.

## Whole-grammar replay

An independent 10,000-universe replay with seed 20261001 reproduces:

| tolerance | mean cells | P(>=1) | P(>=2) | max |
|---|---:|---:|---:|---:|
| 0.1% | 0.4416 | 0.2745 | 0.0961 | 16 |
| 0.25% | 1.1308 | 0.5512 | 0.2781 | 24 |
| 0.5% | 2.2844 | 0.7847 | 0.5462 | 31 |
| 1.0% | 4.5919 | 0.9474 | 0.8387 | 47 |

The fraction whose sharpest cell beats the real (0.17417473%) cell is
**0.4275**. This independently reproduces the manuscript's rounded Table IV
values.

## Corrected conditional replay

A local audit run of the corrected conditional implementation (same formulas
as `conditional_replay.py`) gives, at N=2,000,000 and seed 731994:

- analytic finite-window Koide frequency: about **5.17e-6**
- corrected conditional five-primary-slot joint: about **3.39e-17**
- true six-slot joint including the c-b product strip: about **2.31e-17**
- sixth-slot cost: **joint6/joint5 ~ 0.681**
- corresponding full-chain frequencies: about **1.75e-22** (five-primary)
  and **1.19e-22** (true six-slot)

A five-seed N=1,000,000 diagnostic gave joint5 values spanning roughly
3.06e-17 to 3.84e-17 and joint6 values 2.08e-17 to 2.61e-17. This visible
between-seed spread is a reason to retain multi-seed convergence reporting
rather than quoting a single run's internal Monte Carlo error.

## Corrected per-slot diagnostics

A weighted N=2,000,000 diagnostic under the actual finite-window conditional
lepton law gives approximately:

- (F^2/m_s): **1.13e-3**
- (mu_* m_d/m_s^2): **6.90e-4**
- (P m_d/m_u^2): **1.02e-3**
- (m_c/(3\alpha_K\mu_*)): **4.65e-4**
- (m_b/(\mu_*/(2\alpha_K))): **1.65e-4**
- (m_cm_b/G^2): **2.56e-4**

These remain close to the manuscript's Table V diagnostics. The substantive
change is therefore not that the individual slots vanish; it is the precise
definition of the conditioned joint event.

## Variance-reduced cross-check

A separate exact-support mixture importance sampler targets the same corrected
conditional law while oversampling the analytically known log-`m_d` interval
in which the two strange windows can overlap. A 5-seed diagnostic with
N=500,000 per seed gave:

- (P_K = (5.16993 \pm 0.00216)\times10^{-6}) across seeds;
- five-primary joint (=(3.40074 \pm 0.01307)\times10^{-17});
- true six-slot joint (=(2.31454 \pm 0.00889)\times10^{-17});
- sixth-slot cost (=0.680597);
- full-chain five-primary (≈1.758\times10^{-22});
- full-chain true-six (≈1.196\times10^{-22}).

The importance sampler restores the original log-uniform (m_d) law with an
explicit likelihood ratio and retains a global proposal component, so it has
full support. Its much tighter seed-to-seed stability is the preferred current
cross-check.

## Current interpretation

The historical ~2.83e-17 number is the same order of magnitude, but its exact
value is sensitive to the conditioning measure, residual orientation, and
whether the sixth product constraint is actually imposed. The corrected
baseline calculation should replace it before any precision claim is frozen.

The central qualitative conclusion survives: isolated sub-percent cells are
cheap inside the declared 3740-cell grammar, while the frozen multi-constraint
cascade is much more selective under the chosen procedural generator.
