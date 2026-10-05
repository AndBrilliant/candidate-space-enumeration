# Marcus reference audit results

These are **independent replication diagnostics**, not evidential significances.

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

The deterministic mass algebra independently reproduces the manuscript's
headline heavy residuals:
- charm: +0.255914%
- bottom: -0.081531%
- c*b/G^2: +0.174175%

Interpretation: the historical ~2.83e-17 number is the same order of magnitude,
but the exact value is sensitive to the conditioning measure and event
orientation. The corrected baseline calculation should replace it before any
precision claim is frozen.
