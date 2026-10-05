# Marcus statistical audit notes

## Status labels

- **M1 — conditional-law mismatch:** the historical exact-surface sampler is not
  automatically the conditional law induced by independent log-uniform lepton
  masses. The Marcus baseline targets the finite Koide-window conditional law
  directly.
- **M2 — five-slot vs six-slot naming:** the historical rare-event product used
  five primary intervals while calling the result a six-slot event. The Marcus
  code reports both five-primary and true six-slot joints.
- **M3 — residual orientation:** light-slot interval inequalities are solved from
  the declared cell residual `V/T - 1`, not from reciprocal approximations.
- **M4 — threshold provenance:** the manuscript's 0.6216% strange-anchor
  threshold is explicit in `inputs.py`; it is not silently regenerated from
  rounded table entries.
- **M5 — archived-code drift:** manuscript-facing releases must cite a commit
  containing the exact six-class census and every reported MC driver/output.
- **M6 — trials rhetoric:** nested anchor-menu growth is a robustness experiment,
  not a proof of the historical trials factor. Likewise grid points in a
  continuous scan are not automatically independent trials.

## What survives

Independent deterministic reconstruction reproduces the headline heavy-sector
numbers to displayed precision. The six-class grammar reconstructs to 3740
cells with class counts P=324, T=252, A=728, S=426, D=426, K=1584 and the
reported tolerance ladder [0,1,3,6].

The whole-grammar procedural replay remains the correct control for the claim
that isolated per-mille coincidences are cheap inside the declared grammar.

## Publication recommendation

Do not quote a precision rare-event number until the corrected conditional
estimator has a multi-seed convergence record and an independent implementation
has reproduced it. Retain `P_proc` language and avoid translating it to a
Gaussian discovery significance.
