# candidate-space-enumeration

Complete reproduction archive for the numerical content of
*"Involutive mass geometry: from the Foot cone to split heavy-sector closed
forms"* (working draft `main_revised_7.tex`, archived under
`manuscript/main_revised_7.tex`), A. M. Brilliant, 2026. Rebuilt 2026-10-11
from the two predecessor archives (see Lineage) against the revised
manuscript, after the heavy-sector relation was restated from the mirror
form to the split closed forms

    m_c = 3 aK mu*        m_b = (2 + aK) mu*,      aK = sqrt(3/2) - 1

whose product closes exactly on the grammar's defining constant,
`3 aK (2 + aK) = 3 (2 aK + aK^2) = 3/2` (the quadratic `aK^2 + 2 aK = 1/2`),
i.e. `m_c m_b = (3/2) mu*^2 = (1 + aK)^2 mu*^2`.

**Naming note.** The old mirror-era symbol `G` (= sqrt(3/2) mu*, the
"geometric-mean center") is retired throughout this repository. The anchor is
now carried under its construction, `sqrt32mu = sqrt(3/2) mu*`, and the heavy
relation is stated as the two split forms with the product anchor
`(3/2) mu*^2`. Machine label `sqrt32mu^2` in the census outputs denotes the
squared anchor, exactly `(3/2) mu*^2`.

## Layout

- `manuscript/main_revised_7.tex` — pinned manuscript snapshot (the document
  this archive verifies).
- `inputs/pdg2026_inputs.csv` — every external input with convention and
  provenance, including the PDG 2024 historical rows for the edition table.
  Light quarks enter at the PDG 2 GeV reference and are run to
  mu* = m_e + m_mu + m_tau = 1883.099374 MeV under the four-loop factor
  R_m = 1.01750.
- `census/census_mustar.py` (+ `census_mustar.csv`, `census_summary.txt`) —
  the original 1304-comparison multiplicative census (classes P/T/A), the
  content of the archived pre-rebuild commit `e67211b` (labels renamed per
  the naming note; numbers bit-identical).
- `pipeline/` — the declared six-class grammar (3740 comparisons) and its
  Monte Carlo null treatments. Single source of truth `census_grammar.py`;
  drivers: `extended_census.py`, `lite_mc.py`, `sensitivity_mc.py`,
  `grammar_robustness.py`, `test_b_exact.py`, `test_b_convergence.py`,
  `coefficient_grammar.py` (new), `importance_sampling.py` (new). All
  outputs committed and pinned in `provenance/checksums.txt`.
- `sweep/` — the continuum sweep: 9788 geometric-mean configurations,
  300-point logarithmic grid, 1-400 GeV (`enumerate_crossovers.py`,
  `sweep_table_pdg2026.csv`, `tier1_census.csv`, `crossings.csv`,
  `nearest_misses.csv`).
- `checks/manuscript_checks.py` (+ `checks_summary.txt`) — every standalone
  manuscript number not otherwise reproduced by the artifacts above, as
  PASS/FAIL gates.
- `provenance/conventions_provenance.txt` — conventions and provenance of
  the sweep run (script-generated); `provenance/checksums.txt` — SHA-256
  of every file.

## Validators (self-contained)

Requirements: python3 (+numpy for `test_b_*`, `importance_sampling.py`,
`sensitivity_mc.py`; stdlib otherwise). No network access, nothing fetched
at runtime.

1. `sha256sum -c provenance/checksums.txt` — every file reports OK.
2. `cd census && python3 census_mustar.py` — prints `1304 ok`, sweep
   `0/1/1/2`, unique sub-0.5% hit `+0.1742% c*b vs sqrt32mu^2`, images at
   `+0.7872%` / `+1.1826%`, nearest non-declared `-1.2031%`. Regenerates
   outputs byte-identically.
3. `cd pipeline && python3 extended_census.py` — `3740` comparisons
   (P=324, T=252, A=728, S=426, D=426, K=1584), tolerance ladder
   `0/1/3/6`, disclosed Descartes accidental `-0.2901%` (counted twice via
   the tau mass/anchor duplication), classes S and D empty at every rung.
4. `python3 coefficient_grammar.py` — raw grid 528 points (448 multiplicative
   + 80 additive); the only points within 0.3% of `m_c(m_c)` / `m_b(m_b)` are
   the two declared forms and their quadratic rewritings (`6/2 = 3`;
   `(1/2) aK^-1 = 2 + aK`); at the 0.1% rung only the bottom form remains;
   log-uniform coverage of [0.5, 20] GeV: 0.246 (0.3% windows) / 0.087
   (0.1% windows), the manuscript's "~0.25 / ~0.09".
5. `python3 lite_mc.py` — whole-grammar replay, 10^4 seeded universes:
   means 0.44 / 1.13 / 2.28 / 4.59 cells at 0.1 / 0.25 / 0.5 / 1%; 43% of
   universes beat the real sharpest cell.
6. `python3 sensitivity_mc.py all` — Koide front-end frequency 5.4e-6 under
   lepton nulls (2e8 draws); per-slot cascade frequencies (diagnostics);
   whole-grammar variants (Koide-locus conditioning: 2.24 vs 2.28 baseline).
7. `python3 grammar_robustness.py` — 9/18/27-anchor nested grammars,
   accident means scale linearly 2.28/3.60/4.80.
8. `python3 test_b_exact.py` — Rao-Blackwellized joint, five-cut event
   (matching-or-better on the five mass slots): primary `(2.93 +- 0.05)e-17`;
   product-of-marginals 4.9e-17 (factor 1.7); `P(T <= 1.183%) = 5.57e-15`
   (T-statistic defined on the five mass slots, unaffected by the sixth cut);
   shared-m-s pair validation 3.43e-7 brute vs 3.39e-7 analytic; T-tail curve
   in `test_b_exact_curve.csv`.

   **Sixth cut (post-review).** The product slot
   `|m_c m_b/(1.5 mu*^2) - 1| <= 0.174%` is NOT implied by the charm and
   bottom cuts: same-sign corners of the individual windows give
   `|r_cb| ~ 0.34%`. It is enforced as the exact coupled region
   `|r_c + r_b + r_c r_b| <= t6` in the (c, b) draw plane -- rectangle times
   diagonal strip, closed-form clipped area (`_G`), replacing the independent
   `pc*pb` factor. Six-cut joint6 is reported next to joint5 in all three
   estimators (expected ratio ~0.68). `test_b_convergence.py` — five frozen
   seeds for both events (joint5 `(2.77-2.89)e-17`, combined
   `(2.83 +- 0.02)e-17`; joint6 combined ~1.9e-17), ESS and weight-dominance
   stats for each. `importance_sampling.py` — independent log-normal proposal
   around `d* = aK^4 mu*`, 4 x 5e6 runs at two widths, both events;
   statistically identical to the conditional-integration estimates.
9. `cd sweep && python3 enumerate_crossovers.py` — gates: input-table
   integrity at 91.1876 GeV, formal `c*b = (3/2) mu*^2` crossing at
   mu = 2.068807 GeV, `9788 configurations, 360 crossings`, none at a
   distinguished scale; `92` unphysical continuations. Chance baseline for
   the null: uniform placement of 360 crossings over ln[1, 400] collects
   ~3.0 in the five 5e-3 windows (`checks/manuscript_checks.py` gate).
10. `cd checks && python3 manuscript_checks.py` — ALL GATES PASS covers:
    predictions 1269.65 +- 0.06 / 4189.4 +- 0.2 MeV; residuals +0.256% /
    -0.082% (2026) and +0.264% / -0.153% (2024); grading units 0.354% /
    0.143%; ratio seed 3.2997 vs measured 3.2886 (-0.34%); aK-completion
    shifts 46 keV / 15 keV; light sector F^2 = 95.1158 MeV, images -0.622%
    (archived unrounded), -0.78%, terminal model 2.21585 MeV at -0.81%,
    ratio aK^-2 = 19.798 (-1.5% vs FLAG 20.1); s_eff = 5 (displacement
    2.20e-6); up-type Q 0.888 / 0.914; top probe aK^-3 mu* = 165.88 GeV
    (-3.9% vs direct, +2.1% vs m_t(m_t), schemes 6.2% apart); common-scale
    misses 66% / 37% at M_Z, 45% at 10 GeV, 27% at m_b; transported-product
    residual invariant across the grid at +0.174%.

Scripts in `pipeline/` import each other; keep them in one directory.

## Claim-to-file map (manuscript -> reproduction)

| Manuscript claim | Artifact |
|---|---|
| 3740 six-class census; ladder 0/1/3/6; sub-1% six cells itemized | `pipeline/extended_census.{py,csv}`, `extended_summary.txt` |
| c-b/(3/2)mu*^2 cell at +0.174% (records the grammar; no second empirical success) | `extended_census.csv` class P row |
| Descartes accidental -0.290% counted twice; +0.603% extra at 1%; images +0.622% / +0.787%; terminal +1.183% outside; nearest non-declared -1.203% | `extended_summary.txt` |
| Coefficient grammar: 528 points, no local competitors, coverage 0.25/0.09 | `pipeline/coefficient_grammar.{py,txt}` |
| Whole-grammar replay table | `pipeline/lite_mc.{py,summary}` |
| Six-slot per-slot rates 1e-3..1e-4; Koide front-end 5.4e-6 | `pipeline/sensitivity_*` |
| P_proc primary 2.93e-17; replicated 2.83e-17; marginals 4.9e-17; brute/analytic 3.43/3.39e-7; P(T<=T_obs) 5.6e-15 | `pipeline/test_b_exact_*` |
| Convergence record (5 seeds, ESS, dominance) | `pipeline/test_b_convergence.{py,txt}` |
| Importance-sampling cross-check (two widths) | `pipeline/importance_sampling.{py,txt}` |
| Menu growth 9/18/27 anchors linear | `pipeline/grammar_robustness.py` (stdout) |
| Sweep: 9788 configs, 360 crossings, 92 unphysical, mu(c*b) = 2.069 GeV, no distinguished scale, AM-GM -0.084% nearest anchor-internal, chance baseline ~3 | `sweep/`, `checks/` |
| Heavy forms: predictions, 2024/2026 residuals, grading units, ratio, keV shifts | `checks/manuscript_checks.py` |
| Light sector residuals and closed-form values | `checks/manuscript_checks.py` |
| s_eff = 5; up-type 0.888/0.914; top probe 165.88 GeV; transported-residual invariance | `checks/manuscript_checks.py` |
| Zenodo-only: running pipeline details behind 7.9889/7.9900, Rivero crossing ~280 TeV, pole-mass substitution misses 30%/14%, prescription grid 2.3-17% / 30-73% / -21.9% | supplementary calculations document |

## Lineage

- `candidate-space-enumeration` (commit `e67211b`, pre-rebuild history):
  the original 1304-comparison multiplicative census cited by the
  manuscript's bibliography entry; archived on GitHub as
  `AndBrilliant/candidate-space-enumeration-archived` after this rebuild.
- `foot-mirror-repro` (private): the September/October complete repro; its
  contents were merged into this repository, which was then recreated
  against `main_revised_7.tex`. Bundles of both predecessors are kept at
  `ganymede:/ARCHIVE/GIT_REPOS/` and locally.
- The mirror-era `G`-center presentation was replaced by the split closed
  forms after the 2026-10 review; no numerical result changed (the product
  identity is exact algebra).
