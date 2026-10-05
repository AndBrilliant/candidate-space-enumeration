# Multiplicative Mass Mirrors: candidate-space enumeration and Monte Carlo nulls

Manuscript-facing statistical analysis for *From the Foot Cone to Multiplicative
Mass Mirrors*. This branch (`marcus-reform`) implements the repository reform
required by the independent statistical audit (the "Marcus" review).

## Layout

- `inputs/pdg2026.yaml` — single frozen source of truth for all masses, the
  null generator, tolerances, and seeds. **No threshold is hand-typed anywhere
  downstream** (audit §IX).
- `analysis/` — manuscript-facing lineage:
  - `deterministic_relations.py` — mu*, G, heavy-sector predictions, Q_l, Foot angle
  - `census.py` — six-class 3740-cell grammar + tolerance ladder
  - `whole_grammar_mc.py` — 10^4-universe replay (Table IV)
  - `conditional_joint.py` — corrected finite-window conditional Koide law,
    exact P_K, five-primary (P5) and genuine six-slot (P6) joint frequencies
  - `importance_joint.py` — variance-reduced cross-check (audit §VII)
  - `continuum_sweep.py` — continuous-scale control (audit §XIII-1)
  - `running.py` — light-quark running factor (documentation; see §Caveats)
- `tests/` — automated regression tests (audit §XII)
- `outputs/` — regenerated tables (census.csv, table_iv.csv,
  conditional_joint.json, importance_joint.json, continuum_sweep.csv)
- `marcus_repro/` — independent checker lineage (re-derived from the audit's
  stated numbers; kept conceptually separate per audit §XV)
- `pipeline/`, `candidate_space_enumeration.py` — **historical** code, preserved
  unchanged (commits e67211b, 4fc061c, 5414e31, fec9592). Not manuscript-facing.

## Reproduce

```bash
pip install -r requirements.txt
./reproduce_all.sh
```

## Conventions (binding)

- Residual statistic is `V/T - 1` everywhere (audit §VIII).
- The baseline lepton generator is three independent log-uniform masses
  conditioned on the finite Koide window `|Q - 2/3| <= 2.2e-6`; the exact
  Q=2/3 surface generator is a sensitivity variant, not the baseline
  (audit §IV).
- The "six-slot event" is reported only for the calculation that imposes all
  six conditions; the five-primary event P5 and the genuine six-slot event P6
  are reported separately (audit §V).
- Procedural frequencies are never described as Standard Model p-values,
  global look-elsewhere significances, or HEP discovery significances
  (audit §VI, §XIII).

## Caveats / known provenance notes

- `analysis/running.py` gives a 2-loop estimate R_m ~ 1.003 that does NOT
  reproduce the frozen 1.01750 used by the manuscript; the frozen value is
  authoritative (from the archived CRunDec pipeline) and `running.py` is a
  documentation placeholder pending restoration of the full CRunDec matching
  chain (audit §XI task 8).
- The observed lepton Q deviates from 2/3 by 2.2033e-6, marginally OUTSIDE the
  2.2e-6 conditioning window. The window center/half-width should be reviewed
  against the frozen inputs before final publication (audit §IX provenance).

## Threshold-robustness statement (audit §IX resolution)

All slot thresholds are derived at runtime from `inputs/pdg2026.yaml`
(rounded to 4 decimal places in percent). With these derived thresholds the
five-primary joint frequency is

    P5 = 2.96e-17  (direct, n=500000, 5 seeds)
       = 2.88e-17  (importance-sampling cross-check, agreement ~4%)

The historically quoted value 3.40e-17 corresponds to thresholds rounded from
the manuscript table (0.6216% vs the derived 0.6204% for slot 1). The estimate
is therefore stable at O(10%) under the threshold-rounding choice; we report
this as the dominant methodological uncertainty on P5/P6. No conclusion in the
manuscript depends on it: the joint-frequency claim lives at "P ~ 1e-17",
six orders of magnitude below the look-elsewhere scale of the grammar.

## Koide-window boundary note

The finite conditioning window |Q - 2/3| <= 2.2e-6 is frozen as part of the
null specification; it is NOT re-centered on the observed lepton triple.
The observed deviation 2.2033e-6 sits marginally outside this half-width.
We report P_K = 5.17e-6 with the frozen window as primary and note the
boundary sensitivity explicitly: P_K remains O(1e-6) for any window
center/half-width in this neighborhood (the Koide surface is a measure-zero
line in the log-uniform cube; the exact-window probability scales linearly
with half-width, so a 0.2% shift in the boundary moves P_K by < 1%).
Re-centering the window post hoc would constitute the fitting procedure this
analysis exists to exclude, and is not done.
