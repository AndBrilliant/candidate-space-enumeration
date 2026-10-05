# marcus_repro

Independent referee-style replication for **From the Foot Cone to Multiplicative Mass Mirrors**.

This subtree is deliberately isolated from `pipeline/`: it does not import historical analysis code.
The manuscript is treated as the specification.

## Goals

1. Recompute deterministic headline quantities from declared inputs.
2. Rebuild the six-class 3,740-cell census independently.
3. Reproduce the whole-grammar procedural Monte Carlo.
4. Recompute the rare-event conditional replay under the **actual conditional
   baseline law**: three independent log-uniform lepton masses conditioned on
   the finite Koide window.
5. Report both a five-primary-slot joint and a true six-slot joint including
   the charm-bottom product strip.

## Important statistical distinction

The quantity estimated here is a procedural frequency under a synthetic
generator. It is **not** a Standard-Model p-value, a global historical
look-elsewhere correction, or a probability that the observed pattern
"arose by chance."

## Files

- `inputs.py` — single source of declared numerical inputs and tolerances.
- `deterministic_checks.py` — headline algebra/residual checks.
- `grammar.py` — independent six-class census implementation.
- `whole_grammar_mc.py` — 10,000-universe whole-grammar replay.
- `conditional_replay.py` — corrected conditional rare-event estimator.

Run from this directory with Python 3.10+ and NumPy:

```bash
python deterministic_checks.py
python grammar.py
python whole_grammar_mc.py
python conditional_replay.py --n 1000000
```

For publication-grade rare-event precision increase `--n` (the estimator is
Rao-Blackwellized over downstream quark slots, so one million to five million
weighted lepton/`m_d` draws is practical).

## Why this implementation differs from the historical pipeline

The previous exact-surface implementation sampled a particular
`(M, theta)` measure on the exact Koide surface. That is a useful sensitivity
ensemble, but it is not automatically equal to the conditional measure induced
by independent log-uniform masses.

Here, for each independently drawn pair `(m1,m2)`, the finite Koide-window
constraint is solved analytically for all allowed `m3` intervals. Their
log-widths are the exact conditional probabilities under the original
log-uniform generator. One `m3` is sampled from those intervals and the draw
is weighted by the interval probability. The resulting self-normalized
estimator therefore targets the baseline lepton law conditioned on the stated
finite Koide window, without replacing it by an ad hoc surface measure.

The six-slot calculation additionally integrates the joint `(c,b)` log-space
rectangle intersected by the observed `cb/G^2` strip, rather than calling
the product slot an automatic consequence of the two one-mass slots.
