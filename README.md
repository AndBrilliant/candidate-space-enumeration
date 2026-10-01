# candidate-space-enumeration

Exhaustive enumeration and null-model stress tests of the discrete candidate
space declared in the manuscript "From the Foot Cone to Multiplicative Mass
Mirrors" (working draft, 2026).

## Original multiplicative census (frozen, cited commit e67211b)

- 9 PDG-2026 central charged-fermion masses
- frozen anchor menu {mu_star, 2m_e, G, m_tau, mu+tau, 2tau, mu_star/2,
  2mu_star, 1.5mu_star}; mu_star = 1883.0994 MeV, G = sqrt(3/2) mu_star
- three comparison classes: pair products m_a m_b vs squared anchors A^2;
  center-squared triples m_x^2 = m_a m_b (x not in {a,b});
  anchor-mass triples m_x^2 = A m_a over the full anchor x mass grid
- match: |value - target|/target < 5e-3; algebraic identities excluded
- result: 1304 comparisons, one hit (class P, c*b vs G^2, +0.1742%)

Run: python3 candidate_space_enumeration.py  (Python 3, stdlib only)
The script regenerates candidate_space_enumeration.csv in place; the
regenerated file must be byte-identical to the committed one.

## Extended six-class grammar (pipeline/)

The declared space completed with additive, subtractive, and Descartes
classes, sharing census_grammar.py as the single source of truth:

- classes S/D: sums m_a+m_b and differences |m_a-m_b| vs the 18 targets
  (9 masses + 9 anchors)
- class K: Descartes-Soddy completions k4 = sum(sqrt m) +- 2 sqrt(sum sqrt
  m_i m_j) of all 84 triples, both signs, k4^2 vs the 18 targets
- S/D/K exclude any comparison whose target construction contains an operand
  (structural disjointness; removes hierarchy trivialities by construction)
- result: 3740 comparisons; tolerance ladder 0/1/3/6 at 0.1/0.25/0.5/1%;
  the declared c*b/G^2 candidate (+0.1742%) is the sharpest cell; one
  Descartes-class accidental disclosed (k4(e,mu,s,+)^2 vs m_tau, -0.2901%,
  counted twice via the tau mass/anchor duplication); classes S and D empty
  at every tolerance scanned

Run: python3 pipeline/extended_census.py

## Monte Carlo null treatments (pipeline/)

All drivers are seeded and frozen; outputs are committed and pinned by
checksums.txt.

- lite_mc.py — whole-grammar null: 10^4 mock universes (3 leptons log-U over
  the observed lepton span, 6 quarks log-U over the observed quark span,
  anchors rebuilt per universe from the mock lepton triple). Result: sub-0.5%
  coincidences are generic (mean 2.28/universe; 43% of universes beat the
  real census's sharpest cell). Singles are cheap.
- sensitivity_mc.py — the conditional cascade null plus null-model
  sensitivity variants (APPB A1-A4 discipline):
  * part B: Koide front-end frequency P(|Q-2/3|<2.2e-6) = 5.4e-6 baseline,
    stable across span x2/x0.5 and uniform nulls
  * part D: per-slot cascade frequencies (validated diagnostics only;
    the slots share inputs, so their product is an uncontrolled
    approximation superseded by test_b_exact.py), including a
    rejection-sampled reference ensemble
  * part A: whole-grammar null variants, including lepton triples
    constrained to lie exactly on the Koide locus (accident mean 2.24 vs
    2.28 baseline: conditioning the leptons does not manufacture the
    quark-sector matches)
- grammar_robustness.py — nested grammars with 9/18/27 anchors
  (3740/5717/7694 cells): accident means scale linearly (2.28/3.60/4.80).
  Discrete menu growth costs linearly; a free continuous parameter would
  cost orders of magnitude more.
- test_b_exact.py — conditional analytic integration with Monte Carlo
  averaging (a Rao-Blackwellized estimator) for the joint 5-slot cascade
  frequency: conditional on the anchors and the m_d draw each slot is an
  interval on a single mass, so the per-universe joint is a product of
  analytic interval probabilities. Baseline joint (2.93+-0.05)e-17
  (5e6 universes, seed frozen); interval machinery validated against a
  3e8-draw brute-force count of the shared-m_s slot pair (3.43e-7 vs
  analytic 3.39e-7). Also emits the joint-discrepancy curve
  P(T <= t) for T = max_i |r_i|.
- test_b_convergence.py — estimator convergence documentation: five
  independent frozen seeds x 5e6 universes give per-seed estimates
  (2.77-2.89)e-17, combined (2.83+-0.02)e-17; effective sample size
  ~3.7e3 per run (~75% of contributing universes); no dominance (largest
  single weight 0.03% of total, ten largest 0.3%, max nonzero ~1.5x
  median nonzero).

Run from the pipeline directory: python3 lite_mc.py;
python3 sensitivity_mc.py all; python3 grammar_robustness.py;
python3 test_b_exact.py; python3 test_b_convergence.py
(Python 3 + numpy/scipy; lite_mc is stdlib only). The scripts import each
other and must stay in one directory.

## Integrity

checksums.txt pins sha256 of every file in this repository. Verify:
sha256sum -c checksums.txt

The complete comparison tables are also deposited with the preregistered
bundle on Zenodo (see manuscript citation).
