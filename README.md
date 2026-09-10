# candidate-space-enumeration

Exhaustive enumeration of the discrete candidate space declared in the
manuscript "From the Foot Cone to Multiplicative Mass Mirrors"
(working draft, 2026).

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

The complete comparison table is also deposited with the preregistered
bundle on Zenodo (see manuscript citation).
