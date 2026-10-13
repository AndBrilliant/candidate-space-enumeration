# koide_production_v24 — corrected finite-band Koide conditional law

Independent, portable replacement for the archived signed-amplitude surface
generator (see the "Correction to the conditional Koide generator" paragraph of
the manuscript, Sec. procedural replay).

The archived pipeline (`pipeline/sensitivity_mc.py`, `pipeline/test_b_exact.py`)
generated triples as m_i = M [1 + sqrt(2) cos(theta + 2 pi i/3)]^2 with only an
in-span mass check. Negative bracketed amplitudes were squared into positive
masses, so a fraction of emitted triples do not satisfy the physical Koide
condition. This package instead samples (ln M, theta, Q) over the finite band
|Q - 2/3| < 2.2e-6, retains only positive-amplitude in-span triples, and
weights by the change-of-variables Jacobian, giving the conditional
distribution induced on the band by independently log-uniform lepton masses.

Run (numpy + scipy):
    python3 run_all.py --n 9000 --seeds 4 --nodes 96 --out corrected_results.json

The committed `corrected_results.json` is the four-seed, 96-node regression
whose means are quoted in the manuscript: baseline five/six-cut 3.39e-17 /
2.31e-17, doubled-span 3.33e-18 / 2.27e-18, uniform-in-mass 1.38e-27 /
9.4e-28, halved-span exactly zero by support incompatibility (paired
strange-sector windows force m_d <= 13.88 MeV; the halved-span draw requires
m_d >= 36.79 MeV). Slot-1 threshold is the exact observed F^2 residual
0.6243%.
