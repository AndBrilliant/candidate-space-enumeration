#!/usr/bin/env python3
"""
Candidate-space enumeration accompanying
"From the Foot Cone to Multiplicative Mass Mirrors" (working draft, 2026).

Exhaustive evaluation of the discrete candidate space declared in the
manuscript's "Candidate-space enumeration" subsection.

Frozen inputs (PDG 2026 central values) and frozen anchor menu; no
continuously tunable parameter. Tolerance: |dev|/target < 5e-3 (five per
mille). Output: candidate_space_enumeration.csv with one row per comparison
and a summary printed to stdout.

Comparison classes (all products, matching the manuscript's convention):
  P  unordered mass pairs (m_a, m_b), a < b, tested against A^2 for each
     anchor A  -> m_a*m_b vs A^2
  T  center-squared triples m_x^2 = m_a*m_b, x not in {a,b}
  A  anchor-mass triples m_x^2 = A*m_a, for each anchor A and masses x, a
Algebraic tautologies (value identical to target, e.g. tau^2 = m_tau*tau)
are excluded as trivial.
"""

import csv
import itertools
import math

# --- Frozen mass inputs, PDG 2026 central values (MeV) -----------------------
MASSES = {
    "e": 0.51099895069,
    "mu": 105.6583755,
    "tau": 1776.93,
    "u": 2.16,
    "d": 4.70,
    "s": 93.5,
    "c": 1272.9,
    "b": 4186.0,
    "t": 172690.0,
}

# --- Frozen anchor menu (MeV) ------------------------------------------------
MU_STAR = MASSES["e"] + MASSES["mu"] + MASSES["tau"]          # 1883.0994 MeV
ANCHORS = {
    "mu_star": MU_STAR,
    "2m_e": 2.0 * MASSES["e"],
    "G": math.sqrt(1.5) * MU_STAR,                            # 2306.316 MeV
    "m_tau": MASSES["tau"],
    "mu+tau": MASSES["mu"] + MASSES["tau"],
    "2tau": 2.0 * MASSES["tau"],
    "mu_star/2": MU_STAR / 2.0,
    "2mu_star": 2.0 * MU_STAR,
    "1.5mu_star": 1.5 * MU_STAR,
}

TOL = 5e-3  # five per mille, relative

rows = []
names = list(MASSES)
pairs = list(itertools.combinations(names, 2))

def add(cls, label, anchor, value, target):
    if value == target:
        return  # algebraic tautology, excluded as trivial
    dev = (value - target) / target
    rows.append((cls, label, anchor, value, target, dev, abs(dev) < TOL))

# Class P: m_a*m_b vs anchor^2
for (a, b) in pairs:
    prod = MASSES[a] * MASSES[b]
    for aname, aval in ANCHORS.items():
        add("P", f"{a}*{b}", aname + "^2", prod, aval ** 2)

# Class T: m_x^2 = m_a*m_b triples, x distinct from a,b
for (a, b) in pairs:
    gm2 = MASSES[a] * MASSES[b]
    for x in names:
        if x in (a, b):
            continue
        add("T", f"{x}^2={a}*{b}", x, gm2, MASSES[x] ** 2)

# Class A: m_x^2 = A*m_a triples
for aname, aval in ANCHORS.items():
    for a in names:
        prod = aval * MASSES[a]
        for x in names:
            add("A", f"{x}^2={aname}*{a}", aname, prod, MASSES[x] ** 2)

with open("candidate_space_enumeration.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["class", "comparison", "anchor_or_target", "value_MeV2_or_MeV",
                "target", "relative_deviation", "within_5_per_mille"])
    for r in rows:
        w.writerow([r[0], r[1], r[2], f"{r[3]:.6g}", f"{r[4]:.6g}", f"{r[5]:+.6e}", r[6]])

hits = [r for r in rows if r[6]]
print(f"total comparisons: {len(rows)}")
print(f"hits within {TOL:.0e} tolerance: {len(hits)}")
for r in hits:
    print(f"  class {r[0]}: {r[1]} vs {r[2]}  rel. dev = {r[5]:+.4%}")
