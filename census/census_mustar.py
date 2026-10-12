#!/usr/bin/env python3
"""Candidate-space census at the self-organized scale (chronicle version).

Reproduces the enumeration reported in "From the Foot Cone to Multiplicative
Mass Mirrors" (chronicle branch): 1304 geometric-mean comparisons over the
nine charged-fermion masses and nine frozen anchors, with light quarks run
from the PDG 2 GeV reference to mu* = m_e + m_mu + m_tau under the four-loop
factor R_m = 1.01750 documented in the manuscript appendix.

Classes:
  P: unordered mass pairs m_a*m_b vs squared anchors A^2        (324)
  T: center-squared triples m_x^2 vs m_a*m_b, x not in {a,b}    (252)
  A: anchor-mass products A*m_a vs m_x^2                        (728)
Algebraic identities (value identical to target by construction) are excluded.

Expected result (PDG 2026 central values):
  tolerance sweep 0.1% / 0.25% / 0.5% / 1.0%  ->  0 / 1 / 1 / 2 hits
  unique sub-0.5% hit: c*b vs sqrt32mu^2 at +0.174%  (sqrt32mu^2 = (3/2)mu*^2)
  declared light-mirror images: mu*d vs s^2 at +0.787%; P*d vs u^2 at +1.183%
  nearest non-declared comparison: sqrt32mu*d vs mu_lepton^2 at -1.203%
"""
import csv, itertools, sys

# ---- PDG 2026 inputs (MeV) -------------------------------------------------
M_E, M_MU, M_TAU = 0.51099895069, 105.6583755, 1776.930
U_2GEV, D_2GEV, S_2GEV = 2.16, 4.70, 92.9
R_M = 1.01750  # four-loop 2 GeV -> mu* mass factor, n_f=4 (manuscript App. A)
M_C, M_B, M_T = 1272.9, 4186.0, 172600.0  # c(m_c), b(m_b), t direct

MUSTAR = M_E + M_MU + M_TAU          # 1883.099374
SQRT32MU = (1.5) ** 0.5 * MUSTAR    # sqrt(3/2) mu* = 2306.316...
P = 2 * M_E                          # 1.0219979

masses = {
    "e": M_E, "mu": M_MU, "tau": M_TAU,
    "u": U_2GEV * R_M, "d": D_2GEV * R_M, "s": S_2GEV * R_M,
    "c": M_C, "b": M_B, "t": M_T,
}
anchors = {
    "mustar": MUSTAR, "P_2me": P, "sqrt32mu": SQRT32MU, "tau": M_TAU,
    "mustar_plus_tau": MUSTAR + M_TAU, "two_tau": 2 * M_TAU,
    "mustar_over_2": MUSTAR / 2, "two_mustar": 2 * MUSTAR,
    "one_point_five_mustar": 1.5 * MUSTAR,
}

rows = []
def add(cls, val_name, val, tgt_name, tgt, tag=""):
    if abs(val - tgt) < 1e-12 * max(abs(tgt), 1e-300):
        return  # algebraic identity, excluded as trivial
    rows.append({"class": cls, "value": val_name, "target": tgt_name,
                 "value_num": f"{val:.6g}", "target_num": f"{tgt:.6g}",
                 "residual_pct": f"{(val / tgt - 1) * 100:+.6f}", "note": tag})

names = list(masses)
for a, b in itertools.combinations(names, 2):                      # Class P
    for A, Av in anchors.items():
        add("P", f"{a}*{b}", masses[a] * masses[b], f"{A}^2", Av * Av)
for a, b in itertools.combinations(names, 2):                      # Class T
    for x in names:
        if x in (a, b):
            continue
        add("T", f"{x}^2", masses[x] ** 2, f"{a}*{b}", masses[a] * masses[b])
for A, Av in anchors.items():                                      # Class A
    for a in names:
        for x in names:
            add("A", f"{A}*{a}", Av * masses[a], f"{x}^2", masses[x] ** 2)

DECLARED = {("c*b", "sqrt32mu^2"), ("mustar*d", "s^2"), ("P_2me*d", "u^2")}
for r in rows:
    if (r["value"], r["target"]) in DECLARED:
        r["note"] = "declared relation (measured-mass image)"

assert len(rows) == 1304, f"comparison count drifted: {len(rows)}"

with open("census_mustar.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys(), lineterminator="\n")
    w.writeheader(); w.writerows(rows)

res = sorted(((float(r["residual_pct"]), r) for r in rows), key=lambda t: abs(t[0]))
sweep = {t: sum(1 for x, _ in res if abs(x) < t) for t in (0.1, 0.25, 0.5, 1.0)}
with open("census_summary.txt", "w") as f:
    f.write("mu* census, PDG 2026 central values\n")
    f.write(f"comparisons: {len(rows)} (P=324, T=252, A=728)\n")
    f.write("tolerance sweep 0.1/0.25/0.5/1.0%: "
            + "/".join(str(sweep[t]) for t in (0.1, 0.25, 0.5, 1.0)) + "\n\n")
    f.write("five smallest |residual|:\n")
    for x, r in res[:5]:
        f.write(f"  {x:+.4f}%  {r['value']} vs {r['target']}  {r['note']}\n")
    nd = next((x, r) for x, r in res if not r["note"])
    f.write(f"\nnearest non-declared comparison: {nd[0]:+.4f}%  "
            f"{nd[1]['value']} vs {nd[1]['target']}\n")

print(f"1304 ok; sweep: {sweep}")
for x, r in res[:5]:
    print(f"  {x:+.4f}%  {r['value']} vs {r['target']}  {r['note']}")
