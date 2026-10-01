#!/usr/bin/env python3
"""Extended candidate-space census: additive, subtractive, and Descartes classes.

Extends the mu* census (census_mustar.py, 1304 multiplicative comparisons)
with the remaining dimensionally homogeneous operations on the same frozen
universe: nine charged-fermion masses and nine frozen anchors.

Declared grammar (closed under dimension discipline):
  quadratic classes (identical to census_mustar.py, re-derived):
    P: m_a*m_b vs A^2 (324); T: m_x^2 vs m_a*m_b (252); A: A*m_a vs m_x^2 (728)
  linear classes (targets: 9 masses + 9 anchors = 18):
    S: unordered pair sums        m_a + m_b      vs single scale
    D: unordered pair differences |m_a - m_b|    vs single scale
    K: Descartes-Soddy completions of unordered triples, bends = sqrt(m):
       k4 = k_a+k_b+k_c +/- 2*sqrt(k_a k_b + k_b k_c + k_c k_a),
       compared as k4^2           vs single scale

Exclusions (structural, declared in advance):
  - algebraic identities: value identical to target by construction;
  - hierarchy trivialities (S/D/K only): no target whose construction contains
    an operand. Adding a negligible mass to a large one returns the large one
    trivially (e+t ~ t at 0.0003%); such comparisons are excluded the same way
    class T excludes x in {a,b}. Anchor construction sets: mustar, G,
    mustar+/-tau, two_mustar, 1.5*mustar built from {e,mu,tau}; P_2me from {e};
    two_tau and anchor-tau from {tau}.
Bosons are outside the declared charged-fermion universe by construction.
"""
import csv, itertools

# ---- PDG 2026 inputs (MeV), identical to census_mustar.py -------------------
M_E, M_MU, M_TAU = 0.51099895069, 105.6583755, 1776.930
U_2GEV, D_2GEV, S_2GEV = 2.16, 4.70, 92.9
R_M = 1.01750
M_C, M_B, M_T = 1272.9, 4186.0, 172600.0

MUSTAR = M_E + M_MU + M_TAU
G = (1.5) ** 0.5 * MUSTAR
P = 2 * M_E

masses = {
    "e": M_E, "mu": M_MU, "tau": M_TAU,
    "u": U_2GEV * R_M, "d": D_2GEV * R_M, "s": S_2GEV * R_M,
    "c": M_C, "b": M_B, "t": M_T,
}
anchors = {
    "mustar": MUSTAR, "P_2me": P, "G": G, "tau": M_TAU,
    "mustar_plus_tau": MUSTAR + M_TAU, "two_tau": 2 * M_TAU,
    "mustar_over_2": MUSTAR / 2, "two_mustar": 2 * MUSTAR,
    "one_point_five_mustar": 1.5 * MUSTAR,
}
CONSTRUCTION = {  # which universe members each target is built from
    **{n: {n} for n in masses},
    "A:mustar": {"e", "mu", "tau"}, "A:P_2me": {"e"}, "A:G": {"e", "mu", "tau"},
    "A:tau": {"tau"}, "A:mustar_plus_tau": {"e", "mu", "tau"},
    "A:two_tau": {"tau"}, "A:mustar_over_2": {"e", "mu", "tau"},
    "A:two_mustar": {"e", "mu", "tau"}, "A:one_point_five_mustar": {"e", "mu", "tau"},
}
targets = dict(masses)
targets.update({f"A:{k}": v for k, v in anchors.items()})

rows = []
def add(cls, val_name, val, tgt_name, tgt, operands=None):
    if tgt <= 0:
        return
    if operands is not None and CONSTRUCTION[tgt_name] & set(operands):
        return  # hierarchy-trivial: target contains an operand
    if abs(val - tgt) < 1e-12 * max(abs(tgt), 1e-300):
        return  # algebraic identity, excluded as trivial
    rows.append({"class": cls, "value": val_name, "target": tgt_name,
                 "value_num": f"{val:.6g}", "target_num": f"{tgt:.6g}",
                 "residual_pct": f"{(val / tgt - 1) * 100:+.6f}", "note": ""})

names = list(masses)
tnames = list(targets)

# ---- quadratic classes (census_mustar.py verbatim logic, unchanged) --------
for a, b in itertools.combinations(names, 2):
    for A, Av in anchors.items():
        add("P", f"{a}*{b}", masses[a] * masses[b], f"{A}^2", Av * Av)
for a, b in itertools.combinations(names, 2):
    for x in names:
        if x in (a, b):
            continue
        add("T", f"{x}^2", masses[x] ** 2, f"{a}*{b}", masses[a] * masses[b])
for A, Av in anchors.items():
    for a in names:
        for x in names:
            add("A", f"{A}*{a}", Av * masses[a], f"{x}^2", masses[x] ** 2)

# ---- linear classes (hierarchy-trivial exclusion active) -------------------
for a, b in itertools.combinations(names, 2):
    for tn in tnames:
        add("S", f"{a}+{b}", masses[a] + masses[b], tn, targets[tn], (a, b))
for a, b in itertools.combinations(names, 2):
    for tn in tnames:
        add("D", f"|{a}-{b}|", abs(masses[a] - masses[b]), tn, targets[tn], (a, b))

bends = {n: m ** 0.5 for n, m in masses.items()}
for a, b, c in itertools.combinations(names, 3):
    s = bends[a] + bends[b] + bends[c]
    w = 2 * (bends[a] * bends[b] + bends[b] * bends[c] + bends[c] * bends[a]) ** 0.5
    for k4, sym in ((s + w, "+"), (abs(s - w), "-")):
        for tn in tnames:
            add("K", f"k4({a},{b},{c},{sym})^2", k4 * k4, tn, targets[tn], (a, b, c))

assert sum(1 for r in rows if r["class"] in "PTA") == 1304

DECLARED = {("c*b", "G^2"), ("mustar*d", "s^2"), ("P_2me*d", "u^2"),
            ("k4(e,mu,tau,-)^2", "s")}
for r in rows:
    if (r["value"], r["target"]) in DECLARED:
        r["note"] = "declared relation"

with open("extended_census.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys(), lineterminator="\n")
    w.writeheader(); w.writerows(rows)

res = sorted(((float(r["residual_pct"]), r) for r in rows), key=lambda t: abs(t[0]))
ladder = (0.1, 0.25, 0.5, 1.0)
by = {c: [r for r in rows if r["class"] == c] for c in "PTASDK"}
lines = ["extended mu* census, PDG 2026 central values",
         f"comparisons: {len(rows)} (" + ", ".join(f"{c}={len(by[c])}" for c in "PTASDK") + ")",
         ""]
for scope, sel in (("all classes", rows), *[(f"class {c}", by[c]) for c in "SDK"]):
    sw = {t: sum(1 for r in sel if abs(float(r["residual_pct"])) < t) for t in ladder}
    lines.append(f"{scope}: {len(sel)} comparisons; tolerance sweep "
                 + "/".join(str(sw[t]) for t in ladder))
lines.append("\ntwenty smallest |residual| (all classes):")
for x, r in res[:20]:
    lines.append(f"  {x:+.4f}%  [{r['class']}] {r['value']} vs {r['target']}  {r['note']}")
nd = next((x, r) for x, r in res if not r["note"])
lines.append(f"\nnearest non-declared comparison: {nd[0]:+.4f}%  "
             f"[{nd[1]['class']}] {nd[1]['value']} vs {nd[1]['target']}")

with open("extended_summary.txt", "w") as f:
    f.write("\n".join(lines) + "\n")
print("\n".join(lines))
