#!/usr/bin/env python3
"""Shared census grammar: the declared comparison space as a function of nine
charged-fermion masses (order: e, mu, tau, u, d, s, c, b, t).

Anchors are derived from the lepton triple (indices 0,1,2) by the frozen
rules; targets are the nine masses plus the nine anchors. Classes:
  P: m_a*m_b vs A^2          T: m_x^2 vs m_a*m_b (x not in {a,b})
  A: A*m_a vs m_x^2          S: m_a+m_b vs scale   D: |m_a-m_b| vs scale
  K: Descartes-Soddy completion k4^2 vs scale, bends sqrt(m), both signs
S/D/K exclude targets whose construction contains an operand (hierarchy
trivialities). A excludes its one algebraic identity (tau-anchor * tau vs
tau^2). The comparison structure is identical for any input spectrum.
"""
import itertools

NAMES = ["e", "mu", "tau", "u", "d", "s", "c", "b", "t"]
ANCHOR_KEYS = ["mustar", "P_2me", "G", "tau", "mustar_plus_tau", "two_tau",
               "mustar_over_2", "two_mustar", "one_point_five_mustar"]
LEP = {0, 1, 2}
# construction sets: index targets 0-8 are the masses; 9-17 the anchors
CONSTRUCTION = {i: {i} for i in range(9)}
CONSTRUCTION.update({9: LEP, 10: {0}, 11: LEP, 12: {2}, 13: LEP, 14: {2},
                     15: LEP, 16: LEP, 17: LEP})

def anchors_of(m):
    """m: sequence of 9 masses -> 9 anchor values, frozen construction."""
    mustar = m[0] + m[1] + m[2]
    return [mustar, 2 * m[0], 1.5 ** 0.5 * mustar, m[2],
            mustar + m[2], 2 * m[2], mustar / 2, 2 * mustar, 1.5 * mustar]

def build():
    """Static comparison list: (op, (i,j,k,sign), tgt_index_0..17)."""
    comps = []
    idx = range(9)
    pairs = list(itertools.combinations(idx, 2))
    for a, b in pairs:                                    # P
        for A in range(9):
            comps.append(("P", (a, b, -1, 0), 9 + A))
    for a, b in pairs:                                    # T
        for x in idx:
            if x not in (a, b):
                comps.append(("T", (a, b, x, 0), None))
    for A in range(9):                                    # A
        for a in idx:
            for x in idx:
                if A == 3 and a == 2 and x == 2:
                    continue  # tau*tau vs tau^2 identity
                comps.append(("A", (A, a, x, 0), None))
    for a, b in pairs:                                    # S
        for tn in range(18):
            if not CONSTRUCTION[tn] & {a, b}:
                comps.append(("S", (a, b, -1, 0), tn))
    for a, b in pairs:                                    # D
        for tn in range(18):
            if not CONSTRUCTION[tn] & {a, b}:
                comps.append(("D", (a, b, -1, 0), tn))
    for a, b, c in itertools.combinations(idx, 3):        # K
        for sgn in (1, -1):
            for tn in range(18):
                if not CONSTRUCTION[tn] & {a, b, c}:
                    comps.append(("K", (a, b, c, sgn), tn))
    return comps

COMPS = build()

def targets_of(m):
    return list(m) + anchors_of(m)

def residuals(m):
    """All comparison residuals (val/tgt - 1) for spectrum m (9 masses)."""
    t = targets_of(m)
    sq = [x ** 0.5 for x in m]
    out = []
    for op, (i, j, k, sgn), tn in COMPS:
        if op == "P":
            val, tgt = m[i] * m[j], t[tn] ** 2
        elif op == "T":
            val, tgt = m[k] ** 2, m[i] * m[j]
        elif op == "A":
            val, tgt = t[9 + i] * m[j], m[k] ** 2
        elif op == "S":
            val, tgt = m[i] + m[j], t[tn]
        elif op == "D":
            val, tgt = abs(m[i] - m[j]), t[tn]
        else:  # K
            s = sq[i] + sq[j] + sq[k]
            w = 2 * (sq[i] * sq[j] + sq[j] * sq[k] + sq[k] * sq[i]) ** 0.5
            k4 = s + w if sgn > 0 else abs(s - w)
            val, tgt = k4 * k4, t[tn]
        if tgt > 0:
            out.append(val / tgt - 1.0)
    return out

# frozen PDG 2026 central spectrum (MeV): light quarks at mu* via R_m
R_M = 1.01750
PDG = [0.51099895069, 105.6583755, 1776.930,
       2.16 * R_M, 4.70 * R_M, 92.9 * R_M, 1272.9, 4186.0, 172600.0]

if __name__ == "__main__":
    from collections import Counter
    assert len(COMPS) == 3740, len(COMPS)
    res = residuals(PDG)
    assert len(res) == 3740, len(res)
    ladder = (0.1, 0.25, 0.5, 1.0)
    c = Counter()
    for r in res:
        for t_ in ladder:
            if abs(r) * 100 < t_:
                c[t_] += 1
    print("grammar self-check: 3740 comparisons; PDG ladder:",
          [c[t_] for t_ in ladder], "(expect [0, 1, 3, 6])")
