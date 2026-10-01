#!/usr/bin/env python3
"""Grammar-size robustness of the lite MC null: what if the enumerated space
had been 2x or 3x larger?

Same 10,000 mock universes (same seed, same mass draws) evaluated under
three nested grammars:
  X1: the frozen extended census (9 anchors, 3740 comparisons)
  X2: 18 anchors (base menu + 9 declared filler constructions)
  X3: 27 anchors (base + 18 fillers)
Fillers are integer/half-integer multiples and sums of the lepton-built
scales {mu*, tau, e, G}; construction sets propagate through the same
hierarchy-triviality exclusions. Nested design: identical mock universes,
only the grammar grows, so the accident-rate scaling is measured directly.

Frozen: seed 20261001, N = 10000. stdlib + census_grammar only.
"""
import math, random
import census_grammar as cg

SEED, N = 20261001, 10000
LEP_LO, LEP_HI = 0.51099895069, 1776.930
Q_LO, Q_HI = 2.16 * 1.01750, 172600.0

F2 = ["mus/4", "3*mus", "4*mus", "3*tau", "4*tau", "mus+2*tau", "mus-tau",
      "4*e", "G/2"]
F3 = F2 + ["5*mus", "6*mus", "mus/8", "5*tau", "6*tau", "2*mus+tau",
           "3*mus+tau", "mus+3*tau", "8*e", "2*G"]

def anchor_values(m, extras):
    mu_ = m[0] + m[1] + m[2]
    G = 1.5 ** 0.5 * mu_
    base = cg.anchors_of(m)
    env = {"mus": mu_, "tau": m[2], "e": m[0], "G": G}
    return base + [eval(x, {"__builtins__": {}}, env) for x in extras]

def anchor_construction(extras):
    cons = dict(cg.CONSTRUCTION)
    for k, x in enumerate(extras):
        s = set()
        if "mus" in x or "G" in x: s |= {0, 1, 2}
        if "tau" in x: s.add(2)
        if "e" in x and "mus" not in x and "G" not in x: s = {0}
        cons[18 + k] = s
    return cons

import itertools
def build_comps(n_anchors, extras):
    cons = anchor_construction(extras)
    n_targets = 9 + n_anchors
    comps = []
    idx = range(9)
    pairs = list(itertools.combinations(idx, 2))
    for a, b in pairs:
        for A in range(n_anchors):
            comps.append(("P", (a, b, -1, 0), 9 + A))
    for a, b in pairs:
        for x in idx:
            if x not in (a, b):
                comps.append(("T", (a, b, x, 0), None))
    for A in range(n_anchors):
        for a in idx:
            for x in idx:
                if A == 3 and a == 2 and x == 2:
                    continue
                comps.append(("A", (A, a, x, 0), None))
    for a, b in pairs:
        for tn in range(n_targets):
            if not cons[tn] & {a, b}:
                comps.append(("S", (a, b, -1, 0), tn))
                comps.append(("D", (a, b, -1, 0), tn))
    for a, b, c in itertools.combinations(idx, 3):
        for sgn in (1, -1):
            for tn in range(n_targets):
                if not cons[tn] & {a, b, c}:
                    comps.append(("K", (a, b, c, sgn), tn))
    return comps

def stats(m, comps, anchors):
    t = list(m) + anchors
    sq = [math.sqrt(x) for x in m]
    h01 = h05 = 0
    best = 1e9
    for op, (i, j, k, sgn), tn in comps:
        if op == "P":   val, tgt = m[i] * m[j], t[tn] ** 2
        elif op == "T": val, tgt = m[k] ** 2, m[i] * m[j]
        elif op == "A": val, tgt = t[9 + i] * m[j], m[k] ** 2
        elif op == "S": val, tgt = m[i] + m[j], t[tn]
        elif op == "D": val, tgt = abs(m[i] - m[j]), t[tn]
        else:
            s = sq[i] + sq[j] + sq[k]
            w = 2.0 * math.sqrt(sq[i]*sq[j] + sq[j]*sq[k] + sq[k]*sq[i])
            k4 = s + w if sgn > 0 else abs(s - w)
            val, tgt = k4 * k4, t[tn]
        if tgt <= 0: continue
        r = abs(val / tgt - 1.0)
        if r < best: best = r
        if r < 0.001: h01 += 1
        if r < 0.005: h05 += 1
    return h01, h05, best

GRAMMARS = [("X1 (9 anchors)", 9, []), ("X2 (18 anchors)", 18, F2),
            ("X3 (27 anchors)", 27, F3)]
built = [(name, build_comps(n, ex), ex) for name, n, ex in GRAMMARS]
for name, comps, _ in built:
    print(f"{name}: {len(comps)} comparisons")

rng = random.Random(SEED)
ll = (math.log(LEP_LO), math.log(LEP_HI)); lq = (math.log(Q_LO), math.log(Q_HI))
acc = {name: [[], [], []] for name, _, _ in built}
for u in range(N):
    m = ([math.exp(rng.uniform(*ll)) for _ in range(3)]
         + [math.exp(rng.uniform(*lq)) for _ in range(6)])
    for name, comps, ex in built:
        h01, h05, best = stats(m, comps, anchor_values(m, ex))
        acc[name][0].append(h01); acc[name][1].append(h05); acc[name][2].append(best)

print(f"\nseed {SEED}, N={N} nested universes")
for name, _, _ in built:
    h01, h05, best = acc[name]
    print(f"{name}: sub-0.1% mean {sum(h01)/N:.3f} | sub-0.5% mean {sum(h05)/N:.3f}, "
          f"P(>=1) {sum(1 for x in h05 if x)/N:.4f} | P(best < 0.1742%) "
          f"{sum(1 for x in best if x*100 < 0.1742)/N:.4f}")
