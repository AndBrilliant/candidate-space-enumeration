"""Six-class census over the declared grammar (Marcus §II).

Comparison structure identical to the frozen historical grammar
(pipeline/census_grammar.py), re-derived here as the manuscript-facing
lineage. Residual convention V/T - 1 (Marcus §VIII)."""
import itertools, math, os, sys, csv
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from inputs import load

LEP = {0, 1, 2}
ANCHOR_KEYS = ["mustar", "P_2me", "G", "tau", "mustar_plus_tau", "two_tau",
               "mustar_over_2", "two_mustar", "one_point_five_mustar"]

def anchors_of(m):
    mustar = m[0] + m[1] + m[2]
    return [mustar, 2 * m[0], math.sqrt(1.5) * mustar, m[2],
            mustar + m[2], 2 * m[2], mustar / 2, 2 * mustar, 1.5 * mustar]

CONSTRUCTION = {i: {i} for i in range(9)}
CONSTRUCTION.update({9: LEP, 10: {0}, 11: LEP, 12: {2}, 13: LEP, 14: {2},
                     15: LEP, 16: LEP, 17: LEP})

def build():
    comps = []
    idx = range(9)
    pairs = list(itertools.combinations(idx, 2))
    for a, b in pairs:
        for A in range(9):
            comps.append(("P", (a, b, -1, 0), 9 + A))
    for a, b in pairs:
        for x in idx:
            if x not in (a, b):
                comps.append(("T", (a, b, x, 0), None))
    for A in range(9):
        for a in idx:
            for x in idx:
                if A == 3 and a == 2 and x == 2:
                    continue
                comps.append(("A", (A, a, x, 0), None))
    for a, b in pairs:
        for tn in range(18):
            if not CONSTRUCTION[tn] & {a, b}:
                comps.append(("S", (a, b, -1, 0), tn))
    for a, b in pairs:
        for tn in range(18):
            if not CONSTRUCTION[tn] & {a, b}:
                comps.append(("D", (a, b, -1, 0), tn))
    for a, b, c in itertools.combinations(idx, 3):
        for sgn in (1, -1):
            for tn in range(18):
                if not CONSTRUCTION[tn] & {a, b, c}:
                    comps.append(("K", (a, b, c, sgn), tn))
    return comps

COMPS = build()

def residuals(m):
    t = list(m) + anchors_of(m)
    sq = [x ** 0.5 for x in m]
    out = []
    for op, (i, j, k, sgn), tn in COMPS:
        if op == "P":   val, tgt = m[i] * m[j], t[tn] ** 2
        elif op == "T": val, tgt = m[k] ** 2, m[i] * m[j]
        elif op == "A": val, tgt = t[9 + i] * m[j], m[k] ** 2
        elif op == "S": val, tgt = m[i] + m[j], t[tn]
        elif op == "D": val, tgt = abs(m[i] - m[j]), t[tn]
        else:
            s = sq[i] + sq[j] + sq[k]
            w = 2 * (sq[i] * sq[j] + sq[j] * sq[k] + sq[k] * sq[i]) ** 0.5
            k4 = s + w if sgn > 0 else abs(s - w)
            val, tgt = k4 * k4, t[tn]
        if tgt > 0:
            out.append(val / tgt - 1.0)   # V/T - 1
    return out

def class_counts():
    return Counter(op for op, _, _ in COMPS)

def tolerance_ladder(res, tols):
    c = Counter()
    for r in res:
        for tp in tols:
            if abs(r) * 100 < tp:
                c[tp] += 1
    return [c[tp] for tp in tols]

def main():
    d = load()
    tols = d["mc"]["tolerances_pct"]
    res = residuals(d["masses"])
    cc = class_counts()
    lad = tolerance_ladder(res, tols)
    print("class counts:", {k: cc[k] for k in "PTASDK"}, "total", sum(cc.values()))
    print("PDG tolerance ladder:", dict(zip(tols, lad)))
    out = os.path.join(HERE, "..", "outputs", "census.csv")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["class", "residual_pct"])
        for (op, _, _), r in zip(COMPS, res):
            w.writerow([op, f"{r*100:.6f}"])
    print("wrote", out)
    return cc, lad

HERE = os.path.dirname(os.path.abspath(__file__))
if __name__ == "__main__":
    main()
