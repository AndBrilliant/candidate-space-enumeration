"""Whole-grammar Monte Carlo (Marcus §III): N mock universes, each with
leptons log-U[m_e, m_tau], quarks log-U[m_u(mu*), m_t]; anchors rebuilt per
universe; grammar identical to the census. Reproduces manuscript Table IV.

Residual convention V/T - 1 (Marcus §VIII). Frozen seed from inputs."""
import math, os, sys, csv, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from inputs import load
from census import COMPS, residuals
from deterministic_relations import main as det_main

HERE = os.path.dirname(os.path.abspath(__file__))

def draw_universe(rng, lep_lo, lep_hi, qu_lo, qu_hi):
    # stdlib Mersenne Twister to match the frozen historical seed stream
    ll_lo, ll_hi = math.log(lep_lo), math.log(lep_hi)
    lq_lo, lq_hi = math.log(qu_lo), math.log(qu_hi)
    lep = [math.exp(rng.uniform(ll_lo, ll_hi)) for _ in range(3)]
    qu = [math.exp(rng.uniform(lq_lo, lq_hi)) for _ in range(6)]
    return lep + qu

def run(n=None, seed=None, verbose=True):
    d = load()
    mc = d["mc"]
    n = n or mc["n_universes"]
    seed = seed if seed is not None else mc["seed_whole_grammar"]
    tols = mc["tolerances_pct"]
    rng = random.Random(seed)
    lep_lo, lep_hi = d["lep_lo"], d["lep_hi"]
    qu_lo, qu_hi = d["qu_lo"], d["qu_hi"]

    hits = {tp: [0]*n for tp in tols}
    sharpest = [0.0]*n
    # observed heavy-sector product residual (~0.174175%) from frozen inputs
    m = d["masses"]
    obs_beat = abs(m[6] * m[7] / d["G"] ** 2 - 1)

    for u_ in range(n):
        res = residuals(draw_universe(rng, lep_lo, lep_hi, qu_lo, qu_hi))
        ar = [abs(r) * 100.0 for r in res]
        for tp in tols:
            hits[tp][u_] = sum(1 for a in ar if a < tp)
        sharpest[u_] = min(ar)

    rows = []
    for tp in tols:
        h = hits[tp]
        mn = sum(h)/n; p1 = sum(1 for x in h if x>=1)/n
        p2 = sum(1 for x in h if x>=2)/n; mx = max(h)
        rows.append((tp, mn, p1, p2, mx))
    beat = sum(1 for s in sharpest if s < obs_beat*100)/n

    if verbose:
        print(f"whole-grammar MC, N={n}, seed={seed}")
        print(" tol%   mean    P(>=1)  P(>=2)  max")
        for tp, mn, p1, p2, mx in rows:
            print(f" {tp:5.2f}  {mn:6.4f}  {p1:6.4f}  {p2:6.4f}  {mx:4d}")
        print(f" beat-fraction (sharpest < {obs_beat*100:.4f}%): {beat:.4f}")

    out = os.path.join(HERE, "..", "outputs", "table_iv.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tolerance_pct", "mean_cells", "P_ge1", "P_ge2", "max"])
        for r in rows:
            w.writerow(r)
        w.writerow([])
        w.writerow(["beat_fraction_vs_observed_product", beat])
    return rows, beat

if __name__ == "__main__":
    run()
