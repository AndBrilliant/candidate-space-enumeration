#!/usr/bin/env python3
"""Lite Monte Carlo null for the extended candidate-space census.

Question: in a universe whose nine charged-fermion masses are random, how
many sub-tolerance coincidences does the declared 3740-comparison grammar
manufacture? Each mock universe draws three leptons log-uniformly over the
observed lepton span and six quarks log-uniformly over the observed quark
span; the nine anchors are then rebuilt from the mock lepton triple by the
frozen construction rules (mu* = e+mu+tau, P = 2e, G = sqrt(3/2) mu*, ...),
exactly as in the real census. The grammar, exclusions, and tolerance ladder
are identical to extended_census.py via census_grammar.py.

Frozen: seed 20261001, N = 10000 universes, stdlib only (Mersenne Twister
streams are reproducible across CPython 3 versions for random()/uniform()).

Outputs: lite_mc_universes.csv (one row per universe), lite_mc_summary.txt.
"""
import csv, math, random, time
from census_grammar import COMPS, targets_of

SEED, N = 20261001, 10000
LADDER = (0.1, 0.25, 0.5, 1.0)
# observed spans (MeV), from inputs/pdg2026_inputs.csv
LEP_LO, LEP_HI = 0.51099895069, 1776.930            # m_e .. m_tau
Q_LO, Q_HI = 2.16 * 1.01750, 172600.0               # m_u(mu*) .. m_t direct
SHARP = 0.1742 / 100.0  # residual magnitude of the real top hit, as fraction

# split the static grammar into per-class primitive tuples for the hot loop
BY_OP = {op: [(i, j, k, sgn, tn) for (o, (i, j, k, sgn), tn) in COMPS if o == op]
         for op in "PTASDK"}

def universe_stats(m):
    t = targets_of(m)
    sq = [math.sqrt(x) for x in m]
    hits = [0, 0, 0, 0]
    k05 = pta05 = sd05 = 0
    best = 1e9
    for op, lst in BY_OP.items():
        for i, j, k, sgn, tn in lst:
            if op == "P":
                val, tgt = m[i] * m[j], t[tn] * t[tn]
            elif op == "T":
                val, tgt = m[k] * m[k], m[i] * m[j]
            elif op == "A":
                val, tgt = t[9 + i] * m[j], m[k] * m[k]
            elif op == "S":
                val, tgt = m[i] + m[j], t[tn]
            elif op == "D":
                val, tgt = abs(m[i] - m[j]), t[tn]
            else:
                s = sq[i] + sq[j] + sq[k]
                w = 2.0 * math.sqrt(sq[i] * sq[j] + sq[j] * sq[k] + sq[k] * sq[i])
                k4 = s + w if sgn > 0 else abs(s - w)
                val, tgt = k4 * k4, t[tn]
            if tgt <= 0:
                continue
            r = abs(val / tgt - 1.0)
            if r < best:
                best = r
            for n, thr in enumerate(LADDER):
                if r * 100 < thr:
                    hits[n] += 1
            if r < 0.005:
                if op == "K":
                    k05 += 1
                elif op in "SD":
                    sd05 += 1
                else:
                    pta05 += 1
    sharper = sum(1 for _ in [0] if best < SHARP)
    return hits, k05, pta05, sd05, best, int(best < SHARP)

def main():
    rng = random.Random(SEED)
    ll_lo, ll_hi = math.log(LEP_LO), math.log(LEP_HI)
    lq_lo, lq_hi = math.log(Q_LO), math.log(Q_HI)
    t0 = time.time()
    rows = []
    for u in range(N):
        m = ([math.exp(rng.uniform(ll_lo, ll_hi)) for _ in range(3)]
             + [math.exp(rng.uniform(lq_lo, lq_hi)) for _ in range(6)])
        hits, k05, pta05, sd05, best, sharper = universe_stats(m)
        rows.append({"universe": u,
                     **{f"hits_{str(thr).replace('.','p')}pct": hits[n]
                        for n, thr in enumerate(LADDER)},
                     "k_hits_0p5pct": k05, "pta_hits_0p5pct": pta05, "sd_hits_0p5pct": sd05,
                     "min_abs_residual_pct": f"{best * 100:.6f}",
                     "any_shorter_than_real_top": sharper})
        if (u + 1) % 1000 == 0:
            print(f"  {u + 1}/{N}  ({time.time() - t0:.0f}s)", flush=True)

    with open("lite_mc_universes.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys(), lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    keys = ["hits_0p1pct", "hits_0p25pct", "hits_0p5pct", "hits_1p0pct"]
    lines = [f"lite MC null, seed {SEED}, N={N} universes",
             "null: 3 leptons log-U[m_e,m_tau], 6 quarks log-U[m_u(mu*),m_t];",
             "anchors rebuilt per universe from the mock lepton triple;",
             "grammar identical to extended census (3740 comparisons).", ""]
    for k in keys:
        c = [int(r[k]) for r in rows]
        mean = sum(c) / N
        ge1 = sum(1 for x in c if x >= 1) / N
        ge2 = sum(1 for x in c if x >= 2) / N
        lines.append(f"{k}: mean {mean:.3f}/universe, "
                     f"P(>=1) {ge1:.4f}, P(>=2) {ge2:.4f}, max {max(c)}")
    sharper = sum(int(r["any_shorter_than_real_top"]) for r in rows)
    k05 = [int(r["k_hits_0p5pct"]) for r in rows]
    pta05 = [int(r["pta_hits_0p5pct"]) for r in rows]
    sd05 = [int(r["sd_hits_0p5pct"]) for r in rows]
    lines += ["",
              f"universes whose sharpest hit beats the real top (+0.1742%): "
              f"{sharper}/{N} = {sharper / N:.4f}",
              f"sub-0.5% hits per universe by class: "
              f"P/T/A mean {sum(pta05) / N:.3f}; "
              f"S/D mean {sum(sd05) / N:.3f}; "
              f"K mean {sum(k05) / N:.3f}, P(K>=1) {sum(1 for x in k05 if x) / N:.4f}",
              "",
              "real universe for comparison: sub-0.5% rows = 3 "
              "(+0.174% declared c*b/G^2; -0.290% Descartes accidental, "
              "counted twice via tau mass/anchor duplication); "
              "sub-1% rows = 6."]
    with open("lite_mc_summary.txt", "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))

if __name__ == "__main__":
    main()
