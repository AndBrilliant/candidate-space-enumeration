#!/usr/bin/env python3
"""Null-model sensitivity suite for the chronicle census MCs (APPB A1-A4 style).

EXPLORATORY scratch analysis. Manuscript v5 is frozen on Overleaf and is NOT
touched by anything here. Purpose: check whether the two MC conclusions are
stable under reasonable changes of the null model:

  Part A  whole-grammar null variants (Test A family), N=1e4 universes each:
    A1 baseline      leptons log-U[m_e,m_tau], quarks log-U[m_u(mu*),m_t]
                     (reproduces lite_mc seed 20261001 exactly; cross-checked
                     against the frozen lite_mc_universes.csv)
    A2 koide_surface lepton triples constrained to the Koide surface Q=2/3
                     exactly (m_k = M(1+sqrt2 cos(theta+2pi k/3))^2, M log-U
                     over the feasible range, theta uniform, in-span accept),
                     quarks as A1. Question: do Koide-true leptons alone
                     inflate grammar hits?
    A3 span_x2       both log-spans doubled about their log-centers
    A4 span_x0p5     both log-spans halved about their log-centers
    A5 uniform       uniform-in-mass (not log) over the baseline spans
  Part B  Koide rarity P(|Q-2/3|<2.2e-6) under lepton-null variants:
    B1 baseline log-U (2e8 triples; documented value 4.89e-6, 979 hits)
    B2 uniform-in-mass, B3 span x2, B4 span x0.5 (1e8 each)
  Part D  conditional cascade slot rates under downstream quark-null variants:
    D0 rejection reference: leptons from log-U null accepted at
       |Q-2/3|<2.2e-6 (the original Test B ensemble; documented slot rates
       1.10e-3, 7.17e-4, 1.08e-3, 4.58e-4, 1.35e-4, 2.59e-4)
    D1 surface leptons + baseline quark span
    D2/D3 quark log-span x2 / x0.5
    D4 quarks uniform-in-mass
    Slots (matching-or-better vs the observed PDG residuals):
      s vs alpha_K^2 mu* ; mu* d vs s^2 ; (2e) d vs u^2 ;
      c vs 3 alpha_K mu* ; b vs (mu*/2)/alpha_K ; c b vs 1.5 mu*^2
    Joint products quoted for 5 slots (slot 6 is a near-consequence of
    slots 4&5) and all 6.

Seeds frozen per variant. Run: python3 sensitivity_mc.py {A|B|D|summary|all}
Outputs: sensitivity_A_universes.csv, sensitivity_B.csv, sensitivity_D.csv,
         sensitivity_summary.txt
"""
import argparse, csv, math, os, random, sys, time
import numpy as np
from scipy.stats import norm as _norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from census_grammar import COMPS, PDG
import lite_mc

HERE = os.path.dirname(os.path.abspath(__file__))
AK = math.sqrt(1.5) - 1.0
LEP_LO, LEP_HI = lite_mc.LEP_LO, lite_mc.LEP_HI
Q_LO, Q_HI = lite_mc.Q_LO, lite_mc.Q_HI
LADDER = lite_mc.LADDER
SHARP = lite_mc.SHARP
N = 10000
KOIDE_W = 2.2e-6

# ---------- vectorized grammar evaluator (cross-checked vs lite_mc) ----------

def _class_arrays():
    out = {}
    for op in "PTASDK":
        rows = [(i, j, k, sgn, -1 if tn is None else tn)
                for (o, (i, j, k, sgn), tn) in COMPS if o == op]
        out[op] = [np.array([r[c] for r in rows], dtype=np.int64) for c in range(5)]
    return out

CLS = _class_arrays()

def residuals_matrix(M):
    """M: (n,9) masses -> dict of per-class residual arrays (n, n_class)."""
    sq = np.sqrt(M)
    mustar = M[:, 0] + M[:, 1] + M[:, 2]
    A = np.stack([mustar, 2 * M[:, 0], math.sqrt(1.5) * mustar, M[:, 2],
                  mustar + M[:, 2], 2 * M[:, 2], mustar / 2, 2 * mustar,
                  1.5 * mustar], axis=1)
    T = np.concatenate([M, A], axis=1)
    res = {}
    i, j, k, sgn, tn = CLS["P"]
    res["P"] = M[:, i] * M[:, j] / T[:, tn] ** 2 - 1
    i, j, k, sgn, tn = CLS["T"]
    res["T"] = M[:, k] ** 2 / (M[:, i] * M[:, j]) - 1
    i, j, k, sgn, tn = CLS["A"]
    res["A"] = T[:, 9 + i] * M[:, j] / M[:, k] ** 2 - 1
    i, j, k, sgn, tn = CLS["S"]
    res["S"] = (M[:, i] + M[:, j]) / T[:, tn] - 1
    i, j, k, sgn, tn = CLS["D"]
    res["D"] = np.abs(M[:, i] - M[:, j]) / T[:, tn] - 1
    i, j, k, sgn, tn = CLS["K"]
    s = sq[:, i] + sq[:, j] + sq[:, k]
    w = 2 * np.sqrt(sq[:, i] * sq[:, j] + sq[:, j] * sq[:, k] + sq[:, k] * sq[:, i])
    k4 = np.where(sgn[None, :] > 0, s + w, np.abs(s - w))
    res["K"] = k4 * k4 / T[:, tn] - 1
    return res

def universe_stats_vec(M):
    res = residuals_matrix(M)
    R = np.abs(np.concatenate([res[op] for op in "PTASDK"], axis=1))
    hits = [(R < thr / 100).sum(axis=1) for thr in LADDER]
    best = R.min(axis=1)
    k05 = (np.abs(res["K"]) < 0.005).sum(axis=1)
    sd05 = (np.abs(np.concatenate([res["S"], res["D"]], axis=1)) < 0.005).sum(axis=1)
    pta05 = (np.abs(np.concatenate([res["P"], res["T"], res["A"]], axis=1)) < 0.005).sum(axis=1)
    return hits, best, k05, sd05, pta05

def self_check():
    """Vectorized evaluator must reproduce lite_mc.universe_stats exactly."""
    rng = random.Random(7)
    spectra = [list(PDG)] + [[math.exp(rng.uniform(math.log(0.4), math.log(2e5)))
                              for _ in range(9)] for _ in range(50)]
    M = np.array(spectra)
    hits, best, k05, sd05, pta05 = universe_stats_vec(M)
    for n, m in enumerate(spectra):
        h2, k2, p2, s2, b2, _ = lite_mc.universe_stats(m)
        assert [h[n] for h in hits] == h2, (n, "hits")
        assert (k05[n], pta05[n], sd05[n]) == (k2, p2, s2), (n, "class05")
        assert abs(best[n] - b2) < 1e-12, (n, "best")
    print("self-check: vectorized evaluator == lite_mc.universe_stats on 51 spectra")

# ---------- generators ----------

def logspan(lo, hi, factor):
    c = 0.5 * (math.log(lo) + math.log(hi))
    h = 0.5 * (math.log(hi) - math.log(lo)) * factor
    return math.exp(c - h), math.exp(c + h)

def gen_logu(n, rng, lo, hi):
    return np.exp(rng.uniform(math.log(lo), math.log(hi), n))

def koide_surface_triples(n, rng, lo=LEP_LO, hi=LEP_HI):
    """Exact Q=2/3 triples: m_k = M(1+sqrt2 cos(theta+2pi k/3))^2, M log-U over
    the feasible range [lo*3/6, hi*3/6]... full feasible M range is
    [3*lo/6, 3*hi/6] only if all masses equal; use wide [lo/2, hi/2] and
    accept in-span triples."""
    m_lo, m_hi = lo / 2.0, hi / 2.0
    out, need = [], n
    while need > 0:
        c = max(4 * need, 8192)
        Mv = gen_logu(c, rng, m_lo, m_hi)
        th = rng.uniform(0, 2 * math.pi, c)
        f = 1 + math.sqrt(2) * np.cos(th[:, None] + 2 * math.pi * np.arange(3) / 3)
        trip = Mv[:, None] * f * f
        ok = (trip >= lo).all(axis=1) & (trip <= hi).all(axis=1)
        out.append(trip[ok])
        need -= int(ok.sum())
    return np.concatenate(out)[:n]

# ---------- Part A ----------

def part_A():
    t0 = time.time()
    self_check()
    variants = {}

    # A1 baseline: reproduce lite_mc's exact MT stream
    rng = random.Random(20261001)
    ll_lo, ll_hi = math.log(LEP_LO), math.log(LEP_HI)
    lq_lo, lq_hi = math.log(Q_LO), math.log(Q_HI)
    M1 = np.array([[math.exp(rng.uniform(ll_lo, ll_hi)) for _ in range(3)]
                   + [math.exp(rng.uniform(lq_lo, lq_hi)) for _ in range(6)]
                   for _ in range(N)])
    variants["A1_baseline"] = M1

    rng2 = np.random.default_rng(20261002)
    trips = koide_surface_triples(N, rng2)
    M2 = np.concatenate([trips, gen_logu(N * 6, rng2, Q_LO, Q_HI).reshape(N, 6)], axis=1)
    variants["A2_koide_surface"] = M2

    for tag, fac, seed in [("A3_span_x2", 2.0, 20261003), ("A4_span_x0p5", 0.5, 20261004)]:
        r = np.random.default_rng(seed)
        llo, lhi = logspan(LEP_LO, LEP_HI, fac)
        qlo, qhi = logspan(Q_LO, Q_HI, fac)
        M = np.concatenate([gen_logu(N * 3, r, llo, lhi).reshape(N, 3),
                            gen_logu(N * 6, r, qlo, qhi).reshape(N, 6)], axis=1)
        variants[tag] = M

    r5 = np.random.default_rng(20261005)
    M5 = np.concatenate([r5.uniform(LEP_LO, LEP_HI, N * 3).reshape(N, 3),
                         r5.uniform(Q_LO, Q_HI, N * 6).reshape(N, 6)], axis=1)
    variants["A5_uniform"] = M5

    fields = ["variant", "universe"] + [f"hits_{str(t).replace('.','p')}pct" for t in LADDER] + \
             ["k_hits_0p5pct", "pta_hits_0p5pct", "sd_hits_0p5pct",
              "min_abs_residual_pct", "any_shorter_than_real_top"]
    with open(os.path.join(HERE, "sensitivity_A_universes.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        for tag, M in variants.items():
            hits, best, k05, sd05, pta05 = universe_stats_vec(M)
            for u in range(N):
                w.writerow({"variant": tag, "universe": u,
                            **{f"hits_{str(t).replace('.','p')}pct": int(hits[n][u])
                               for n, t in enumerate(LADDER)},
                            "k_hits_0p5pct": int(k05[u]), "pta_hits_0p5pct": int(pta05[u]),
                            "sd_hits_0p5pct": int(sd05[u]),
                            "min_abs_residual_pct": f"{best[u]*100:.6f}",
                            "any_shorter_than_real_top": int(best[u] < SHARP)})
            print(f"  {tag}: done ({time.time()-t0:.0f}s)", flush=True)

    # cross-check A1 against the frozen lite_mc CSV
    with open(os.path.join(HERE, "lite_mc_universes.csv")) as f:
        rows = list(csv.DictReader(f))
    hits1, best1, _, _, _ = universe_stats_vec(variants["A1_baseline"])
    for u, row in enumerate(rows):
        for n, t in enumerate(LADDER):
            assert int(row[f"hits_{str(t).replace('.','p')}pct"]) == int(hits1[n][u]), (u, t)
        assert abs(float(row["min_abs_residual_pct"]) - best1[u] * 100) < 1e-5, u
    print("cross-check: A1 reproduces frozen lite_mc_universes.csv 10000/10000")

# ---------- Part B ----------

def koide_hits(n_total, rng, gen, chunk=5_000_000):
    hits = got = 0
    while got < n_total:
        c = min(chunk, n_total - got)
        m = gen(c, rng)
        q = m.sum(axis=1) / np.sqrt(m).sum(axis=1) ** 2
        hits += int((np.abs(q - 2.0 / 3.0) < KOIDE_W).sum())
        got += c
    return hits, got

def sigma(p):
    # one-sided Gaussian equivalent, stable for tiny p (1-p underflows otherwise)
    if not 0 < p < 0.5:
        return float("nan")
    return float(_norm.isf(p))  # inverse survival function, stable for tiny p

def part_B():
    t0 = time.time()
    rows = []
    cfgs = [
        ("B1_baseline_logU", 200_000_000, 20261011,
         lambda n, r: gen_logu(n * 3, r, LEP_LO, LEP_HI).reshape(n, 3)),
        ("B2_uniform", 100_000_000, 20261012,
         lambda n, r: r.uniform(LEP_LO, LEP_HI, n * 3).reshape(n, 3)),
        ("B3_span_x2", 100_000_000, 20261013,
         lambda n, r: gen_logu(n * 3, r, *logspan(LEP_LO, LEP_HI, 2.0)).reshape(n, 3)),
        ("B4_span_x0p5", 100_000_000, 20261014,
         lambda n, r: gen_logu(n * 3, r, *logspan(LEP_LO, LEP_HI, 0.5)).reshape(n, 3)),
    ]
    for tag, n_tot, seed, gen in cfgs:
        rng = np.random.default_rng(seed)
        h, tot = koide_hits(n_tot, rng, gen)
        p = h / tot
        rows.append({"variant": tag, "triples": tot, "koide_hits": h,
                     "p": f"{p:.4e}", "sigma_equiv": f"{sigma(p):.2f}"})
        print(f"  {tag}: {h}/{tot} = {p:.4e} ({sigma(p):.2f} sigma) "
              f"({time.time()-t0:.0f}s)", flush=True)
    with open(os.path.join(HERE, "sensitivity_B.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys(), lineterminator="\n")
        w.writeheader(); w.writerows(rows)

# ---------- Part D ----------

# observed slot thresholds from the frozen PDG spectrum
_MS = PDG[0] + PDG[1] + PDG[2]
_OBS = [abs(PDG[5] / (AK ** 2 * _MS) - 1),          # s vs alpha_K^2 mu*
        abs(_MS * PDG[4] / PDG[5] ** 2 - 1),        # mu* d vs s^2
        abs(2 * PDG[0] * PDG[4] / PDG[3] ** 2 - 1),  # (2e) d vs u^2
        abs(PDG[6] / (3 * AK * _MS) - 1),            # c vs 3 alpha_K mu*
        abs(PDG[7] / (0.5 * _MS / AK) - 1),          # b vs (mu*/2)/alpha_K
        abs(PDG[6] * PDG[7] / (1.5 * _MS ** 2) - 1)]  # c b vs 1.5 mu*^2
SLOT_NAMES = ["s/F^2", "mu*d/s^2", "P*d/u^2", "c/(3aK mu*)", "b/(0.5 mu*/aK)", "cb/G^2"]

def slot_rates(lep, quarks):
    """lep: (T,3) leptons; quarks: (T,D,6) u,d,s,c,b,t. -> per-slot pooled rates."""
    mu_s = lep.sum(axis=1)[:, None]          # (T,1)
    p2e = 2 * lep[:, 0][:, None]
    u, d, s = quarks[:, :, 0], quarks[:, :, 1], quarks[:, :, 2]
    cq, bq = quarks[:, :, 3], quarks[:, :, 4]
    r = [np.abs(s / (AK ** 2 * mu_s) - 1),
         np.abs(mu_s * d / s ** 2 - 1),
         np.abs(p2e * d / u ** 2 - 1),
         np.abs(cq / (3 * AK * mu_s) - 1),
         np.abs(bq / (0.5 * mu_s / AK) - 1),
         np.abs(cq * bq / (1.5 * mu_s ** 2) - 1)]
    return [float((x <= t).mean()) for x, t in zip(r, _OBS)], r

def part_D():
    t0 = time.time()
    n_tr, n_dr = 2000, 1000
    rows = []

    # D0 reference: rejection-sampled leptons from the log-U null (original Test B)
    rng0 = np.random.default_rng(20261019)
    acc, drawn = [], 0
    while sum(len(a) for a in acc) < 1000 and drawn < 300_000_000:
        c = 5_000_000
        m = gen_logu(c * 3, rng0, LEP_LO, LEP_HI).reshape(c, 3)
        drawn += c
        q = m.sum(axis=1) / np.sqrt(m).sum(axis=1) ** 2
        acc.append(m[np.abs(q - 2.0 / 3.0) < KOIDE_W])
    lep0 = np.concatenate(acc)[:1000]
    print(f"  D0: {len(lep0)} rejection leptons from {drawn:.2e} draws "
          f"({time.time()-t0:.0f}s)", flush=True)
    q0 = gen_logu(len(lep0) * n_dr * 6, rng0, Q_LO, Q_HI).reshape(len(lep0), n_dr, 6)
    rates0, _ = slot_rates(lep0, q0)
    rows.append({"variant": "D0_rejection_ref", "draws": len(lep0) * n_dr,
                 **{n_: f"{r:.4e}" for n_, r in zip(SLOT_NAMES, rates0)}})
    print(f"  D0 rates: {['%.2e' % r for r in rates0]}", flush=True)

    # D1-D4: surface leptons + quark-null variants
    rng = np.random.default_rng(20261020)
    lep = koide_surface_triples(n_tr, rng)
    qgens = [("D1_baseline", Q_LO, Q_HI, "log"),
             ("D2_span_x2", *logspan(Q_LO, Q_HI, 2.0), "log"),
             ("D3_span_x0p5", *logspan(Q_LO, Q_HI, 0.5), "log"),
             ("D4_uniform", Q_LO, Q_HI, "uni")]
    for tag, lo, hi, kind in qgens:
        if kind == "log":
            qk = gen_logu(n_tr * n_dr * 6, rng, lo, hi).reshape(n_tr, n_dr, 6)
        else:
            qk = rng.uniform(lo, hi, n_tr * n_dr * 6).reshape(n_tr, n_dr, 6)
        rates, r = slot_rates(lep, qk)
        j5 = float(np.prod(rates[:5])); j6 = j5 * rates[5]
        rows.append({"variant": tag, "draws": n_tr * n_dr,
                     **{n_: f"{x:.4e}" for n_, x in zip(SLOT_NAMES, rates)},
                     "joint5": f"{j5:.3e}", "joint6": f"{j6:.3e}",
                     "sigma5": f"{sigma(j5):.2f}"})
        # fixed-tolerance tab for D1
        if tag == "D1_baseline":
            for tol in (0.005, 0.01):
                fx = [float((x <= tol).mean()) for x in r]
                rows.append({"variant": f"D1_tol_{tol}", "draws": n_tr * n_dr,
                             **{n_: f"{x:.4e}" for n_, x in zip(SLOT_NAMES, fx)}})
        print(f"  {tag}: joint5 {j5:.3e} ({sigma(j5):.2f} sigma) "
              f"({time.time()-t0:.0f}s)", flush=True)

    fields = ["variant", "draws"] + SLOT_NAMES + ["joint5", "joint6", "sigma5"]
    with open(os.path.join(HERE, "sensitivity_D.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fields})

# ---------- summary ----------

def part_summary():
    lines = ["null-model sensitivity suite (exploratory; manuscript v5 untouched)",
             f"grammar: 3740 comparisons, anchors rebuilt per universe; N={N} universes (A);",
             f"Koide window |Q-2/3|<{KOIDE_W}; slot criterion: matching-or-better vs PDG residual.",
             "slot thresholds from PDG: " + ", ".join(f"{n} {t*100:.3f}%" for n, t in zip(SLOT_NAMES, _OBS)),
             "", "== Part A: whole-grammar null variants =="]
    import collections
    data = collections.defaultdict(list)
    with open(os.path.join(HERE, "sensitivity_A_universes.csv")) as f:
        for row in csv.DictReader(f):
            data[row["variant"]].append(row)
    for tag, rows in data.items():
        def col(k):
            return [int(r[k]) for r in rows]
        c05 = col("hits_0p5pct")
        c1 = col("hits_1p0pct")
        best = [float(r["min_abs_residual_pct"]) for r in rows]
        sharper = sum(int(r["any_shorter_than_real_top"]) for r in rows)
        lines.append(f"{tag}: sub-0.5% mean {sum(c05)/len(c05):.3f}/univ, "
                     f"P(>=1) {sum(x>=1 for x in c05)/len(c05):.4f}; "
                     f"sub-1% mean {sum(c1)/len(c1):.3f}; "
                     f"P(best<0.1742%) {sharper/len(rows):.4f}; "
                     f"K05 mean {sum(col('k_hits_0p5pct'))/len(rows):.3f}")
    lines += ["", "== Part B: Koide rarity under lepton nulls =="]
    with open(os.path.join(HERE, "sensitivity_B.csv")) as f:
        for row in csv.DictReader(f):
            lines.append(f"{row['variant']}: {row['koide_hits']}/{row['triples']} = "
                         f"{row['p']} ({row['sigma_equiv']} sigma)")
    lines += ["documented 2026-09 value: 4.89e-6 (979/2e8, ~4.4 sigma).", "",
              "== Part D: conditional cascade slot rates =="]
    with open(os.path.join(HERE, "sensitivity_D.csv")) as f:
        for row in csv.DictReader(f):
            seg = " ".join(f"{n}={row[n]}" for n in SLOT_NAMES)
            lines.append(f"{row['variant']}: {seg}" +
                         (f" | joint5 {row.get('joint5','')} ({row.get('sigma5','')} sigma)"
                          if row.get("joint5") else ""))
    lines += ["documented 2026-09 slot rates: 1.10e-3, 7.17e-4, 1.08e-3, 4.58e-4, 1.35e-4, 2.59e-4;",
              "documented joint5 ~5.3e-17 (~8 sigma); full chain with P_Koide ~1e-22."]
    with open(os.path.join(HERE, "sensitivity_summary.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("part", choices=["A", "B", "D", "summary", "all"])
    a = ap.parse_args()
    if a.part in ("A", "all"):
        part_A()
    if a.part in ("B", "all"):
        part_B()
    if a.part in ("D", "all"):
        part_D()
    if a.part in ("summary", "all"):
        part_summary()
