#!/usr/bin/env python3
"""Importance-sampling cross-check of the frozen six-slot joint frequency.

Independent confirmation of the Rao-Blackwellized conditional-integration
estimate (test_b_exact.py) with a different sampling law, as described in
the manuscript: draw m_d from a log-normal proposal centered per universe on
the slot-overlap locus d* = alpha_K^4 mu*, and restore the null measure
through likelihood-ratio weights (ordinary importance sampling for a
rare-event frequency, Buckwew 2004).

  null measure (log space, span [ln Q_LO, ln Q_HI]): uniform, density 1/L
  proposal: ln d = ln(aK^4 mu*) + sigma * Z, Z ~ N(0,1)
  weight: w(d) = (sigma*sqrt(2pi)/L) * exp(((ln d - ln d*)^2) / (2 sigma^2)),
          w = 0 outside the span
  estimator: mean over universes of w * P_joint(anchors, d), where P_joint is
          the SAME conditional-analytic per-universe joint as test_b_exact
          (lepton triple law also identical: Koide-surface draws).

Four frozen seeds of 5e6 universes at two proposal widths (sigma = 0.5, 1.0
in ln units). Outputs importance_sampling.txt. Requires numpy.

Run from pipeline/:  python3 importance_sampling.py
"""
import math, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from census_grammar import PDG
from sensitivity_mc import LEP_LO, LEP_HI
from test_b_exact import (koide_surface_triples, joint_components, AK,
                          Q_LO, Q_HI)

N = 5_000_000
SEEDS = {"sigma0p5_a": (20264001, 0.5), "sigma0p5_b": (20264002, 0.5),
         "sigma1p0_a": (20264003, 1.0), "sigma1p0_b": (20264004, 1.0)}

def run(seed, sigma, n=N):
    rng = np.random.default_rng(seed)
    lep = koide_surface_triples(n, rng)
    mu_s = lep.sum(axis=1)
    p2e = 2 * lep[:, 0]
    lo, hi = math.log(Q_LO), math.log(Q_HI)
    L = hi - lo
    center = np.log(AK ** 4 * mu_s)
    lnd = center + sigma * rng.standard_normal(n)
    inb = (lnd >= lo) & (lnd <= hi)
    w = np.where(inb,
                 sigma * math.sqrt(2 * math.pi) / L
                 * np.exp(np.where(inb, (lnd - center) ** 2, np.inf)
                          / (2 * sigma ** 2)),
                 0.0)
    # np.exp(inf) overflow guard: use clipped exponent
    expo = np.clip((lnd - center) ** 2 / (2 * sigma ** 2), None, 700)
    w = np.where(inb, sigma * math.sqrt(2 * math.pi) / L * np.exp(expo), 0.0)
    joint, _ = joint_components(mu_s, p2e, np.exp(lnd), lo, hi)
    est_w = w * joint
    est = float(est_w.mean())
    err = float(est_w.std() / math.sqrt(n))
    s1, s2 = float(est_w.sum()), float((est_w ** 2).sum())
    ess = s1 * s1 / s2 if s2 > 0 else 0.0
    maxshare = float(est_w.max() / s1) if s1 > 0 else 0.0
    return est, err, ess, maxshare

def main():
    lines = [f"importance-sampling cross-check: N={N:,} per run, "
             f"span log-uniform [{Q_LO:.4g}, {Q_HI:.4g}] MeV, "
             f"proposal log-normal around d*=aK^4 mu*", ""]
    per, allw, allww, tot = [], 0.0, 0.0, []
    for tag, (seed, sigma) in SEEDS.items():
        est, err, ess, maxshare = run(seed, sigma)
        per.append((tag, sigma, est, err, ess, maxshare))
        lines.append(f"{tag}: sigma={sigma}, est {est:.3e} +- {err:.1e}; "
                     f"ESS {ess:.3e}; max weight share {maxshare:.2e}")
        tot.append(est)
        allw += est * N          # track simple means for pooled estimate
        allww += (err * N) ** 2  # variance of summed estimator
    pooled = allw / (N * len(tot))
    pooled_err = math.sqrt(allww) / (N * len(tot))
    ess_pooled = sum(p[4] for p in per)
    lines.append("")
    lines.append(f"combined (4 runs): {pooled:.3e} +- {pooled_err:.1e}")
    lines.append(f"pooled ESS: {ess_pooled:.3e}")
    for sig in (0.5, 1.0):
        vals = [p[2] for p in per if p[1] == sig]
        lines.append(f"width sigma={sig}: " + ", ".join(f"{v:.3e}" for v in vals))
    lines.append("reference (test_b_exact primary): 2.93e-17; "
                 "replicated: 2.83e-17 +- 0.02e-17")
    text = "\n".join(lines) + "\n"
    with open("importance_sampling.txt", "w") as f:
        f.write(text)
    print(text)

if __name__ == "__main__":
    main()
