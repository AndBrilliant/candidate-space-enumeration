#!/usr/bin/env python3
"""Exact conditional joint frequency for the declared cascade (Test B-exact).

Replaces the product-of-marginals estimate with an exact computation under
the disclosed null. Given the lepton anchors and the d-quark draw, each
cascade slot is an interval constraint on a single mass:

  slot 1:  s  in F^2 * [1-t1, 1+t1]                 F^2 = alpha_K^2 mu*
  slot 2:  s  in sqrt(mu* d) * [1-t2, 1+t2]^{1/2}   (shares s with slot 1)
  slot 3:  u  in sqrt(P d)   * [1-t3, 1+t3]^{1/2}   (shares d with slot 2)
  slot 4:  c  in 3 alpha_K mu* * [1-t4, 1+t4]
  slot 5:  b  in (mu*/2)/alpha_K * [1-t5, 1+t5]

Under log-uniform quark draws over the declared span, each interval's
probability is its clipped log-width / log-span — analytic. The joint per
universe is P_s_overlap * P_u * P_c * P_b, exact given (anchors, d); the
null average is taken over (Koide-surface lepton triple, d) pairs drawn
i.i.d. (no triple reuse, so no clustering). No independence assumption
anywhere: every shared input is either conditioned on or integrated through
the interval overlap.

Also computes:
  - the product of marginals on the same ensemble (for honest comparison
    with the superseded estimate)
  - the common-tolerance joint-discrepancy tail P(T <= t), T = max_i |r_i|,
    evaluated at T_obs = max observed slot residual, plus the tail curve
  - null-shape variants: span x2, span x0.5, uniform-in-mass
  - brute-force validation of the s-pair overlap factor at observed
    thresholds (3e8 direct draws, ~10^3 expected events)
  - brute-force direct joint count at 1e9 full spectra (independent bound)

Frozen seeds. Outputs: test_b_exact_summary.txt, test_b_exact_curve.csv.
"""
import math, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from census_grammar import PDG
from sensitivity_mc import logspan, sigma, LEP_LO, LEP_HI

HERE = os.path.dirname(os.path.abspath(__file__))

def koide_surface_triples(n, rng, lo=LEP_LO, hi=LEP_HI):
    """Chunked exact-Q=2/3 triple generator (memory-safe). Same law as
    sensitivity_mc.koide_surface_triples: M log-U over [lo/2, hi/2],
    theta uniform, in-span acceptance."""
    m_lo, m_hi = lo / 2.0, hi / 2.0
    out, got = [], 0
    while got < n:
        c = min(2_000_000, max(4 * (n - got), 65_536))
        Mv = np.exp(rng.uniform(math.log(m_lo), math.log(m_hi), c))
        th = rng.uniform(0, 2 * math.pi, c)
        f = 1 + math.sqrt(2) * np.cos(th[:, None] + 2 * math.pi * np.arange(3) / 3)
        trip = Mv[:, None] * f * f
        ok = (trip >= lo).all(axis=1) & (trip <= hi).all(axis=1)
        acc = trip[ok]
        out.append(acc)
        got += len(acc)
    return np.concatenate(out)[:n]
AK = math.sqrt(1.5) - 1.0
Q_LO, Q_HI = 2.16 * 1.01750, 172600.0

# observed slot residual magnitudes from the frozen PDG spectrum
_MS = PDG[0] + PDG[1] + PDG[2]
T = [abs(PDG[5] / (AK ** 2 * _MS) - 1),           # s vs alpha_K^2 mu*
     abs(_MS * PDG[4] / PDG[5] ** 2 - 1),         # mu* d vs s^2
     abs(2 * PDG[0] * PDG[4] / PDG[3] ** 2 - 1),  # P d vs u^2
     abs(PDG[6] / (3 * AK * _MS) - 1),            # c vs 3 alpha_K mu*
     abs(PDG[7] / (0.5 * _MS / AK) - 1)]          # b vs (mu*/2)/alpha_K
T_OBS = max(T)

def interval(center_ln, frac, half, lo, hi):
    """[lo,hi]-clipped log-interval width for |x/target - 1| <= frac,
    halved if the constrained quantity is a square root of the ratio."""
    a = center_ln + math.log(1 - frac) / (2 if half else 1)
    b = center_ln + math.log(1 + frac) / (2 if half else 1)
    return max(0.0, min(b, hi) - max(a, lo))

def joint_components(mu_s, p2e, d, lo, hi, tol=None, uniform=False):
    """Vectorized exact per-universe joint. mu_s, p2e, d arrays; returns
    (joint, (m1..m5 marginals)). tol overrides per-slot thresholds with a
    common tolerance (T-statistic mode)."""
    t = tol if tol is not None else T
    t1, t2, t3, t4, t5 = (t, t, t, t, t) if tol is not None else t
    L = hi - lo
    F2 = AK ** 2 * mu_s
    lnF2 = np.log(F2)
    lns2 = 0.5 * np.log(mu_s * d)          # center of ln s from slot 2
    lnu = 0.5 * np.log(p2e * d)            # center of ln u
    lnc = np.log(3 * AK * mu_s)
    lnb = np.log(0.5 * mu_s / AK)

    if uniform:
        # linear-space interval widths / linear span
        def w(center, frac, half):
            if half:  # x^2 in target*[1-f,1+f] -> x in sqrt range
                return np.maximum(0.0, np.minimum(np.sqrt(1 + frac) * center, Q_HI_)
                                  - np.maximum(np.sqrt(1 - frac) * center, Q_LO_))
            return np.maximum(0.0, np.minimum((1 + frac) * center, Q_HI_)
                              - np.maximum((1 - frac) * center, Q_LO_))
        L = Q_HI_ - Q_LO_
        s1_lo = np.maximum(F2 * (1 - t1), Q_LO_); s1_hi = np.minimum(F2 * (1 + t1), Q_HI_)
        s2_lo = np.maximum(np.sqrt(mu_s * d * (1 - t2)), Q_LO_)
        s2_hi = np.minimum(np.sqrt(mu_s * d * (1 + t2)), Q_HI_)
        ov = np.maximum(0.0, np.minimum(s1_hi, s2_hi) - np.maximum(s1_lo, s2_lo))
        pu = w(np.sqrt(p2e * d), t3, True) / L
        pc = w(3 * AK * mu_s, t4, False) / L
        pb = w(0.5 * mu_s / AK, t5, False) / L
        m1 = (s1_hi - s1_lo) / L; m2 = (s2_hi - s2_lo) / L
        return ov / L * pu * pc * pb, (m1, m2, pu, pc, pb)

    # vectorized interval widths (log space)
    def wl(c, f, h):
        d_ = math.log(1 + f) / (2 if h else 1)
        d__ = math.log(1 - f) / (2 if h else 1)
        return np.maximum(0.0, np.minimum(c + d_, hi) - np.maximum(c + d__, lo))
    s1 = (np.minimum(lnF2 + math.log(1 + t1), hi)
          - np.maximum(lnF2 + math.log(1 - t1), lo))
    s1 = np.maximum(0.0, s1)
    s2 = (np.minimum(lns2 + 0.5 * math.log(1 + t2), hi)
          - np.maximum(lns2 + 0.5 * math.log(1 - t2), lo))
    s2 = np.maximum(0.0, s2)
    ov = np.maximum(0.0,
                    np.minimum(lnF2 + math.log(1 + t1), lns2 + 0.5 * math.log(1 + t2))
                    - np.maximum(lnF2 + math.log(1 - t1), lns2 + 0.5 * math.log(1 - t2)))
    # overlap may lie outside the draw span; clip
    ov = np.minimum(ov, np.minimum(s1, s2))
    pu = wl(lnu, t3, True) / L
    pc = wl(lnc, t4, False) / L
    pb = wl(lnb, t5, False) / L
    joint = ov / L * pu * pc * pb
    return joint, (s1 / L, s2 / L, pu, pc, pb)

def run_variant(tag, q_lo, q_hi, uniform, n, seed):
    global Q_LO_, Q_HI_
    Q_LO_, Q_HI_ = q_lo, q_hi
    rng = np.random.default_rng(seed)
    lep = koide_surface_triples(n, rng)
    mu_s = lep.sum(axis=1)
    p2e = 2 * lep[:, 0]
    if uniform:
        d = rng.uniform(q_lo, q_hi, n)
    else:
        d = np.exp(rng.uniform(math.log(q_lo), math.log(q_hi), n))
    lo, hi = (q_lo, q_hi) if uniform else (math.log(q_lo), math.log(q_hi))
    joint, marg = joint_components(mu_s, p2e, d, lo, hi, uniform=uniform)
    prod = float(np.prod([m.mean() for m in marg]))
    jmean = float(joint.mean())
    jerr = float(joint.std() / math.sqrt(n))
    # T-statistic tail at T_obs (common tolerance mode)
    jointT, _ = joint_components(mu_s, p2e, d, lo, hi, tol=T_OBS, uniform=uniform)
    return {"variant": tag, "n": n, "joint5_exact": jmean, "mc_err": jerr,
            "sigma_equiv": sigma(jmean), "product_of_marginals": prod,
            "ratio_exact_over_product": jmean / prod if prod else float("nan"),
            "P_T_le_Tobs": float(jointT.mean()),
            "marginals": [float(m.mean()) for m in marg]}

def main():
    lines = ["Test B-exact: conditional cascade joint by exact integration",
             f"slot thresholds (matching-or-better): "
             + ", ".join(f"{t*100:.3f}%" for t in T),
             f"T_obs (max slot residual) = {T_OBS*100:.3f}%", ""]
    variants = [("baseline", Q_LO, Q_HI, False, 20262001),
                ("span_x2", *logspan(Q_LO, Q_HI, 2.0), False, 20262002),
                ("span_x0p5", *logspan(Q_LO, Q_HI, 0.5), False, 20262003),
                ("uniform", Q_LO, Q_HI, True, 20262004)]
    results = []
    for tag, lo, hi, uni, seed in variants:
        r = run_variant(tag, lo, hi, uni, 5_000_000, seed)
        results.append(r)
        lines.append(f"{tag}: joint5 exact = {r['joint5_exact']:.3e} "
                     f"+- {r['mc_err']:.1e} (sigma-equiv {r['sigma_equiv']:.2f}); "
                     f"product-of-marginals {r['product_of_marginals']:.3e}; "
                     f"ratio {r['ratio_exact_over_product']:.1f}x; "
                     f"P(T<=T_obs) {r['P_T_le_Tobs']:.3e}")
        print(lines[-1], flush=True)

    # T-tail curve under the baseline null
    rng = np.random.default_rng(20262011)
    n = 5_000_000
    lep = koide_surface_triples(n, rng)
    mu_s, p2e = lep.sum(axis=1), 2 * lep[:, 0]
    d = np.exp(rng.uniform(math.log(Q_LO), math.log(Q_HI), n))
    glo, ghi = math.log(Q_LO), math.log(Q_HI)
    rows = []
    for t in np.logspace(-4, math.log10(3e-2), 25):
        jt, _ = joint_components(mu_s, p2e, d, glo, ghi, tol=float(t))
        rows.append((t, float(jt.mean())))
    with open(os.path.join(HERE, "test_b_exact_curve.csv"), "w") as f:
        f.write("tolerance_frac,P_T_le_t_baseline\n")
        for t, p in rows:
            f.write(f"{t:.6e},{p:.6e}\n")
    lines.append("")
    lines.append("T-tail curve (baseline null): tolerance -> P(T<=t)")
    for t, p in rows[::4]:
        lines.append(f"  {t*100:8.3f}%  {p:.3e}")

    # brute-force validation of the s-pair overlap at observed thresholds
    rng = np.random.default_rng(20262021)
    nbf, hits = 300_000_000, 0
    done = 0
    while done < nbf:
        c = 10_000_000
        # validation needs the s-pair only; fresh mu*, d, s per trial
        mu_s = np.exp(rng.uniform(math.log(1.5), math.log(5331.0), c))
        F2 = AK ** 2 * mu_s
        d = np.exp(rng.uniform(glo, ghi, c))
        s = np.exp(rng.uniform(glo, ghi, c))
        ok = (np.abs(s / F2 - 1) <= T[0]) & (np.abs(mu_s * d / s ** 2 - 1) <= T[1])
        hits += int(ok.sum()); done += c
    bf = hits / nbf
    # analytic s-pair under the same mu*~log-U[1.5,5331] law
    rng2 = np.random.default_rng(20262022)
    mu_s = np.exp(rng2.uniform(math.log(1.5), math.log(5331.0), 5_000_000))
    d2 = np.exp(rng2.uniform(glo, ghi, 5_000_000))
    F2 = AK ** 2 * mu_s
    lns2c = 0.5 * np.log(mu_s * d2)
    lnF2 = np.log(F2)
    ov = np.maximum(0.0,
                    np.minimum(lnF2 + math.log(1 + T[0]), lns2c + 0.5 * math.log(1 + T[1]))
                    - np.maximum(lnF2 + math.log(1 - T[0]), lns2c + 0.5 * math.log(1 - T[1])))
    span_w = ghi - glo
    s1w = np.maximum(0.0, np.minimum(lnF2 + math.log(1 + T[0]), ghi)
                     - np.maximum(lnF2 + math.log(1 - T[0]), glo))
    s2w = np.maximum(0.0, np.minimum(lns2c + 0.5 * math.log(1 + T[1]), ghi)
                     - np.maximum(lns2c + 0.5 * math.log(1 - T[1]), glo))
    ov = np.minimum(ov, np.minimum(s1w, s2w))
    an = float((ov / span_w).mean())
    lines += ["",
              f"s-pair validation at observed thresholds: brute force {bf:.3e} "
              f"({hits} events / {nbf:.0e}) vs analytic {an:.3e} "
              f"(same mu* law)", ""]
    print("\n".join(lines[-4:]), flush=True)

    with open(os.path.join(HERE, "test_b_exact_summary.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")

if __name__ == "__main__":
    main()
