"""Variance-reduced cross-check of the corrected conditional law (Marcus §VII).

Importance sampling oversamples the analytically known d-region where the two
strange-sector windows can overlap, with a global proposal component for full
support. Likelihood-ratio weights restore the log-uniform d target.

Five seeds x N per seed; reports across-seed MC standard errors (diagnostic
only, Marcus §VII)."""
import math, os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from inputs import load
from conditional_joint import koide_window_prob

HERE = os.path.dirname(os.path.abspath(__file__))

def surface_leptons(n, rng, lep_lo, lep_hi):
    m_lo, m_hi = lep_lo / 2, lep_hi / 2
    out, got = [], 0
    while got < n:
        c = min(1_000_000, max(4 * (n - got), 65536))
        Mv = np.exp(rng.uniform(math.log(m_lo), math.log(m_hi), c))
        th = rng.uniform(0, 2 * math.pi, c)
        f = 1 + math.sqrt(2) * np.cos(th[:, None] + 2 * math.pi * np.arange(3) / 3)
        tr = Mv[:, None] * f * f
        ok = (tr >= lep_lo).all(1) & (tr <= lep_hi).all(1)
        out.append(tr[ok]); got += len(tr[ok])
    return np.concatenate(out)[:n]

def joint5_array(mu_s, p2e, d, slots, glo, ghi):
    ak = math.sqrt(1.5) - 1.0
    L = ghi - glo
    t1, t2, t3, t4, t5 = slots
    lnF2 = np.log(ak ** 2 * mu_s); lns2 = 0.5 * np.log(mu_s * d)
    lnu = 0.5 * np.log(p2e * d); lnc = np.log(3 * ak * mu_s); lnb = np.log(0.5 * mu_s / ak)
    def wl(c, f, h):
        return np.maximum(0.0, np.minimum(c + math.log(1 + f) / (2 if h else 1), ghi)
                          - np.maximum(c + math.log(1 - f) / (2 if h else 1), glo))
    s1 = wl(lnF2, t1, False); s2 = wl(lns2, t2, True)
    ov = np.maximum(0.0, np.minimum(lnF2 + math.log(1 + t1), lns2 + 0.5 * math.log(1 + t2))
                    - np.maximum(lnF2 + math.log(1 - t1), lns2 + 0.5 * math.log(1 - t2)))
    ov = np.minimum(ov, np.minimum(s1, s2))
    pu = wl(lnu, t3, True) / L
    pc = wl(lnc, t4, False) / L
    pb = wl(lnb, t5, False) / L
    return ov / L * pu * pc * pb

def run_seed(seed, n, d, slots):
    rng = np.random.default_rng(seed)
    glo, ghi = math.log(d["qu_lo"]), math.log(d["qu_hi"])
    L = ghi - glo
    lep = surface_leptons(n, rng, d["lep_lo"], d["lep_hi"])
    mu_s = lep.sum(1); p2e = 2 * lep[:, 0]
    # The s-window overlap support concentrates in a narrow band at the
    # bottom of the quark span (ln d ~ [glo, glo+~1.5]; see diagnostic).
    # Oversample that band + keep a global component for full support.
    band_lo, band_hi = glo, min(glo + 2.0, ghi)
    half = rng.random(n) < 0.5
    lprop = np.where(half, rng.uniform(band_lo, band_hi, n),
                     rng.uniform(glo, ghi, n))
    dq = np.exp(lprop)
    inband = (lprop >= band_lo) & (lprop <= band_hi)
    p_prop = 0.5 * inband / (band_hi - band_lo) + 0.5 / L
    p_targ = 1.0 / L
    w = p_targ / np.maximum(p_prop, 1e-300)
    j5 = joint5_array(mu_s, p2e, dq, slots, glo, ghi)
    est = float((w * j5).mean())
    return est

def main():
    d = load()
    m = d["masses"]; ak = d["alpha_K"]; ms = d["mustar"]
    e, mu, tau, u, dw, s, c, b, t = m
    slots = [abs(s / (ak**2*ms) - 1), abs(ms*dw/s**2 - 1), abs(2*e*dw/u**2 - 1),
             abs(c/(3*ak*ms) - 1), abs(b/(0.5*ms/ak) - 1)]
    seeds = d["mc"]["seeds_importance"]; n = d["mc"]["n_per_seed_importance"]
    ests = []
    for sd in seeds:
        ests.append(run_seed(sd, n, d, slots))
    ests = np.array(ests)
    print(f"importance P5 per seed: {[f'{x:.4e}' for x in ests]}")
    print(f"P5 = {ests.mean():.4e} +- {ests.std(ddof=1)/math.sqrt(len(ests)):.1e}"
          f"  (across-seed MC se, diagnostic only)")
    out = dict(per_seed=ests.tolist(), mean=float(ests.mean()),
               across_seed_se=float(ests.std(ddof=1)/math.sqrt(len(ests))))
    json.dump(out, open(os.path.join(HERE, "..", "outputs", "importance_joint.json"), "w"), indent=2)
    return out

if __name__ == "__main__":
    main()
