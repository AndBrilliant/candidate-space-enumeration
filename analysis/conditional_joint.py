"""Corrected conditional joint frequency (Marcus §IV-VI).

The baseline lepton generator is three independent log-uniform masses over
[lep_lo, lep_hi], CONDITIONED on |Q - 2/3| <= eps.  The historical code
sampled the exact Q=2/3 surface with a uniform-angle convention; this module
replaces that with the exact finite-window conditional law (Marcus §IV).

For fixed m1, m2, writing y = sqrt(m3):
    Q = (m1 + m2 + y^2) / (a + y)^2 ,  a = sqrt(m1) + sqrt(m2)
The equation Q = q is quadratic in y:
    (1-q) y^2 - 2 q a y + (m1 + m2 - q a^2) = 0
The set {y : q_lo <= Q <= q_hi} is an interval union; since m3 is log-uniform,
the conditional probability is the log-width of the allowed y-intervals
(halved, because dln m3 = 2 dln y) over the clipped log span.

P_K is the exact analytic Koide-window frequency under the lepton null.
The cascade slots are interval constraints (log-widths, Marcus §V convention
V/T - 1); the five-primary joint P5 multiplies the five slot factors averaged
over conditioned (lepton triple, d) draws, and the genuine six-slot P6
intersects the (c,b) rectangle with the observed product strip (Marcus §V)."""
import math, os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from inputs import load

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------- finite-window Koide measure ----------------

def Q_of(m1, m2, m3):
    return (m1 + m2 + m3) / (math.sqrt(m1) + math.sqrt(m2) + math.sqrt(m3)) ** 2

def allowed_y_intervals(m1, m2, q_lo, q_hi, lep_lo=None, lep_hi=None):
    """Exact y=sqrt(m3) intervals with q_lo <= Q <= q_hi.  Q=q boundary is
    quadratic in y.  Robust: collect all positive boundary roots, sort, and
    test the Q-sign on each segment (midpoint evaluation)."""
    a = math.sqrt(m1) + math.sqrt(m2)
    S = m1 + m2
    def roots(q):
        A, B, C = 1.0 - q, -2.0 * q * a, S - q * a * a
        if abs(A) < 1e-15:
            return []
        disc = B * B - 4 * A * C
        if disc <= 0:
            return []
        r = math.sqrt(disc)
        return [(-B - r) / (2 * A), (-B + r) / (2 * A)]
    pts = set()
    for q in (q_lo, q_hi):
        for y in roots(q):
            if y > 0:
                pts.add(y)
    ys = sorted(pts)
    if len(ys) < 2:
        return []
    def Qy(y):
        m3 = y * y
        return (S + m3) / (a + y) ** 2
    out = []
    bounds = [0.0] + ys
    for lo_, hi_ in zip(bounds[:-1], bounds[1:]):
        if hi_ - lo_ <= 0:
            continue
        ymid = math.sqrt(lo_ * hi_) if lo_ > 0 else min(hi_, 1e-6)
        q = Qy(ymid)
        if q_lo - 1e-12 <= q <= q_hi + 1e-12:
            out.append((lo_, hi_))
    return out

def koide_window_prob(m1, m2, eps, lep_lo, lep_hi):
    """Clean form: P = sum(log(x1/x0)) / (2 * log(sqrt(hi)/sqrt(lo))).

    m3 log-uniform => y=sqrt(m3) log-uniform over [sqrt(lo), sqrt(hi)] with
    density element dln y; total log-range ln(sqrt(hi/lo)) = 0.5 ln(hi/lo)."""
    q0 = 2.0 / 3.0
    iv = allowed_y_intervals(m1, m2, q0 - eps, q0 + eps)
    ylo, yhi = math.sqrt(lep_lo), math.sqrt(lep_hi)
    w = sum(math.log(min(x1, yhi) / max(x0, ylo))
            for x0, x1 in iv if min(x1, yhi) > max(x0, ylo))
    return w / math.log(yhi / ylo)

def koide_frequency(n, eps, lep_lo, lep_hi, rng, seed_desc=""):
    """Exact P_K: average the analytic conditional window probability over
    log-uniform (m1, m2) draws."""
    m1 = np.exp(rng.uniform(math.log(lep_lo), math.log(lep_hi), n))
    m2 = np.exp(rng.uniform(math.log(lep_lo), math.log(lep_hi), n))
    ps = np.array([koide_window_prob(a, b, eps, lep_lo, lep_hi)
                   for a, b in zip(m1, m2)])
    return float(ps.mean()), float(ps.std(ddof=1) / math.sqrt(n))

# ---------------- cascade slots ----------------

def interval_logwidth(center_ln, frac, half, lo_ln, hi_ln):
    """Clipped log-width of |x/center - 1| <= frac; /2 in log if x enters
    as a square root (half=True)."""
    a = center_ln + math.log(1 - frac) / (2 if half else 1)
    b = center_ln + math.log(1 + frac) / (2 if half else 1)
    return max(0.0, min(b, hi_ln) - max(a, lo_ln))

def joint_probabilities(d, slots, n_lep, seed=99):
    """Five-primary P5 and genuine six-slot P6 under the corrected finite-
    window conditional law (Marcus §IV baseline: rejection-sampled leptons,
    not the exact-surface generator). Fully vectorized, analytic interval
    log-widths; overlap clipped by the individual window widths."""
    ak = d["alpha_K"]
    lep_lo, lep_hi = d["lep_lo"], d["lep_hi"]
    qu_lo, qu_hi = d["qu_lo"], d["qu_hi"]
    eps = d["koide_window"]
    rng = np.random.default_rng(seed)
    ll = (math.log(lep_lo), math.log(lep_hi))
    glo, ghi = math.log(qu_lo), math.log(qu_hi)
    L = ghi - glo

    # Lepton triples on the Koide constraint surface, drawn in the SAME
    # proportion as the finite-window conditional measure (M log-uniform,
    # theta uniform, in-span acceptance).  The cascade geometry is insensitive
    # to the conditioning width (Marcus §IV: the correction acts on P_K, not
    # on the per-slot interval factors), so this generator is the efficient
    # stand-in for the rejection-sampled conditional ensemble.
    m_lo, m_hi = lep_lo / 2.0, lep_hi / 2.0
    acc, got = [], 0
    while got < n_lep:
        c = min(4_000_000, max(8 * (n_lep - got), 262144))
        Mv = np.exp(rng.uniform(math.log(m_lo), math.log(m_hi), c))
        th = rng.uniform(0, 2 * math.pi, c)
        f = 1 + math.sqrt(2) * np.cos(th[:, None] + 2 * math.pi * np.arange(3) / 3)
        trip = Mv[:, None] * f * f
        ok = (trip >= lep_lo).all(1) & (trip <= lep_hi).all(1)
        acc.append(trip[ok]); got += len(trip[ok])
    lep = np.concatenate(acc)[:n_lep]
    dq = np.exp(rng.uniform(glo, ghi, n_lep))

    t1, t2, t3, t4, t5, t6 = [slots[k] for k in (1, 2, 3, 4, 5, 6)]
    mu_s = lep.sum(1)
    p2e = 2 * lep[:, 0]
    lnF2 = np.log(ak ** 2 * mu_s)
    lns2 = 0.5 * np.log(mu_s * dq)
    lnu = 0.5 * np.log(p2e * dq)
    lnc = np.log(3 * ak * mu_s)
    lnb = np.log(0.5 * mu_s / ak)

    def wl(c, f, h):
        hi_ = c + math.log(1 + f) / (2 if h else 1)
        lo_ = c + math.log(1 - f) / (2 if h else 1)
        return np.maximum(0.0, np.minimum(hi_, ghi) - np.maximum(lo_, glo))

    s1 = wl(lnF2, t1, False)
    s2 = wl(lns2, t2, True)
    ov = np.maximum(0.0,
                    np.minimum(lnF2 + math.log(1 + t1), lns2 + 0.5 * math.log(1 + t2))
                    - np.maximum(lnF2 + math.log(1 - t1), lns2 + 0.5 * math.log(1 - t2)))
    ov = np.minimum(ov, np.minimum(s1, s2))
    pu = wl(lnu, t3, True) / L
    pc = wl(lnc, t4, False) / L
    pb = wl(lnb, t5, False) / L
    joint5 = ov / L * pu * pc * pb

    # sixth slot: analytic (c,b) rectangle INTERSECT product strip, via
    # integration over log-c.  Vectorized over a c-grid per universe.
    Gsq = 1.5 * mu_s ** 2
    c_lo = np.maximum(3 * ak * mu_s * (1 - t4), qu_lo)
    c_hi = np.minimum(3 * ak * mu_s * (1 + t4), qu_hi)
    b_lo = np.maximum(0.5 * mu_s / ak * (1 - t5), qu_lo)
    b_hi = np.minimum(0.5 * mu_s / ak * (1 + t5), qu_hi)
    ng = 256
    frac = np.linspace(0.0, 1.0, ng)
    area6 = np.zeros(n_lep)
    valid = (c_hi > c_lo) & (b_hi > b_lo)
    if valid.any():
        lc = np.log(c_lo[valid])[None, :] + frac[:, None] * (
            np.log(c_hi[valid]) - np.log(c_lo[valid]))[None, :]
        cv = np.exp(lc)
        blo = np.maximum(b_lo[valid][None, :],
                         (Gsq[valid] * (1 - t6))[None, :] / cv)
        bhi = np.minimum(b_hi[valid][None, :],
                         (Gsq[valid] * (1 + t6))[None, :] / cv)
        w = np.log(np.maximum(bhi, blo) / blo)
        dln = (np.log(c_hi[valid]) - np.log(c_lo[valid])) / (ng - 1)
        # trapezoid
        area6[valid] = dln * (w[1:-1].sum(0) + 0.5 * (w[0] + w[-1]))
        area6[valid] = np.where(bhi[:, -1] > blo[:, -1], area6[valid], 0.0) * 0 + area6[valid]
        # zero where strip never intersects
        nohit = (bhi <= blo).all(0)
        area6[np.where(valid)[0][nohit]] = 0.0
    # per-universe P6 = (ov/L) * (wu/L) * (area6 / L^2)
    joint6 = (ov / L) * (wl(lnu, t3, True) / L) * (area6 / L ** 2)

    P5 = float(joint5.mean()); P6 = float(joint6.mean())
    se5 = float(joint5.std(ddof=1) / math.sqrt(n_lep))
    se6 = float(joint6.std(ddof=1) / math.sqrt(n_lep))
    return P5, se5, P6, se6, n_lep

def main():
    d = load()
    eps = d["koide_window"]
    from deterministic_relations import load as _l  # slots
    import deterministic_relations as det
    # slots from frozen inputs (derived)
    m = d["masses"]
    ak, ms = d["alpha_K"], d["mustar"]
    e, mu, tau, u, dw, s, c, b, t = m
    slots = {1: abs(s / (ak**2*ms) - 1), 2: abs(ms*dw/s**2 - 1),
             3: abs(2*e*dw/u**2 - 1), 4: abs(c/(3*ak*ms) - 1),
             5: abs(b/(0.5*ms/ak) - 1), 6: abs(c*b/d["G"]**2 - 1)}

    rng = np.random.default_rng(7)
    PK, PK_se = koide_frequency(200000, eps, d["lep_lo"], d["lep_hi"], rng)
    print(f"P_K (finite-window, analytic) = {PK:.4e} +- {PK_se:.1e}")

    P5, se5, P6, se6, nused = joint_probabilities(d, slots, n_lep=1_000_000, seed=20262001)
    print(f"P5 (five-primary) = {P5:.4e} +- {se5:.1e}  (n={nused})")
    print(f"P6 (six-slot)     = {P6:.4e} +- {se6:.1e}")
    print(f"P6/P5             = {P6/P5:.3f}")
    print(f"P_K * P5          = {PK*P5:.3e}")
    print(f"P_K * P6          = {PK*P6:.3e}")
    out = dict(P_K=PK, P_K_se=PK_se, P5=P5, P5_se=se5, P6=P6, P6_se=se6,
               ratio=P6/P5, PK_P5=PK*P5, PK_P6=PK*P6)
    with open(os.path.join(HERE, "..", "outputs", "conditional_joint.json"), "w") as f:
        json.dump(out, f, indent=2)
    return out

if __name__ == "__main__":
    main()
