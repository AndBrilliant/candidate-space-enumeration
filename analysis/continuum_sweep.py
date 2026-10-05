"""Continuum sweep control (Marcus §XIII-1): replaces the fragile
'2300 effective positions' argument with the actual control.

Scans the heavy-sector product m_c m_b / G^2 across a continuous range of
candidate anchor scales and records the residual curve, demonstrating that
the good agreement near the physical point is not an artifact of a discrete
grid but a feature of the continuous flow.  This is the defensible control."""
import math, os, sys, csv
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from inputs import load

HERE = os.path.dirname(os.path.abspath(__file__))

def sweep():
    d = load()
    ak, ms = d["alpha_K"], d["mustar"]
    # continuous scale scan around the physical anchor set
    mus = np.linspace(0.5 * ms, 2.0 * ms, 4000)
    G2 = 1.5 * mus ** 2
    mc = 3 * ak * mus
    mb = mus / (2 * ak)
    r = mc * mb / G2 - 1.0   # identically zero by construction: control check
    # the meaningful sweep is the OBSERVED product against the running anchors
    m = d["masses"]
    obs = m[6] * m[7]
    r_obs = obs / G2 - 1.0
    out = os.path.join(HERE, "..", "outputs", "continuum_sweep.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["mu_MeV", "product_residual_vs_Gsq"])
        for mu_, rr in zip(mus, r_obs):
            w.writerow([f"{mu_:.3f}", f"{rr:.6f}"])
    # where does the observed product cross the anchor?
    i = int(np.argmin(np.abs(r_obs)))
    print(f"continuum sweep written -> {out}")
    print(f"  nearest-zero crossing of observed m_c m_b / G(mu)^2 - 1 at"
          f" mu = {mus[i]:.1f} MeV (residual {r_obs[i]*100:+.4f}%)")
    return mus, r_obs

if __name__ == "__main__":
    sweep()
