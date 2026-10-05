"""Light-quark running factor R_m (2 GeV -> mu*), restored (Marcus §XI t8).

2-loop MSbar with standard conventions:
  da/dt  = -b0 a^2 - b1 a^3,   a = alpha_s/pi,   t = ln(mu^2)
  d ln m / dt = -(g0 a + g1 a^2)
with
  b0 = (11 - 2 nf/3) / 4,  b1 = (102 - 38 nf/3) / 16,
  g0 = 1,                  g1 = (202/3 - 20 nf/9) / 16.
PDG alpha_s(mZ) = 0.1179.  nf = 5 above m_b, 4 for m_b..m_c region.  The
mu* = 1.883 GeV step sits in the nf=4 regime (mu* > m_c = 1.273 GeV)."""
import math

def b0(nf): return (11 - 2 * nf / 3) / 4.0
def b1(nf): return (102 - 38 * nf / 3) / 16.0
def g1(nf): return (202.0 / 3 - 20 * nf / 9) / 16.0

def _rk(a, lnm, t0, t1, nf, steps=4000):
    B0, B1, G1 = b0(nf), b1(nf), g1(nf)
    h = (t1 - t0) / steps
    t = t0
    def fa(a): return -B0 * a * a - B1 * a ** 3
    def fm(a): return -(a + G1 * a * a)
    for _ in range(steps):
        k1a, k1m = fa(a), fm(a)
        k2a, k2m = fa(a + 0.5 * h * k1a), fm(a + 0.5 * h * k1a)
        k3a, k3m = fa(a + 0.5 * h * k2a), fm(a + 0.5 * h * k2a)
        k4a, k4m = fa(a + h * k3a), fm(a + h * k3a)
        a += h * (k1a + 2 * k2a + 2 * k3a + k4a) / 6
        lnm += h * (k1m + 2 * k2m + 2 * k3m + k4m) / 6
        t += h
    return a, lnm

def alpha_s_over_pi(mu, nf, a_mz=0.1179 / math.pi, mz=91.1876):
    a, _ = _rk(a_mz, 0.0, math.log(mz ** 2), math.log(mu ** 2), nf)
    return a

def R_m(mu_star_mev=1883.09937445, mb=4186.0):
    mu0, mu1 = 2000.0, mu_star_mev
    # alpha_s/pi at m_b with nf=5, then continue with nf=4
    a_mb = alpha_s_over_pi(mb, nf=5)
    a_mu0, _ = _rk(a_mb, 0.0, math.log(mb ** 2), math.log(mu0 ** 2), 4)
    _, lnm = _rk(a_mu0, 0.0, math.log(mu0 ** 2), math.log(mu1 ** 2), 4)
    return math.exp(lnm)

if __name__ == "__main__":
    r = R_m()
    print(f"alpha_s(2 GeV) = {alpha_s_over_pi(2000.0, 4) * math.pi:.4f}")
    print(f"2-loop running m(2 GeV) -> m(mu*): R_m = {r:.5f}  (frozen: 1.01750)")
