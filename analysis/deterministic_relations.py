"""Deterministic heavy-sector + lepton relations (Marcus §I reconstruction).

Residual convention is V/T - 1 everywhere (Marcus §VIII)."""
import math
from inputs import load

def koide_Q(a, b, c):
    return (a + b + c) / (math.sqrt(a) + math.sqrt(b) + math.sqrt(c)) ** 2

def foot_angle_deg(a, b, c):
    """Foot cone angle: angle between the sqrt-mass vector and the
    democratic direction (1,1,1)/sqrt(3).  Q = 1/(3 cos^2 alpha), so
    Q = 2/3 iff alpha = 45 deg exactly."""
    Q = koide_Q(a, b, c)
    return math.degrees(math.acos(1.0 / math.sqrt(3.0 * Q)))

def main():
    d = load()
    m = d["masses"]
    ak, ms, G = d["alpha_K"], d["mustar"], d["G"]
    e, mu, tau, u, dw, s, c, b, t = m

    mc_pred = 3 * ak * ms
    mb_pred = ms / (2 * ak)
    r_c = c / mc_pred - 1
    r_b = b / mb_pred - 1
    r_cb = c * b / G ** 2 - 1
    Q = koide_Q(e, mu, tau)
    ang = foot_angle_deg(e, mu, tau)

    # slot thresholds (derived, Marcus §IX)
    P2e = 2 * e
    slots = {
        1: abs(s / (ak ** 2 * ms) - 1),
        2: abs(ms * dw / s ** 2 - 1),
        3: abs(P2e * dw / u ** 2 - 1),
        4: abs(c / mc_pred - 1),
        5: abs(b / mb_pred - 1),
        6: abs(r_cb),
    }
    print(f"mustar            = {ms:.8f} MeV")
    print(f"G                 = {G:.6f} MeV")
    print(f"m_c_pred          = {mc_pred:.6f} MeV   r_c = {r_c*100:+.6f}%")
    print(f"m_b_pred          = {mb_pred:.6f} MeV   r_b = {r_b*100:+.6f}%")
    print(f"m_c m_b / G^2 - 1 = {r_cb*100:+.6f}%")
    print(f"Q_leptons         = {Q:.10f}  (2/3 - Q = {2/3 - Q:.3e})")
    print(f"Foot angle        = {ang:.7f} deg")
    print("slot thresholds (%):",
          {k: round(v * 100, 4) for k, v in slots.items()})
    return slots

if __name__ == "__main__":
    main()
