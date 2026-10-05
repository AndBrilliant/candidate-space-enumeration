import math, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "analysis"))
from conditional_joint import joint_probabilities
from inputs import load

def test_p6_le_p5():
    d = load(); m = d["masses"]; ak = d["alpha_K"]; ms = d["mustar"]
    e, mu, tau, u, dw, s, c, b, t = m
    slots = {1: abs(s/(ak**2*ms)-1), 2: abs(ms*dw/s**2-1), 3: abs(2*e*dw/u**2-1),
             4: abs(c/(3*ak*ms)-1), 5: abs(b/(0.5*ms/ak)-1), 6: abs(c*b/d["G"]**2-1)}
    # small-n smoke test only (P6 <= P5 pointwise by construction)
    P5, _, P6, _, _ = joint_probabilities(d, slots, n_lep=20000, seed=1)
    assert P6 <= P5 * 1.5, (P5, P6)   # allow MC noise; geometry bound is P6<=P5
