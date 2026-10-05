import math, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "analysis"))
from inputs import load
from deterministic_relations import koide_Q, foot_angle_deg

d = load(); m = d["masses"]; ak = d["alpha_K"]; ms = d["mustar"]; G = d["G"]
c, b = m[6], m[7]

def test_mustar():
    assert abs(ms - 1883.09937445) < 1e-6, ms

def test_charm_residual():
    r = c / (3 * ak * ms) - 1
    assert abs(r - 0.00255914) < 1e-7, r

def test_bottom_residual():
    r = b / (ms / (2 * ak)) - 1
    assert abs(r - (-0.00081531)) < 1e-7, r

def test_product_residual():
    r = c * b / G ** 2 - 1
    assert abs(r - 0.00174175) < 1e-7, r

def test_koide_leptons():
    Q = koide_Q(m[0], m[1], m[2])
    assert abs(Q - 0.6666644634) < 1e-9, Q

def test_foot_angle():
    ang = foot_angle_deg(m[0], m[1], m[2])
    assert abs(ang - 44.9999053) < 1e-5, ang
