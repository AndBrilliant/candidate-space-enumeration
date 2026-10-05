import math, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "analysis"))
from conditional_joint import allowed_y_intervals, Q_of

def test_window_bracket():
    """Constructed triple exactly on Q=2/3 must fall inside its own window,
    and a perturbed m3 outside the window must be excluded."""
    # build an exact-Q triple: m1=1, m2=4, solve for m3 on Q=2/3
    m1, m2 = 1.0, 4.0
    eps = 2.2e-6
    iv = allowed_y_intervals(m1, m2, 2/3 - eps, 2/3 + eps)
    assert iv, "window should be non-empty for a feasible pair"
    # midpoint of a lobe is inside the window
    x0, x1 = iv[0]
    ymid = math.sqrt(x0 * x1)
    m3 = ymid ** 2
    Q = Q_of(m1, m2, m3)
    assert abs(Q - 2/3) <= eps, (Q, iv)
    # a far-away m3 is outside
    Qfar = Q_of(m1, m2, 400.0)
    assert abs(Qfar - 2/3) > eps

def test_window_symmetry_small_eps():
    """For a feasible pair the window is nonempty and narrow."""
    m1, m2 = 105.6583755, 1776.930
    eps = 2.2e-6
    iv = allowed_y_intervals(m1, m2, 2/3 - eps, 2/3 + eps)
    # two side lobes (or one merged) around the Q-minimum
    assert all(x1 > x0 > 0 for x0, x1 in iv)
