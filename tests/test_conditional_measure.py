import math, sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "analysis"))
from conditional_joint import koide_frequency
from inputs import load

def test_koide_frequency_scale():
    d = load()
    rng = np.random.default_rng(7)
    PK, se = koide_frequency(50000, d["koide_window"], d["lep_lo"], d["lep_hi"], rng)
    # Marcus: P_K ~ 5.17e-6 (order-of-magnitude + 10% window)
    assert 4.5e-6 < PK < 6.0e-6, PK
