import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "analysis"))
from whole_grammar_mc import run

def test_frozen_seed_table_iv():
    rows, beat = run(verbose=False)
    exp = {0.1: (0.4416, 0.2745, 0.0961, 16), 0.25: (1.1308, 0.5512, 0.2781, 24),
           0.5: (2.2844, 0.7847, 0.5462, 31), 1.0: (4.5919, 0.9474, 0.8387, 47)}
    for tp, mn, p1, p2, mx in rows:
        em, e1, e2, ex = exp[tp]
        assert abs(mn - em) < 1e-3 and abs(p1 - e1) < 1e-3 and abs(p2 - e2) < 1e-3 and mx == ex, (tp, mn, p1, p2, mx)
    assert abs(beat - 0.4275) < 1e-3, beat
