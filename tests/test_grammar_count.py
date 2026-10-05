import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "analysis"))
from census import class_counts, residuals, tolerance_ladder
from inputs import load

def test_counts():
    cc = class_counts()
    assert cc["P"] == 324 and cc["T"] == 252 and cc["A"] == 728
    assert cc["S"] == 426 and cc["D"] == 426 and cc["K"] == 1584
    assert sum(cc.values()) == 3740

def test_ladder():
    d = load()
    lad = tolerance_ladder(residuals(d["masses"]), [0.1, 0.25, 0.5, 1.0])
    assert lad == [0, 1, 3, 6], lad
