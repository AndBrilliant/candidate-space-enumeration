"""Load the frozen input file. Everything downstream imports from here;
nothing hard-codes a mass or threshold (Marcus §IX)."""
import math, os
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
_YAML = os.path.join(HERE, "..", "inputs", "pdg2026.yaml")

def load(path=None):
    with open(path or _YAML) as f:
        cfg = yaml.safe_load(f)
    sp = cfg["spectrum"]
    R = sp["light_quark_scale_factor"]
    raw = sp["masses_mev"]
    # light quarks u,d,s (idx 3,4,5) evaluated at mu* via R_m
    masses = raw[:3] + [x * R for x in raw[3:6]] + raw[6:]
    ak = math.sqrt(1.5) - 1.0
    mustar = sum(masses[:3])
    G = math.sqrt(1.5) * mustar
    ng = cfg["null_generator"]
    return dict(
        cfg=cfg, names=sp["names"], masses=masses, R_m=R,
        alpha_K=ak, mustar=mustar, G=G,
        lep_lo=ng["leptons"]["lo_mev"], lep_hi=ng["leptons"]["hi_mev"],
        qu_lo=raw[3]*R, qu_hi=raw[8],  # m_u(mu*) .. m_t
        koide_window=ng["koide_window"],
        mc=cfg["mc"],
    )

if __name__ == "__main__":
    d = load()
    print("masses:", d["masses"])
    print(f"mustar={d['mustar']:.8f}  G={d['G']:.6f}  alpha_K={d['alpha_K']:.8f}")
