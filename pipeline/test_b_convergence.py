#!/usr/bin/env python3
"""Convergence documentation for the Test B Rao-Blackwellized estimator.

Referee round-5 request: effective sample size of the weighted calculation,
independent-seed replication, and dominance assessment (does a small number
of unusually favorable anchor configurations dominate the estimate?).

Baseline null only, five independent seeds x 5e6 universes each. Per seed:
joint estimate, ESS = (sum w)^2 / sum w^2, contributor count, top-1 and
top-1% weight shares. Combined estimate with seed-spread error.

Frozen seeds. Output: test_b_convergence.txt.
"""
import math, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from test_b_exact import koide_surface_triples, joint_components, Q_LO, Q_HI

HERE = os.path.dirname(os.path.abspath(__file__))
N = 5_000_000
SEEDS = [20263001, 20263002, 20263003, 20263004, 20263005]

def one_seed(seed):
    rng = np.random.default_rng(seed)
    lep = koide_surface_triples(N, rng)
    mu_s, p2e = lep.sum(axis=1), 2 * lep[:, 0]
    d = np.exp(rng.uniform(math.log(Q_LO), math.log(Q_HI), N))
    w, _ = joint_components(mu_s, p2e, d, math.log(Q_LO), math.log(Q_HI))
    w = np.asarray(w)
    est = w.mean()
    ess = w.sum() ** 2 / (w ** 2).sum()
    contrib = int((w > 0).sum())
    srt = np.sort(w)[::-1]
    nz = srt[srt > 0]
    top1 = srt[0] / w.sum()
    top10 = srt[:10].sum() / w.sum()
    k = max(1, len(nz) // 100)
    top1pct = nz[:k].sum() / w.sum()
    maxmed = float(srt[0] / np.median(nz))
    return est, ess, contrib, float(top1), float(top10), float(top1pct), maxmed

def main():
    lines = [f"Test B estimator convergence, baseline null, {N} universes/seed", ""]
    ests = []
    for s in SEEDS:
        est, ess, contrib, top1, top10, top1pct, maxmed = one_seed(s)
        ests.append(est)
        lines.append(f"seed {s}: joint5 {est:.4e}; ESS {ess:.0f} "
                     f"({contrib} contributors); top weight {top1:.2e}; "
                     f"top-10 share {top10:.3f}; top-1%-of-contributors share "
                     f"{top1pct:.3f}; max/median nonzero weight {maxmed:.1f}")
        print(lines[-1], flush=True)
    ests = np.array(ests)
    comb = ests.mean()
    spread = ests.std(ddof=1) / math.sqrt(len(ests))
    lines += ["",
              f"combined (25e6 universes): {comb:.4e} +- {spread:.1e} "
              f"(seed spread); per-seed range {ests.min():.3e}..{ests.max():.3e}",
              f"relative seed spread of the mean: {100*spread/comb:.1f}%",
              "",
              "reference: v7 manuscript value (2.93+-0.05)e-17 from a single",
              "5e6-universe run (seed 20262001)."]
    with open(os.path.join(HERE, "test_b_convergence.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines[-5:]))

if __name__ == "__main__":
    main()
