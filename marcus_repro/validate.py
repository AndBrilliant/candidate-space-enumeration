#!/usr/bin/env python3
"""Fast invariant tests for the independent conditional replay."""
import math
import numpy as np
from inputs import *
from conditional_replay import conditioned_leptons, downstream_joint

def qvalue(a,b,c):
    x=np.stack([a,b,c],axis=1)
    return x.sum(axis=1)/np.sqrt(x).sum(axis=1)**2

rng=np.random.default_rng(99173)
a,b,c,w=conditioned_leptons(250_000,rng)
mask=w>0
q=qvalue(a[mask],b[mask],c[mask])
assert np.all(np.abs(q-2/3) <= KOIDE_EPS*(1+1e-10))
assert np.all((c[mask]>=LEP_LO)&(c[mask]<=LEP_HI))
assert np.all(w>=0)
# The exact two-dimensional integration is much lower-variance than a raw
# three-dimensional rare-event count. Broad guardrail, not a fitted target.
p=float(w.mean())
assert 4.8e-6 < p < 5.5e-6, p

d=np.exp(rng.uniform(math.log(Q_LO),math.log(Q_HI),len(w)))
j5,j6,t=downstream_joint(a,b,c,d)
assert np.all(j5>=0) and np.all(j6>=0)
assert np.all(j6 <= j5*(1+1e-12))
ratio=float(np.sum(w*j6)/np.sum(w*j5))
assert 0.60 < ratio < 0.76, ratio

print(f"conditional invariants: PASS; P_Koide={p:.7e}; joint6/joint5={ratio:.4f}")
