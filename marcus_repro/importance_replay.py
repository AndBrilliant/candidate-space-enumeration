#!/usr/bin/env python3
"""Variance-reduced corrected conditional replay.

Uses the same exact finite-window lepton conditioning as conditional_replay.py.
For m_d, instead of always drawing from the full log-uniform quark span, use a
mixture proposal:
  - local component: uniform in the exact log-d interval where the two strange
    windows can overlap;
  - global component: the original full log-uniform law.

The global component guarantees full support. Exact likelihood-ratio weights
restore the target log-uniform m_d law. No SciPy dependency is required.
"""
import argparse, math
import numpy as np
from inputs import *
from conditional_replay import conditioned_leptons, downstream_joint

T2=abs(MU_STAR*M_D/M_S**2-1.0)

def local_d_interval(mus):
    """Exact log-d interval where slot-1 and slot-2 s-windows can overlap."""
    L0,L1=math.log(Q_LO),math.log(Q_HI)
    F2=ALPHA_K**2*mus
    s1lo=np.log(F2)-math.log(1+STRANGE_ANCHOR_THRESHOLD)
    s1hi=np.log(F2)-math.log(1-STRANGE_ANCHOR_THRESHOLD)
    zlo=2*s1lo-np.log(mus)+math.log(1-T2)
    zhi=2*s1hi-np.log(mus)+math.log(1+T2)
    return np.maximum(zlo,L0),np.minimum(zhi,L1)

def one(seed,n,rho):
    rng=np.random.default_rng(seed)
    a,b,c,wlep=conditioned_leptons(n,rng)
    mus=a+b+c
    lo,hi=local_d_interval(mus)
    width=np.maximum(0.0,hi-lo)
    valid=width>0

    L0,L1=math.log(Q_LO),math.log(Q_HI)
    span=L1-L0
    choose=(rng.random(n)<rho)&valid
    z=rng.uniform(L0,L1,n)
    z[choose]=lo[choose]+rng.random(int(choose.sum()))*width[choose]

    # Proposal density q(z|lep). When the local interval is absent the proposal
    # falls back to the target global uniform exactly.
    q=np.full(n,(1-rho)/span)
    inlocal=valid&(z>=lo)&(z<=hi)
    q += np.where(inlocal,rho/np.maximum(width,1e-300),0.0)
    q[~valid]=1.0/span
    lr=(1.0/span)/q

    j5,j6,_=downstream_joint(a,b,c,np.exp(z))
    W=wlep*lr
    sw=W.sum()
    e5=float(np.sum(W*j5)/sw)
    e6=float(np.sum(W*j6)/sw)
    ess=float(sw*sw/np.sum(W*W))
    return {
        "P_K":float(wlep.mean()),"joint5":e5,"joint6":e6,
        "ratio":e6/e5,"ESS":ess,
        "nonzero":float(np.mean(j5>0))
    }

def main(n,seeds,rho):
    rows=[]
    for s in seeds:
        r=one(s,n,rho);rows.append(r)
        print(f"seed {s}: P_K={r['P_K']:.8e} joint5={r['joint5']:.8e} "
              f"joint6={r['joint6']:.8e} ratio={r['ratio']:.6f} "
              f"ESS={r['ESS']:.0f} nonzero={r['nonzero']:.3f}")
    for key in ("P_K","joint5","joint6"):
        x=np.array([r[key] for r in rows])
        print(f"{key} mean +/- seed-SEM = {x.mean():.8e} +/- "
              f"{x.std(ddof=1)/math.sqrt(len(x)):.2e}")
    m5=np.mean([r["joint5"] for r in rows])
    m6=np.mean([r["joint6"] for r in rows])
    pk=np.mean([r["P_K"] for r in rows])
    print(f"slot6 cost = {m6/m5:.6f}")
    print(f"full-chain five-primary = {pk*m5:.8e}")
    print(f"full-chain true-six = {pk*m6:.8e}")

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--n",type=int,default=500_000)
    ap.add_argument("--rho",type=float,default=0.95)
    ap.add_argument("--seeds",type=int,nargs="+",default=[9001,9002,9003,9004,9005])
    a=ap.parse_args()
    if not 0<a.rho<1: raise SystemExit("--rho must lie strictly between 0 and 1")
    main(a.n,a.seeds,a.rho)
