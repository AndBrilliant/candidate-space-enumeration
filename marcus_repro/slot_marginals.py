#!/usr/bin/env python3
"""Corrected conditional per-slot diagnostics under the actual finite-window law."""
import argparse, math
import numpy as np
from inputs import *
from conditional_replay import conditioned_leptons

def main(n,seed):
    rng=np.random.default_rng(seed)
    e,mu,tau,w=conditioned_leptons(n,rng)
    mus=e+mu+tau
    q=np.exp(rng.uniform(math.log(Q_LO),math.log(Q_HI),(n,5)))
    u,d,s,c,b=q.T
    t=[
        STRANGE_ANCHOR_THRESHOLD,
        abs(MU_STAR*M_D/M_S**2-1.0),
        abs(P*M_D/M_U**2-1.0),
        abs(M_C_SELF/(3*ALPHA_K*MU_STAR)-1.0),
        abs(M_B_SELF/(0.5*MU_STAR/ALPHA_K)-1.0),
        abs(M_C_SELF*M_B_SELF/G**2-1.0),
    ]
    events=[
        np.abs((ALPHA_K**2*mus)/s-1)<=t[0],
        np.abs(mus*d/s**2-1)<=t[1],
        np.abs((2*e)*d/u**2-1)<=t[2],
        np.abs(c/(3*ALPHA_K*mus)-1)<=t[3],
        np.abs(b/(0.5*mus/ALPHA_K)-1)<=t[4],
        np.abs(c*b/(1.5*mus**2)-1)<=t[5],
    ]
    names=["F2/s","mu*d/s2","P*d/u2","c/(3aKmu)","b/(0.5mu/aK)","cb/G2"]
    sw=w.sum()
    print(f"N={n} seed={seed}")
    for name,thr,ev in zip(names,t,events):
        rate=float(np.sum(w*ev)/sw)
        print(f"{name:12s} threshold={100*thr:.6f}%  conditional rate={rate:.8e}")

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--n",type=int,default=2_000_000)
    ap.add_argument("--seed",type=int,default=12345)
    a=ap.parse_args();main(a.n,a.seed)
