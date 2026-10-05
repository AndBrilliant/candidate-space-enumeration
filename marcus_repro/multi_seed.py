#!/usr/bin/env python3
"""Independent multi-seed convergence driver for corrected conditional replay."""
import argparse, math
import numpy as np
from inputs import *
from conditional_replay import conditioned_leptons, downstream_joint

def one(seed,n):
    rng=np.random.default_rng(seed)
    a,b,c,w=conditioned_leptons(n,rng)
    d=np.exp(rng.uniform(math.log(Q_LO),math.log(Q_HI),n))
    j5,j6,_=downstream_joint(a,b,c,d)
    sw=w.sum()
    return float(w.mean()),float(np.sum(w*j5)/sw),float(np.sum(w*j6)/sw)

def main(n,seeds):
    rows=[]
    for s in seeds:
        r=one(s,n); rows.append(r)
        print(f"seed {s}: P_K={r[0]:.8e} joint5={r[1]:.8e} joint6={r[2]:.8e}")
    A=np.array(rows)
    mean=A.mean(axis=0); sem=A.std(axis=0,ddof=1)/math.sqrt(len(rows))
    print("across-seed mean +/- SEM:")
    print(f"P_K    {mean[0]:.8e} +/- {sem[0]:.2e}")
    print(f"joint5 {mean[1]:.8e} +/- {sem[1]:.2e}")
    print(f"joint6 {mean[2]:.8e} +/- {sem[2]:.2e}")
    print(f"full5  {mean[0]*mean[1]:.8e}")
    print(f"full6  {mean[0]*mean[2]:.8e}")
    print(f"cost of slot6 = {mean[2]/mean[1]:.6f}")

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--n",type=int,default=1_000_000)
    ap.add_argument("--seeds",type=int,nargs="+",default=[731994,731995,731996,731997,731998])
    a=ap.parse_args(); main(a.n,a.seeds)
