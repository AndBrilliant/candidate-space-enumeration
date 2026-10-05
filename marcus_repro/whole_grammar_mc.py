#!/usr/bin/env python3
"""Independent whole-grammar procedural replay."""
import math, random, argparse
from grammar import evaluated
from inputs import LEP_LO,LEP_HI,Q_LO,Q_HI

def main(n,seed):
    rng=random.Random(seed)
    ll=(math.log(LEP_LO),math.log(LEP_HI))
    lq=(math.log(Q_LO),math.log(Q_HI))
    thrs=(.001,.0025,.005,.01)
    sums=[0]*4; ge1=[0]*4; ge2=[0]*4; sharper=0
    sharp=0.0017417473
    maxima=[0]*4
    for z in range(n):
        m=[math.exp(rng.uniform(*ll)) for _ in range(3)]
        m += [math.exp(rng.uniform(*lq)) for _ in range(6)]
        rr=[abs(x[3]) for x in evaluated(m)]
        best=min(rr)
        sharper += best < sharp
        for j,t in enumerate(thrs):
            c=sum(x<t for x in rr)
            sums[j]+=c; ge1[j]+=c>=1; ge2[j]+=c>=2; maxima[j]=max(maxima[j],c)
    print(f"seed={seed} N={n}")
    for t,s,a,b,mx in zip(thrs,sums,ge1,ge2,maxima):
        print(f"{100*t:.2g}%: mean={s/n:.4f} P>=1={a/n:.4f} P>=2={b/n:.4f} max={mx}")
    print(f"P(best < real 0.17417473%)={sharper/n:.4f}")

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--n",type=int,default=10000)
    ap.add_argument("--seed",type=int,default=20261001)
    a=ap.parse_args(); main(a.n,a.seed)
