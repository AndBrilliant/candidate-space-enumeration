#!/usr/bin/env python3
"""Corrected conditional rare-event replay.

Targets the manuscript's baseline lepton generator exactly:
three independent log-uniform masses on [m_e,m_tau], conditioned on the
finite event |Q-2/3| < KOIDE_EPS.

For each sampled (m1,m2), the allowed m3 intervals are solved analytically.
Their total log-width is P(Koide window | m1,m2). One m3 is sampled from
that conditional interval union, and the downstream statistic is averaged
with those exact conditional weights.

The downstream quark slots are Rao-Blackwellized analytically. Both:
  joint5 = five primary slots
  joint6 = joint5 plus the cb/G^2 product strip
are reported.
"""
import argparse, math
import numpy as np
from inputs import *

def _roots_for_q(m1,m2,q):
    A=np.sqrt(m1)+np.sqrt(m2)
    B=m1+m2
    aa=1.0-q; bb=-2.0*q*A; cc=B-q*A*A
    disc=bb*bb-4*aa*cc
    good=disc>=0
    sd=np.sqrt(np.maximum(disc,0.0))
    r1=(-bb-sd)/(2*aa); r2=(-bb+sd)/(2*aa)
    return np.where(good,r1,np.nan),np.where(good,r2,np.nan)

def conditioned_leptons(n,rng):
    """Return (m1,m2,m3,weights), where weights are exact
    P(Koide-window | m1,m2) under the original log-U m3 law."""
    lm0,lm1=math.log(LEP_LO),math.log(LEP_HI)
    m1=np.exp(rng.uniform(lm0,lm1,n))
    m2=np.exp(rng.uniform(lm0,lm1,n))
    y0,y1=math.sqrt(LEP_LO),math.sqrt(LEP_HI)
    q0,q1=2/3-KOIDE_EPS,2/3+KOIDE_EPS
    a0,b0=_roots_for_q(m1,m2,q0)
    a1,b1=_roots_for_q(m1,m2,q1)
    roots=[np.nan_to_num(x,nan=y0,posinf=y1,neginf=y0) for x in (a0,b0,a1,b1)]
    pts=np.stack([np.full(n,y0),np.full(n,y1),*roots],axis=1)
    pts=np.clip(pts,y0,y1)
    pts.sort(axis=1)
    left,right=pts[:,:-1],pts[:,1:]
    mid=0.5*(left+right)
    mm=mid*mid
    qmid=(m1[:,None]+m2[:,None]+mm)/(np.sqrt(m1)[:,None]+np.sqrt(m2)[:,None]+mid)**2
    inside=(qmid>=q0)&(qmid<=q1)&(right>left)
    segw=np.where(inside,2.0*(np.log(right)-np.log(left)),0.0) # log-mass width
    total=segw.sum(axis=1)
    weights=total/(lm1-lm0)

    # Draw one m3 from the allowed union, uniformly in log(m3).
    u=rng.random(n)*total
    cum=np.cumsum(segw,axis=1)
    idx=(cum>=u[:,None]).argmax(axis=1)
    row=np.arange(n)
    prev=np.where(idx>0,cum[row,np.maximum(idx-1,0)],0.0)
    lsel=left[row,idx]
    logm3=2*np.log(np.maximum(lsel,y0))+(u-prev)
    m3=np.exp(logm3)
    m3=np.where(total>0,m3,LEP_LO)
    return m1,m2,m3,weights

def _width(lo,hi,span_lo,span_hi):
    return np.maximum(0.0,np.minimum(hi,span_hi)-np.maximum(lo,span_lo))

def _band_area(a,b,c,d,z):
    """Area in rectangle [a,b]x[c,d] satisfying u+v <= z."""
    pos=lambda x: np.maximum(x,0.0)
    return 0.5*(pos(z-a-c)**2-pos(z-b-c)**2-pos(z-a-d)**2+pos(z-b-d)**2)

def downstream_joint(m1,m2,m3,d):
    mus=m1+m2+m3
    p2e=2*m1
    L0,L1=math.log(Q_LO),math.log(Q_HI)
    span=L1-L0

    # Observed absolute thresholds in the paper's cell orientation V/T - 1.
    t1=STRANGE_ANCHOR_THRESHOLD
    t2=abs(MU_STAR*M_D/M_S**2-1.0)
    t3=abs(P*M_D/M_U**2-1.0)
    t4=abs(M_C_SELF/(3*ALPHA_K*MU_STAR)-1.0)
    t5=abs(M_B_SELF/(0.5*MU_STAR/ALPHA_K)-1.0)
    t6=abs(M_C_SELF*M_B_SELF/G**2-1.0)

    # slot 1: F^2 / s - 1
    F2=ALPHA_K**2*mus
    s1lo=np.log(F2)-math.log(1+t1)
    s1hi=np.log(F2)-math.log(1-t1)

    # slot 2: mus*d / s^2 - 1
    cs=0.5*np.log(mus*d)
    s2lo=cs-0.5*math.log(1+t2)
    s2hi=cs-0.5*math.log(1-t2)
    slo=np.maximum(np.maximum(s1lo,s2lo),L0)
    shi=np.minimum(np.minimum(s1hi,s2hi),L1)
    ps=np.maximum(0.0,shi-slo)/span

    # slot 3: P*d / u^2 - 1
    cu=0.5*np.log(p2e*d)
    ulo=cu-0.5*math.log(1+t3)
    uhi=cu-0.5*math.log(1-t3)
    pu=_width(ulo,uhi,L0,L1)/span

    # slots 4/5: direct observed/predicted residuals
    c0=3*ALPHA_K*mus
    b0=0.5*mus/ALPHA_K
    lc0=np.log(c0); lb0=np.log(b0)
    clo=np.maximum(lc0+math.log(1-t4),L0)
    chi=np.minimum(lc0+math.log(1+t4),L1)
    blo=np.maximum(lb0+math.log(1-t5),L0)
    bhi=np.minimum(lb0+math.log(1+t5),L1)
    pc=np.maximum(0.0,chi-clo)/span
    pb=np.maximum(0.0,bhi-blo)/span
    joint5=ps*pu*pc*pb

    # slot 6: cb/G(mus)^2 - 1. Work in offsets from c0,b0, whose
    # product is identically 1.5*mus^2.
    a=clo-lc0; b=chi-lc0; c=blo-lb0; dd=bhi-lb0
    zlo=math.log(1-t6); zhi=math.log(1+t6)
    area=_band_area(a,b,c,dd,zhi)-_band_area(a,b,c,dd,zlo)
    pcb=np.maximum(0.0,area)/(span*span)
    joint6=ps*pu*pcb
    return joint5,joint6,(t1,t2,t3,t4,t5,t6)

def _ratio_and_batch_se(w,y,batches=50):
    est=float(np.sum(w*y)/np.sum(w))
    n=len(w); b=max(2,min(batches,n//1000))
    vals=[]
    for ix in np.array_split(np.arange(n),b):
        sw=w[ix].sum()
        if sw>0: vals.append(float(np.sum(w[ix]*y[ix])/sw))
    se=float(np.std(vals,ddof=1)/math.sqrt(len(vals))) if len(vals)>1 else float("nan")
    return est,se

def main(n,seed):
    rng=np.random.default_rng(seed)
    m1,m2,m3,w=conditioned_leptons(n,rng)
    d=np.exp(rng.uniform(math.log(Q_LO),math.log(Q_HI),n))
    j5,j6,t=downstream_joint(m1,m2,m3,d)

    pkoide=float(w.mean())
    pkoide_se=float(w.std(ddof=1)/math.sqrt(n))
    e5,se5=_ratio_and_batch_se(w,j5)
    e6,se6=_ratio_and_batch_se(w,j6)

    print(f"N={n} seed={seed}")
    print("thresholds (%):",", ".join(f"{100*x:.6f}" for x in t))
    print(f"P(Koide window) = {pkoide:.8e} +/- {pkoide_se:.2e}")
    print(f"conditional joint5 = {e5:.8e} +/- {se5:.2e}")
    print(f"conditional joint6 = {e6:.8e} +/- {se6:.2e}")
    print(f"full-chain joint5 = {pkoide*e5:.8e}")
    print(f"full-chain joint6 = {pkoide*e6:.8e}")
    print(f"joint6/joint5 = {e6/e5:.6f}")

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--n",type=int,default=1_000_000)
    ap.add_argument("--seed",type=int,default=731994)
    a=ap.parse_args(); main(a.n,a.seed)
