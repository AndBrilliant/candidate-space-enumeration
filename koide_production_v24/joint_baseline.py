#!/usr/bin/env python3
"""Corrected finite-band positive-amplitude Koide law and five/six-cut tables.
Separate independent code; no repository imports. Quad in log d, analytic log c,b strip.
"""
import json,math,sys
import numpy as np
from scipy.special import roots_legendre
from koide_conditioning import finite_band_proposal

def proposal(n,seed):
    return finite_band_proposal(n,np.random.default_rng(seed),b.LEPLO,b.LEPHI)
import independent_rare_event_audit as b

ORIG=(b.QL,b.QH,b.L)
def span_bounds(f):
    lo,hi=ORIG[:2];mid=(math.log(lo)+math.log(hi))/2;h=(math.log(hi)-math.log(lo))*f/2
    return math.exp(mid-h),math.exp(mid+h)
def evaluate(m,nodes=128):
    star=m.sum(axis=1);E=2*m[:,0];L=b.L
    s0=np.log(b.A**2*star)
    sL=s0+math.log1p(-b.T1);sH=s0+math.log1p(b.T1)
    dL=np.maximum(2*sL-np.log(star)-math.log1p(b.T2),math.log(b.QL))
    dH=np.minimum(2*sH-np.log(star)-math.log1p(-b.T2),math.log(b.QH))
    W=np.maximum(0,dH-dL)
    z,weights=roots_legendre(nodes)
    # memory-safe batches should be imposed by caller
    logd=(dL+dH)[:,None]/2+W[:,None]*z[None,:]/2
    low2=.5*(np.log(star)[:,None]+logd+math.log1p(-b.T2))
    high2=.5*(np.log(star)[:,None]+logd+math.log1p(b.T2))
    sw=np.maximum(0,np.minimum(sH[:,None],high2)-np.maximum(sL[:,None],low2))
    ulo=.5*(np.log(E)[:,None]+logd+math.log1p(-b.T3))
    uhi=.5*(np.log(E)[:,None]+logd+math.log1p(b.T3))
    uw=np.maximum(0,np.minimum(uhi,math.log(b.QH))-np.maximum(ulo,math.log(b.QL)))
    dintegral=((sw*uw)@weights)*W/2
    cl=np.maximum(np.log(3*b.A*star)+math.log1p(-b.T4),math.log(b.QL))
    ch=np.minimum(np.log(3*b.A*star)+math.log1p(b.T4),math.log(b.QH))
    bl=np.maximum(np.log((2+b.A)*star)+math.log1p(-b.T5),math.log(b.QL))
    bh=np.minimum(np.log((2+b.A)*star)+math.log1p(b.T5),math.log(b.QH))
    aw=np.maximum(0,ch-cl);bw=np.maximum(0,bh-bl)
    five=dintegral*aw*bw/L**5
    # use reference's independently tested positive-part integration primitive
    upper=np.log(1.5*star**2)+math.log1p(b.T6)-cl-bl
    lower=np.log(1.5*star**2)+math.log1p(-b.T6)-cl-bl
    area=np.maximum(0,np.minimum(aw*bw,b.posint(upper,aw,bw)-b.posint(lower,aw,bw)))
    six=dintegral*area/L**5
    return five,six

def run(tag,factor,seed,n=90000,order=128):
    b.QL,b.QH=span_bounds(factor);b.L=math.log(b.QH/b.QL)
    m,w=proposal(n,seed)
    sumw=sw2=0.;sy5=sy6=sv5=sv6=0.;yy5=yy6=0.;nz=0
    for start in range(0,n,4000):
        end=min(start+4000,n);w0=w[start:end]
        f5,f6=evaluate(m[start:end],order)
        sumw+=w0.sum();sw2+=(w0*w0).sum()
        sy5+=w0@f5;sy6+=w0@f6
        sv5+=(w0*f5)@(w0*f5);sv6+=(w0*f6)@(w0*f6)
        yy5+=(w0*f5*f5).sum(); yy6+=(w0*f6*f6).sum()
        nz+=np.count_nonzero(f6)
    p5=sy5/sumw;p6=sy6/sumw
    # self-normalized IS delta-method (mean of weighted residuals)
    # derive from weighted second moments, includes mean normalization uncertainty
    # variance of normalized importance estimator: sum(w*(y-p))^2/sum(w)^2
    se5=math.sqrt(max(0,sv5-2*p5*sy5+p5*p5*sw2))/sumw
    se6=math.sqrt(max(0,sv6-2*p6*sy6+p6*p6*sw2))/sumw
    return {'model':tag,'seed':seed,'n':n,'quad_nodes':order,'five':p5,'six':p6,'five_se':se5,'six_se':se6,'ratio':p6/p5 if p5 else None,'proposal_ess':sumw**2/sw2,'six_event_ess':sy6**2/sv6 if sv6 else 0,'nonzero_six':int(nz),'support':[b.QL,b.QH]}

if __name__=='__main__':
    out=[]
    for tag,f in [('baseline',1.),('double_logspan',2.)]:
        for seed in [953201,953202,953203]:
            x=run(tag,f,seed)
            print(tag,seed,f"5={x['five']:.8g} 6={x['six']:.8g} se6={x['six_se']:.2g} ESS={x['six_event_ess']:.0f}",flush=True)
            out.append(x)
    b.QL,b.QH,b.L=ORIG
    with open(str(__import__('pathlib').Path(__file__).resolve().parent/'corrected_joint_results.json'),'w') as h:json.dump({'runs':out},h,indent=2)
