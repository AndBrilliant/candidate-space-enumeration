#!/usr/bin/env python3
"""Independent mass-space six-cut probability, and analytic half-log-span support bound.
Finite-band positive-root Jacobian conditional lepton ensemble. Compare with old surface
law by using unweighted vs weighted expectations from the same proposal.
"""
import math,json,argparse
import numpy as np
from scipy.special import roots_legendre
import independent_rare_event_audit as b
from koide_conditioning import finite_band_proposal

def proposal(n,seed):
    return finite_band_proposal(n,np.random.default_rng(seed),b.LEPLO,b.LEPHI)
def uniform_conditional(m,nD=48,nC=24):
    n=len(m); star=m.sum(axis=1); e=2*m[:,0]
    low,high=b.QL,b.QH;scale=high-low
    slo=np.maximum(b.A*b.A*star*(1-b.T1),low)
    shi=np.minimum(b.A*b.A*star*(1+b.T1),high)
    # d range from s-pair overlap (and native draw span)
    dlow=np.maximum(low, slo*slo/(star*(1+b.T2)))
    dhi=np.minimum(high, shi*shi/(star*(1-b.T2)))
    # d must also enable u in draw range
    dlow=np.maximum(dlow,low*low*(1-b.T3)/e)
    width=np.maximum(0,dhi-dlow)
    nodes,weights=roots_legendre(nD)
    dc=(dhi+dlow)/2
    d=dc[:,None]+width[:,None]*nodes/2
    s2lo=np.sqrt(star[:,None]*d*(1-b.T2))
    s2hi=np.sqrt(star[:,None]*d*(1+b.T2))
    sov=np.maximum(0,np.minimum(shi[:,None],s2hi)-np.maximum(slo[:,None],s2lo))
    ulo=np.sqrt(e[:,None]*d*(1-b.T3))
    uhi=np.sqrt(e[:,None]*d*(1+b.T3))
    uw=np.maximum(0,np.minimum(uhi,high)-np.maximum(ulo,low))
    # coupled c,b mass-space area via 1-dimensional integral in c.
    cl=np.maximum(low,3*b.A*star*(1-b.T4))
    ch=np.minimum(high,3*b.A*star*(1+b.T4))
    bl=np.maximum(low,(2+b.A)*star*(1-b.T5))
    bh=np.minimum(high,(2+b.A)*star*(1+b.T5))
    cw=np.maximum(0,ch-cl)
    nc,wc=roots_legendre(nC)
    cc=(ch+cl)[:,None]/2+cw[:,None]*nc[None,:]/2
    prodlo=(1-b.T6)*1.5*star*star
    prodhi=(1+b.T6)*1.5*star*star
    # integrate clipped B interval w.r.t. actual charm mass
    bw=np.maximum(0,np.minimum(bh[:,None],prodhi[:,None]/cc)-np.maximum(bl[:,None],prodlo[:,None]/cc))
    area=np.sum(bw*wc[None,:],axis=1)*cw/2
    return np.sum(sov*uw*weights[None,:],axis=1)*width/2*area/scale**5
def one(seed,n,nd,nc):
    m,w=proposal(n,seed)
    y=uniform_conditional(m,nd,nc)
    v=np.average(y,weights=w)
    se=(np.sum((w*(y-v))**2)**.5)/np.sum(w)
    return {"seed":seed,"n":n,"d_nodes":nd,"c_nodes":nc,"unweighted":float(np.mean(y)),
            "weighted":float(v),"weighted_se":float(se),"weight_ess":float(w.sum()**2/np.dot(w,w)),
            "event_ess":float((w@y)**2/np.dot(w*y,w*y)) if np.any(y) else 0,
            "nonzero":int(np.count_nonzero(y))}
def half_bound():
    lo=math.sqrt(b.QL*b.QH)/(b.QH/b.QL)**.25
    hi=math.sqrt(b.QL*b.QH)*(b.QH/b.QL)**.25
    # Max lepton star, allowing all three masses independently <= upper
    star_max=3*b.LEPHI
    d_max=(b.A**4)*star_max*(1+b.T1)**2/(1-b.T2)
    return {"half_span_lower_quark_MeV":lo,"half_span_upper_quark_MeV":hi,
            "max_d_from_s_overlap_MeV":d_max,
            "no_overlap_proven":bool(d_max<lo),"starmass_upper_MeV":star_max}
