#!/usr/bin/env python3
"""Six individual slot probabilities for finite-band Jacobian-weighted Koide lepton law.
Independent integration over conditional quark log-uniform law; 3 independent seeds.
"""
import math,json,sys
import numpy as np
from scipy.special import roots_legendre
from koide_conditioning import finite_band_proposal
import independent_rare_event_audit as ref
A=ref.A; T1,T2,T3,T4,T5,T6=ref.T
lo,hi=ref.QL,ref.QH; lnlo,lnhi=np.log(lo),np.log(hi); span=lnhi-lnlo

def runs(seed,n=140000,nodes=96):
 rng=np.random.default_rng(seed)
 m,w=finite_band_proposal(n,rng,ref.LEPLO,ref.LEPHI)
 star=m.sum(axis=1); E=2*m[:,0]
 def logwidth(center,f,half=False):
  mult=.5 if half else 1
  return np.maximum(0,np.minimum(center+mult*np.log1p(f),lnhi)-np.maximum(center+mult*np.log1p(-f),lnlo))/span
 p1=logwidth(np.log(A*A*star),T1)
 p4=logwidth(np.log(3*A*star),T4)
 p5=logwidth(np.log((2+A)*star),T5)
 # p2,p3 average over d ~ log-uniform. Fixed high-order Legendre integration.
 xi,wi=roots_legendre(nodes)
 logd=(lnhi+lnlo)/2+span*xi/2
 p2=np.zeros(n);p3=np.zeros(n)
 for x,ww in zip(logd,wi):
  p2+=ww*logwidth(.5*(np.log(star)+x),T2,True)*.5
  p3+=ww*logwidth(.5*(np.log(E)+x),T3,True)*.5
 # p6 product-only: log mass product at 1.5star^2. Integrate c exactly in
 # log space via two uniform log draws using clipped triangular CDF.
 def G(z,l):
  return .5*(np.maximum(z,0)**2 -2*np.maximum(z-l,0)**2 + np.maximum(z-2*l,0)**2)
 # C and B are independent uniform on [lnlo,lnhi]. Sum shift is triangular.
 target=np.log(1.5*star*star)-2*lnlo
 p6=(G(target+math.log1p(T6),span)-G(target+math.log1p(-T6),span))/span**2
 p6=np.clip(p6,0,1)
 mat=np.stack([p1,p2,p3,p4,p5,p6],axis=1)
 sm=w.sum();means=(w[:,None]*mat).sum(axis=0)/sm
 centered=(mat-means)*w[:,None];se=np.sqrt(np.sum(centered**2,axis=0))/sm
 ess=((w[:,None]*mat).sum(axis=0))**2/np.sum((w[:,None]*mat)**2,axis=0)
 return {'seed':seed,'n':n,'marginal':means.tolist(),'se':se.tolist(),'ess':ess.tolist(),'lepton_ess':float(sm**2/np.dot(w,w)),'nodes':nodes}
if __name__=='__main__':
 out={'thresholds':ref.T,'runs':[runs(k) for k in (742301,742302,742303)]}
 a=np.array([r['marginal'] for r in out['runs']]);s=np.array([r['se'] for r in out['runs']])
 out['pooled']={'mean':a.mean(axis=0).tolist(),'within_se':(np.sqrt((s*s).sum(axis=0))/3).tolist(),'between_se':(a.std(axis=0,ddof=1)/math.sqrt(3)).tolist()}
 with open(str(__import__('pathlib').Path(__file__).resolve().parent/'corrected_slot_results.json'),'w') as f:json.dump(out,f,indent=2)
 print(json.dumps(out['pooled'],indent=2))
