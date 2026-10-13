#!/usr/bin/env python3
"""Reproducible, relocatable corrected Koide statistical regression.
Run: python run_all.py --n 8000 --seeds 3 --nodes 64 --out results.json
"""
import argparse,json,math
from pathlib import Path
import numpy as np
import independent_rare_event_audit as b
import joint_baseline as j
import slot_marginals as sm
import uniform_mass as um
from koide_conditioning import finite_band_proposal,physical_Q
REF_MARGINALS=[.001106404,.0006984004,.001045038,.0004541013,.0001446706,.0002653203]
REF={'baseline':(3.421e-17,2.328e-17),'double':(3.340e-18,2.273e-18),'uniform':(1.390e-27,9.455e-28)}
def aggregate(rows,key):
    y=np.array([x[key] for x in rows],float)
    return {'mean':float(y.mean()),'between_seed_se':float(y.std(ddof=1)/np.sqrt(len(y))) if len(y)>1 else None,
            'within_seed_se':float(np.sqrt(sum(x.get(key+'_se',0)**2 for x in rows))/len(rows))}
def run(n=8000,seeds=3,nodes=64):
    if seeds<2:raise ValueError('at least 2 seeds required')
    result={'input':{'n_per_seed':n,'seeds':seeds,'nodes':nodes},'laws':{},'checks':{}}
    try:
      for tag,factor in [('baseline',1.),('double',2.)]:
        rr=[j.run(tag,factor,29000+101*i,n=n,order=nodes) for i in range(seeds)]
        result['laws'][tag]={'runs':rr,'five':aggregate(rr,'five'),'six':aggregate(rr,'six')}
        assert all(0<=x['six']<=x['five'] for x in rr)
      b.QL,b.QH,b.L=j.ORIG
      margin=[sm.runs(32000+101*i,n=n,nodes=nodes) for i in range(seeds)]
      marr=np.array([x['marginal'] for x in margin]);errs=np.array([x['se'] for x in margin])
      result['marginals']={'runs':margin,'mean':marr.mean(axis=0).tolist(),'between_seed_se':(marr.std(axis=0,ddof=1)/np.sqrt(seeds)).tolist(), 'within_seed_se':(np.sqrt(np.sum(errs**2,axis=0))/seeds).tolist()}
      uniform=[]
      for i in range(seeds):
        m,w=um.proposal(n,35000+101*i)
        y6=um.uniform_conditional(m,max(32,nodes),max(24,nodes//2))
        # uniform five-cut independently, using same d quadrature but no coupled product cut:
        # recompute charm and bottom individual widths via 1D down-draw quadrature of s*u.
        from scipy.special import roots_legendre
        star=m.sum(axis=1);E=2*m[:,0];low=b.QL;high=b.QH;span=high-low
        sL=np.maximum(b.A*b.A*star*(1-b.T1),low);sH=np.minimum(b.A*b.A*star*(1+b.T1),high)
        dL=np.maximum(low,sL*sL/(star*(1+b.T2)));dH=np.minimum(high,sH*sH/(star*(1-b.T2)))
        dL=np.maximum(dL,low*low*(1-b.T3)/E);dw=np.maximum(0,dH-dL)
        z,zw=roots_legendre(max(32,nodes));d=(dL+dH)[:,None]/2+dw[:,None]*z/2
        sw=np.maximum(0,np.minimum(sH[:,None],np.sqrt(star[:,None]*d*(1+b.T2)))-np.maximum(sL[:,None],np.sqrt(star[:,None]*d*(1-b.T2))))
        uw=np.maximum(0,np.minimum(high,np.sqrt(E[:,None]*d*(1+b.T3)))-np.maximum(low,np.sqrt(E[:,None]*d*(1-b.T3))))
        prod=(sw*uw)@zw*dw/2
        cw=np.maximum(0,np.minimum(high,3*b.A*star*(1+b.T4))-np.maximum(low,3*b.A*star*(1-b.T4)))
        bw=np.maximum(0,np.minimum(high,(2+b.A)*star*(1+b.T5))-np.maximum(low,(2+b.A)*star*(1-b.T5)))
        y5=prod*cw*bw/span**5
        def weighted(y):
            avg=np.average(y,weights=w);se=np.sqrt(np.sum((w*(y-avg))**2))/w.sum()
            return float(avg),float(se),float((w@y)**2/np.dot(w*y,w*y)) if np.any(y) else 0.
        five,fs,ess5=weighted(y5);six,ss,ess6=weighted(y6)
        uniform.append({'seed':35000+101*i,'n':n,'five':five,'six':six,'five_se':fs,'six_se':ss,'six_event_ess':ess6,'weight_ess':float(w.sum()**2/np.dot(w,w))})
      result['laws']['uniform']={'runs':uniform,'five':aggregate(uniform,'five'),'six':aggregate(uniform,'six')}
      bound=um.half_bound();result['laws']['half']={'five':{'mean':0.},'six':{'mean':0.},'support_proof':bound}
      assert bound['no_overlap_proven']
      m,w=finite_band_proposal(3000,np.random.default_rng(45000),b.LEPLO,b.LEPHI)
      err=np.max(np.abs(physical_Q(m)-2/3))
      assert err<2.2e-6 and np.all(np.isfinite(w)) and np.all(w>0)
      result['checks']['physical_koide_max_error']=float(err)
      for tag in ('baseline','double','uniform'):
        for key,v in zip(('five','six'),REF[tag]):
          got=result['laws'][tag][key]['mean']
          # deliberately loose smoke-test tolerance at n=8k; reproducibility is supported by multi-seed SE instead
          assert abs(got/v-1)<.25,(tag,key,got,v)
      assert all(abs(v/r-1)<.08 for v,r in zip(result['marginals']['mean'],REF_MARGINALS))
      result['checks']['all_regression_checks_pass']=True
      return result
    finally:b.QL,b.QH,b.L=j.ORIG
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=8000);p.add_argument('--seeds',type=int,default=3);p.add_argument('--nodes',type=int,default=64);p.add_argument('--out',default='results.json');a=p.parse_args()
 result=run(a.n,a.seeds,a.nodes)
 Path(a.out).write_text(json.dumps(result,indent=2));print(json.dumps({'checks':result['checks'],'summary':{k:[v['five']['mean'],v['six']['mean']] for k,v in result['laws'].items()}},indent=2))
