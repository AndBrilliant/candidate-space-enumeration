"""Independent six-slot conditional importance sampler; no repository imports."""
import numpy as np, math, json
from datetime import datetime,timezone
A=math.sqrt(1.5)-1
LEPLO=.51099895069; LEPHI=1776.930
QL=2.16*1.0175; QH=172600.; L=math.log(QH/QL)
e,mu,tau=LEPLO,105.6583755,LEPHI
S=e+mu+tau; u,d,s,c,b=2.16*1.0175,4.70*1.0175,92.9*1.0175,1272.9,4186.
# printed unrounded-archive threshold for slot 1, as disclosed in manuscript
T1=0.006242571  # exact observed F^2 vs m_s residual (F2 audit fix)
T2=abs(S*d/s**2-1);T3=abs(2*e*d/u**2-1);T4=abs(c/(3*A*S)-1);T5=abs(b/((2+A)*S)-1);T6=abs(c*b/(1.5*S**2)-1)
T=[T1,T2,T3,T4,T5,T6]

def generate_balanced(n,rng):
    out=np.empty((n,3)); filled=0; accepted=0; attempted=0
    while filled<n:
        batch=min(2_000_000,max(50000,4*(n-filled)))
        m=np.exp(rng.uniform(math.log(LEPLO/2),math.log(LEPHI/2),batch))
        theta=rng.uniform(0,2*math.pi,batch)
        v=(1+math.sqrt(2)*np.cos(theta[:,None]+2*math.pi*np.arange(3)/3))**2*m[:,None]
        good=v[(v>=LEPLO).all(axis=1)&(v<=LEPHI).all(axis=1)]
        k=min(len(good),n-filled);out[filled:filled+k]=good[:k]
        filled+=k;accepted+=len(good);attempted+=batch
    return out,accepted/attempted

def posint(x,a,b):
    # rectangular CDF for sum of two uniform coordinate lengths
    return .5*(np.maximum(x,0)**2-np.maximum(x-a,0)**2-np.maximum(x-b,0)**2+np.maximum(x-a-b,0)**2)

def run(seed,n,sigma):
    rng=np.random.default_rng(seed)
    lep,accept=generate_balanced(n,rng)
    star=lep.sum(axis=1); E=2*lep[:,0]
    center=np.log(A**4*star)
    xd=center+sigma*rng.normal(size=n)
    dval=np.exp(xd)
    inside=(xd>=math.log(QL))&(xd<=math.log(QH))
    lw=np.zeros(n)
    lw[inside]=sigma*math.sqrt(2*math.pi)/L*np.exp(((xd[inside]-center[inside])**2)/(2*sigma**2))
    # construct exact log windows independently of source code
    s1lo=np.log(A*A*star)+math.log1p(-T1);s1hi=np.log(A*A*star)+math.log1p(T1)
    s2lo=.5*(np.log(star)+xd+math.log1p(-T2));s2hi=.5*(np.log(star)+xd+math.log1p(T2))
    overlap=np.maximum(0,np.minimum(np.minimum(s1hi,s2hi),math.log(QH))-np.maximum(np.maximum(s1lo,s2lo),math.log(QL)))
    ulo=.5*(np.log(E)+xd+math.log1p(-T3));uhi=.5*(np.log(E)+xd+math.log1p(T3))
    uw=np.maximum(0,np.minimum(uhi,math.log(QH))-np.maximum(ulo,math.log(QL)))
    clo=np.log(3*A*star)+math.log1p(-T4);chi=np.log(3*A*star)+math.log1p(T4)
    blo=np.log((2+A)*star)+math.log1p(-T5);bhi=np.log((2+A)*star)+math.log1p(T5)
    cl=np.maximum(clo,math.log(QL)); ch=np.minimum(chi,math.log(QH))
    bl=np.maximum(blo,math.log(QL)); bh=np.minimum(bhi,math.log(QH))
    widthc=np.maximum(ch-cl,0);widthb=np.maximum(bh-bl,0)
    # shift to lower-left of log-space rectangle, intersection with log product strip
    baseline=np.log(1.5*star*star)
    high=baseline+math.log1p(T6)-cl-bl
    low=baseline+math.log1p(-T6)-cl-bl
    area=np.maximum(0,np.minimum(widthc*widthb,posint(high,widthc,widthb)-posint(low,widthc,widthb)))
    common=lw*overlap*uw/L**2
    y5=common*widthc*widthb/L**2
    y6=common*area/L**2
    ans={'seed':seed,'n':n,'sigma':sigma,'acceptance_koide_surface_proposal':accept}
    for name,y in [('joint5',y5),('joint6',y6)]:
        sm=y.sum(); ess=sm**2/np.dot(y,y) if sm>0 else 0
        ans[name]={'mean':float(y.mean()),'se':float(y.std(ddof=1)/math.sqrt(n)),'ess':float(ess),'contributors':int((y>0).sum()),'max_weight_fraction':float(y.max()/sm) if sm else 0,'top10_weight_fraction':float(np.partition(y,-10)[-10:].sum()/sm) if sm else 0}
    ans['ratio']=ans['joint6']['mean']/ans['joint5']['mean']
    return ans

if __name__=='__main__':
    import sys
    n=int(sys.argv[1]) if len(sys.argv)>1 else 1000000
    cases=[run(48151,n,.5),run(48152,n,.5),run(48153,n,1.),run(48154,n,1.)]
    result={'timestamp_utc':datetime.now(timezone.utc).isoformat(),'thresholds':T,'quark_log_span':[QL,QH],'runs':cases}
    for tag in ['joint5','joint6']:
        means=[x[tag]['mean'] for x in cases]
        result[tag+'_pooled']={'mean':float(np.mean(means)), 'within_run_se':float(np.sqrt(sum(x[tag]['se']**2 for x in cases))/4),'seed_sd':float(np.std(means,ddof=1)),'pooled_ess':float(sum(x[tag]['ess'] for x in cases))}
    result['ratio_pooled']=result['joint6_pooled']['mean']/result['joint5_pooled']['mean']
    open('independent_rare_event_results.json','w').write(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
