import numpy as np, math, time, json
rng=np.random.default_rng(20262001)
AK=math.sqrt(1.5)-1.0
LEP_LO,LEP_HI=0.51099895069,1776.930
Q_LO,Q_HI=2.16*1.01750,172600.0
T=[0.006216,0.007872,0.011826,0.002559,0.000815]
T6=0.00174175
def surface(n):
    out,got=[],0
    m_lo,m_hi=LEP_LO/2,LEP_HI/2
    while got<n:
        c=min(1_000_000,max(4*(n-got),65536))
        Mv=np.exp(rng.uniform(math.log(m_lo),math.log(m_hi),c))
        th=rng.uniform(0,2*math.pi,c)
        f=1+math.sqrt(2)*np.cos(th[:,None]+2*math.pi*np.arange(3)/3)
        tr=Mv[:,None]*f*f
        ok=(tr>=LEP_LO).all(1)&(tr<=LEP_HI).all(1)
        out.append(tr[ok]); got+=len(tr[ok])
    return np.concatenate(out)[:n]
lep=surface(500_000)
mu_s=lep.sum(1); p2e=2*lep[:,0]; n=len(lep)
d=np.exp(rng.uniform(math.log(Q_LO),math.log(Q_HI),n))
glo,ghi=math.log(Q_LO),math.log(Q_HI); L=ghi-glo
t1,t2,t3,t4,t5=T
lnF2=np.log(AK**2*mu_s); lns2=0.5*np.log(mu_s*d)
lnu=0.5*np.log(p2e*d); lnc=np.log(3*AK*mu_s); lnb=np.log(0.5*mu_s/AK)
def wl(c,f,h):
    return np.maximum(0.0,np.minimum(c+math.log(1+f)/(2 if h else 1),ghi)-np.maximum(c+math.log(1-f)/(2 if h else 1),glo))
s1=wl(lnF2,t1,False); s2=wl(lns2,t2,True)
ov=np.maximum(0.0,np.minimum(lnF2+math.log(1+t1),lns2+0.5*math.log(1+t2))-np.maximum(lnF2+math.log(1-t1),lns2+0.5*math.log(1-t2)))
ov=np.minimum(ov,np.minimum(s1,s2))
pu=wl(lnu,t3,True)/L; pc=wl(lnc,t4,False)/L; pb=wl(lnb,t5,False)/L
j5=ov/L*pu*pc*pb
# P6: (c,b) rectangle ^ product strip, log-area via 256-pt c-grid
Gsq=1.5*mu_s**2
c_lo=np.maximum(3*AK*mu_s*(1-t4),Q_LO); c_hi=np.minimum(3*AK*mu_s*(1+t4),Q_HI)
b_lo=np.maximum(0.5*mu_s/AK*(1-t5),Q_LO); b_hi=np.minimum(0.5*mu_s/AK*(1+t5),Q_HI)
ng=256; frac=np.linspace(0,1,ng); area6=np.zeros(n)
valid=(c_hi>c_lo)&(b_hi>b_lo)
idx=np.where(valid)[0]
if len(idx):
    lc=np.log(c_lo[idx])[None,:]+frac[:,None]*(np.log(c_hi[idx])-np.log(c_lo[idx]))[None,:]
    cv=np.exp(lc)
    blo=np.maximum(b_lo[idx][None,:],(Gsq[idx]*(1-T6))[None,:]/cv)
    bhi=np.minimum(b_hi[idx][None,:],(Gsq[idx]*(1+T6))[None,:]/cv)
    w=np.log(np.maximum(bhi,blo)/blo)
    dln=(np.log(c_hi[idx])-np.log(c_lo[idx]))/(ng-1)
    a=dln*(w[1:-1].sum(0)+0.5*(w[0]+w[-1]))
    nohit=(bhi<=blo).all(0)
    a[nohit]=0.0
    area6[idx]=a
wu=wl(lnu,t3,True)
j6=(ov/L)*(wu/L)*(area6/L**2)
res=dict(n=n,P5=float(j5.mean()),P5_se=float(j5.std(ddof=1)/math.sqrt(n)),
         P6=float(j6.mean()),P6_se=float(j6.std(ddof=1)/math.sqrt(n)),
         ratio=float(j6.mean()/j5.mean()))
print(res)
json.dump(res,open('/workspace/repo/outputs/conditional_joint.json','w'),indent=2)
