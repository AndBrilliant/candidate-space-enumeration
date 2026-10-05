import numpy as np, math, sys, time
sys.path.insert(0,'.')
import conditional_joint as cj
from inputs import load
d=load(); m=d["masses"]; ak=d["alpha_K"]; ms=d["mustar"]
e,mu,tau,u,dw,s,c,b,t=m
slots={1:abs(s/(ak**2*ms)-1),2:abs(ms*dw/s**2-1),3:abs(2*e*dw/u**2-1),4:abs(c/(3*ak*ms)-1),5:abs(b/(0.5*ms/ak)-1),6:abs(c*b/d["G"]**2-1)}
t0=time.time()
P5,se5,P6,se6,n=cj.joint_probabilities(d,slots,1_000_000,seed=20262001)
print(f"({time.time()-t0:.1f}s) P5={P5:.4e} +-{se5:.1e}  P6={P6:.4e} +-{se6:.1e}  ratio={P6/P5:.3f} n={n}")
