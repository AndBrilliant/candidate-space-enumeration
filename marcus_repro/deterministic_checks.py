#!/usr/bin/env python3
import math
from inputs import *

def pct(x): return 100.0*x

q = (M_E+M_MU+M_TAU)/(math.sqrt(M_E)+math.sqrt(M_MU)+math.sqrt(M_TAU))**2
d=(1/math.sqrt(3),)*3
x=(math.sqrt(M_E),math.sqrt(M_MU),math.sqrt(M_TAU))
cos=sum(a*b for a,b in zip(x,d))/math.sqrt(sum(a*a for a in x))
theta=math.degrees(math.acos(cos))

mc_pred=3*ALPHA_K*MU_STAR
mb_pred=0.5*MU_STAR/ALPHA_K
rc=M_C_SELF/mc_pred-1
rb=M_B_SELF/mb_pred-1
rp=M_C_SELF*M_B_SELF/G**2-1

print(f"mu_star = {MU_STAR:.12f} MeV")
print(f"alpha_K = {ALPHA_K:.15f}")
print(f"Q_l = {q:.10f}")
print(f"theta = {theta:.8f} deg")
print(f"G = {G:.9f} MeV")
print(f"mc_pred = {mc_pred:.9f} MeV ; residual = {pct(rc):+.6f}%")
print(f"mb_pred = {mb_pred:.9f} MeV ; residual = {pct(rb):+.6f}%")
print(f"cb/G^2 residual = {pct(rp):+.6f}%")

assert abs(MU_STAR-1883.09937445069) < 1e-9
assert abs(q-0.6666644634) < 5e-11
assert abs(theta-44.99990532) < 5e-8
assert abs(pct(rc)-0.25591445) < 1e-6
assert abs(pct(rb)+0.08153107) < 1e-6
assert abs(pct(rp)-0.17417473) < 1e-6
print("deterministic checks: PASS")
