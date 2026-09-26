#!/usr/bin/env python3
import math, sys
import rundec

crd = rundec.CRunDec()
print("rundec_module", getattr(rundec, "__version__", "unknown"))

AS_MZ=0.1180
MZ=91.1879
MB=4.186
MC=1.2729
MU=0.00216
MT=162.5
MU_STAR=(0.51099895069+105.6583755+1776.93)/1000.0

def mass_ratio_fixed_nf(mu0, mu1, nf, loops, as0):
    as1=crd.AlphasExact(as0,mu0,mu1,nf,loops)
    return crd.mMS2mMS(1.0,as0,as1,nf,loops), as1

for loops in (3,4,5):
    as5_mb=crd.AlphasExact(AS_MZ,MZ,MB,5,loops)
    as4_mb=crd.DecAsDownSI(as5_mb,MB,MB,4,loops)
    as4_2=crd.AlphasExact(as4_mb,MB,2.0,4,loops)
    as4_mus=crd.AlphasExact(as4_mb,MB,MU_STAR,4,loops)
    rm=crd.mMS2mMS(1.0,as4_2,as4_mus,4,loops)
    print(f"LOOP {loops} as5_mb={as5_mb:.12f} as4_mb={as4_mb:.12f} as4_2={as4_2:.12f} as4_mustar={as4_mus:.12f} Rm={rm:.12f}")

# Diagnostic: same 4-loop with legacy MZ value, to identify source of historical number.
legacy_mz=91.1876
loops=4
a5=crd.AlphasExact(AS_MZ,legacy_mz,MB,5,loops)
a4=crd.DecAsDownSI(a5,MB,MB,4,loops)
a2=crd.AlphasExact(a4,MB,2.0,4,loops)
amus=crd.AlphasExact(a4,MB,MU_STAR,4,loops)
print(f"LEGACY_MZ Rm={crd.mMS2mMS(1.0,a2,amus,4,loops):.12f}")

# Common-scale up-type diagnostic from PDG Live masses only.
# Run u from 2 GeV and c from mc to mt, crossing b at mb. Use mL2mH.
bottom=rundec.TriplenfMmu()
bottom.nf=5
bottom.Mth=MB
bottom.muth=MB
decb=rundec.TriplenfMmuArray(4)
decb[0]=bottom

as4_2=crd.AlphasExact(crd.DecAsDownSI(crd.AlphasExact(AS_MZ,MZ,MB,5,4),MB,MB,4,4),MB,2.0,4,4)
try:
    u_mt=crd.mL2mH(MU,as4_2,2.0,decb,MT,4)
    as4_mc=crd.AlphasExact(crd.DecAsDownSI(crd.AlphasExact(AS_MZ,MZ,MB,5,4),MB,MB,4,4),MB,MC,4,4)
    c_mt=crd.mL2mH(MC,as4_mc,MC,decb,MT,4)
    def q9(vals):
        return 9*sum(vals)/(sum(math.sqrt(x) for x in vals)**2)
    print(f"COMMON_TOP u_mt={u_mt:.12g} c_mt={c_mt:.12g} mt={MT:.12g} 9Q={q9([u_mt,c_mt,MT]):.12f}")
except Exception as e:
    print("COMMON_TOP_ERROR",repr(e))

# Also calculate exact mt needed for 9Q=8 given u,c at top from the above convention, via bisection.
try:
    def q9_mt(mt):
        return 9*(u_mt+c_mt+mt)/(math.sqrt(u_mt)+math.sqrt(c_mt)+math.sqrt(mt))**2
    lo,hi=150.,180.
    for _ in range(100):
        mid=(lo+hi)/2
        if q9_mt(mid)<8: lo=mid
        else: hi=mid
    print(f"MT_REQUIRED_COMMON_TOP={(lo+hi)/2:.12f} Q={q9_mt((lo+hi)/2):.12f}")
except Exception as e:
    print("MT_REQUIRED_ERROR",repr(e))

print("DIRECT_NF4_NO_THRESHOLD")
for loops in (3,4,5):
    a2=crd.AlphasExact(AS_MZ,MZ,2.0,4,loops)
    amus=crd.AlphasExact(AS_MZ,MZ,MU_STAR,4,loops)
    rm=crd.mMS2mMS(1.0,a2,amus,4,loops)
    print(f"DIRECT_NF4 LOOP {loops} a2={a2:.12f} amus={amus:.12f} Rm={rm:.12f}")
