#!/usr/bin/env python3
"""Independent reconstruction of the manuscript's six-class 3740-cell grammar."""
import itertools, math
from inputs import PDG, NAMES, MU_STAR

LEP={0,1,2}
ANCHOR_KEYS=["mustar","P_2me","G","tau","mustar_plus_tau","two_tau",
             "mustar_over_2","two_mustar","one_point_five_mustar"]
CONSTRUCTION={i:{i} for i in range(9)}
CONSTRUCTION.update({9:LEP,10:{0},11:LEP,12:{2},13:LEP,14:{2},
                     15:LEP,16:LEP,17:LEP})

def anchors(m):
    ms=m[0]+m[1]+m[2]
    return [ms,2*m[0],math.sqrt(1.5)*ms,m[2],ms+m[2],2*m[2],
            ms/2,2*ms,1.5*ms]

def build():
    out=[]; idx=range(9); pairs=list(itertools.combinations(idx,2))
    for a,b in pairs:
        for A in range(9): out.append(("P",(a,b,-1,0),9+A))
    for a,b in pairs:
        for x in idx:
            if x not in (a,b): out.append(("T",(a,b,x,0),None))
    for A in range(9):
        for a in idx:
            for x in idx:
                if A==3 and a==2 and x==2: continue
                out.append(("A",(A,a,x,0),None))
    for op in ("S","D"):
        for a,b in pairs:
            for tn in range(18):
                if not CONSTRUCTION[tn] & {a,b}:
                    out.append((op,(a,b,-1,0),tn))
    for a,b,c in itertools.combinations(idx,3):
        for sign in (1,-1):
            for tn in range(18):
                if not CONSTRUCTION[tn] & {a,b,c}:
                    out.append(("K",(a,b,c,sign),tn))
    return out

COMPS=build()

def evaluated(m):
    t=list(m)+anchors(m); sq=[math.sqrt(x) for x in m]; rows=[]
    for op,(i,j,k,sgn),tn in COMPS:
        if op=="P": val,tgt=m[i]*m[j],t[tn]**2
        elif op=="T": val,tgt=m[k]**2,m[i]*m[j]
        elif op=="A": val,tgt=t[9+i]*m[j],m[k]**2
        elif op=="S": val,tgt=m[i]+m[j],t[tn]
        elif op=="D": val,tgt=abs(m[i]-m[j]),t[tn]
        else:
            s=sq[i]+sq[j]+sq[k]
            w=2*math.sqrt(sq[i]*sq[j]+sq[j]*sq[k]+sq[k]*sq[i])
            k4=s+w if sgn>0 else abs(s-w)
            val,tgt=k4*k4,t[tn]
        rows.append((op,(i,j,k,sgn),tn,val/tgt-1))
    return rows

if __name__=="__main__":
    assert len(COMPS)==3740
    rows=evaluated(PDG)
    ladder=[.001,.0025,.005,.01]
    counts=[sum(abs(r[3])<x for r in rows) for x in ladder]
    from collections import Counter
    classes=Counter(r[0] for r in rows)
    print("comparisons:",len(rows),dict(classes))
    print("tolerance ladder:",counts,"expected [0,1,3,6]")
    assert counts==[0,1,3,6]
    print("grammar checks: PASS")
