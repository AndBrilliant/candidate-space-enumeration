"""Physical finite-width Koide-conditioned lepton proposals.

The Jacobian weights transform a uniform (ln M, theta, Q) proposal into the
conditional distribution of three independently log-uniform *labeled* masses,
restricted to |Q-2/3| < tolerance. Mass limits are in common units.

IMPORTANT: these are WEIGHTED proposals. Averages under the conditioned null
must use sum(weight*f)/sum(weight), not the unweighted sample mean.
"""
import math
import numpy as np

def finite_band_proposal(n, rng, lo, hi, tolerance=2.2e-6):
    out, weights, got = [], [], 0
    while got < n:
        k = min(500000, max(4*(n-got), 8192))
        M = np.exp(rng.uniform(math.log(lo/2), math.log(hi/2), k))
        th = rng.uniform(0, 2*math.pi, k)
        Q = rng.uniform(2/3-tolerance, 2/3+tolerance, k)
        r = np.sqrt(6*Q-2)
        phase = th[:, None]+2*math.pi*np.arange(3)/3
        amp = 1+r[:,None]*np.cos(phase)
        masses = M[:,None]*amp**2
        ok=(amp>0).all(axis=1)&(masses>=lo).all(axis=1)&(masses<=hi).all(axis=1)
        ix=np.flatnonzero(ok)[:n-got]
        if not len(ix): continue
        f, rr, ang = amp[ix], r[ix], phase[ix]
        dtheta=-2*rr[:,None]*np.sin(ang)/f
        dq=6*np.cos(ang)/(rr[:,None]*f)
        j=np.abs(np.linalg.det(np.stack((np.ones_like(f),dtheta,dq),axis=-1)))
        out.append(masses[ix]);weights.append(j);got+=len(ix)
    return np.concatenate(out),np.concatenate(weights)

def weighted_mean(values, weights):
    return float(np.average(values,weights=weights))

def physical_Q(masses):
    return masses.sum(axis=-1)/np.sqrt(masses).sum(axis=-1)**2

def validated_exact_surface(n,rng,lo,hi):
    """Optional exact-surface procedural generator; NOT the rejection-conditioned law."""
    out=[];got=0
    while got<n:
        k=min(500000,max(4*(n-got),8192))
        M=np.exp(rng.uniform(math.log(lo/2),math.log(hi/2),k))
        th=rng.uniform(0,2*math.pi,k)
        amp=1+math.sqrt(2)*np.cos(th[:,None]+2*math.pi*np.arange(3)/3)
        m=M[:,None]*amp**2
        ok=(amp>0).all(axis=1)&(m>=lo).all(axis=1)&(m<=hi).all(axis=1)
        arr=m[ok][:n-got];out.append(arr);got+=len(arr)
    return np.concatenate(out)
