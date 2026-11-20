#!/usr/bin/env python3
"""Exact image/rank and symmetric-channel controls for K=C (Y Y†)^T.

The algebra concerns the selected source coefficient, not physical dipoles.
Standalone, pinned SymPy only, no repository imports or Python assertions.
"""
from pathlib import Path
import json
import sympy as s

counts={}
negatives={}


def check(name,a,b):
    d=a-b
    ok=all(s.expand(v)==0 for v in d) if isinstance(d,s.MatrixBase) else s.simplify(d)==0
    if not ok:raise RuntimeError(name+': '+str(d))
    counts[name]=counts.get(name,0)+1


def reject(name,a,b):
    d=a-b
    detected=any(s.expand(v)!=0 for v in d) if isinstance(d,s.MatrixBase) else s.simplify(d)!=0
    if not detected:raise RuntimeError('Undetected negative control: '+name)
    negatives[name]=negatives.get(name,0)+1


def nonnegative(name,a):
    if s.simplify(a).is_nonnegative is not True:raise RuntimeError(name+': '+str(a))
    counts[name]=counts.get(name,0)+1


def symbasis(n):
    out=[]
    for i in range(n):
        for j in range(i,n):
            e=s.zeros(n);e[i,j]=e[j,i]=1;out.append(e)
    return out


def flatten(m):return s.Matrix(list(m))
def stack(columns):return s.Matrix.hstack(*columns)
def snorm(a):return s.expand(sum(s.conjugate(v)*v for v in a))


def symmetric(n,seed):
    out=s.zeros(n)
    for i in range(n):
        for j in range(i,n):
            out[i,j]=out[j,i]=(i+1)*(j+2)+seed+s.I*(seed*(i+j+1)-i*j)
    return out


def rotation(n):
    out=s.eye(n)
    for j in range(n-1):
        u=s.eye(n);u[j,j]=u[j+1,j+1]=s.Rational(3,5)
        u[j,j+1]=u[j+1,j]=4*s.I/5
        out=u*out
    return out


def reconstructed_a(sm,h):
    n=len(h)
    return s.Matrix(n,n,lambda i,j:s.Rational(h[j]-h[i],h[i]+h[j])*sm[i,j] if h[i]+h[j] else 0)


def run():
    cases=set()
    for n in range(1,6):
        for r in range(n+1):
            cases.add(tuple([0]*(n-r)+list(range(1,r+1))))
            cases.add(tuple([0]*(n-r)+[2]*r))
            cases.add(tuple([0]*(n-r)+[1+2*(j%2) for j in range(r)]))
    ranks=[]
    for h in sorted(cases,key=lambda x:(len(x),x)):
        n=len(h);d=s.diag(*h);p0=s.diag(*(int(v==0) for v in h));sb=symbasis(n)
        cols=[e*d for e in sb]
        m=stack([flatten(v) for v in cols])
        ms=stack([flatten((v+v.T)/2) for v in cols])
        ma=stack([flatten((v-v.T)/2) for v in cols])
        r=sum(v>0 for v in h);expected=r*(2*n-r+1)//2
        expected_a=sum(h[i]!=h[j] for i in range(n) for j in range(i+1,n))
        check('full seed-image rank',m.rank(),expected)
        check('symmetric-channel rank',ms.rank(),expected)
        check('inherited antisymmetric spectral rank',ma.rank(),expected_a)
        # Test the stated image characterization on the *whole* unrestricted
        # n-by-n coefficient space, including its kernel-kernel block.
        constraints=[]
        for i in range(n):
            for j in range(n):
                e=s.zeros(n);e[i,j]=1
                constraints.append(flatten(d*e-e.T*d).col_join(flatten(e*p0)))
        t=stack(constraints)
        check('image constraints annihilate source map',t*m,s.zeros(2*n*n,len(sb)))
        check('image constraint dimension is sufficient',n*n-t.rank(),expected)
        c=symmetric(n,1);k=c*d;sm=(k+k.T)/2;am=(k-k.T)/2
        check('symmetric-channel reconstruction',reconstructed_a(sm,h),am)
        check('Sylvester source equation',d*am+am*d,sm*d-d*sm)
        check('null-block source support',k*p0,s.zeros(n))
        # Conversely start with an arbitrary symmetric channel with supported
        # slots, then recover a symmetric seed representative.
        ss=symmetric(n,3)
        for i in range(n):
            for j in range(n):
                if not h[i] and not h[j]:ss[i,j]=0
        aa=reconstructed_a(ss,h);kk=ss+aa
        cc=s.Matrix(n,n,lambda i,j:kk[i,j]/h[j] if h[j] else (kk[j,i]/h[i] if h[i] else 0))
        check('arbitrary supported symmetric channel has symmetric seed',cc,cc.T)
        check('supported symmetric-channel surjectivity',cc*d,kk)
        check('supported channel reproduces S',(kk+kk.T)/2,ss)
        # Exact complex unitary congruence tests the required transpose,
        # projector placement, and basis-invariant Frobenius statement.
        u=rotation(n);check('complex unitary rotation',u.conjugate().T*u,s.eye(n))
        ch=u*c*u.T;hh=u*d*u.conjugate().T;kh=ch*hh.T
        sh=(kh+kh.T)/2;ah=(kh-kh.T)/2;pp=u*p0*u.conjugate().T
        check('source congruence',kh,u*k*u.T)
        check('global flavor image equation',hh*kh,kh.T*hh.T)
        check('global null-projector constraint',kh*pp.T,s.zeros(n))
        check('global Sylvester reconstruction',hh*ah+ah*hh.T,sh*hh.T-hh*sh)
        check('symmetric Frobenius congruence',snorm(sh),snorm(sm))
        check('antisymmetric Frobenius congruence',snorm(ah),snorm(am))
        eta=s.Rational(max(h)-min(h),max(h)+min(h)) if max(h) else s.Integer(0)
        nonnegative('relative source-channel norm bound',eta**2*snorm(sh)-snorm(ah))
        nonnegative('unit source-channel norm bound',snorm(sh)-snorm(ah))
        if n==2 and h==(1,2):
            reject('transpose omitted from complex H',ch*hh,kh)
            reject('Hermitian instead of transpose source projection',(kh+kh.conjugate().T)/2,sh)
            reject('antisymmetric weak source forced to zero',ah,s.zeros(n))
            reject('symmetric weak source discarded',sh,s.zeros(n))
            reject('missing projection half',kh+kh.T,sh)
        if n==2 and h==(0,1):
            reject('missing null-projector transpose',kh*pp,s.zeros(n))
        ranks.append({'n':n,'eigenvalues':list(h),'rank_H':r,'seed_image_rank':expected,'antisymmetric_image_rank':expected_a})
    # A zero eigenspace requires the projector condition in addition to HK=K^T H^T.
    h=s.diag(0,2);bad=s.diag(1,0);p0=s.diag(1,0)
    check('singular-image incomplete condition witness',h*bad,bad.T*h.T)
    reject('omit null-projector image constraint',bad*p0.T,s.zeros(2))
    reject('rank-deficient source treated as full symmetric seed',3,6)
    # Positivity is essential: an indefinite Hermitian matrix can have a
    # nonzero antisymmetric source while its symmetric channel is zero.
    indefinite=s.diag(1,-1);off=s.Matrix([[0,1],[1,0]])*indefinite
    check('indefinite premise witness symmetric channel',off+off.T,s.zeros(2))
    reject('extend PSD reconstruction to indefinite H',off-off.T,s.zeros(2))
    # The bound is sharp: an off-diagonal seed on the extremal eigenspaces.
    for h in ((1,5),(0,3)):
        c=s.Matrix([[0,1+s.I],[1+s.I,0]]);k=c*s.diag(*h)
        ss=(k+k.T)/2;aa=(k-k.T)/2
        eta=s.Rational(max(h)-min(h),sum(h))
        check('sharp relative source-channel bound',snorm(aa),eta**2*snorm(ss))
    # Check the actual conjugate-Yukawa ordering, not only an abstract PSD H.
    y=s.Matrix([[1+s.I,2],[3-s.I,s.I],[-2,1-2*s.I]])
    c=symmetric(3,2);f=c*y.conjugate();h=y*y.conjugate().T;k=f*y.T
    check('actual source Yukawa ordering',k,c*h.T)
    check('actual source Hermitian PSD matrix',h,h.conjugate().T)
    check('actual source image equation',h*k,k.T*h.T)
    reject('unconjugated Yukawa source',c*y*y.T,k)
    v=rotation(2);yr=y*v.conjugate().T
    check('sterile unitary invariance of source',(c*yr.conjugate())*yr.T,k)
    result={'status':'PASS','input_commit':'34acf849172231282a0a6e7a07a81157b2d48799',
        'scope':'Exact selected seed-source image, supported symmetric weak channel, and relative algebraic Frobenius bound. Not a physical dipole bound, anomalous dimension, or full EFT basis.',
        'positive_checks':sum(counts.values()),'positive_categories':counts,
        'negative_checks':sum(negatives.values()),'negative_categories':negatives,
        'spectra_checked':ranks,'repository_imports':False,'assertions_required':False}
    Path(__file__).with_name('SEED_FLAVOR_IMAGE_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':run()
