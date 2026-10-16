#!/usr/bin/env python3
"""Exhaustive algebraic basis controls for 2L 2phi strength tensors.

Standard library only; exact complex rational arithmetic, polynomial weak
covariance, and complete chiral Grassmann monomial coefficient expansion.
No sampled fields, equation-of-motion reduction, or repository imports.
"""
from collections import Counter
from fractions import Fraction as R
from itertools import combinations, product
import json
from pathlib import Path
import sys


class Q:
    __slots__=('a','b')
    def __init__(self,a=0,b=0):
        if isinstance(a,Q):self.a,self.b=a.a,a.b
        else:self.a,self.b=R(a),R(b)
    def __add__(self,o):o=Q(o);return Q(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self):return Q(-self.a,-self.b)
    def __sub__(self,o):return self+-Q(o)
    def __rsub__(self,o):return Q(o)+-self
    def __mul__(self,o):
        o=Q(o);return Q(self.a*o.a-self.b*o.b,self.a*o.b+self.b*o.a)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=Q(o);d=o.a*o.a+o.b*o.b
        return Q((self.a*o.a+self.b*o.b)/d,(self.b*o.a-self.a*o.b)/d)
    def __eq__(self,o):o=Q(o);return self.a==o.a and self.b==o.b
    def __bool__(self):return bool(self.a or self.b)
    def conj(self):return Q(self.a,-self.b)
    def serial(self):return [str(self.a),str(self.b)]
    def __repr__(self):return 'Q(%s,%s)'%(self.a,self.b)


I=Q(0,1);eta=(1,-1,-1,-1)
def matrix(n,m,f):return tuple(tuple(Q(f(a,b)) for b in range(m)) for a in range(n))
def identity(n):return matrix(n,n,lambda a,b:a==b)
def transpose(a):return tuple(zip(*a))
def scale(c,a):return tuple(tuple(c*x for x in row) for row in a)
def add(a,b):return tuple(tuple(x+y for x,y in zip(ar,br)) for ar,br in zip(a,b))
def sub(a,b):return add(a,scale(-1,b))
def mul(a,b):return tuple(tuple(sum((a[i][k]*b[k][j] for k in range(len(b))),Q()) for j in range(len(b[0]))) for i in range(len(a)))
def flat(a):return tuple(x for row in a for x in row)
def eps(a,b,c):
    if len({a,b,c})<3:return 0
    return 1 if (a,b,c) in ((0,1,2),(1,2,0),(2,0,1)) else -1
pauli=(matrix(2,2,lambda a,b:a!=b),((Q(),-I),(I,Q())),((Q(1),Q()),(Q(),Q(-1))))
T=tuple(scale(R(1,2),p) for p in pauli)
g0=matrix(4,4,lambda a,b:(1 if a<2 else -1) if a==b else 0)
def spatial(p):return tuple(tuple(p[a][b-2] if a<2 and b>=2 else -p[a-2][b] if a>=2 and b<2 else Q() for b in range(4)) for a in range(4))
gamma=(g0,)+tuple(spatial(p) for p in pauli)
sigma=tuple(tuple(scale(I/2,sub(mul(gamma[m],gamma[n]),mul(gamma[n],gamma[m]))) for n in range(4)) for m in range(4))
g5=scale(I,mul(mul(mul(gamma[0],gamma[1]),gamma[2]),gamma[3]))
PR=scale(R(1,2),add(identity(4),g5))
C=scale(I,mul(gamma[2],gamma[0]))
# A full basis of the right-chiral row/column space; its normalization is
# irrelevant to rank and used identically for every coefficient map.
U=((Q(1),Q()),(Q(),Q(1)),(Q(1),Q()),(Q(),Q(1)))
pairs=tuple(combinations(range(4),2))
spin=tuple(mul(mul(transpose(U),mul(sigma[m][n],C)),U) for m,n in pairs)
massspin=tuple(mul(mul(transpose(U),sigma[m][n]),U) for m,n in pairs)


counts=Counter();mutations=set()
def check(label,a,b):
    if a!=b:raise RuntimeError('check failed: '+label+'; '+repr(a)+' != '+repr(b))
    counts[label]+=1
def reject(label,a,b):
    if a==b:raise RuntimeError('mutation survived: '+label)
    mutations.add(label)
def rref(rows,ncols=None):
    if ncols is None:ncols=len(rows[0]) if rows else 0
    a=[list(row) for row in rows];piv=[];r=0
    for c in range(ncols):
        pick=next((i for i in range(r,len(a)) if a[i][c]),None)
        if pick is None:continue
        a[r],a[pick]=a[pick],a[r];v=a[r][c]
        a[r]=[x/v for x in a[r]]
        for i in range(len(a)):
            if i!=r and a[i][c]:
                f=a[i][c];a[i]=[x-f*y for x,y in zip(a[i],a[r])]
        piv.append(c);r+=1
        if r==len(a):break
    return a,piv
def rank(rows,ncols=None):return len(rref(rows,ncols)[1])
def nullspace(rows,ncols):
    a,piv=rref(rows,ncols);out=[]
    for free in range(ncols):
        if free in piv:continue
        v=[Q() for j in range(ncols)];v[free]=Q(1)
        for r,c in enumerate(piv):v[c]=-a[r][free]
        out.append(tuple(v))
    return out
def mvec(a,v):return tuple(sum((x*y for x,y in zip(row,v)),Q()) for row in a)


def covariance_system(degree,legs,weak,adjoint_sign=1):
    mons=((1,0),(0,1)) if degree==1 else ((2,0),(1,1),(0,2))
    strengths=range(3) if weak else range(1)
    outputs=tuple(product(range(2),repeat=legs))
    keys=tuple((j,ab,q) for j in strengths for ab in outputs for q in range(len(mons)))
    index={k:i for i,k in enumerate(keys)}
    rowkeys=tuple((a,j,ab,q) for a in range(3) for j in strengths for ab in outputs for q in range(len(mons)))
    rows=[[Q() for k in keys] for row in rowkeys];ri={k:i for i,k in enumerate(rowkeys)}
    for c,(j,ab,q) in enumerate(keys):
        for a in range(3):
            mon=mons[q]
            for r in range(2):
                if mon[r]:
                    for s in range(2):
                        mm=list(mon);mm[r]-=1;mm[s]+=1;qq=mons.index(tuple(mm))
                        rows[ri[a,j,ab,qq]][c]+=mon[r]*T[a][r][s]
            for slot in range(legs):
                for out in range(2):
                    new=list(ab);new[slot]=out
                    rows[ri[a,j,tuple(new),q]][c]-=T[a][out][ab[slot]]
            if weak:
                for jout in range(3):
                    rows[ri[a,jout,ab,q]][c]+=adjoint_sign*I*eps(a,jout,j)
    return mons,keys,tuple(tuple(row) for row in rows)


def polynomial_monomial(a,b=None):
    mon=[0,0];mon[a]+=1
    if b is not None:mon[b]+=1
    return tuple(mon)
def candidate(degree,legs,weak,kind='raw'):
    mons,keys,rows=covariance_system(degree,legs,weak)
    out=[]
    for j,ab,q in keys:
        if legs==1:
            v=(Q(int(mons[q]==polynomial_monomial(ab[0]))) if not weak else sum((T[j][ab[0]][c]*int(mons[q]==polynomial_monomial(c)) for c in range(2)),Q()))
        elif not weak:v=Q(int(mons[q]==polynomial_monomial(ab[0],ab[1])))
        else:
            a,b=ab
            p=sum((T[j][a][c]*int(mons[q]==polynomial_monomial(c,b)) for c in range(2)),Q())
            pt=sum((T[j][b][c]*int(mons[q]==polynomial_monomial(c,a)) for c in range(2)),Q())
            v=p if kind=='raw' else (p+pt)/2 if kind=='plus' else (p-pt)/2
        out.append(v)
    return tuple(out)


weak_receipts=[]
for degree,legs,name,expected in ((2,2,'active_B',1),(2,2,'active_W',2),(1,1,'mass_B',1),(1,1,'mass_W',1)):
    weak=name.endswith('W');mons,keys,rows=covariance_system(degree,legs,weak)
    ker=nullspace(rows,len(keys));check('complete weak invariant nullity',len(ker),expected)
    candidates=[candidate(degree,legs,weak,k) for k in ('plus','minus')] if degree==2 and weak else [candidate(degree,legs,weak)]
    for v in candidates:check('explicit weak invariant solves every covariance coefficient',mvec(rows,v),tuple(Q() for row in rows))
    check('explicit weak invariants span full nullspace',rank(candidates,len(keys)),expected)
    weak_receipts.append(dict(sector=name,unknown_coefficients=len(keys),covariance_equations=len(rows),rank=len(keys)-len(ker),nullity=len(ker),nullspace_basis=[[[str(x.a),str(x.b)] for x in v] for v in ker]))
    if name=='active_W':
        wrong=covariance_system(degree,legs,weak,-1)[2]
        reject('opposite adjoint transformation sign',mvec(wrong,candidates[0]),tuple(Q() for row in wrong))
        reject('omit the second weak invariant',rank(candidates[:1],len(keys)),expected)


check('right chiral embedding',mul(PR,U),U)
for ss in spin:check('symmetric chiral sigma C',transpose(ss),ss)
check('chiral strength spin map has rank three',rank([flat(s) for s in spin],4),3)
def eps4(a,b,c,d):
    if len({a,b,c,d})!=4:return 0
    v=(a,b,c,d);return (-1)**sum(v[i]>v[j] for i in range(4) for j in range(i+1,4))
dualspin=[]
for r,s in pairs:
    # F_rs=1,F_sr=-1.  The dual lower-index tensor is
    # eta_m eta_n epsilon^{mnrs}; sum all m,n gives the factor two.
    dual=matrix(2,2,lambda a,b:sum((2*eta[m]*eta[n]*eps4(m,n,r,s)*spin[k][a][b] for k,(m,n) in enumerate(pairs)),Q()))
    ordinary=scale(2,spin[pairs.index((r,s))])
    check('dual chiral strength is minus i times ordinary',dual,scale(-I,ordinary))
    dualspin.append(dual)
check('dual strengths add no independent Lorentz tensor',rank([flat(s) for s in spin]+[flat(s) for s in dualspin],4),3)
massdual=[]
for r,s in pairs:
    dual=matrix(2,2,lambda a,b:sum((2*eta[m]*eta[n]*eps4(m,n,r,s)*massspin[k][a][b] for k,(m,n) in enumerate(pairs)),Q()))
    check('retained-mass dual strength is minus i times ordinary',dual,scale(-2*I,massspin[pairs.index((r,s))]))
    massdual.append(dual)
check('retained-mass chiral strength spin map rank',rank([flat(s) for s in massspin],4),3)
check('retained-mass dual strengths add no Lorentz direction',rank([flat(s) for s in massspin]+[flat(s) for s in massdual],4),3)


def clean(d):return {k:v for k,v in d.items() if v}
def padd(a,b):
    d=dict(a)
    for k,v in b.items():d[k]=d.get(k,Q())+v
    return clean(d)
def pscale(c,a):return clean({k:c*v for k,v in a.items()})
def psum(xs):
    d={}
    for x in xs:d=padd(d,x)
    return d
def expansion(n,al,be,weak,kind='raw'):
    mons,keys,rows=covariance_system(2,2,weak);vec=candidate(2,2,weak,kind);out={}
    for (j,(a,b),q),v in zip(keys,vec):
        if not v:continue
        for k,ss in enumerate(spin):
            for s,t in product(range(2),repeat=2):
                i=4*al+2*a+s;ii=4*be+2*b+t
                if i==ii:continue
                lo,hi=sorted((i,ii));sign=1 if i<ii else -1
                key=('W' if weak else 'B',j,k,q,lo,hi)
                out[key]=out.get(key,Q())+2*sign*v*ss[s][t]
    return clean(out)
def sparse_rank(columns):
    keys=sorted(set().union(*(set(c) for c in columns))) if columns else []
    return rank([[col.get(k,Q()) for col in columns] for k in keys],len(columns))


flavor_receipts=[]
for n in (1,2,3):
    bops={(a,b):expansion(n,a,b,False) for a,b in product(range(n),repeat=2)}
    wops={(a,b):expansion(n,a,b,True) for a,b in product(range(n),repeat=2)}
    wplus={(a,b):expansion(n,a,b,True,'plus') for a,b in product(range(n),repeat=2)}
    wminus={(a,b):expansion(n,a,b,True,'minus') for a,b in product(range(n),repeat=2)}
    independent_b=[bops[a,b] for a in range(n) for b in range(a+1,n)]
    independent_w=[wplus[a,b] for a in range(n) for b in range(a+1,n)]+[wminus[a,b] for a in range(n) for b in range(a,n)]
    brank=sparse_rank(list(bops.values()));wrank=sparse_rank(list(wops.values()))
    check('exhaustive active B Grassmann rank',brank,n*(n-1)//2)
    check('exhaustive active W Grassmann rank',wrank,n*n)
    check('independent B pair basis rank',sparse_rank(independent_b),n*(n-1)//2)
    check('independent W symmetric and antisymmetric pair basis rank',sparse_rank(independent_w),n*n)
    check('combined B and W basis rank',sparse_rank(independent_b+independent_w),n*(n-1)//2+n*n)
    for a,b in product(range(n),repeat=2):
        check('all B flavor exchange identities',bops[b,a],pscale(-1,bops[a,b]))
        check('all W plus flavor exchange identities',wplus[b,a],pscale(-1,wplus[a,b]))
        check('all W minus flavor exchange identities',wminus[b,a],wminus[a,b])
        check('raw W splits into two weak tensors',wops[a,b],padd(wplus[a,b],wminus[a,b]))
    # Matrix units are a complete coefficient-space basis. These identities
    # establish the source map for every complex K by linearity.
    for ka,kb in product(range(n),repeat=2):
        k=lambda a,b:Q(int(a==ka and b==kb))
        mapped_b=psum(pscale(k(a,b)-k(b,a),bops[a,b]) for a in range(n) for b in range(a+1,n))
        mapped_w=psum([pscale(k(a,b)-k(b,a),wplus[a,b]) for a in range(n) for b in range(a+1,n)]+[pscale(k(a,b)+k(b,a),wminus[a,b]) for a in range(n) for b in range(a+1,n)]+[pscale(k(a,a),wminus[a,a]) for a in range(n)])
        check('exhaustive independent-pair B source map',pscale(R(-1,12),bops[ka,kb]),pscale(R(-1,12),mapped_b))
        check('exhaustive independent-pair W source map',pscale(R(-1,12),wops[ka,kb]),pscale(R(-1,12),mapped_w))
    # Exhibit the complete B kernel, not just some vanishing combinations.
    symvectors=[]
    for a in range(n):
        for b in range(a,n):
            v=[Q(int((i==a and j==b) or (i==b and j==a))) for i,j in product(range(n),repeat=2)]
            symvectors.append(v)
            check('symmetric flavor coefficient lies in full B kernel',psum(pscale(v[i*n+j],bops[i,j]) for i,j in product(range(n),repeat=2)),{})
    check('B symmetric kernel is complete',rank(symvectors,n*n),n*n-brank)
    flavor_receipts.append(dict(active_flavors=n,ordered_B_columns=n*n,B_rank=brank,B_kernel_dimension=n*n-brank,ordered_W_columns=n*n,W_rank=wrank,W_kernel_dimension=n*n-wrank,independent_pair_W_plus=n*(n-1)//2,independent_pair_W_minus=n*(n+1)//2))
    if n==1:
        reject('flavor-only antisymmetry removes charged one-flavor W',wops[0,0],{})
        reject('double the diagonal symmetric-flavor term',wops[0,0],pscale(2,wminus[0,0]))
        reject('Pauli normalization silently replaces T',wops[0,0],pscale(2,wops[0,0]))
    if n==2:
        reject('omit independent B pair factor two',pscale(2,bops[0,1]),bops[0,1])
        reject('omit independent W offdiagonal symmetric factor two',padd(wops[0,1],wops[1,0]),wminus[0,1])
        reject('transpose ordered source flavor matrix',wops[0,1],wops[1,0])


def mass_expansion(n,ns,al,sterile,weak):
    mons,keys,rows=covariance_system(1,1,weak);vec=candidate(1,1,weak);out={}
    for (j,(a,),q),v in zip(keys,vec):
        for k,ss in enumerate(massspin):
            for s,t in product(range(2),repeat=2):
                left=4*al+2*a+s;right=4*n+2*sterile+t
                key=('W' if weak else 'B',j,k,q,left,right)
                out[key]=out.get(key,Q())+2*v*ss[s][t]
    return clean(out)
mass_receipts=[]
for n,ns in product((1,2,3),(1,2)):
    bb=[mass_expansion(n,ns,a,i,False) for a in range(n) for i in range(ns)]
    ww=[mass_expansion(n,ns,a,i,True) for a in range(n) for i in range(ns)]
    check('retained-mass B sector independent rank',sparse_rank(bb),n*ns)
    check('retained-mass W sector independent rank',sparse_rank(ww),n*ns)
    check('retained-mass combined sector rank',sparse_rank(bb+ww),2*n*ns)
    mass_receipts.append(dict(active_flavors=n,sterile_flavors=ns,B_rank=n*ns,W_rank=n*ns,combined_rank=2*n*ns))


result={
    'status':'passed','python':sys.version.split()[0],
    'input_candidate':'34acf849172231282a0a6e7a07a81157b2d48799',
    'arithmetic':'exact Gaussian rational arithmetic; standard library only',
    'method':'Complete polynomial weak covariance nullspaces and chiral Grassmann monomial coefficient ranks; no sampled field values.',
    'positive_check_count':sum(counts.values()),'checks_by_category':dict(sorted(counts.items())),
    'negative_control_count':len(mutations),'negative_controls_rejected':sorted(mutations),
    'weak_invariant_systems':weak_receipts,'active_flavor_rank_results':flavor_receipts,'retained_mass_rank_results':mass_receipts,
    'scope':'Algebraic independence and coefficient transport inside the stated no-derivative field-content sectors; no full EFT/EOM basis closure, physical Wilson prediction, finite amplitude or beta function.'
}
Path(__file__).with_name('ACTIVE_TENSOR_BASIS_RESULTS.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='weak_invariant_systems'},indent=2,sort_keys=True))
