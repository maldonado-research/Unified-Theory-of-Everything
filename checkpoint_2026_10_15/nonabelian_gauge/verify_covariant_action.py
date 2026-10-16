#!/usr/bin/env python3
"""Independent exact Fourier differentiation of the D-plus matrix action.

No project imports, loop formulas, symbolic algebra package, or floating point.
Coefficients of distinct formal gauge-mode variables are literal functional
derivatives: every field assignment occurs in polynomial multiplication.
"""
from fractions import Fraction as F
from itertools import permutations
import json
from pathlib import Path


class Q:
    __slots__ = ('a', 'b')
    def __init__(self, a=0, b=0):
        if isinstance(a, Q): self.a, self.b = a.a, a.b
        else: self.a, self.b = F(a), F(b)
    def __add__(self, o):
        o = Q(o); return Q(self.a+o.a, self.b+o.b)
    __radd__ = __add__
    def __neg__(self): return Q(-self.a, -self.b)
    def __sub__(self, o): return self + -Q(o)
    def __rsub__(self, o): return Q(o) + -self
    def __mul__(self, o):
        o = Q(o); return Q(self.a*o.a-self.b*o.b, self.a*o.b+self.b*o.a)
    __rmul__ = __mul__
    def __eq__(self, o):
        o = Q(o); return self.a == o.a and self.b == o.b
    def __bool__(self): return bool(self.a or self.b)
    def __repr__(self): return 'Q(%s,%s)' % (self.a,self.b)


I = Q(0, 1)
eta = (1, -1, -1, -1)
def eye(n): return tuple(tuple(Q(i == j) for j in range(n)) for i in range(n))
def zero(n): return tuple(tuple(Q() for j in range(n)) for i in range(n))
def add(a,b): return tuple(tuple(x+y for x,y in zip(ar,br)) for ar,br in zip(a,b))
def neg(a): return scale(-1,a)
def sub(a,b): return add(a,neg(b))
def scale(c,a): return tuple(tuple(c*x for x in ar) for ar in a)
def mul(a,b):
    n=len(a)
    return tuple(tuple(sum((a[i][k]*b[k][j] for k in range(n) if a[i][k] and b[k][j]), Q()) for j in range(n)) for i in range(n))
def comm(a,b): return sub(mul(a,b),mul(b,a))
def anti(a,b): return add(mul(a,b),mul(b,a))
def msum(terms,n):
    out=zero(n)
    for t in terms: out=add(out,t)
    return out
def kron(a,b):
    return tuple(tuple(a[i//len(b)][j//len(b)]*b[i%len(b)][j%len(b)] for j in range(len(a)*len(b))) for i in range(len(a)*len(b)))
def vecadd(*vs): return tuple(sum(v[i] for v in vs) for i in range(4))
def vecscale(c,v): return tuple(c*x for x in v)
def dot(a,b): return sum(eta[i]*a[i]*b[i] for i in range(4))

id2, id4, id8 = eye(2), eye(4), eye(8)
p1=((Q(),Q(1)),(Q(1),Q()))
p2=((Q(),-I),(I,Q()))
p3=((Q(1),Q()),(Q(),Q(-1)))
pauli=(p1,p2,p3)
g0=tuple(tuple(Q((1 if i<2 else -1) if i==j else 0) for j in range(4)) for i in range(4))
def spatial(s):
    return tuple(tuple(s[i][j-2] if i<2 and j>=2 else -s[i-2][j] if i>=2 and j<2 else Q() for j in range(4)) for i in range(4))
gamma=(g0,)+(tuple(spatial(s) for s in pauli))
def slash(v): return msum((scale(eta[i]*v[i],gamma[i]) for i in range(4)),4)
sigma=tuple(tuple(scale(I*F(1,2),comm(gamma[i],gamma[j])) for j in range(4)) for i in range(4))
g5=scale(I,mul(mul(mul(gamma[0],gamma[1]),gamma[2]),gamma[3]))
pl=scale(F(1,2),sub(id4,g5))


# A polynomial maps a subset of labelled gauge modes to a weak-space matrix.
# Products keep square-free coefficients only; discarded repeated modes cannot
# contribute to a functional derivative that uses each external gauge once.
def padd(a,b):
    out=dict(a)
    for k,v in b.items(): out[k]=add(out.get(k,zero(2)),v)
    return {k:v for k,v in out.items() if v!=zero(2)}
def pscale(c,a): return {k:scale(c,v) for k,v in a.items()}
def psub(a,b): return padd(a,pscale(-1,b))
def pmul(a,b):
    out={}
    for ka,va in a.items():
        for kb,vb in b.items():
            if ka&kb: continue
            key=ka|kb
            out[key]=add(out.get(key,zero(2)),mul(va,vb))
    return out
def pd(a,base,ks,mu):
    return {mask:scale(I*eta[mu]*(base[mu]+sum(k[mu] for j,k in enumerate(ks) if mask&(1<<j))),v) for mask,v in a.items()}
def component(a,mask): return a.get(mask,zero(2))


def direct(r,h,ks,es,xs,ms=F(2,3), curvature_sign=1,
           row_order=False, adjoint_sign=1):
    """Differentiate the action itself, before simplifying Clifford/group algebra."""
    n=len(ks); mask=(1<<n)-1; p=vecadd(r,h,*ks); z=(0,0,0,0)
    aa=[{1<<j:scale(eta[mu]*es[j][mu],xs[j]) for j in range(n)} for mu in range(4)]
    def row(a,mu): return psub(pd(a,r,ks,mu),pscale(I,pmul(a,aa[mu])))
    def scalar(a,mu): return padd(pd(a,h,ks,mu),pscale(I,pmul(aa[mu],a)))
    unit={0:id2}; row1=[row(unit,mu) for mu in range(4)]
    minus_row2={}; minus_scalar2={}
    for mu in range(4):
        minus_row2=padd(minus_row2,pscale(-eta[mu],row(row1[mu],mu)))
        minus_scalar2=padd(minus_scalar2,pscale(-eta[mu],scalar(scalar(unit,mu),mu)))
    # Literal six placements of O_J, without importing its gauge kernels.
    oj=zero(8)
    if mask==0: oj=add(oj,kron(scale(F(5,6)*dot(p,p)+ms/2,id2),slash(p)))
    oj=add(oj,kron(scale(F(1,3),component(minus_row2,mask)),slash(p)))
    oj=add(oj,kron(scale(F(1,6),component(minus_scalar2,mask)),slash(p)))
    for rho in range(4):
        krow=scale(-I,component(row1[rho],mask))
        oj=add(oj,kron(scale(-dot(p,p)/6+ms/2,krow),gamma[rho]))
    # F = dA + i[A,A]; both index orders are retained explicitly.
    ff=[[padd(psub(pd(aa[nu],z,ks,mu),pd(aa[mu],z,ks,nu)),
                  pscale(I*curvature_sign,psub(pmul(aa[mu],aa[nu]),pmul(aa[nu],aa[mu]))))
         for nu in range(4)] for mu in range(4)]
    sig=zero(8); div=zero(8)
    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                prod=pmul(ff[mu][nu],row1[rho]) if row_order else pmul(row1[rho],ff[mu][nu])
                sig=add(sig,kron(scale(I*F(1,4),component(prod,mask)),mul(sigma[mu][nu],gamma[rho])))
    for mu in range(4):
        for nu in range(4):
            cov=padd(pd(ff[nu][mu],z,ks,nu),pscale(I*adjoint_sign,psub(pmul(aa[nu],ff[nu][mu]),pmul(ff[nu][mu],aa[nu]))))
            div=add(div,kron(scale(eta[nu]*F(1,3),component(cov,mask)),gamma[mu]))
    return add(add(oj,sig),div), (oj,sig,div)


def r1(r,h,k,e,ms=F(2,3)):
    p=vecadd(r,h,k)
    j=add(scale(-dot(vecadd(vecscale(2,r),k),e)/3+dot(vecadd(vecscale(2,h),k),e)/6,slash(p)),
          scale(dot(p,p)/6-ms/2,slash(e)))
    s=scale(F(1,4),mul(comm(slash(k),slash(e)),slash(r)))
    d=scale(F(1,3),sub(scale(dot(k,e),slash(k)),scale(dot(k,k),slash(e))))
    return add(add(j,s),d)


def k2(r,h,ks,es,xs):
    k,l=ks; e,f=es; x,y=xs; p=vecadd(r,h,k,l)
    j=kron(scale(dot(e,f)/2,anti(x,y)),slash(p))
    s=add(kron(comm(x,y),mul(comm(slash(e),slash(f)),slash(r))),
          neg(add(kron(mul(y,x),mul(comm(slash(k),slash(e)),slash(f))),
                  kron(mul(x,y),mul(comm(slash(l),slash(f)),slash(e))))))
    v=vecadd(vecscale(dot(vecadd(k,vecscale(2,l)),e),f),
             vecscale(-dot(vecadd(vecscale(2,k),l),f),e),
             vecscale(dot(e,f),vecadd(k,vecscale(-1,l))))
    return add(add(j,scale(F(1,4),s)),kron(scale(F(-1,3),comm(x,y)),slash(v)))


def k3(es,xs):
    out=zero(8)
    for i,j,l in ((0,1,2),(0,2,1),(1,2,0)):
        ei,ej,el=es[i],es[j],es[l]; xi,xj,xl=xs[i],xs[j],xs[l]
        out=add(out,kron(scale(F(-1,4),mul(xl,comm(xi,xj))),mul(comm(slash(ei),slash(ej)),slash(el))))
        v=vecadd(vecscale(dot(el,ei),ej),vecscale(-dot(el,ej),ei))
        out=add(out,kron(scale(F(-1,3),comm(xl,comm(xi,xj))),slash(v)))
    return out


checks=0; rejected={}
def check(label,a,b):
    global checks
    if a!=b: raise AssertionError(label)
    checks+=1
def bad(label,a,b):
    global checks
    if a==b: raise AssertionError('negative control survived: '+label)
    rejected[label]=rejected.get(label,0)+1
    checks+=1
def rndvec(s,shift):
    return tuple(F(((s+2)*(i+3)+shift*(i*i+2))%13-6, ((s+i+shift)%3)+1) for i in range(4))


def run():
    for mu in range(4):
        for nu in range(4):
            check('Clifford',anti(gamma[mu],gamma[nu]),scale(2*eta[mu] if mu==nu else 0,id4))
    check('chirality projector',mul(pl,pl),pl)
    for a in range(3):
        for b in range(3):
            check('Pauli anticommutator',anti(pauli[a],pauli[b]),scale(2 if a==b else 0,id2))
    # A square-free two-mode action has momentum degree one. Exhaust all
    # 16 momentum-component coefficients and 16 polarization-bilinear
    # coefficients; E12,E21 give independent ordered products E11,E22.
    basis=tuple(tuple(F(i==j) for i in range(4)) for j in range(4))
    z=(F(0),)*4
    e12=((Q(),Q(1)),(Q(),Q()))
    e21=((Q(),Q()),(Q(1),Q()))
    exhaustive_two=0
    for slot in range(4):
        for mom in basis:
            mm=[z,z,z,z]; mm[slot]=mom
            for e in basis:
                for f in basis:
                    got,_=direct(mm[0],mm[1],tuple(mm[2:]),(e,f),(e12,e21))
                    check('exhaustive two-gauge coefficient',got,k2(mm[0],mm[1],tuple(mm[2:]),(e,f),(e12,e21)))
                    exhaustive_two+=1
    # Cubic kernels have momentum degree zero, and are trilinear in the
    # polarization vectors. Non-orthogonal generators retain the nested
    # commutators in the adjoint-curvature term.
    cubic_x=(add(p1,p2),add(p2,p3),add(p3,p1))
    exhaustive_three=0
    for e in basis:
        for f in basis:
            for g in basis:
                got,_=direct(z,z,(z,z,z),(e,f,g),cubic_x)
                check('exhaustive three-gauge polarization',got,k3((e,f,g),cubic_x))
                exhaustive_three+=1
    cases=12
    selected_negative=0
    for s in range(cases):
        r=rndvec(s,1); h=rndvec(s,2)
        ks=tuple(rndvec(s,3+j) for j in range(3))
        es=tuple(rndvec(s,6+j) for j in range(3))
        xs=tuple(scale(F(1,2),add(pauli[j],scale(F(s+1,j+2),pauli[(j+1)%3]))) for j in range(3))
        if s%3==0: xs=(scale(F(-1,2),id2),xs[1],xs[2])
        # First, derive all action kernels from polynomial fields.
        one,_=direct(r,h,ks[:1],es[:1],xs[:1])
        check('literal one-gauge Fourier',one,kron(xs[0],r1(r,h,ks[0],es[0])))
        two,parts2=direct(r,h,ks[:2],es[:2],xs[:2])
        check('literal two-gauge Fourier',two,k2(r,h,ks[:2],es[:2],xs[:2]))
        three,parts3=direct(r,h,ks,es,xs)
        check('literal three-gauge Fourier',three,k3(es,xs))
        check('O_J has no cubic',parts3[0],zero(8))
        # All label permutations, and the Abelian and mixed specializations.
        for perm in permutations(range(3)):
            check('three-gauge Bose',k3(tuple(es[i] for i in perm),tuple(xs[i] for i in perm)),three)
        check('two-gauge Bose',k2(r,h,ks[:2][::-1],es[:2][::-1],xs[:2][::-1]),two)
        ax=scale(F(-1,2),id2)
        ae=dot(es[0],es[1])
        ar2=sub(scale(ae,slash(vecadd(r,h,*ks[:2]))),
                scale(F(1,4),add(mul(comm(slash(ks[0]),slash(es[0])),slash(es[1])),mul(comm(slash(ks[1]),slash(es[1])),slash(es[0])))))
        check('Abelian reduction',k2(r,h,ks[:2],es[:2],(ax,ax)),kron(mul(ax,ax),ar2))
        check('Abelian cubic',k3(es,(ax,ax,ax)),zero(8))
        mixed=scale(F(-1,4),kron(mul(ax,comm(xs[1],xs[2])),mul(comm(slash(es[1]),slash(es[2])),slash(es[0]))))
        check('BWW cubic',k3(es,(ax,xs[1],xs[2])),mixed)
        # Functional Ward recursion uses fixed incoming p on every term.
        for leg in range(2):
            other=1-leg; k=ks[leg]; x=xs[leg]; y=xs[other]; e=es[other]; q=ks[other]
            ws=list(es[:2]); ws[leg]=k
            lhs=k2(r,h,ks[:2],tuple(ws),xs[:2])
            rhs=add(sub(kron(mul(y,x),r1(r,vecadd(h,k),q,e)),kron(mul(x,y),r1(vecadd(r,k),h,q,e))),
                    kron(comm(x,y),r1(r,h,vecadd(k,q),e)))
            check('two-gauge non-Abelian Ward',lhs,rhs)
        for leg in range(3):
            others=[j for j in range(3) if j!=leg]; k=ks[leg]; x=xs[leg]
            q=tuple(ks[j] for j in others); ee=tuple(es[j] for j in others); xx=tuple(xs[j] for j in others)
            ws=list(es); ws[leg]=k
            lhs=k3(tuple(ws),xs)
            rhs=sub(mul(k2(r,vecadd(h,k),q,ee,xx),kron(x,id4)),mul(kron(x,id4),k2(vecadd(r,k),h,q,ee,xx)))
            for pos in range(2):
                q2=list(q); q2[pos]=vecadd(k,q[pos]); xx2=list(xx); xx2[pos]=comm(x,xx[pos])
                rhs=add(rhs,k2(r,h,tuple(q2),ee,tuple(xx2)))
            check('three-gauge non-Abelian Ward',lhs,rhs)
        check('Jacobi',add(add(comm(xs[0],comm(xs[1],xs[2])),comm(xs[1],comm(xs[2],xs[0]))),comm(xs[2],comm(xs[0],xs[1]))),zero(2))
        # Mutations are tested only where generic noncommutativity makes them visible.
        if s%3!=0:
            wrong,_=direct(r,h,ks[:2],es[:2],xs[:2],curvature_sign=-1)
            bad('wrong nonlinear curvature sign',wrong,two)
            wrong,_=direct(r,h,ks[:2],es[:2],xs[:2],row_order=True)
            bad('reversed row-curvature order',wrong,two)
            wrong,_=direct(r,h,ks[:2],es[:2],xs[:2],adjoint_sign=-1)
            bad('wrong adjoint derivative sign',wrong,two)
            # Incorrectly Abelianizing the group kills the commutator completion.
            sym=scale(F(1,2),add(two,k2(r,h,ks[:2],es[:2],xs[:2][::-1])))
            bad('dropped antisymmetric group sector',sym,two)
            bad('dropped cubic',zero(8),three)
            wrong3,_=direct(r,h,ks,es,xs,curvature_sign=-1)
            bad('wrong cubic curvature sign',wrong3,three)
            contracted=k2(r,h,ks[:2],(ks[0],es[1]),xs[:2])
            abelian_rhs=sub(kron(mul(xs[1],xs[0]),r1(r,vecadd(h,ks[0]),ks[1],es[1])),
                            kron(mul(xs[0],xs[1]),r1(vecadd(r,ks[0]),h,ks[1],es[1])))
            bad('omitted lower-point Ward commutator',contracted,abelian_rhs)
            selected_negative+=1
    result={'status':'PASS','checks':checks,'kinematic_points':cases,
            'exhaustive_two_gauge_momentum_polarization_coefficients':exhaustive_two,
            'exhaustive_three_gauge_polarization_coefficients':exhaustive_three,
            'negative_controls_rejected':sum(rejected.values()),'negative_control_counts':rejected,
            'method':'literal square-free Fourier polynomial differentiation, exact Q(i) matrices; independent explicit kernels and two-/three-gauge Ward recursion',
            'scope':'classical local action Fourier/group/Ward checks; no loop integration or physical-running calculation'}
    path=Path(__file__).with_name('COVARIANT_ACTION_CHECKS.json')
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__': run()
