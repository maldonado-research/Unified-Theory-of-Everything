#!/usr/bin/env python3
"""Exact local Fourier-operator substitution for the actual selected packet.

Standard library only. Matrix coefficients are Gaussian rationals; finite
Fourier sums are arbitrary off-shell test functions, not on-shell particles.
The fifth key component is a formal parameter, truncated only when requested.
No calculation modules or previously saved coefficients are imported.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import time


class Q:
    __slots__ = ('r', 'i')
    def __init__(self, r=0, i=0):
        if isinstance(r, Q): self.r, self.i = r.r, r.i
        else: self.r, self.i = F(r), F(i)
    def __add__(self, o):
        o=Q(o); return Q(self.r+o.r,self.i+o.i)
    __radd__=__add__
    def __neg__(self): return Q(-self.r,-self.i)
    def __sub__(self,o): return self+-Q(o)
    def __rsub__(self,o): return Q(o)+-self
    def __mul__(self,o):
        o=Q(o); return Q(self.r*o.r-self.i*o.i,self.r*o.i+self.i*o.r)
    __rmul__=__mul__
    def __bool__(self): return bool(self.r or self.i)
    def __eq__(self,o):
        o=Q(o); return self.r==o.r and self.i==o.i
    def star(self): return Q(self.r,-self.i)


I=Q(0,1)
ETA=(1,-1,-1,-1)
ORIGIN=(0,0,0,0,0)
LIMIT=1


class M:
    """Sparse rectangular matrix, with exact scalar coefficients."""
    def __init__(self,n,m,a=None):
        self.n,self.m=n,m
        self.a={k:Q(v) for k,v in (a or {}).items() if Q(v)}
    def __eq__(self,b): return self.n==b.n and self.m==b.m and self.a==b.a
    def __bool__(self): return bool(self.a)
    def add(self,b):
        if (self.n,self.m)!=(b.n,b.m): raise ValueError('matrix sum shape')
        out=dict(self.a)
        for k,v in b.a.items(): out[k]=out.get(k,Q())+v
        return M(self.n,self.m,out)
    def scale(self,c): return M(self.n,self.m,{k:c*v for k,v in self.a.items()})
    def mul(self,b):
        if self.m!=b.n: raise ValueError('matrix product shape')
        out={}; rows={}
        for (k,j),v in b.a.items(): rows.setdefault(k,[]).append((j,v))
        for (i,k),v in self.a.items():
            for j,w in rows.get(k,[]): out[i,j]=out.get((i,j),Q())+v*w
        return M(self.n,b.m,out)
    def transpose(self): return M(self.m,self.n,{(j,i):v for (i,j),v in self.a.items()})
    def star(self): return M(self.n,self.m,{k:v.star() for k,v in self.a.items()})
    def dagger(self): return self.star().transpose()


def eye(n): return M(n,n,{(i,i):1 for i in range(n)})
def kron(a,b):
    return M(a.n*b.n,a.m*b.m,{(i*b.n+k,j*b.m+l):v*w
              for (i,j),v in a.a.items() for (k,l),w in b.a.items()})
def dense(a): return M(len(a),len(a[0]),{(i,j):v for i,row in enumerate(a) for j,v in enumerate(row)})


class P:
    """Fourier functions valued in rectangular matrices."""
    def __init__(self,n,m,a=None):
        self.n,self.m=n,m
        self.a={k:v for k,v in (a or {}).items() if v and k[4]<=LIMIT}
    def __eq__(self,b): return self.n==b.n and self.m==b.m and self.a==b.a
    def __bool__(self): return bool(self.a)
    def add(self,b):
        if (self.n,self.m)!=(b.n,b.m): raise ValueError('Fourier sum shape')
        out=dict(self.a)
        for k,v in b.a.items(): out[k]=out.get(k,M(self.n,self.m)).add(v)
        return P(self.n,self.m,out)
    def scale(self,c): return P(self.n,self.m,{k:v.scale(c) for k,v in self.a.items()})
    def mul(self,b):
        if self.m!=b.n: raise ValueError('Fourier product shape')
        out={}
        for ka,va in self.a.items():
            for kb,vb in b.a.items():
                key=tuple(x+y for x,y in zip(ka,kb))
                if key[4]>LIMIT: continue
                out[key]=out.get(key,M(self.n,b.m)).add(va.mul(vb))
        return P(self.n,b.m,out)
    def deriv(self,mu): return P(self.n,self.m,{k:v.scale(I*ETA[mu]*k[mu]) for k,v in self.a.items()})
    def transpose(self): return P(self.m,self.n,{k:v.transpose() for k,v in self.a.items()})
    def star(self): return P(self.n,self.m,{tuple(-x for x in k[:4])+(k[4],):v.star() for k,v in self.a.items()})
    def dagger(self): return self.star().transpose()
    def degree(self,j): return P(self.n,self.m,{k[:4]+(0,):v for k,v in self.a.items() if k[4]==j})
    def parameter(self,j=1): return P(self.n,self.m,{k[:4]+(k[4]+j,):v for k,v in self.a.items()})


def const(m): return P(m.n,m.m,{ORIGIN:m})
def mode(k,m): return P(m.n,m.m,{tuple(k)+(0,):m})
def psum(terms,n,m):
    out=P(n,m)
    for t in terms: out=out.add(t)
    return out
def p0(a): return psum((a.deriv(mu).deriv(mu).scale(-ETA[mu]) for mu in range(4)),a.n,a.m)
def pp(a): return p0(a).add(a.scale(MH2))
def comm(a,b): return a.mul(b).add(b.mul(a).scale(-1))
def integral(a): return a.a.get(ORIGIN,M(a.n,a.m))
def integral_pair(a,b):
    """Zero Fourier coefficient directly, without irrelevant product modes."""
    if a.m!=b.n: raise ValueError('integrated product shape')
    out=M(a.n,b.m)
    for k,v in a.a.items():
        if k[4]: continue
        opposite=tuple(-x for x in k[:4])+(0,)
        if opposite in b.a: out=out.add(v.mul(b.a[opposite]))
    return out


pauli=(dense(((0,1),(1,0))),dense(((0,-I),(I,0))),dense(((1,0),(0,-1))))
g0=dense(((1,0,0,0),(0,1,0,0),(0,0,-1,0),(0,0,0,-1)))
def gi(s):
    return M(4,4,{(i,j+2):v for (i,j),v in s.a.items()}|
                  {(i+2,j):-v for (i,j),v in s.a.items()})
gamma=(g0,)+(tuple(gi(s) for s in pauli))
GA=tuple(const(kron(eye(4),g)) for g in gamma)
GN=tuple(const(kron(eye(2),g)) for g in gamma)
SA=tuple(tuple(comm(GA[mu],GA[nu]).scale(I*F(1,2)) for nu in range(4)) for mu in range(4))
BAR_A=const(kron(eye(4),g0)); BAR_N=const(kron(eye(2),g0))
MH2=F(7,3)
FF=dense(((Q(1,2),Q(-2,1)),(Q(3,-1),Q(1,-3))))
YY=dense(((Q(2,-1),Q(1,3)),(Q(-1,2),Q(4,1))))
MM=dense(((Q(5,-1),Q(2,1)),(Q(2,1),Q(7,2))))
MNS=const(kron(MM.star(),eye(4)))


def lift_phi(phi,flavor):
    """phi is a weak column; result maps sterile to weak x active flavors."""
    return P(16,8,{k:kron(kron(v,flavor),eye(4)) for k,v in phi.a.items()})
def lift_weak(w): return P(16,16,{k:kron(kron(v,eye(2)),eye(4)) for k,v in w.a.items()})
def cov(a,mu,aa): return a.deriv(mu).add(aa[mu].mul(a).scale(I))
def lap_cov(a,aa): return psum((cov(cov(a,mu,aa),mu,aa).scale(-ETA[mu]) for mu in range(4)),a.n,a.m)
def ka(a,aa): return psum((GA[mu].mul(cov(a,mu,aa)).scale(I) for mu in range(4)),a.n,a.m)
def kn(a): return psum((GN[mu].mul(a.deriv(mu)).scale(I) for mu in range(4)),a.n,a.m)


def normalization(phi,n,aa):
    return ka(phi.mul(pp(n)),aa).add(phi.mul(pp(kn(n)))).scale(F(1,2))
def original_j(phi,n,aa):
    """The six literal placements after integration by parts on the row."""
    an=kn(n); pn=p0(n)
    return psum((phi.mul(kn(pn)).scale(F(5,6)),
                 lap_cov(phi.mul(an),aa).scale(F(1,3)),
                 lap_cov(phi,aa).mul(an).scale(F(1,6)),
                 ka(phi.mul(pn),aa).scale(-F(1,6)),
                 phi.mul(an).scale(MH2/2),ka(phi.mul(n),aa).scale(MH2/2)),16,1)
def remainder_j(phi,n,aa,gradient_sign=1,drop_gradient=False,drop_box=False):
    an=kn(n); pn=p0(n)
    terms=[lap_cov(phi,aa).mul(an).scale(F(1,2))]
    for mu in range(4):
        if not drop_gradient:
            terms.append(cov(phi,mu,aa).mul(an.deriv(mu)).scale(-gradient_sign*F(2,3)*ETA[mu]))
        if not drop_box: terms.append(GA[mu].mul(cov(phi,mu,aa)).mul(pn).scale(-I*F(2,3)))
    return psum(terms,16,1)


def curvature(aa):
    return tuple(tuple(aa[nu].deriv(mu).add(aa[mu].deriv(nu).scale(-1)).add(comm(aa[mu],aa[nu]).scale(I))
                       for nu in range(4)) for mu in range(4))
def adj(a,mu,aa): return a.deriv(mu).add(comm(aa[mu],a).scale(I))
def curvature_packet(phi,n,aa,expanded=False,drop_phi_derivative=False):
    ff=curvature(aa); terms=[]; phin=phi.mul(n)
    phidn=tuple(phi.mul(n.deriv(rho)) for rho in range(4))
    dphin=tuple(cov(phi,rho,aa).mul(n) for rho in range(4)) if expanded else ()
    for mu in range(4):
        for nu in range(4):
            fphin=ff[mu][nu].mul(phin)
            for rho in range(4):
                if expanded:
                    prod=adj(ff[mu][nu],rho,aa).mul(phin).add(ff[mu][nu].mul(phidn[rho]))
                    if not drop_phi_derivative: prod=prod.add(ff[mu][nu].mul(dphin[rho]))
                else: prod=cov(fphin,rho,aa)
                terms.append(SA[mu][nu].mul(GA[rho]).mul(prod).scale(-I*F(1,4)))
            terms.append(GA[mu].mul(adj(ff[nu][mu],nu,aa)).mul(phin).scale(ETA[nu]*F(1,3)))
    return psum(terms,16,1)


def inner(a,b,charged=True): return integral_pair(a.dagger(),(BAR_A if charged else BAR_N).mul(b))
def plain(a,b): return integral_pair(a.transpose(),b)
def z(q,phi):
    l,n=q; return phi.mul(pp(n)),pp(phi.dagger().mul(l))
def zt(q,phi):
    l,n=q; return phi.star().mul(pp(n)),pp(phi.transpose().mul(l))
def mass(q,b):
    l,n=q; return b.mul(n),b.transpose().mul(l).add(MNS.mul(n))
def pairadd(q,r): return tuple(a.add(b) for a,b in zip(q,r))
def pairscale(c,q): return tuple(a.scale(c) for a in q)
def dm(q,phi,b,wrong_adjoint=False):
    left=z(mass(q,b),phi) if wrong_adjoint else zt(mass(q,b),phi)
    return pairscale(-F(1,2),pairadd(left,mass(z(q,phi),b)))
def pair_plain(a,b): return plain(a[0],b[0]).add(plain(a[1],b[1]))


checks={}; mutations={}
def check(category,label,a,b):
    if a!=b: raise RuntimeError('FAILED '+category+': '+label)
    checks[category]=checks.get(category,0)+1
def reject(label,a,b):
    if a==b: raise RuntimeError('MUTATION SURVIVED '+label)
    mutations[label]=mutations.get(label,0)+1


def field(n,seed,momenta):
    out=P(n,1)
    for j,k in enumerate(momenta):
        vals=M(n,1,{(i,0):Q(((i+2)*(seed+j+1))%7-3,((i+seed+2*j)%5)-2) for i in range(n)})
        out=out.add(mode(k,vals))
    return out


def setup(seed):
    p=(seed+1,-1,2,1); h=(1,seed+1,-1,2); h2=(-1,1,seed+2,0)
    phi=mode(h,dense(((Q(1,seed),),(Q(2,-1),)))).add(mode(h2,dense(((Q(-1,2),),(Q(1,seed+1),)))))
    aa=[]
    for mu in range(4):
        k=(1,mu+1,-1,seed+1)
        t=pauli[mu%3].scale(F(seed+1,mu+2)).add(eye(2).scale(F(mu+1,5)))
        t2=pauli[(mu+1)%3].scale(F(1,seed+2))
        w=mode(k,t.add(t2.scale(I))).add(mode(tuple(-x for x in k),t.add(t2.scale(-I))))
        aa.append(lift_weak(w))
    n=field(8,seed,((-p[0],-p[1],-p[2],-p[3]),(2,-1,seed+1,-2)))
    # Pairing test fields cover the complete Fourier support when needed.
    return phi,lift_phi(phi,FF),tuple(aa),n


def cover(v,seed): return field(v.n,seed,[k[:4] for k in v.a])
def cover_plain(v,seed): return field(v.n,seed,[tuple(-x for x in k[:4]) for k in v.a])


def run():
    global LIMIT
    started=time.perf_counter()
    LIMIT=1
    for mu in range(4):
        for nu in range(4):
            check('clifford','gamma anticommutator',
                  gamma[mu].mul(gamma[nu]).add(gamma[nu].mul(gamma[mu])),
                  eye(4).scale(2*ETA[mu] if mu==nu else 0))
    for seed in (1,2,3):
        trial_started=time.perf_counter()
        print('Starting exact nonconstant-background trial '+str(seed),flush=True)
        raw,phi,aa,n=setup(seed)
        oj=original_j(phi,n,aa); norm=normalization(phi,n,aa); rem=remainder_j(phi,n,aa)
        check('actual_six_placement_split','generic off-shell backgrounds',oj,norm.add(rem))
        reject('omit_covariant_gradient',oj,norm.add(remainder_j(phi,n,aa,drop_gradient=True)))
        reject('reverse_covariant_gradient_sign',oj,norm.add(remainder_j(phi,n,aa,gradient_sign=-1)))
        reject('omit_gradient_times_box',oj,norm.add(remainder_j(phi,n,aa,drop_box=True)))
        cd=curvature_packet(phi,n,aa)
        check('curvature_product_rule','all derivatives retained',cd,curvature_packet(phi,n,aa,expanded=True))
        reject('omit_curvature_times_Dphi',cd,curvature_packet(phi,n,aa,expanded=True,drop_phi_derivative=True))
        # Arbitrary infinitesimal local weak transformation; two noncommuting modes.
        th1=mode((1,-1,2,seed),pauli[0]); th2=mode((2,1,-seed,1),pauli[2].scale(F(2,3)))
        theta=lift_weak(th1.add(th1.dagger()).add(th2).add(th2.dagger()))
        dphi=theta.mul(phi).scale(-I)
        daa=tuple(theta.deriv(mu).add(comm(aa[mu],theta).scale(I)) for mu in range(4))
        ap=tuple(a.add(da.parameter()) for a,da in zip(aa,daa)); phip=phi.add(dphi.parameter())
        for label,fun,v in (('six_placements',original_j,oj),('normalization',normalization,norm),('remainder',remainder_j,rem),('curvature',curvature_packet,cd)):
            check('local_weak_covariance',label,fun(phip,n,ap).degree(1),theta.mul(v).scale(-I))
        lap=cover(phi.mul(n),seed+10)
        lap_p=lap.add(theta.mul(lap).scale(-I).parameter())
        neutral=pp(phi.dagger().mul(lap))
        check('neutral_composite_covariance','operator on phi dagger L',pp(phip.dagger().mul(lap_p)).degree(1),P(8,1))
        # Local current derivative, using a direction independent of theta.
        vv=tuple(lift_weak(mode((mu+1,seed,-2,1),pauli[(mu+2)%3])) for mu in range(4))
        av=tuple(a.add(v.parameter()) for a,v in zip(aa,vv))
        current=psum((GA[mu].mul(vv[mu]).mul(phi).mul(pp(n)).scale(-F(1,2)) for mu in range(4)),16,1)
        check('current_transport','variation of normalization contact',normalization(phi,n,av).degree(1),current)
        shifted_kinetic=psum((GA[mu].mul(vv[mu]).mul(phi).mul(pp(n)).scale(F(1,2)) for mu in range(4)),16,1)
        check('current_transport','bare weak current substitution cancellation',current.add(shifted_kinetic),P(16,1))
        check('current_transport','retained residual current',original_j(phi,n,av).degree(1).add(current.scale(-1)),remainder_j(phi,n,av).degree(1))
        reject('omit_current_transport',current,P(16,1))
        reject('naive_ordinary_Dphi_remainder',rem,remainder_j(phi,n,tuple(P(16,16) for _ in range(4))))
        # Integrated derivative adjoints are tested on fields covering their support.
        hn=phi.mul(pp(n)); l=cover(hn,seed+4)
        check('integrated_operator_adjoints','H versus P H dagger',inner(l,hn),inner(pp(phi.dagger().mul(l)),n,False))
        klhn=ka(hn,aa); lk=cover(klhn,seed+8)
        check('integrated_operator_adjoints','covariant Dirac adjoint',inner(lk,klhn),inner(ka(lk,aa),hn))
        # Distinct row sources; source derivatives must transport through P.
        src=cover_plain(pp(phi.dagger().mul(l)),seed+11).transpose()
        lhs=integral(src.mul(pp(phi.dagger().mul(l))))
        rhs=integral(pp(src).mul(phi.dagger()).mul(l))
        check('source_derivative_transport','full local source row',lhs,rhs)
        reject('drop_source_derivatives',lhs,integral(src.scale(MH2).mul(phi.dagger()).mul(l)))
        upper_source=cover_plain(hn,seed+16).transpose()
        upper_lhs=integral(upper_source.mul(phi).mul(pp(n)))
        upper_rhs=integral(pp(upper_source.mul(phi)).mul(n))
        check('source_derivative_transport','upper source and scalar composite',upper_lhs,upper_rhs)
        reject('drop_upper_source_scalar_derivatives',upper_lhs,integral(pp(upper_source).mul(phi).mul(n)))
        # Retained Majorana mass, nonconstant Higgs Yukawa blocks, differential transpose.
        b=lift_phi(raw.star(),YY.star())
        q=(field(16,seed+1,((1,2,3,1),(-2,1,-1,seed))),n)
        dv=dm(q,phi,b)
        u=(cover_plain(dv[0],seed+21),cover_plain(dv[1],seed+31))
        check('majorana_local_mass_congruence','formal differential transpose symmetry',pair_plain(u,dv),pair_plain(dm(u,phi,b),q))
        direct=pair_plain(z(u,phi),mass(q,b)).add(pair_plain(u,mass(z(q,phi),b))).scale(-F(1,2))
        check('majorana_local_mass_congruence','actual substitution in both slots',pair_plain(u,dv),direct)
        reject('mass_adjoint_instead_of_transpose',pair_plain(u,dv),pair_plain(u,dm(q,phi,b,wrong_adjoint=True)))
        # Powers two and three of arbitrary full noncommuting gauge background.
        LIMIT=3
        at=tuple(a.parameter() for a in aa)
        nt=normalization(phi,n,at); ct=curvature_packet(phi,n,at)
        ot=original_j(phi,n,at).add(ct)
        rt=remainder_j(phi,n,at).add(ct)
        for deg in (2,3):
            check('higher_gauge_contacts_retained','normalization has no degree '+str(deg),nt.degree(deg),P(16,1))
            check('higher_gauge_contacts_retained','actual contact degree '+str(deg),ot.degree(deg),rt.degree(deg))
            reject('discard_higher_gauge_contact_'+str(deg),ot.degree(deg),P(16,1))
        LIMIT=1
        # Exact previous neutral background, without imposing scalar stationarity.
        rawc=const(dense(((F(3,2),),(0,))))
        phic=lift_phi(rawc,FF); az=tuple(P(16,16) for _ in range(4))
        check('neutral_projection_reproduction','previous derivative normalization',original_j(phic,n,az),phic.mul(kn(pp(n))))
        check('neutral_projection_reproduction','interaction remainder vanishes',remainder_j(phic,n,az),P(16,1))
        check('neutral_projection_reproduction','curvature vanishes',curvature_packet(phic,n,az),P(16,1))
        print('Completed trial %s in %.2f seconds' % (seed,time.perf_counter()-trial_started),flush=True)
    result={'status':'PASS','calculation':'actual selected six-placement packet: local differential normalization and retained interactions',
            'arithmetic':'exact Gaussian rational finite Fourier sums; Python standard library only',
            'tests_by_category':checks,'positive_checks':sum(checks.values()),
            'mutation_rejections':mutations,'negative_checks':sum(mutations.values()),
            'total_checks':sum(checks.values())+sum(mutations.values()),
            'nonconstant_background_trials':3,'active_flavors':2,'sterile_flavors':2,'weak_components':2,
            'scope':['first order in selected counterterm','full local differential operator substitution, generic off-shell test fields',
                     'full covariant six-term O_J and curvature product rule','retained symmetric complex Majorana mass',
                     'selected normalization has at most one explicit gauge connection',
                     'not a complete basis, physical beta function, finite amplitude, or loop re-integration']}
    target=Path(__file__).with_name('INTERACTING_CANONICAL_RESULTS.json')
    target.write_text(json.dumps(result,indent=2)+'\n')
    print('Completed exact operator replay in %.2f seconds' % (time.perf_counter()-started),flush=True)
    print(json.dumps(result,indent=2))


if __name__=='__main__': run()
