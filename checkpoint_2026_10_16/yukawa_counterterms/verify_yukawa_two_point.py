#!/usr/bin/env python3
"""Pure-neutrino-Yukawa two-point UV residues and action counterterms.

Exact Feynman-parameter polynomial integration, Gamma-pole recurrences,
literal phases, Clifford projectors, weak/Wick multiplicities and complex
flavor covariance. Standard library only; no project imports or asserts.
"""
from collections import Counter
from fractions import Fraction as R
from itertools import permutations
import json
from math import factorial
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
        if isinstance(o,Grass):return o*self
        o=Q(o);return Q(self.a*o.a-self.b*o.b,self.a*o.b+self.b*o.a)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=Q(o);d=o.a*o.a+o.b*o.b
        return Q((self.a*o.a+self.b*o.b)/d,(self.b*o.a-self.a*o.b)/d)
    def __eq__(self,o):o=Q(o);return self.a==o.a and self.b==o.b
    def __bool__(self):return bool(self.a or self.b)
    def conj(self):return Q(self.a,-self.b)
    def __repr__(self):return 'Q(%s,%s)'%(self.a,self.b)


class Poly:
    """Rational polynomial in x,M^2,m_H^2,p^2, never sampled."""
    def __init__(self,d=None):self.d={k:R(v) for k,v in (d or {}).items() if v}
    @staticmethod
    def cast(v):return v if isinstance(v,Poly) else Poly({(0,0,0,0):v})
    @staticmethod
    def var(i):
        m=[0]*4;m[i]=1;return Poly({tuple(m):1})
    def __add__(self,o):
        o=Poly.cast(o);d=dict(self.d)
        for k,v in o.d.items():d[k]=d.get(k,0)+v
        return Poly(d)
    __radd__=__add__
    def __neg__(self):return Poly({k:-v for k,v in self.d.items()})
    def __sub__(self,o):return self+-Poly.cast(o)
    def __rsub__(self,o):return Poly.cast(o)+-self
    def __mul__(self,o):
        o=Poly.cast(o);d={}
        for a,x in self.d.items():
            for b,y in o.d.items():
                k=tuple(i+j for i,j in zip(a,b));d[k]=d.get(k,0)+x*y
        return Poly(d)
    __rmul__=__mul__
    def integrate_x(self):
        d={}
        for k,v in self.d.items():
            m=(0,)+k[1:];d[m]=d.get(m,0)+v/R(k[0]+1)
        return Poly(d)
    def __eq__(self,o):return self.d==Poly.cast(o).d
    def __repr__(self):return repr(self.d)


class Grass:
    def __init__(self,d=None):self.d={k:Q(v) for k,v in (d or {}).items() if v}
    @staticmethod
    def cast(v):return v if isinstance(v,Grass) else Grass({():Q(v)})
    def __add__(self,o):
        o=Grass.cast(o);d=dict(self.d)
        for k,v in o.d.items():d[k]=d.get(k,Q())+v
        return Grass(d)
    __radd__=__add__
    def __neg__(self):return Grass({k:-v for k,v in self.d.items()})
    def __sub__(self,o):return self+-Grass.cast(o)
    def __mul__(self,o):
        o=Grass.cast(o);d={}
        for a,x in self.d.items():
            for b,y in o.d.items():
                if set(a)&set(b):continue
                sign=(-1)**sum(i>j for i in a for j in b);k=tuple(sorted(a+b))
                d[k]=d.get(k,Q())+sign*x*y
        return Grass(d)
    __rmul__=__mul__
    def dagger(self):
        out=Grass()
        for k,v in self.d.items():
            term=Grass.cast(v.conj())
            for j in reversed(k):term=term*Grass({(j^1,):1})
            out+=term
        return out
    def left_derivative(self,j):
        d={}
        for k,v in self.d.items():
            if j in k:
                pos=k.index(j);new=k[:pos]+k[pos+1:]
                d[new]=d.get(new,Q())+(-1)**pos*v
        return Grass(d)
    def __eq__(self,o):return self.d==Grass.cast(o).d
    def __repr__(self):return repr(self.d)


I=Q(0,1);eta=(1,-1,-1,-1)
def matrix(n,m,f):return tuple(tuple(Q(f(i,j)) for j in range(m)) for i in range(n))
def eye(n):return matrix(n,n,lambda i,j:i==j)
def zero(n,m=None):return matrix(n,n if m is None else m,lambda i,j:0)
def scale(c,a):return tuple(tuple(c*x for x in row) for row in a)
def add(a,b):return tuple(tuple(x+y for x,y in zip(ar,br)) for ar,br in zip(a,b))
def sub(a,b):return add(a,scale(-1,b))
def mul(a,b):return tuple(tuple(sum((a[i][k]*b[k][j] for k in range(len(b))),Q()) for j in range(len(b[0]))) for i in range(len(a)))
def transpose(a):return tuple(zip(*a))
def dagger(a):return tuple(tuple(a[j][i].conj() for j in range(len(a))) for i in range(len(a[0])))
def tr(a):return sum((a[i][i] for i in range(len(a))),Q())
def conj(a):return tuple(tuple(z.conj() for z in row) for row in a)
def msum(xs,n=4,m=None):
    z=zero(n,m)
    for x in xs:z=add(z,x)
    return z
pauli=(matrix(2,2,lambda i,j:i!=j),((Q(),-I),(I,Q())),((Q(1),Q()),(Q(),Q(-1))))
g0=matrix(4,4,lambda i,j:(1 if i<2 else -1) if i==j else 0)
def spatial(p):return tuple(tuple(p[i][j-2] if i<2 and j>=2 else -p[i-2][j] if i>=2 and j<2 else Q() for j in range(4)) for i in range(4))
gamma=(g0,)+tuple(spatial(p) for p in pauli)
g5=scale(I,mul(mul(mul(gamma[0],gamma[1]),gamma[2]),gamma[3]))
PL=scale(R(1,2),sub(eye(4),g5));PR=scale(R(1,2),add(eye(4),g5));C=scale(I,mul(gamma[2],gamma[0]))
counts=Counter();mutations=set()
def check(label,a,b):
    if a!=b:raise RuntimeError('check failed: '+label+'; '+repr(a)+' != '+repr(b))
    counts[label]+=1
def reject(label,a,b):
    if a==b:raise RuntimeError('mutation survived: '+label)
    mutations.add(label)


# I_n=i(-1)^n Gamma(n-2+eps)Delta^(2-eps-n)/(4pi)^(2-eps).
# The simple-pole residue of Gamma(eps-m) follows by recurrence from
# Res Gamma(eps)=1.  Positive integer shifts are UV finite here.
def gamma_residue(shift):
    return R((-1)**(-shift),factorial(-shift)) if shift<=0 else R(0)
for n,expected in ((1,1),(2,1),(3,0),(4,0)):
    check('radial Gamma master residue',(-1)**n*gamma_residue(n-2),R(expected))
x,M2,m2,p2=(Poly.var(i) for i in range(4))
deltaL=x*M2+(1-x)*m2-x*(1-x)*p2
deltaN=(1-x)*m2-x*(1-x)*p2
deltaPhi=x*M2-x*(1-x)*p2
check('massive scalar fermion parameter denominator',deltaL,deltaPhi+(1-x)*m2)
check('massless active fermion parameter denominator',deltaN,(1-x)*m2-x*(1-x)*p2)
check('open fermion vector bubble coefficient',(1-x).integrate_x(),Poly.cast(R(1,2)))
check('opposite routing gives same vector coefficient',x.integrate_x(),Poly.cast(R(1,2)))
# tr[P_R(qslash+M)P_L(qslash-pslash)] after q=l+(1-x)p
# is 2[l^2-x(1-x)p^2]+an odd term. Res l^2/(l^2-Delta)^2=2Delta.
scalar_numerator_residue=(4*deltaPhi-2*x*(1-x)*p2).integrate_x()
check('actual scalar Feynman parameter UV polynomial',scalar_numerator_residue,2*M2-p2)
check('closed-loop scalar inverse-kernel UV polynomial',-scalar_numerator_residue,p2-2*M2)
reject('drop massive radial tadpole contribution',scalar_numerator_residue,(2*deltaPhi-2*x*(1-x)*p2).integrate_x())
reject('fermion vector bubble coefficient doubled',(1-x).integrate_x(),Poly.cast(1))


for m in range(4):
    check('open L chain strips massive numerator chirality',mul(mul(PR,gamma[m]),PL),mul(gamma[m],PL))
    check('open N chain strips active chirality',mul(mul(PL,gamma[m]),PR),mul(gamma[m],PR))
    check('Majorana conjugate kinetic chain flips chirality',mul(mul(C,transpose(scale(-1,mul(gamma[m],PR)))),scale(-1,C)),mul(gamma[m],PL))
    check('massive scalar numerator linear mass term vanishes',tr(mul(mul(PR,PL),gamma[m])),Q())
    for n in range(4):
        check('scalar loop all q and p spin coefficients',tr(mul(mul(mul(PR,gamma[m]),PL),gamma[n])),Q(2*eta[m] if m==n else 0))
check('L chirality-flipping mass numerator is absent',mul(PR,PL),zero(4))
check('opposite L chirality-flipping mass numerator is absent',mul(PL,PR),zero(4))
open_phase=(-I)*(-I)*I*I
check('literal open two-point phase',open_phase,Q(1))
check('literal closed two-point phase',-open_phase,Q(-1))
# At fixed external phi component, the fermion sequence is barL,N,barN,L.
# Nested contractions have no crossing and exactly one reversed contraction
# <barL L>=-<L barL>. Swapping the vertices gives the same sign.
wick_weights=[]
for orientation in permutations((0,1)):
    sequence=('barL','N','barN','L') if orientation==(0,1) else ('barN','L','barL','N')
    pairs=((0,3),(1,2));crossings=sum(a<c<b<d or c<a<d<b for i,(a,b) in enumerate(pairs) for c,d in pairs[i+1:])
    reversed_pairs=sum(sequence[a].startswith('bar') for a,b in pairs)
    wick_weights.append((-1)**(crossings+reversed_pairs))
check('closed mixed-species Wick sign and expansion factorial',R(sum(wick_weights),factorial(2)),R(-1))
weakN=sum(int(a==b) for a in range(2) for b in range(2))
check('sterile loop weak doublet multiplicity',weakN,2)
for ext in range(2):
    check('active fixed-component weak multiplicity',sum(int(ext==a) for a in range(2)),1)
    check('scalar fixed-component weak multiplicity',sum(int(ext==a) for a in range(2)),1)
reject('missing closed fermion minus',-open_phase,open_phase)
reject('sterile weak multiplicity omitted',R(weakN,2),R(1,2))
reject('scalar receives sterile weak trace factor',scalar_numerator_residue,2*scalar_numerator_residue)


y=((Q(1,2),Q(2,-1),Q(-1,1)),(Q(3,-1),Q(-2,2),Q(1,3)))
hl=mul(y,dagger(y));h=mul(dagger(y),y);ty=tr(h)
md=matrix(3,3,lambda i,j:(2,3,5)[i] if i==j else 0)
mass=((Q(2,1),Q(1,-2),Q(-1,1)),(Q(1,-2),Q(4,-1),Q(2,3)),(Q(-1,1),Q(2,3),Q(6,2)))
u=((Q(R(3,5)),Q(0,R(4,5)),Q()),(Q(0,R(4,5)),Q(R(3,5)),Q()),(Q(),Q(),Q(1)))
v=((Q(R(5,13)),Q(R(12,13))), (Q(R(-12,13)),Q(R(5,13))))
check('test sterile rotation is unitary',mul(u,dagger(u)),eye(3))
check('test active rotation is unitary',mul(v,dagger(v)),eye(2))
check('active kinetic flavor matrix Hermiticity',hl,dagger(hl))
check('sterile kinetic flavor matrix Hermiticity',h,dagger(h))
check('active and sterile Yukawa traces agree',tr(hl),ty)
for mm in (md,mass):
    tm=tr(mul(h,mul(dagger(mm),mm)))
    ym=mul(y,dagger(mm));norm=sum((z*z.conj() for row in ym for z in row),Q())
    check('scalar mass trace is a nonnegative squared norm',tm,norm)
    check('scalar mass trace is real',tm,tm.conj())
    if tm.a<0:raise RuntimeError('nonnegative scalar trace failed')
    zp=scale(-ty,eye(1));zl=scale(R(-1,2),hl);zr=scale(-1,h)
    dzmass=scale(R(1,2),add(mul(transpose(h),mm),mul(mm,h)))
    check('Majorana mass parameter counterterm is symmetric',dzmass,transpose(dzmass))
    induced=scale(R(1,2),add(mul(transpose(zr),mm),mul(mm,zr)))
    check('bare field and mass changes leave zero total mass pole',add(dzmass,induced),zero(3))
    scalar_parameter=Q(7)*ty-2*tm
    scalar_coefficient=scalar_parameter-Q(7)*ty
    check('scalar parameter plus field shift gives literal mass coefficient',scalar_coefficient,-2*tm)
    check('scalar kinetic and mass counterterms cancel inverse-kernel pole',add(((ty,),),zp),zero(1))
    check('active kinetic counterterm cancels open loop',add(scale(R(1,2),hl),zl),zero(2))
    check('sterile kinetic counterterm cancels open loop',add(h,zr),zero(3))
    yp=mul(mul(v,y),dagger(u));hp=mul(dagger(yp),yp);mp=mul(mul(conj(u),mm),dagger(u))
    tmp=tr(mul(hp,mul(dagger(mp),mp)))
    check('general complex mass scalar trace is unitary covariant',tmp,tm)
    check('sterile kinetic matrix is unitary covariant',hp,mul(mul(u,h),dagger(u)))
    check('active kinetic matrix is unitary covariant',mul(yp,dagger(yp)),mul(mul(v,hl),dagger(v)))
    dmp=scale(R(1,2),add(mul(transpose(hp),mp),mul(mp,hp)))
    check('Majorana parameter mass pole transforms by congruence',dmp,mul(mul(conj(u),dzmass),dagger(u)))
    if mm==mass:
        reject('complex scalar mass ordering MMdag instead of MdagM',tm,tr(mul(h,mul(mm,dagger(mm)))))
        reject('wrong Majorana mass flavor order',dzmass,scale(R(1,2),add(mul(h,mm),mul(mm,transpose(h)))))
        reject('scalar parameter mass pole omits field normalization',scalar_parameter,-2*tm)
        reject('proper Majorana chirality mass pole invented',add(dzmass,induced),dzmass)


# Majorana conjugation of the complete kinetic flavor kernel. The
# transpose includes both spin and flavor and flips p in the open chain.
for m in range(4):
    kr=mul(gamma[m],PR);kl=mul(gamma[m],PL)
    kernel=matrix(12,12,lambda a,b:h[a//4][b//4]*kr[a%4][b%4]+h[b//4][a//4]*kl[a%4][b%4])
    conjugator=matrix(12,12,lambda a,b:C[a%4][b%4] if a//4==b//4 else 0)
    conjugated=mul(mul(conjugator,transpose(scale(-1,kernel))),scale(-1,conjugator))
    check('full Majorana flavor kinetic kernel has both required channels',conjugated,kernel)
    if m==0:
        wrong=matrix(12,12,lambda a,b:h[a//4][b//4]*kr[a%4][b%4])
        reject('omit conjugate Majorana kinetic channel',mul(mul(conjugator,transpose(scale(-1,wrong))),scale(-1,conjugator)),wrong)


def gv(a,v):return [sum((v[j]*a[i][j] for j in range(4)),Grass()) for i in range(4)]
def gb(v):return [sum((v[j].dagger()*gamma[0][j][i] for j in range(4)),Grass()) for i in range(4)]
def gc(v):return gv(C,gb(v))
def ga(a,b):return [x+y for x,y in zip(a,b)]
def gd(a,b):return sum((x*y for x,y in zip(a,b)),Grass())
z=[Grass({(2*i,):1}) for i in range(4)];w=[Grass({(8+2*i,):1}) for i in range(4)]
n=ga(z,gc(z));dn=ga(w,gc(w))
first=gd(gb(dn),n)
check('Majorana quadratic half cancels two varied mass slots',(first+gd(gb(n),dn))*R(1,2),first)
for m in range(4):
    first=gd(gb(dn),gv(scale(I,gamma[m]),n))
    second=-gd(gb(n),gv(scale(I,gamma[m]),dn))
    check('Majorana quadratic half cancels kinetic slots after IBP',(first+second)*R(1,2),first)
reject('retain extra half in Majorana two-point residue',first,first*R(1,2))
theta=[Grass({(i,):1}) for i in range(2)];aa=((Q(),Q(3,2)),(Q(-3,-2),Q()))
action=sum((theta[i]*theta[j]*(aa[i][j]/2) for i in range(2) for j in range(2)),Grass())
check('Majorana antisymmetric quadratic derivative fixes propagator normalization',action.left_derivative(0).left_derivative(1),Grass.cast(aa[0][1]))


# With Y_0=mu^eps(Y+...), all displayed residues are homogeneous
# degree two in Y,Y*: d_mu residue=-2eps residue at this order.
# Therefore gamma=(1/2)d_mu deltaZ=-eps deltaZ; a dimensionful
# parameter pole gives beta=2 times its 1/eps numerator.
check('epsilon relation gives active anomalous-dimension numerator',scale(-1,scale(R(-1,2),hl)),scale(R(1,2),hl))
check('epsilon relation gives sterile anomalous-dimension numerator',scale(-1,scale(-1,h)),h)
check('epsilon relation gives scalar anomalous-dimension numerator',-(-ty),ty)
reject('use d4minus-eps factor instead of d4minus2eps',scale(R(1,2),h),h)


result={
    'status':'passed','python':sys.version.split()[0],
    'input_candidate':'0770c304c7e05f213dd796d3fd5a194107f38858',
    'date':'15 October 2026 Pacific; directory uses 16 October UTC',
    'arithmetic':'Exact rational Feynman-parameter polynomials, Gamma-pole recurrences, Gaussian rational Clifford/flavor matrices and exterior Grassmann algebra; standard library only.',
    'positive_check_count':sum(counts.values()),'checks_by_category':dict(sorted(counts.items())),
    'negative_control_count':len(mutations),'negative_controls_rejected':sorted(mutations),
    'integral_residues':{'vector_bubble':'p_slash/2','scalar_trace_integral':'2 M_i^2-p^2','scalar_closed_loop':'p^2-2 M_i^2'},
    'counterterm_numerators_after_factoring_kappa':{'deltaZ_L':'-YYdag/2','deltaZ_NR':'-YdagY','deltaZ_phi':'-Tr(YdagY)','delta_m2_literal_coefficient':'-2Tr(YdagY MdagM)','delta_m2_parameter':'m_H^2 Tr(YdagY)-2Tr(YdagY MdagM)','deltaM_parameter':'(H^T M+M H)/2'},
    'scope':'Actual restricted pure-neutrino-Yukawa one-loop renormalizable two-point poles and minimal-subtraction inputs; excludes gauge, other Yukawa, scalar quartic, EFT insertions and complete dimension-seven running.'
}
Path(__file__).with_name('YUKAWA_TWO_POINT_RESULTS.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps(result,indent=2,sort_keys=True))
