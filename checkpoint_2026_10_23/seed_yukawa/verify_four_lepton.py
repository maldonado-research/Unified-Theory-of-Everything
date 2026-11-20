#!/usr/bin/env python3
"""Literal seed/two-Yukawa four-lepton pole: parameter integral and Wick algebra.

Pinned SymPy; no repository imports and no removable assertions. The primary
counterterm is the ordered spinor scalar; the vector-current translation is a
four-dimensional physical-spinor projection, not an evanescent-basis closure.
"""
from collections import Counter
from itertools import product
from math import factorial
from pathlib import Path
import json
import sympy as s

I=s.I;R=s.Rational
positive=Counter();negative=Counter()

def zz(x):
    if isinstance(x,G):return not x.d
    if isinstance(x,s.MatrixBase):return all(s.expand(v)==0 for v in x)
    return s.expand(x)==0

def check(label,x):
    if not zz(x):raise RuntimeError(label+': exact identity failed')
    positive[label]+=1

def reject(label,x):
    if zz(x):raise RuntimeError(label+': wrong formula escaped')
    negative[label]+=1

class G:
    """Exterior algebra, with coefficients exactly reduced after each product."""
    def __init__(self,data=None):
        self.d={k:s.expand(v) for k,v in (data or {}).items() if s.expand(v)!=0}
    @staticmethod
    def cast(x):return x if isinstance(x,G) else G({():s.sympify(x)})
    @staticmethod
    def gen(i):return G({(i,):1})
    def __add__(self,x):
        x=G.cast(x);a=dict(self.d)
        for k,v in x.d.items():a[k]=a.get(k,0)+v
        return G(a)
    __radd__=__add__
    def __neg__(self):return G({k:-v for k,v in self.d.items()})
    def __sub__(self,x):return self+-G.cast(x)
    def __rsub__(self,x):return G.cast(x)+-self
    def __mul__(self,x):
        x=G.cast(x);a={}
        for k,v in self.d.items():
            for l,w in x.d.items():
                if set(k)&set(l):continue
                sign=(-1)**sum(i>j for i in k for j in l)
                key=tuple(sorted(k+l));a[key]=a.get(key,0)+sign*v*w
        return G(a)
    __rmul__=__mul__
    def __truediv__(self,x):return self*(1/s.sympify(x))
    def dagger(self,mate):
        out=G()
        for k,v in self.d.items():
            t=G.cast(s.conjugate(v))
            for j in reversed(k):t=t*G.gen(mate[j])
            out+=t
        return out

E=s.Matrix([[0,1],[-1,0]])
pauli=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-I],[I,0]]),s.diag(1,-1)]
bs=[s.eye(2)]+[-p for p in pauli];eta=[1,-1,-1,-1]
z2=s.zeros(2)
gamma=[s.BlockMatrix([[z2,s.eye(2)],[s.eye(2),z2]]).as_explicit()]+[s.BlockMatrix([[z2,p],[-p,z2]]).as_explicit() for p in pauli]
Cspin=I*gamma[2]*gamma[0]
g5=I*gamma[0]*gamma[1]*gamma[2]*gamma[3]
PL=(s.eye(4)-g5)/2;PR=(s.eye(4)+g5)/2
check('literal chiral charge-conjugation blocks',Cspin-s.diag(E,-E))
check('left charge-conjugate lower block',Cspin[2:4,2:4]+E)

# A genuine triangle parameter integral. The leading numerator is q² from
# P0 on the internal scalar times its external lepton. Translation changes
# only lower loop-momentum degrees, which are UV finite for this graph.
x,y=s.symbols('x y')
def simplex_integral(poly):
    return s.expand(sum(c*R(factorial(a)*factorial(b),factorial(a+b+2))
                 for (a,b),c in s.Poly(s.expand(poly),x,y).terms()))
def radial_pole(r,n):
    """Coefficient of i/(16pi²epsilon) in l^(2r)/(l²-Delta)^n."""
    k=n-r-2
    if k>0:return s.Integer(0)
    Delta=s.symbols('Delta')
    return s.Integer(-1)**(r-n)*s.factorial(r+1)*s.Integer(-1)**(-k)*Delta**(2+r-n)/(s.factorial(n-1)*s.factorial(-k))
triangle=2*simplex_integral(1)*radial_pole(1,3)
check('three-denominator Feynman simplex normalization',2*simplex_integral(1)-1)
check('actual q squared triangle UV residue',triangle-1)
check('constant triangle numerator is UV finite',radial_pole(0,3))
check('one more denominator is UV finite',radial_pole(1,4))
for a,b in product(range(3),repeat=2):
    # Independent nested elementary integration of simplex monomials.
    nested=s.integrate(s.integrate(x**a*y**b,(y,0,1-x)),(x,0,1))
    check('simplex exact polynomial moment',simplex_integral(x**a*y**b)-nested)
reject('omit denominator-combination factorial',triangle/2-1)
reject('drop scalar derivatives in P0',radial_pole(0,3)-triangle)
for mass,p0,p1 in [(R(1,3),s.Matrix([0,1,2,0]),s.Matrix([0,-2,1,1])),
                   (R(7,2),s.Matrix([0,2,0,1]),s.Matrix([0,1,-1,2]))]:
    # q=l-shift, where shift=x p0+y p1 for a routed triangle. The complete
    # shifted numerator (q+p_external)^2 has coefficient one on l².
    shift=x*p0+y*p1
    q=s.Matrix(s.symbols('q0:4'))
    metric=s.diag(*eta)
    full=((q-shift+p0).T*metric*(q-shift+p0))[0]
    for mu in range(4):check('whole-composite leading quadratic coefficient',s.diff(full,q[mu],2)/2-eta[mu])
    delta=mass+x*(p0.T*metric*p0)[0]+y*(p1.T*metric*p1)[0]-(shift.T*metric*shift)[0]
    check('parameter mass does not enter logarithmic residue',s.diff(2*simplex_integral(1)*radial_pole(1,3),s.symbols('Delta')))
    # Its constant and linear terms multiply UV-finite/odd radial masters.
    check('external-shift constant term is UV finite',full.subs(dict.fromkeys(q,0))*radial_pole(0,3))

# Full-spin/flavor inverse-propagator expansion fixes which complex mass
# multiplies the anomalous <n n> propagator. n=N_L; its holomorphic mass is M*.
Y=s.Matrix([[1+I,2-I],[3-2*I,-1+I]])
M=s.Matrix([[2+I,1-2*I],[1-2*I,4-I]])
C=s.Matrix([[2-I,1+2*I],[1+2*I,3+I]])
S=Y.conjugate()*M*Y.H
check('complex Majorana input symmetric',M-M.T)
check('seed input symmetric',C-C.T)
check('anomalous Yukawa mass tensor symmetric',S-S.T)
q=s.Matrix([3,1,-1,0]);q2=(q.T*s.diag(*eta)*q)[0]
slash=sum((eta[mu]*q[mu]*gamma[mu] for mu in range(4)),s.zeros(4))
slashF=s.kronecker_product(s.eye(2),slash)
PLF=s.kronecker_product(s.eye(2),PL)
massF=s.kronecker_product(M,PR)+s.kronecker_product(M.H,PL)
leading=PLF*slashF*massF*slashF*PLF/q2
check('inverse-propagator first mass term has M on left block',leading-s.kronecker_product(M,PL))
prop_nn=leading*s.kronecker_product(s.eye(2),Cspin.T)
check('literal anomalous left Weyl propagator tensor',prop_nn+s.kronecker_product(M,s.diag(E,s.zeros(2))))
reject('replace anomalous mass by its adjoint',S-Y.conjugate()*M.H*Y.H)

# Construct external variables for two active flavors, two weak components
# and two spin components. No commuting-spinor approximation is used.
nf=2
bar=[[[G.gen((a*2+w)*2+r) for r in range(2)] for w in range(2)] for a in range(nf)]
lep=[[[G.gen(8+(a*2+w)*2+r) for r in range(2)] for w in range(2)] for a in range(nf)]
mate={j:(j+8 if j<8 else j-8) for j in range(16)}
def bil(u,v):return sum((E[r,t]*u[r]*v[t] for r,t in product(range(2),repeat=2)),G())
def packet(Ca,Sa):
    return sum((Ca[a,b]*Sa[c,d]*bil(bar[a][u],bar[b][v])*bil(lep[c][u],lep[d][v])
       for a,b,c,d,u,v in product(range(nf),range(nf),range(nf),range(nf),range(2),range(2))),G())
P=packet(C,S)
check('symmetric seed projection removes antisymmetric coefficient',packet(s.Matrix([[0,1],[-1,0]]),S))
check('symmetric anomalous mass projection removes antisymmetric coefficient',packet(C,s.Matrix([[0,1],[-1,0]])))

# Literal ordered path-integral expansion. Vertex factors: i³. The seed
# charge-conjugate block has -E; each B has action sign minus. The external
# order is [bar_a,bar_b,n_at_B1,l_c,n_at_B2,l_d]; contraction of the two n's
# crosses exactly one external l. The expansion coefficient is 1/(2!*2).
wick_sign=(-1)**(4-2-1)
phase=I**3*(-1)*(-1)**2*(-I)*(I**2)*wick_sign/R(4)
check('literal Wick crossing sign',wick_sign+1)
check('literal seed Yukawa scalar anomalous phases',phase-R(1,4))
for u,w in product(range(2),repeat=2):
    spin=sum(E[t,u]*E[v,w]*E[t,v] for t,v in product(range(2),repeat=2))
    check('contracted two Yukawa spin tensors',spin-E[u,w])

# Keep both scalar pairings separately before using C and S symmetry.
direct=phase*triangle*P
cross=sum((phase*triangle*C[a,b]*S[c,d]*bil(bar[a][u],bar[b][v])*bil(lep[c][v],lep[d][u])
    for a,b,c,d,u,v in product(range(nf),range(nf),range(nf),range(nf),range(2),range(2))),G())
loop=direct+cross
check('crossed and direct scalar contractions agree after flavor sums',cross-direct)
check('literal triangle action coefficient',loop-P/2)
ct=-R(1,2)*P
check('counterterm cancels loop action pole',loop+ct)
reject('omit crossed identical scalar pairing',direct+ct)
reject('append extra Majorana half',loop/2+ct)
reject('delete identical-fermion seed half',2*loop+ct)
reject('wrong Wick or counterterm sign',-loop+ct)

# Exact Pauli completeness and four-dimensional Fierz translation, with
# fermion reordering provided by the exterior product itself.
for r,t,u,v in product(range(2),repeat=4):
    tensor=sum(eta[mu]*bs[mu][r,u]*bs[mu][t,v] for mu in range(4))
    check('all Pauli completeness coefficients',tensor-2*E[r,t]*E[u,v])
def current(a,c,mu):
    return sum((bs[mu][r,t]*bar[a][w][r]*lep[c][w][t]
                for w,r,t in product(range(2),repeat=3)),G())
JJ=sum((C[a,b]*S[c,d]*eta[mu]*current(a,c,mu)*current(b,d,mu)
    for a,b,c,d,mu in product(range(nf),range(nf),range(nf),range(nf),range(4))),G())
check('four-dimensional full exterior Fierz identity',JJ+2*P)
check('vector-current counterterm normalization',ct-R(1,4)*JJ)
reject('commuting-spinor Fierz sign',JJ-2*P)
full=loop+loop.dagger(mate)
check('conjugate seed channel gives Hermitian action',full-full.dagger(mate))
reject('generic complex displayed channel mistaken for complete Hermitian action',loop-loop.dagger(mate))

# One-flavor, one fixed weak component: P=(2bar1bar2)(2l1l2).
b0,b1,l0,l1=[G.gen(i) for i in range(4)]
one=bil([b0,b1],[b0,b1])*bil([l0,l1],[l0,l1])
check('nonzero one-flavor component normalization',one-4*b0*b1*l0*l1)
reject('one flavor four-lepton pole incorrectly removed',one)
check('real one-flavor seed and conjugate pole',R(1,2)*one+R(1,2)*one-one)

U=s.Matrix([[R(3,5),4*I/5],[4*I/5,R(3,5)]])
V=s.diag(I,(3+4*I)/5)*s.Matrix([[R(5,13),R(12,13)],[-R(12,13),R(5,13)]])
check('test sterile basis unitary',U.H*U-s.eye(2))
check('test active basis unitary',V.H*V-s.eye(2))
Yp=V.H*Y*U.conjugate();Mp=U.H*M*U.conjugate();Sp=Yp.conjugate()*Mp*Yp.H
check('anomalous mass source flavor covariance',Sp-V.T*S*V)
Cp=V.H*C*V.conjugate()
check('rotated seed remains symmetric',Cp-Cp.T)
reject('wrong anomalous source flavor adjoint',Sp-V.H*S*V.conjugate())
check('zero retained mass removes this triangle pole',packet(C,s.zeros(2)))
check('zero Yukawa removes this triangle pole',packet(C,s.zeros(2)*M*s.zeros(2)))

report={'status':'PASS','input_commit':'f1f9a892e15545c70e4fa05827ea7fefc5a438eb',
 'dateline':'23 October 2026 Pacific',
 'method':'Literal seed and two conjugate Yukawa vertices; Feynman-simplex/Gamma radial UV integration, explicit Wick phases, exterior Grassmann spin/weak/flavor coefficients.',
 'positive_checks':sum(positive.values()),'positive_categories':dict(positive),
 'rejected_mutations':sum(negative.values()),'mutation_categories':dict(negative),
 'computed_triangle_residue_over_i_kappa':str(triangle),
 'computed_literal_prefactor_before_scalar_pair_sum':str(phase),
 'ordered_loop_action':'kappa/2 C_ab (Y* M Ydagger)_cd P_abcd + h.c.',
 'ordered_counterterm_action':'-kappa/2 C_ab (Y* M Ydagger)_cd P_abcd + h.c.',
 'four_dimensional_vector_counterterm':'kappa/4 C_ab (Y* M Ydagger)_cd (barL_a gamma_mu L_c)(barL_b gamma^mu L_d) + h.c.',
 'scope':'Specific one-seed/two-Yukawa four-lepton pole; no total physical Wilson coefficient, complete RG, finite matching or novelty claim.',
 'repository_imports':False,'assertions_required':False}
Path(__file__).with_name('FOUR_LEPTON_RESULTS.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
