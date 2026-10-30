#!/usr/bin/env python3
"""Independent zero-retained-N projection: exact Fourier/exterior algebra.

SymPy 1.14.0. No project imports or asserts. Prints JSON under python -I -B -O.
Zero-N means literal field-content projection, never integration or an EOM.
"""
import itertools
import json
import sys
import sympy as s
I,R=s.I,s.Rational
eta=[1,-1,-1,-1]
ZERO=(0,0,0,0)
conj={}
checks,negatives=[],[]

class P:
    def __init__(self,d=None):self.d={k:s.expand(v) for k,v in (d or {}).items() if s.expand(v)!=0}
    @staticmethod
    def num(x):return P({(ZERO,()):x})
    def __add__(self,q):
        if not isinstance(q,P):q=P.num(q)
        out=self.d.copy()
        for k,v in q.d.items():out[k]=out.get(k,0)+v
        return P(out)
    __radd__=__add__
    def __neg__(self):return P({k:-v for k,v in self.d.items()})
    def __sub__(self,q):return self+-q
    def __rsub__(self,q):return -self+q
    def __mul__(self,q):
        if not isinstance(q,P):q=P.num(q)
        out={}
        for (pa,ga),va in self.d.items():
            for (pb,gb),vb in q.d.items():
                if set(ga)&set(gb):continue
                sign=(-1)**sum(a>b for a in ga for b in gb)
                key=(tuple(a+b for a,b in zip(pa,pb)),tuple(sorted(ga+gb)))
                out[key]=out.get(key,0)+sign*va*vb
        return P(out)
    __rmul__=__mul__
    def __truediv__(self,c):return self*(1/s.sympify(c))
    def dag(self):
        out=P()
        for (p,g),v in self.d.items():
            gs=tuple(conj[a] for a in reversed(g))
            sign=(-1)**sum(gs[a]>gs[b] for a in range(len(gs)) for b in range(a+1,len(gs)))
            out=out+P({(tuple(-x for x in p),tuple(sorted(gs))):sign*s.conjugate(v)})
        return out
    def der(self,mu):return P({(p,g):-I*eta[mu]*p[mu]*v for (p,g),v in self.d.items()})
    def p0(self):return P({(p,g):sum(eta[i]*p[i]*p[i] for i in range(4))*v for (p,g),v in self.d.items()})
    def integral(self):return P({(p,g):v for (p,g),v in self.d.items() if p==ZERO})

def check(name,p):
    bad=p.d if isinstance(p,P) else (any(s.expand(v)!=0 for v in p) if isinstance(p,s.MatrixBase) else s.expand(p)!=0)
    if bad:raise RuntimeError(name+': '+str(bad))
    checks.append(name)
def reject(name,p):
    bad=p.d if isinstance(p,P) else (any(s.expand(v)!=0 for v in p) if isinstance(p,s.MatrixBase) else s.expand(p)!=0)
    if not bad:raise RuntimeError(name+': mutation escaped')
    negatives.append(name)
def addv(a,b):return [x+y for x,y in zip(a,b)]
def sv(a,c):return [p*c for p in a]
def sumv(vs):
    out=[P() for _ in range(4)]
    for v in vs:out=addv(out,v)
    return out
def mv(m,v):return [sum((v[j]*m[i,j] for j in range(len(v))),P()) for i in range(m.rows)]
def dot(a,b):return sum((x*y for x,y in zip(a,b)),P())
def flav(m,v):return [sumv([sv(v[j],m[i,j]) for j in range(m.cols)]) for i in range(m.rows)]
def bil(rows,m,cols):return sum((dot(a,b)*m[i,j] for i,a in enumerate(rows) for j,b in enumerate(cols)),P())
z=s.zeros(2)
pauli=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-I],[I,0]]),s.diag(1,-1)]
T=[m/2 for m in pauli]
gamma=[s.BlockMatrix([[z,s.eye(2)],[s.eye(2),z]]).as_explicit()]+[s.BlockMatrix([[z,p],[-p,z]]).as_explicit() for p in pauli]
C=I*gamma[2]*gamma[0]
sigma=[[I*(gamma[m]*gamma[n]-gamma[n]*gamma[m])/2 for n in range(4)] for m in range(4)]
def bar(v):return [sum((v[j].dag()*gamma[0][j,i] for j in range(4)),P()) for i in range(4)]
def cc(v):return mv(C,bar(v))
def A(v):return sumv([sv(mv(gamma[mu],[p.der(mu) for p in v]),I) for mu in range(4)])
def p0(v):return [p.p0() for p in v]
def kr(row):return [sum((p.der(mu)*(-I)*gamma[mu][j,i] for mu in range(4) for j,p in enumerate(row)),P()) for i in range(4)]
def hc(p):return p+p.dag()

n=2;ns=2
mom=(3,1,-1,1);wave=(1,1,0,0)
L=[]
for alpha in range(n):
    weak=[]
    for a in range(2):
        spin=[P() for _ in range(4)]
        for h in range(2):
            for sign in [-1,1]:
                g=len(conj);conj[g]=g+1;conj[g+1]=g
                spin[h]=spin[h]+P({(tuple(sign*x for x in mom),(g,)):1})
        weak.append(spin)
    L.append(weak)
phi=[P({(ZERO,()):R(2,3),(wave,()):1+I,(tuple(-x for x in wave),()):R(1,2)}),
     P({(ZERO,()):1-I,(wave,()):R(1,3),(tuple(-x for x in wave),()):-I})]
rho=sum((p.dag()*p for p in phi),P())
X=[sumv([sv(L[a][b],phi[b].dag()) for b in range(2)]) for a in range(n)]
bX=[bar(x) for x in X];Xc=[cc(x) for x in X]
Y=s.Matrix([[1+I,2-I],[3,I]])
seed=s.Matrix([[2,1+I],[1+I,3-I]])
F=seed*Y.conjugate();M=s.Matrix([[2,1],[1,3]])
m2=R(7,3);lam=R(2,5);gp=R(2,3);g=R(3,2);y=R(-1,2)
K=F*Y.T;K2=F*M*M*Y.T;V=F*M*Y.conjugate().T
KS=(K+K.T)/2;K2S=(K2+K2.T)/2
N=[[P() for _ in range(4)] for _ in range(ns)]
NR=[cc(v) for v in N]
SL=[[sumv([sv(NR[i],Y[a,i]*phi[b]) for i in range(ns)]) for b in range(2)] for a in range(n)]
bSL=[[bar(v) for v in weak] for weak in SL]
Sphi=[-sum((dot(bar(NR[i]),L[a][b])*Y[a,i].conjugate() for a in range(n) for i in range(ns)),P()) for b in range(2)]
SNL=flav(Y.conjugate().T,X);SR=[cc(v) for v in SNL]
RR=[addv(a,b) for a,b in zip(flav(M,NR),SR)]
TT=[sumv(v) for v in zip(flav(M*M,N),flav(M,SNL),[A(v) for v in SR])]
UU=[sumv(v) for v in zip(flav(M*M*M,NR),flav(M*M,SR),flav(M,[A(v) for v in SNL]),[p0(v) for v in SR])]

# Every line of P_red is evaluated after the literal zero-N substitution.
def placed(rows,cols,scalar=phi):
    return sum((dot(rows[a][b],cols[i])*F[a,i]*scalar[b] for a in range(n) for i in range(ns) for b in range(2)),P())
bL=[[bar(v) for v in weak] for weak in L]
terms={
 'sterile cubic source':bil(bX,F,UU)*R(5,6),
 'scalar coefficient source':bil(bX,F,RR)*(P.num(2*m2/3)+rho*lam/3),
 'row source derivative':placed([[kr(v) for v in weak] for weak in bSL],RR)/3,
 'scalar source':-placed(bL,RR,Sphi)/6,
 'row source quadratic sterile':-placed(bSL,TT)/6,
 'row source retained N':placed(bSL,N)*m2/2,
}
for name in list(terms)[2:]:check('vanishing full packet line: '+name,terms[name])
for label,rows in [('S_L',SL),('bar S_L',bSL)]:
    for a in range(n):
        for b in range(2):
            for h in range(4):check(f'projection {label} slot {a}{b}{h}',rows[a][b][h])
for a in range(2):check(f'projection scalar source {a}',Sphi[a])

# Generic Hermitian weak curvatures; the projection itself invokes no Bianchi
# or curvature EOM and therefore may test arbitrary antisymmetric jets.
strength=[[s.zeros(2) for _ in range(4)] for _ in range(4)]
for mu in range(4):
    for nu in range(mu+1,4):
        b=mu+2*nu+1
        w=[mu-nu+1,mu+nu+2,2*mu-nu]
        strength[mu][nu]=gp*y*b*s.eye(2)+g*sum((T[i]*w[i] for i in range(3)),s.zeros(2))
        strength[nu][mu]=-strength[mu][nu]
def curve(matrix,cols):
    return sum((dot(bL[a][b],mv(sigma[mu][nu],cols[j]))*matrix[a,j]*strength[mu][nu][b,c]*phi[c]
                for a in range(n) for j in range(len(cols)) for b in range(2) for c in range(2) for mu in range(4) for nu in range(4)),P())
terms['curvature source']=-curve(F,RR)/12
terms['Higgs derivative curvature']=curve(F,N) # actual spin/derivative factors are irrelevant once the explicit N slot is zero
check('zero explicit N curvature-derivative family',terms['Higgs derivative curvature'])
current=sum((dot(bL[a][b],mv(gamma[mu],N[i]))*F[a,i]*phi[b] for a in range(n) for b in range(2) for i in range(ns) for mu in range(4)),P())
check('zero explicit N gauge-current family',current)
terms['gauge currents']=current
full=sum(terms.values(),P())
scalar=lambda mat:bil(bX,mat,Xc)
box=lambda mat:bil(bX,mat,[p0(v) for v in Xc])
lnc=lambda mat:bil(bX,mat,[A(v) for v in X])
pred=scalar(K2)*R(5,6)+lnc(V)*R(5,6)+box(K)*R(5,6)+scalar(K)*(P.num(2*m2/3)+rho*lam/3)-curve(K,Xc)/12
check('complete active projection from actual R T U, all lines exhausted',full-pred)
check('scalar source symmetric flavor pointwise',scalar(K)-scalar(KS))
check('M squared scalar source symmetric flavor pointwise',scalar(K2)-scalar(K2S))
check('Box source symmetric flavor after integration',(box(K)-box(KS)).integral())
check('Hermitian LNC composite after integration',(hc(lnc(V))-lnc(V+V.conjugate().T)).integral())
reject('Box coefficient replaced by isolated packet value',(box(K)*R(5,6)-box(K)*R(1,4)).integral())
reject('M squared source removed',scalar(K2).integral())
reject('LNC mass source omitted',hc(lnc(V)).integral())
reject('quartic source omitted',(rho*scalar(K)).integral())
reject('weak tensor incorrectly flavor-antisymmetrized',(curve(K,Xc)-curve((K-K.T)/2,Xc)).integral())
reject('retained mass order reversed',(scalar(K2)-scalar(M*M*F*Y.T)).integral())

# Independent coefficient-space curvature projection: spin exchange is
# symmetric for sigma C, so antisymmetrize combined flavor/weak slots.
ph=s.Matrix(s.symbols('ph0 ph1'))
spin=(sigma[0][1]*C)[2:4,2:4]
KA=(K-K.T)/2
check('curvature chiral spin exchange sign',spin-spin.T)
for label,weakgen in [('B',s.eye(2)),('W1',T[0]),('W2',T[1]),('W3',T[2])]:
    wp=(weakgen*ph)*ph.T
    coeff=s.kronecker_product(K,wp,spin)
    target=s.kronecker_product(KA,(wp+wp.T)/2,spin)+s.kronecker_product(KS,(wp-wp.T)/2,spin)
    check(label+' combined weak flavor curvature projection',(coeff-coeff.T)/2-target)
    if label=='B':check('B curvature flavor-antisymmetric only',target-s.kronecker_product(KA,wp,spin))
    else:reject(label+' symmetric flavor weak channel discarded',target-s.kronecker_product(KA,(wp+wp.T)/2,spin))

# A Box on the whole neutral composite is not a Box on the lepton alone.
wrongXc=[cc(sumv([sv(p0(L[a][b]),phi[b].dag()) for b in range(2)])) for a in range(n)]
reject('scalar derivatives lost from composite Box',(box(K)-bil(bX,K,wrongXc)).integral())

# Independent flavor-pair factors in the integrated derivative and pointwise
# scalar sectors, including diagonal same-flavor Grassmann monomials.
def scalar_entry(a,b,derivative=False):return dot(bX[a],p0(Xc[b]) if derivative else Xc[b])
for mat,func,label,integ in [(K2S,scalar,'mass',False),(KS,box,'Box',True)]:
    explicit=sum((scalar_entry(a,a,integ)*mat[a,a] for a in range(n)),P())
    explicit+=sum((scalar_entry(a,b,integ)*2*mat[a,b] for a in range(n) for b in range(a+1,n)),P())
    diff=func(mat)-explicit
    check('independent pair normalization '+label,diff.integral() if integ else diff)
    reject('same-flavor scalar incorrectly discarded '+label,scalar_entry(0,0,integ).integral())

# Actual residual covectors at N=0, not a generic chain-rule placeholder.
W=[]
for i in range(ns):W.append([sum((bX[a][h]*F[a,i] for a in range(n)),P()) for h in range(4)])
bx=[]
for i in range(ns):
    row=sumv([sv(W[j],R(5,6)*(M*M)[j,i]) for j in range(ns)])
    row=addv(row,sv([p.p0() for p in W[i]],R(5,6)))
    row=addv(row,sv(W[i],P.num(2*m2/3)+rho*lam/3))
    cur=[sum((bL[a][b][h]*sigma[mu][nu][h,k]*strength[mu][nu][b,c]*phi[c]*F[a,i]
              for a in range(n) for b in range(2) for c in range(2) for mu in range(4) for nu in range(4) for h in range(4)),P()) for k in range(4)]
    bx.append(addv(row,sv(cur,-R(1,12))))
by=[sumv([sv(kr(W[j]),R(5,6)*M[j,i]) for j in range(ns)]) for i in range(ns)]
certificate=sum((dot(bx[i],SR[i])+dot(by[i],SNL[i]) for i in range(ns)),P())
check('active projection equals actual sterile residual covectors',(full-certificate).integral())
reject('original zero-N slice equated to reduced zero-N slice',hc(full).integral())

# The additional scalar-current EOM transport and projection do not commute.
d=(g*g+gp*gp)/24
extra=-rho*scalar(K)*d
eR=[sv(v,-1) for v in SR]
residual=-rho*bil(bX,F,eR)*d
check('zero-N scalar-current residual certificate',-extra-residual)
# Its actual additional sterile shift is positive d rho F†X plus conjugate.
zshift=flav(F.conjugate().T,X)
deltaNL=[sv(v,rho*d) for v in zshift]
deltaNR=[cc(v) for v in deltaNL]
delta_yuk=-sum((dot(bL[a][b],deltaNR[i])*Y[a,i]*phi[b] for a in range(n) for b in range(2) for i in range(ns)),P())
check('extra active scalar from actual Majorana Yukawa variation',(hc(delta_yuk)-hc(extra)).integral())
reject('projection and scalar-current EOM reduction forced to commute',hc(extra).integral())
reject('extra scalar source sign reversed',(hc(delta_yuk)+hc(extra)).integral())

result={'status':'pass','python':sys.version.split()[0],'sympy':s.__version__,
 'scope':'Complete literal zero-retained-N field-content projection of the declared selected sourced representative; not integration, EOM elimination, matching or a physical coefficient.',
 'exact_checks':len(checks),'rejected_mutations':len(negatives),'grassmann_generators':len(conj),
 'scalar_modes_per_component':3,'active_projection_density_monomials':len(full.d),
 'active_projection_integrated_monomials':len(full.integral().d),'exhausted_packet_lines':list(terms),
 'checks':checks,'mutations':negatives}
print(json.dumps(result,indent=2))
