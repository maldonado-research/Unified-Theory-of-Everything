#!/usr/bin/env python3
"""Independent exact controls for the action-level selected-packet shifts.
No repository imports, loop integrals, floating-point arithmetic, or network use.
"""
import json
import random
import sys
from collections import Counter
from pathlib import Path
import sympy as s

I=s.I
R=s.Rational
eta=(1,-1,-1,-1)
z2=s.zeros(2)
pauli=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-I],[I,0]]),s.diag(1,-1)]
g=[s.diag(1,1,-1,-1)]+[z2.row_join(p).col_join((-p).row_join(z2)) for p in pauli]
g5=I*g[0]*g[1]*g[2]*g[3]
PL=(s.eye(4)-g5)/2
PR=(s.eye(4)+g5)/2
C=I*g[2]*g[0]
counts=Counter()
controls=[]
rng=random.Random(15102026)

def zero(q):
    if isinstance(q,s.MatrixBase):
        return all(s.expand(x)==0 for x in q)
    if "Grass" in globals() and isinstance(q,Grass):
        return not q.d
    return s.expand(q)==0

def check(name,q):
    if not zero(q):
        raise RuntimeError('check failed: '+name+'; '+str(q))
    counts[name]+=1

def reject(name,q):
    if zero(q):
        raise RuntimeError('negative control failed to reject: '+name)
    controls.append(name)

def scalar():
    return R(rng.randrange(-5,6),rng.randrange(1,5))+I*R(rng.randrange(-4,5),rng.randrange(1,4))

def real():
    return R(rng.randrange(-5,6),rng.randrange(1,5))

def row(P):
    return s.Matrix([[scalar() for _ in range(4)]])*P

def total(xs,rows=1,cols=4):
    return sum(xs,s.zeros(rows,cols))

def adj(row_or_col):
    if row_or_col.rows==1:
        return g[0]*row_or_col.conjugate().T
    return row_or_col.conjugate().T*g[0]

def slash(v):
    # v has covariant components, so no extra metric factor here.
    return total((g[m]*v[m] for m in range(4)),4,4)

for m in range(4):
    for n in range(4):
        check('Clifford convention',g[m]*g[n]+g[n]*g[m]-(2*eta[m]*s.eye(4) if m==n else s.zeros(4)))
    check('charge-conjugation convention',C*g[m].T*C.inv()+g[m])

# Functional differentiation of the charged scalar action. Its source is
# denoted q=bar N_R Y^dagger L and not preassigned an equation-of-motion sign.
coords=s.symbols('x0:4',real=True)
f=s.Function('phi')(*coords)
fc=s.Function('phi_conj')(*coords)
b=[s.Function('B'+str(m))(*coords) for m in range(4)]
q=s.Function('q')(*coords)
qc=s.Function('q_conj')(*coords)
gy,m2,lam=s.symbols('g_y m_H_squared lambda',real=True)
kin=sum(eta[m]*(s.diff(fc,coords[m])-I*gy*b[m]*fc)*(s.diff(f,coords[m])+I*gy*b[m]*f) for m in range(4))
lag=kin-m2*fc*f-lam*(fc*f)**2-q*fc-qc*f
EL=s.diff(lag,fc)-sum(s.diff(s.diff(lag,s.diff(fc,coords[m])),coords[m]) for m in range(4))
D2f=sum(eta[m]*(s.diff(f,coords[m],2)+I*gy*s.diff(b[m],coords[m])*f+2*I*gy*b[m]*s.diff(f,coords[m])-gy**2*b[m]**2*f) for m in range(4))
check('literal scalar Euler source sign',EL-(-D2f-m2*f-2*lam*fc*f*f-q))
reject('positive rather than negative scalar Yukawa source',EL-(-D2f-m2*f-2*lam*fc*f*f+q))

# Maxwell variation: derivative with respect to a lower-index B_mu gives
# the upper-index residual partial_nu B^{nu mu}.
strength=[[s.diff(b[n],coords[m])-s.diff(b[m],coords[n]) for n in range(4)] for m in range(4)]
maxwell=-R(1,4)*sum(eta[m]*eta[n]*strength[m][n]**2 for m in range(4) for n in range(4))
for m in range(4):
    eb=s.diff(maxwell,b[m])-sum(s.diff(s.diff(maxwell,s.diff(b[m],coords[n])),coords[n]) for n in range(4))
    expected=sum(eta[n]*eta[m]*s.diff(strength[n][m],coords[n]) for n in range(4))
    check('Maxwell Euler index and sign',eb-expected)

# Fourier verification of the light-field integrated adjoint, with a
# nonzero U(1) connection. All momentum components are covariant.
for trial in range(8):
    h=[real() for _ in range(4)]
    p=[real() for _ in range(4)]
    er=[-h[m]-p[m] for m in range(4)]
    conn=[real() for _ in range(4)]
    charge=R(-3,7)
    eb=row(PL)
    n=adj(row(PR))  # left-chiral column
    ph=scalar()
    mn=real()
    an=-slash(p)*n
    a2n=slash(p)**2*n
    ker=eb*slash([er[m]-charge*conn[m] for m in range(4)])
    original=R(1,3)*ker*ph*an-R(1,6)*eb*ph*a2n+mn*eb*ph*n/2
    dph=[I*(h[m]+charge*conn[m])*ph for m in range(4)]
    ell=I*total((g[m]*dph[m]*an for m in range(4)),4,1)/3+ph*a2n/6+mn*ph*n/2
    check('lepton integrated adjoint with connection',original-eb*ell)
    check('lepton chirality',PL*ell-ell)
    if trial==0:
        reject('unintegrated minus one-sixth kept in lepton shift',original-eb*(ell-ph*a2n/3))
        reject('wrong covariant connection in integrated adjoint',original-eb*(ell+R(2,3)*charge*ph*slash(conn)*an))

# Jet tests for the sterile rows: compare direct integration by parts,
# ordinary derivatives of neutral composites, expanded covariant shifts,
# and their actual Dirac adjoints. The jets are arbitrary and nonzero;
# in particular they are not solutions of the reference equations.
for trial in range(6):
    na,ns,nw=2,3,2
    F=s.Matrix(na,ns,lambda a,i:scalar())
    M=s.diag(2,3,5)
    if trial%2:
        M=s.Matrix([[2,1,-1],[1,4,2],[-1,2,6]])
    FM=F*M
    FM2=F*M*M
    charge=R(-2,5)
    conn=[real() for _ in range(4)]
    divB=real()
    B2=sum(eta[m]*conn[m]**2 for m in range(4))
    ph=[scalar() for a in range(nw)]
    dph=[[scalar() for m in range(4)] for a in range(nw)]
    boxph=[scalar() for a in range(nw)]
    dphcov=[[dph[a][m]+I*charge*conn[m]*ph[a] for m in range(4)] for a in range(nw)]
    boxphcov=[boxph[a]+I*charge*divB*ph[a]+2*I*charge*sum(eta[m]*conn[m]*dph[a][m] for m in range(4))-charge**2*B2*ph[a] for a in range(nw)]
    br=[[row(PR) for a in range(nw)] for al in range(na)]
    db=[[[row(PR) for m in range(4)] for a in range(nw)] for al in range(na)]
    boxbr=[[row(PR) for a in range(nw)] for al in range(na)]
    sr=[[row(PL) for a in range(nw)] for al in range(na)]
    ds=[[[row(PL) for m in range(4)] for a in range(nw)] for al in range(na)]
    sph=[scalar() for a in range(nw)]
    dbcov=[[[db[al][a][m]-I*charge*conn[m]*br[al][a] for m in range(4)] for a in range(nw)] for al in range(na)]
    dscov=[[[ds[al][a][m]-I*charge*conn[m]*sr[al][a] for m in range(4)] for a in range(nw)] for al in range(na)]
    boxbrcov=[[boxbr[al][a]-I*charge*divB*br[al][a]-2*I*charge*total((eta[m]*conn[m]*db[al][a][m] for m in range(4)))-charge**2*B2*br[al][a] for a in range(nw)] for al in range(na)]
    tensor=s.zeros(4)
    for m in range(4):
        for n in range(m+1,4):
            # Both antisymmetric slots are included in sigma^{mu nu} B_mu nu.
            tensor+=I*(g[m]*g[n]-g[n]*g[m])*real()
    cv=real()
    W=[total((F[al,i]*br[al][a]*ph[a] for al in range(na) for a in range(nw))) for i in range(ns)]
    V=[total((F[al,i]*sr[al][a]*ph[a] for al in range(na) for a in range(nw))) for i in range(ns)]
    KW=[-I*total((F[al,i]*(db[al][a][m]*ph[a]+br[al][a]*dph[a][m])*g[m] for al in range(na) for a in range(nw) for m in range(4))) for i in range(ns)]
    K2W=[-total((F[al,i]*(boxbr[al][a]*ph[a]+2*total((eta[m]*db[al][a][m]*dph[a][m] for m in range(4)))+br[al][a]*boxph[a]) for al in range(na) for a in range(nw))) for i in range(ns)]
    KV=[-I*total((F[al,i]*(ds[al][a][m]*ph[a]+sr[al][a]*dph[a][m])*g[m] for al in range(na) for a in range(nw) for m in range(4))) for i in range(ns)]
    KS=[-I*total((F[al,i]*dscov[al][a][m]*g[m]*ph[a] for al in range(na) for a in range(nw) for m in range(4))) for i in range(ns)]
    Qphi=[total((F[al,i]*br[al][a]*sph[a] for al in range(na) for a in range(nw))) for i in range(ns)]
    for i in range(ns):
        compactx=R(5,6)*(total((W[j]*(M*M)[j,i] for j in range(ns)))+K2W[i])+cv*W[i]+KS[i]/3-Qphi[i]/6-KV[i]/6-charge*W[i]*tensor/12
        compacty=R(5,6)*total((KW[j]*M[j,i] for j in range(ns)))-total((V[j]*M[j,i] for j in range(ns)))/6
        expandedx=R(5,6)*total((FM2[al,i]*br[al][a]*ph[a] for al in range(na) for a in range(nw)))+cv*W[i]
        expandedx-=R(5,6)*total((F[al,i]*(boxbrcov[al][a]*ph[a]+2*total((eta[m]*dbcov[al][a][m]*dphcov[a][m] for m in range(4)))+br[al][a]*boxphcov[a]) for al in range(na) for a in range(nw)))
        expandedx-=I*total((F[al,i]*dscov[al][a][m]*g[m]*ph[a] for al in range(na) for a in range(nw) for m in range(4)))/6
        expandedx+=I*total((F[al,i]*sr[al][a]*g[m]*dphcov[a][m] for al in range(na) for a in range(nw) for m in range(4)))/6
        expandedx-=Qphi[i]/6+charge*W[i]*tensor/12
        expandedy=-5*I*total((FM[al,i]*(dbcov[al][a][m]*ph[a]+br[al][a]*dphcov[a][m])*g[m] for al in range(na) for a in range(nw) for m in range(4)))/6-total((FM[al,i]*sr[al][a]*ph[a] for al in range(na) for a in range(nw)))/6
        check('expanded sterile x row',compactx-expandedx)
        check('expanded sterile y row',compacty-expandedy)
        # Independently assemble the column formula (12).
        xcol=R(5,6)*total((s.conjugate(FM2[al,i])*s.conjugate(ph[a])*adj(br[al][a]) for al in range(na) for a in range(nw)),4,1)+cv*adj(W[i])
        xcol-=R(5,6)*total((s.conjugate(F[al,i])*(s.conjugate(ph[a])*adj(boxbrcov[al][a])+2*total((eta[m]*s.conjugate(dphcov[a][m])*adj(dbcov[al][a][m]) for m in range(4)),4,1)+s.conjugate(boxphcov[a])*adj(br[al][a])) for al in range(na) for a in range(nw)),4,1)
        xcol+=I*total((s.conjugate(F[al,i])*s.conjugate(ph[a])*g[m]*adj(dscov[al][a][m]) for al in range(na) for a in range(nw) for m in range(4)),4,1)/6
        xcol-=I*total((s.conjugate(F[al,i])*s.conjugate(dphcov[a][m])*g[m]*adj(sr[al][a]) for al in range(na) for a in range(nw) for m in range(4)),4,1)/6
        xcol-=total((s.conjugate(F[al,i])*s.conjugate(sph[a])*adj(br[al][a]) for al in range(na) for a in range(nw)),4,1)/6
        xcol-=charge*total((s.conjugate(F[al,i])*s.conjugate(ph[a])*tensor*adj(br[al][a]) for al in range(na) for a in range(nw)),4,1)/12
        ycol=5*I*total((s.conjugate(FM[al,i])*g[m]*(s.conjugate(ph[a])*adj(dbcov[al][a][m])+s.conjugate(dphcov[a][m])*adj(br[al][a])) for al in range(na) for a in range(nw) for m in range(4)),4,1)/6
        ycol-=total((s.conjugate(FM[al,i])*s.conjugate(ph[a])*adj(sr[al][a]) for al in range(na) for a in range(nw)),4,1)/6
        check('explicit x column adjoint',xcol-adj(compactx))
        check('explicit y column adjoint',ycol-adj(compacty))
        check('x left chirality',PL*xcol-xcol)
        check('y right chirality',PR*ycol-ycol)
        if trial==0 and i==0:
            reject('wrong sterile integrated first-derivative sign',compacty-(-R(5,6)*total((KW[j]*M[j,i] for j in range(ns)))-total((V[j]*M[j,i] for j in range(ns)))/6))
            reject('missing differentiated lepton source',compactx-(compactx+KV[i]/6))
            reject('missing scalar source term',compactx-(compactx+Qphi[i]/6))
            reject('wrong sterile source-derivative sign',compactx-(compactx+KV[i]/3))
        if trial==0 and i==1:
            reject('wrong mass flavor slot',compacty-(R(5,6)*KW[i]*M[0,0]-V[i]*M[0,0]/6))

# Independently start with the raw sterile terms of Eq(5), let all row
# composites carry momentum q and residual columns momentum -q, and only
# then collect residual columns. This tests both source slots, both sterile
# chiralities, two derivatives, and ordered mass contractions.
for trial in range(6):
    ns=3
    M=s.Matrix([[2,1,-1],[1,4,2],[-1,2,6]])
    mom=[real() for _ in range(4)]
    aq=slash(mom)
    W=[row(PR) for _ in range(ns)]
    V=[row(PL) for _ in range(ns)]
    KS=[row(PR) for _ in range(ns)]
    QS=[row(PR) for _ in range(ns)]
    eR=[adj(row(PL)) for _ in range(ns)]
    eL=[adj(row(PR)) for _ in range(ns)]
    cv=real()
    charge=R(-2,5)
    tensor=I*(g[0]*g[1]-g[1]*g[0])*R(2,3)
    original=s.zeros(1)
    for i in range(ns):
        m2er=total(((M*M)[i,j]*eR[j] for j in range(ns)),4,1)
        mael=total((M[i,j]*aq*eL[j] for j in range(ns)),4,1)
        mel=total((M[i,j]*eL[j] for j in range(ns)),4,1)
        original+=R(5,6)*W[i]*(m2er+mael+aq**2*eR[i])+cv*W[i]*eR[i]
        original+=KS[i]*eR[i]/3-QS[i]*eR[i]/6-V[i]*(mel+aq*eR[i])/6-charge*W[i]*tensor*eR[i]/12
    factored=s.zeros(1)
    wrong=s.zeros(1)
    for i in range(ns):
        xr=R(5,6)*(total((W[j]*(M*M)[j,i] for j in range(ns)))+W[i]*aq**2)+cv*W[i]+KS[i]/3-QS[i]/6-V[i]*aq/6-charge*W[i]*tensor/12
        yl=R(5,6)*total((W[j]*aq*M[j,i] for j in range(ns)))-total((V[j]*M[j,i] for j in range(ns)))/6
        factored+=xr*eR[i]+yl*eL[i]
        wrong+=xr*eR[i]-yl*eL[i]
    check('raw sterile residual Fourier factorization',original-factored)
    if trial==0:
        reject('wrong sterile chiral residual coefficient',original-wrong)

# Simple finite exterior algebra: dagger reverses products and pairs each
# complex Grassmann generator with its adjoint. This avoids treating
# fermions as commuting vectors in the Majorana-factor checks.
class Grass:
    def __init__(self,d=None):
        self.d={k:s.expand(v) for k,v in (d or {}).items() if s.expand(v)!=0}
    @staticmethod
    def cast(x):
        return x if isinstance(x,Grass) else Grass({():s.sympify(x)})
    def __add__(self,other):
        other=Grass.cast(other)
        d=dict(self.d)
        for k,v in other.d.items():d[k]=d.get(k,0)+v
        return Grass(d)
    __radd__=__add__
    def __neg__(self):return Grass({k:-v for k,v in self.d.items()})
    def __sub__(self,other):return self+-Grass.cast(other)
    def __rsub__(self,other):return Grass.cast(other)+-self
    def __mul__(self,other):
        other=Grass.cast(other)
        d={}
        for a,x in self.d.items():
            for b,y in other.d.items():
                if set(a)&set(b):continue
                sign=(-1)**sum(i>j for i in a for j in b)
                k=tuple(sorted(a+b))
                d[k]=d.get(k,0)+sign*x*y
        return Grass(d)
    __rmul__=__mul__
    def dagger(self):
        out=Grass()
        for k,v in self.d.items():
            term=Grass.cast(s.conjugate(v))
            for j in reversed(k):term=term*Grass({(j^1,):1})
            out+=term
        return out
    def __repr__(self):return str(self.d)

def matvec(mat,v):
    return [sum((v[j]*mat[i,j] for j in range(mat.cols)),Grass()) for i in range(mat.rows)]

def bar(v):
    vd=[x.dagger() for x in v]
    return [sum((vd[j]*g[0][j,i] for j in range(4)),Grass()) for i in range(4)]

def conjc(v):return matvec(C,bar(v))
def addvec(a,b):return [x+y for x,y in zip(a,b)]
def dot(a,b):return sum((x*y for x,y in zip(a,b)),Grass())
def vdiff(a,b):
    for i in range(4):check('Majorana reality of completed shift',a[i]-b[i])

w=[Grass({(2*i,):1}) for i in range(4)]
z=[Grass({(8+2*i,):1}) for i in range(4)]
E=addvec(z,conjc(z))
wreal=addvec(w,conjc(w))
vdiff(wreal,conjc(wreal))
packet=dot(bar(w),E)
realfact=dot(bar(wreal),E)
check('Grassmann Majorana residual plus conjugate',packet+packet.dagger()-realfact)
reject('half-sized Majorana shift',packet+packet.dagger()-realfact*R(1,2))
reject('omitted Majorana conjugate completion',w[0]-conjc(w)[0])
# Mass-slot variation of 1/2 bar N N, with independent Majorana variation.
mass_half=(dot(bar(wreal),E)+dot(bar(E),wreal))*R(1,2)
check('Majorana mass half cancels two slots',mass_half-dot(bar(wreal),E))
for mu in range(4):
    AdN=matvec(I*g[mu],E)
    first=dot(bar(wreal),AdN)
    # The second kinetic-slot variation is integrated by parts before
    # invoking the Grassmann charge-conjugation identity.
    second=-dot(bar(E),matvec(I*g[mu],wreal))
    check('Majorana kinetic half cancels after IBP',(first+second)*R(1,2)-first)

# Check that the light/gauge conjugate completion cancels a complex packet.
a,bv,ev=s.symbols('a b e',real=True)
j=a+I*bv
res=-R(1,6)*gy*j*ev
shift=gy*(j+s.conjugate(j))/6
check('real gauge shift cancels both conjugate channels',res+s.conjugate(res)+ev*shift)
reject('wrong sign gauge shift',res+s.conjugate(res)-ev*shift)

# Homogeneous neutral-background projection of Eq(12), compared to the
# independently compact form in Eq(17). A acts as -slash(p) on columns.
for trial in range(4):
    F=s.Matrix(2,3,lambda al,i:scalar())
    Y=s.Matrix(2,3,lambda al,i:scalar())
    M=s.diag(2,3,5)
    u=R(3,2)
    cv=R(5,7)
    p=[real() for _ in range(4)]
    aa=-slash(p)
    nu=[adj(row(PR)) for al in range(2)]
    nr=[adj(row(PL)) for i in range(3)]
    src=[u*total((Y[al,j]*nr[j] for j in range(3)),4,1) for al in range(2)]
    for i in range(3):
        # All derivatives of the scalar background vanish; S_phi contains
        # two fermions, so it cannot contribute to x at linear degree.
        xdirect=total((s.conjugate((F*M*M)[al,i])*u*nu[al]*R(5,6)+s.conjugate(F[al,i])*u*(R(5,6)*aa**2+cv*s.eye(4))*nu[al]+s.conjugate(F[al,i])*u*aa*src[al]/6 for al in range(2)),4,1)
        xexpected=u*total((s.conjugate((F*M*M)[al,i])*nu[al]*R(5,6)+s.conjugate(F[al,i])*(R(5,6)*aa**2+cv*s.eye(4))*nu[al] for al in range(2)),4,1)+u**2*total((s.conjugate(F[al,i])*Y[al,j]*aa*nr[j]/6 for al in range(2) for j in range(3)),4,1)
        ydirect=total((s.conjugate((F*M)[al,i])*(5*u*aa*nu[al]/6-u*src[al]/6) for al in range(2)),4,1)
        yexpected=5*u*total((s.conjugate((F*M)[al,i])*aa*nu[al] for al in range(2)),4,1)/6-u**2*total((s.conjugate((F*M)[al,i])*Y[al,j]*nr[j] for al in range(2) for j in range(3)),4,1)/6
        check('neutral background x linearization',xdirect-xexpected)
        check('neutral background y linearization',ydirect-yexpected)
        if trial==0 and i==1:
            wrong=total((s.conjugate((F*M)[al,i])*(5*u*aa*nu[al]/6+u*src[al]/6) for al in range(2)),4,1)
            reject('background source sign reversed',ydirect-wrong)

# Direct coefficient audit of the unreduced neutral quadratic packet.
p2=s.symbols('p_squared',real=True)
coefficient=R(5,6)*p2+R(1,3)*p2-R(1,6)*p2+m2/2+m2/2
check('original neutral quadratic packet',coefficient-(p2+m2))
# Off-stationary Hessian warning: the scalar Euler derivative at constant
# phi is exactly the tadpole multiplying the quadratic scalar shift.
u=s.symbols('u',real=True)
check('constant scalar tadpole',(-m2*u-2*lam*u**3)-(-(m2+2*lam*u**2)*u))
reject('dropping off-stationary scalar tadpole',-(m2+2*lam*u**2)*u)

result={
    'status':'passed',
    'scope':'Exact action algebra for selected residual packet; no new loop integration or physical matching.',
    'source_candidate':'89f5867453bc7ed706d249dc82b393245b28d7ea',
    'python':sys.version.split()[0],
    'sympy':s.__version__,
    'checks':dict(sorted(counts.items())),
    'total_checks':sum(counts.values()),
    'negative_controls_rejected':controls,
    'negative_control_count':len(controls),
    'trial_configuration':{'fourier_trials':8,'jet_trials':6,'active_flavors':2,'sterile_flavors':3,'weak_components':2,'neutral_background_trials':4,'raw_sterile_residual_fourier_trials':6},
    'limits':['Does not derive Eq(5) from the loop integrals.','Does not calculate the full EFT counterterm vector.','Does not replace a complete gauge-fixed quantum change-of-variables analysis.']
}
out=Path(__file__).with_name('FIELD_REDEFINITION_RESULTS.json')
out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps(result,indent=2,sort_keys=True))
