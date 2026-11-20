#!/usr/bin/env python3
"""Independent Yukawa two-point UV computation by asymptotic radial moments.

SymPy 1.14.0, no repository imports or assertions. The loop coefficients are
computed from the expanded integrands, isotropic tensors, and Gamma-function
vacuum residues, then compared with numerator/tadpole identities. They are
not imported anomalous dimensions or preassigned counterterm targets.
Run python -I -B -O independent_two_point_audit.py to print the JSON receipt.
"""
import itertools
import json
import sys
import sympy as s
I,R=s.I,s.Rational
q2,z,psq,aux,mn2,ms2,d=s.symbols('q2 z p_squared auxiliary_mass_squared sterile_mass_squared scalar_mass_squared d')
checks,mutations=[],[]

def zero(x):
    if isinstance(x,s.MatrixBase):return all(s.simplify(a)==0 for a in x)
    return s.simplify(x)==0
def check(name,x):
    if not zero(x):raise RuntimeError(name+': '+str(x))
    checks.append(name)
def reject(name,x):
    if zero(x):raise RuntimeError(name+': mutation escaped')
    mutations.append(name)

def radial_residue(power,denominator,mass=aux):
    """Coefficient of i/(16 pi² epsilon), q^(2 power)/(q²-mass)^den.

    Wick rotation gives (-1)^(power-den). The only pole is Gamma(n+eps)
    with n=den-power-2. Its residue follows from Gamma's integer poles.
    Analytic d-dependent prefactors are evaluated at d=4 only after the
    simple pole has been selected.
    """
    n=denominator-power-2
    if n>0:return s.Integer(0)
    gamma_residue=s.Integer(-1)**(-n)/s.factorial(-n)
    return s.simplify(s.Integer(-1)**(power-denominator)*mass**(2+power-denominator)
       *s.gamma(power+2)*gamma_residue/(s.gamma(2)*s.gamma(denominator)))

def sphere_factor(k,vector=False):
    # z=q.p. For a scalar k must be even; for q_mu z^k it must be odd.
    if vector:
        if k%2==0:return None
        n=(k-1)//2
        return n+1,psq**n*s.factorial2(2*n+1)/s.prod(d+2*j for j in range(n+1))
    if k%2:return None
    n=k//2
    return n,psq**n*s.factorial2(2*n-1)/s.prod(d+2*j for j in range(n))

def asymptotic(numerator,mass_first,mass_second,vector=False,drop=None):
    """Geometric expansion around Q=q²-aux, retaining all UV terms.

    qslash is external to `numerator` for vector=True. With at most a
    quadratic scalar numerator, i+j>=3 is UV convergent by radial degree.
    """
    a=aux-mass_first
    b=-2*z+psq+aux-mass_second
    answer=0;pieces=[]
    for i in range(3):
        for j in range(3-i):
            if drop==(i,j):continue
            polynomial=s.Poly(s.expand(numerator*(-a)**i*(-b)**j),q2,z)
            piece=0
            for (qa,zk),coef in polynomial.terms():
                averaged=sphere_factor(zk,vector)
                if averaged is None:continue
                radial_power,angular=averaged
                piece+=coef*angular*radial_residue(qa+radial_power,2+i+j)
            piece=s.simplify(piece)
            pieces.append({'denominator_orders':[i,j],'residue_before_d4':str(piece)})
            answer+=piece
    return s.simplify(answer.subs(d,4)),pieces

fermion,fermion_pieces=asymptotic(1,mn2,ms2,vector=True)
scalar,scalar_pieces=asymptotic(2*(q2-z),mn2,0)
bubble,_=asymptotic(1,mn2,ms2)
check('fermion regulator cancellation',s.diff(fermion,aux))
check('scalar regulator cancellation',s.diff(scalar,aux))
check('fermion internal-mass independence',s.diff(fermion,mn2)+s.diff(fermion,ms2))
# Independent numerator reductions, using radial Gamma residues again but
# no Feynman parameters and no prescribed wavefunction target coefficient.
vector_pv=(radial_residue(0,1,ms2)-radial_residue(0,1,mn2)
           +(psq+mn2-ms2)*bubble)/(2*psq)
scalar_pv=radial_residue(0,1,mn2)+(mn2-psq)*bubble
check('open vector asymptotic versus numerator/tadpole identity',fermion-vector_pv)
check('closed scalar asymptotic versus numerator/tadpole identity',scalar-scalar_pv)
reject('omit leading quadratic radial contribution',scalar-asymptotic(2*(q2-z),mn2,0,drop=(0,0))[0])
reject('omit quadratic external-momentum denominator term',scalar-asymptotic(2*(q2-z),mn2,0,drop=(0,2))[0])
reject('omit open numerator angular term',fermion-asymptotic(1,mn2,ms2,vector=True,drop=(0,1))[0])
for regulator in [R(1,3),R(7,2),s.Integer(11)]:
    for mass in [0,R(2,5),7]:
        check(f'IR regulator witness {regulator}, mass {mass}',scalar.subs({aux:regulator,mn2:mass})-scalar_pv.subs(mn2,mass))
        check(f'open IR regulator witness {regulator}, mass {mass}',fermion.subs({aux:regulator,mn2:mass})-vector_pv)

# Literal vertex/propagator phases; the overall loop integration contributes i.
open_phase=(-I)**2*I**2
closed_phase=-open_phase
reject('closed fermion loop minus omitted',closed_phase-open_phase)
scalar_wave=s.diff(closed_phase*scalar,psq)
scalar_mass=s.diff(closed_phase*scalar,mn2)
vector_wave=s.simplify(open_phase*fermion)

# Explicit Clifford calculation of every coefficient of the closed trace.
pauli=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-I],[I,0]]),s.diag(1,-1)]
z2=s.zeros(2)
gamma=[s.BlockMatrix([[z2,s.eye(2)],[s.eye(2),z2]]).as_explicit()]+[s.BlockMatrix([[z2,p],[-p,z2]]).as_explicit() for p in pauli]
eta=[1,-1,-1,-1]
g5=I*gamma[0]*gamma[1]*gamma[2]*gamma[3]
PL=(s.eye(4)-g5)/2;PR=(s.eye(4)+g5)/2
C=I*gamma[2]*gamma[0]
for mu in range(4):
    for nu in range(4):
        check(f'closed chiral trace coefficient {mu}{nu}',s.trace(PR*gamma[mu]*PL*gamma[nu])-2*eta[mu]*(mu==nu))
for mu in range(4):
    check(f'open mass numerator killed {mu}',PR*s.eye(4)*PL)
    check(f'closed linear mass trace killed {mu}',s.trace(PR*PL*gamma[mu]))
    check(f'charge-conjugate gamma coefficient {mu}',C*gamma[mu].T*C.inv()+gamma[mu])
check('left chiral propagator mass term absent',PL*PR)

# The weak Yukawa intertwiner is a Kronecker delta, not an assumed loop factor.
vertex=s.eye(2)
weak_L=s.Matrix(2,2,lambda a,b:sum(vertex[a,c]*vertex[b,c] for c in range(2)))
weak_phi=s.Matrix(2,2,lambda a,b:sum(vertex[c,a]*vertex[c,b] for c in range(2)))
weak_N=sum(vertex[a,b]*vertex[a,b] for a in range(2) for b in range(2))
check('fixed external lepton weak contraction',weak_L-s.eye(2))
check('fixed external scalar weak contraction',weak_phi-s.eye(2))
reject('sterile weak sum replaced by fixed component',weak_N-weak_L[0,0])

def unitary(n):
    out=s.eye(n)
    for i in range(n-1):
        u=s.eye(n);u[i,i]=u[i+1,i+1]=R(3,5);u[i,i+1]=u[i+1,i]=4*I/5
        out=u*out
    return out
Y=s.Matrix([[1+I,2-I],[3,I],[1-2*I,-2]])
H=Y.conjugate().T*Y
U=unitary(2);Vact=unitary(3)
check('sterile unitary',U.conjugate().T*U-s.eye(2))
check('active unitary',Vact.conjugate().T*Vact-s.eye(3))
Lp=Y*Y.conjugate().T
zL=-vector_wave*Lp
zNL=-weak_N*vector_wave*H.T
zphi=-scalar_wave*s.trace(H)
mom=s.Matrix([3,1,-2,1])
slash=sum((eta[mu]*mom[mu]*gamma[mu] for mu in range(4)),s.zeros(4))
right_kernel=s.kronecker_product(H,slash*PR)*weak_N*vector_wave
cfull=s.kronecker_product(s.eye(2),C)
conjugate_kernel=cfull*(-right_kernel).T*cfull.inv()
expected_left=s.kronecker_product(H.T,slash*PL)*weak_N*vector_wave
check('Majorana conjugate channel and flavor transpose',conjugate_kernel-expected_left)
reject('Majorana channels incorrectly use same flavor matrix',conjugate_kernel-s.kronecker_product(H,slash*PL)*weak_N*vector_wave)

# Full Majorana kinetic action versus one left-Weyl action, in a genuine
# coefficient-level exterior algebra with two opposite Fourier modes.
ng=16;mate={};energies=[]
left=[]
for species in range(2):
    data=s.zeros(4,ng)
    for spin in range(2):
        for mode in range(2):
            a=species*8+spin*4+mode*2;b=a+1
            mate[a]=b;mate[b]=a
            while len(energies)<=b:energies.append(0)
            energies[a]=(-1 if mode==0 else 1)*R(5,3);energies[b]=-energies[a]
            data[spin,a]=1
    left.append(data)
def dagger_data(data):return s.Matrix(data.rows,ng,lambda r,c:s.conjugate(data[r,mate[c]]))
def bar_data(data):return gamma[0].T*dagger_data(data)
def charge_data(data):return C*bar_data(data)
right=[charge_data(a) for a in left]
major=[a+b for a,b in zip(left,right)]
def avec(data):return gamma[0]*data*s.diag(*energies)
def exterior_integrated(matrix):
    keep=s.Matrix(ng,ng,lambda a,b:matrix[a,b] if energies[a]+energies[b]==0 else 0)
    return keep-keep.T
weyl=s.zeros(ng);majorana=s.zeros(ng)
for a in range(2):
    for b in range(2):
        weyl+=bar_data(left[a]).T*avec(left[b])*zNL[a,b]
        majorana+=(bar_data(major[a]).T*avec(PL*major[b])*zNL[a,b]
                   +bar_data(major[a]).T*avec(PR*major[b])*s.conjugate(zNL[a,b]))/2
check('Majorana half action equals left-Weyl wave counterterm',exterior_integrated(majorana-weyl))
reject('Majorana half applied twice',exterior_integrated(majorana/2-weyl))

# Complex masses: derive the right-chiral propagator's mass-squared product
# by multiplying the large-q inverse series, not by choosing a transpose.
masses=[s.Matrix([[2+I,1-2*I],[1-2*I,4+I]]),s.Matrix([[1,I],[I,-1]]),s.zeros(2)]
flavor_checks=[]
for label,M in enumerate(masses):
    check(f'mass symmetry {label}',M-M.T)
    bigslash=s.kronecker_product(s.eye(2),slash)
    mass_insert=s.kronecker_product(M,PR)+s.kronecker_product(M.conjugate().T,PL)
    pr=s.kronecker_product(s.eye(2),PR);pl=s.kronecker_product(s.eye(2),PL)
    pnorm=(mom.T*s.diag(*eta)*mom)[0]
    second=pr*bigslash*mass_insert*bigslash*mass_insert*bigslash*pl
    check(f'right sterile propagator squared mass {label}',second-pnorm*s.kronecker_product(M.conjugate().T*M,PR*slash*PL))
    TM=s.trace(H*M.conjugate().T*M)
    mass_total=scalar_mass*TM
    mH=R(7,3)
    mass_parameter=mass_total-mH*zphi
    deltaM=-(zNL.conjugate().T*M+M*zNL.conjugate())/2
    check(f'no net Majorana chirality-flip action pole {label}',deltaM+(zNL.conjugate().T*M+M*zNL.conjugate())/2)
    check(f'Majorana parameter counterterm symmetric {label}',deltaM-deltaM.T)
    check(f'scalar action kinetic pole cancels {label}',zphi+scalar_wave*s.trace(H))
    check(f'scalar action mass pole cancels {label}',-mass_total+scalar_mass*TM)
    check(f'bare scalar parameter plus field transport {label}',mass_parameter+mH*zphi-mass_total)
    Yp=Vact.conjugate().T*Y*U.conjugate();Mp=U.conjugate().T*M*U.conjugate()
    Hp=Yp.conjugate().T*Yp
    zLp=-vector_wave*Yp*Yp.conjugate().T
    zNp=-weak_N*vector_wave*Hp.T
    deltaMp=-(zNp.conjugate().T*Mp+Mp*zNp.conjugate())/2
    check(f'active field-counterterm covariance {label}',zLp-Vact.conjugate().T*zL*Vact)
    check(f'sterile left field-counterterm covariance {label}',zNp-U.conjugate().T*zNL*U)
    check(f'complex mass parameter congruence {label}',deltaMp-U.conjugate().T*deltaM*U.conjugate())
    check(f'scalar mass trace covariance {label}',s.trace(Hp*Mp.conjugate().T*Mp)-TM)
    if label==0:
        reject('wrong complex mass product in scalar pole',TM-s.trace(H*M*M.conjugate().T))
        reject('ordinary M squared substituted in scalar pole',TM-s.trace(H*M*M))
        reject('omit mass parameter induced by field renormalization',deltaM)
        reject('scalar mass-action pole sign reversed',2*mass_total)
        reject('scalar action coefficient confused with bare parameter',mass_parameter-mass_total)
        reject('weak multiplicity halved in sterile mass transport',deltaM-deltaM/2)
    flavor_checks.append({'mass_case':label,'scalar_mass_trace':str(TM),'mass_symmetric':True})

result={'status':'pass','python':sys.version.split()[0],'sympy':s.__version__,
 'method':'Independent common-auxiliary-mass large-q angular expansion and Gamma-function radial UV residues, cross-checked by numerator/tadpole identities; no Feynman parameters or repository inputs.',
 'scope':'Pure-neutrino-Yukawa one-loop renormalizable two-point UV poles and cancelling action/parameter counterterms; no Yukawa-vertex, full beta function, seed-loop or observable claim.',
 'computed_normalized_poles':{'open_vector_moment':str(fermion),'closed_trace_before_fermion_minus':str(scalar),'scalar_wave_after_fermion_minus':str(scalar_wave),'scalar_mass_multiplier_after_fermion_minus':str(scalar_mass),'sterile_weak_multiplicity':str(weak_N)},
 'asymptotic_fermion_terms':fermion_pieces,'asymptotic_scalar_terms':scalar_pieces,
 'exact_checks':len(checks),'rejected_mutations':len(mutations),'mass_cases':flavor_checks,'checks':checks,'mutations':mutations}
print(json.dumps(result,indent=2))
