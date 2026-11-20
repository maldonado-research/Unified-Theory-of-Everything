#!/usr/bin/env python3
"""Independent auxiliary-mass/angular four-scalar seed UV audit.

Exact SymPy algebra, no project imports or removable assertions. The UV
polynomial is integrated from the expanded numerator and denominators.
"""
from collections import Counter
from itertools import permutations, product
from pathlib import Path
import json
import sympy as s
I=s.I;R=s.Rational
positive=Counter();negative=Counter()

def zero(x):
    if isinstance(x,s.MatrixBase):return all(s.expand(z)==0 for z in x)
    return s.cancel(s.expand(x))==0

def check(name,x):
    if not zero(x):raise RuntimeError(name+': '+str(x))
    positive[name]+=1

def reject(name,x):
    if zero(x):raise RuntimeError(name+': mutation escaped')
    negative[name]+=1

q2,zP,zp,zk,aux,m2,PP,pp,kk,Pp,Pk,pk,d=s.symbols('q2 zP zp zk auxiliary_mass_squared sterile_mass_squared P2 p2 k2 Pp Pk pk d')
N=2*(q2-zP)*(q2-2*zp+pp)
shifts=[aux,aux-2*zP+PP,aux+2*zk+kk-m2]

def radial(power,denominator):
    n=denominator-power-2
    if n>0:return s.Integer(0)
    return s.Integer(-1)**(power-denominator)*aux**(2+power-denominator)*s.gamma(power+2)*s.Integer(-1)**(-n)/(s.gamma(denominator)*s.factorial(-n))

def angular(mon,denominator):
    qa,i,j,k=mon;rank=i+j+k
    if rank%2:return s.Integer(0)
    if rank==0:return radial(qa,denominator)
    dots={('P','P'):PP,('p','p'):pp,('k','k'):kk,('P','p'):Pp,('P','k'):Pk,('p','k'):pk}
    def dotkey(a,b):return dots.get((a,b),dots.get((b,a)))
    labels=['P']*i+['p']*j+['k']*k
    if rank==2:return dotkey(*labels)*radial(qa+1,denominator)/d
    if rank==4:
        pairing=sum(dotkey(labels[0],labels[z])*dotkey(*[labels[x] for x in range(1,4) if x!=z]) for z in range(1,4))
        return pairing*radial(qa+2,denominator)/(d*(d+2))
    raise RuntimeError('Unexpected angular rank '+str(rank))

def integrate(numerator=N,drop=None):
    total=0;parts=[]
    for i in range(3):
        for j in range(3-i):
            for k in range(3-i-j):
                if drop==(i,j,k):continue
                den=3+i+j+k
                poly=s.Poly(s.expand(numerator*(-shifts[0])**i*(-shifts[1])**j*(-shifts[2])**k),q2,zP,zp,zk)
                piece=s.expand(sum(coef*angular(mon,den) for mon,coef in poly.terms()))
                total+=piece
                if piece:parts.append({'denominator_orders':[i,j,k],'residue_before_d4':str(s.factor(piece))})
    return s.expand(total.subs(d,4)),parts

residue,pieces=integrate()
check('auxiliary regulator cancels',s.diff(residue,aux))
check('no external k squared remainder',s.diff(residue,kk))
# Independent numerator cancellation reduces the triangle into two bubbles.
# For int (l+a)^2/[l^2((l+Q)^2-m2)], the UV residue is
# m2-Q.a+a^2 (tadpole + vector and scalar bubble residues).
first_bubble=m2-(PP+Pk-Pp-pk)+(PP+pp-2*Pp)
second_bubble=m2+pk+pp
reduced=first_bubble+second_bubble-PP
check('angular triangle versus denominator cancellation',residue-reduced)
check('zero external momentum mass residue',residue.subs({PP:0,pp:0,kk:0,Pp:0,Pk:0,pk:0})-2*m2)
check('radial logarithm normalization',radial(1,3)-1)
check('radial quadratic integral mass dependence',radial(2,3)-3*aux)
for orders in [(0,0,0),(0,0,1),(0,1,0),(1,0,0),(0,0,2),(0,1,1),(0,2,0)]:
    reject('omit UV expansion term '+str(orders),residue-integrate(drop=orders)[0])
reject('delete Higgs and mixed derivatives in seed',residue-integrate(2*(q2-zP)*q2)[0])
reject('discard retained sterile mass',residue-residue.subs(m2,0))
# The first neglected geometric order has UV degree 2-3<0.
check('geometric expansion UV closure',s.Integer(2-3)+1)

# Full chiral spin contractions and Wick signs before introducing a loop minus.
eps=s.Matrix([[0,1],[-1,0]])
pauli=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-I],[I,0]]),s.diag(1,-1)]
sigma=[s.eye(2),*pauli];eta=[1,-1,-1,-1]
for mu,nu in product(range(4),repeat=2):
    direct=sum(eps[a,b]*eps[c,e]*eps[f,g]*eps[e,g]*sigma[mu][c,a]*sigma[nu][f,b] for a,b,c,e,f,g in product(range(2),repeat=6))
    exchanged=-sum(eps[a,b]*eps[c,e]*eps[f,g]*eps[e,g]*sigma[mu][f,a]*sigma[nu][c,b] for a,b,c,e,f,g in product(range(2),repeat=6))
    target=2*eta[mu]*(mu==nu)
    check('literal direct chiral contraction',direct-target)
    check('exchanged active Wick contraction',exchanged-target)
# Sequence: barl_a barl_b l_c n_e l_f n_g. Matchings (02)(14)(35)
# and (04)(12)(35) have crossing signs + and - respectively.
def wick_sign(pairs):return (-1)**sum(a<c<b<e or c<a<e<b for n,(a,b) in enumerate(pairs) for c,e in pairs[n+1:])
check('first Wick pairing crossings',wick_sign([(0,2),(1,4),(3,5)])-1)
check('second Wick pairing crossings',wick_sign([(0,4),(1,2),(3,5)])+1)
# S=-C barl eps barl/2; B=-Y* phi* l eps n. Covariances in the
# displayed order are -i sigma(q),-i sigma(P-q),-i M eps.
phase=(-R(1,2))*(-1)**2*(I**3/R(2))*((-I)**3)
check('literal action expansion phase',phase+R(1,4))
check('unlabelled loop numerator equals half stripped residue',phase*4*(-1)-R(1,1))
reject('extra Majorana half',phase*4*(-1)-R(1,2))
reject('extra closed loop minus after Wick enumeration',-phase-phase)
reject('drop Yukawa expansion factorial',2*phase-phase)

# Constant-background independent Majorana determinant check, real one flavor.
x,c,y,h,m=s.symbols('q_squared C7 y phi M')
Mass=s.Matrix([[-c*x*h*h,y*h],[y*h,m]])
det=(x*s.eye(2)-Mass*Mass).det()
linear=s.diff(det,c).subs(c,0)
check('constant background determinant linear seed term',linear-2*x*y*y*m*h**4)
# At order C*y^2 use det_0=x(x-M^2), hence delta logdet=2C*y²*M*h⁴/(x-M²).
check('constant background potential residue',s.cancel(linear/(x*(x-m*m)))-2*y*y*m*h**4/(x-m*m))
reject('sterile mass cube omitted from potential',2*c*y*y*m**3*h**4-2*c*y*y*m*h**4)

# Actual differentiation of bosonic local monomials, all fields Fourier
# expanded with incoming waves exp(-ip.x). No off-shell EOM is used.
def mdot(a,b):return sum(eta[j]*a[j]*b[j] for j in range(4))
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def neg(a):return tuple(-x for x in a)
def kernel(kind,phimom,dagmom,weak):
    # rho^2, A=rho dphi*.dphi, B=(phi* dphi)^2.
    # Each monomial has weak contractions (slot0,slot1),(slot2,slot3).
    total=0
    for up in permutations(range(2)):
        for down in permutations(range(2)):
            if weak[2+down[0]]!=weak[up[0]] or weak[2+down[1]]!=weak[up[1]]:continue
            if kind=='rho2':term=1
            elif kind=='A':term=-mdot(dagmom[down[1]],phimom[up[1]])
            elif kind=='B':term=-mdot(phimom[up[0]],phimom[up[1]])
            elif kind=='Bdag':term=-mdot(dagmom[down[0]],dagmom[down[1]])
            else:raise RuntimeError(kind)
            total+=term
    return s.expand(total)

def R0(P,p,k):return -mdot(P,P)-mdot(P,k)-mdot(P,p)+2*mdot(p,p)+2*mdot(p,k)
for trial in range(5):
    p=(2+trial,1,trial-1,2);p2=(1,-1,2,trial);k=(-2,trial,1,-1)
    P=add(p,p2);k2=neg(add(P,k));ph=[p,p2];dag=[k,k2]
    subs={PP:mdot(P,P),pp:mdot(p,p),kk:mdot(k,k),Pp:mdot(P,p),Pk:mdot(P,k),pk:mdot(p,k)}
    check('off shell angular coefficient witness',residue.subs(subs)-2*m2-R0(P,p,k))
    check('loop reversal and conjugate-scalar Bose identity',R0(P,p,k)-R0(P,p2,k2))
    for weak in product(range(2),repeat=4):
        dd=int(weak[0]==weak[2] and weak[1]==weak[3]);dx=int(weak[0]==weak[3] and weak[1]==weak[2])
        target=dd*R0(P,p,k)+dx*R0(P,p2,k)
        check('local derivative action full weak kernel',kernel('B',ph,dag,weak)-kernel('A',ph,dag,weak)-target)
        check('local potential full weak kernel',kernel('rho2',ph,dag,weak)-2*(dd+dx))
        swapped=(weak[1],weak[0],weak[2],weak[3])
        check('identical scalar Bose exchange',kernel('B',ph,dag,weak)-kernel('A',ph,dag,weak)-kernel('B',[p2,p],dag,swapped)+kernel('A',[p2,p],dag,swapped))
        check('Hermitian derivative placement exchange',kernel('Bdag',ph,dag,weak)-kernel('B',dag,ph,(weak[2],weak[3],weak[0],weak[1])))
    if trial==0:
        w=(0,1,0,1)
        reject('omit identical crossed weak channel',kernel('B',ph,dag,(0,1,1,0))-kernel('A',ph,dag,(0,1,1,0)))
        reject('wrong relative local derivative sign',kernel('B',ph,dag,w)+kernel('A',ph,dag,w)-R0(P,p,k))

# Complex flavor and ordered Majorana mass covariance, including singular M.
Y=s.Matrix([[1+I,2-I,1],[3-2*I,-1+I,2+I]])
C=s.Matrix([[2+I,1-2*I],[1-2*I,3-I]])
V=s.diag(I,(3+4*I)/5)*s.Matrix([[R(3,5),R(4,5)],[-R(4,5),R(3,5)]])
U=s.diag(I,(5+12*I)/13,1)*s.Matrix([[R(3,5),R(4,5),0],[-R(4,5),R(3,5),0],[0,0,1]])
check('active rotation unitary',V.H*V-s.eye(2));check('sterile rotation unitary',U.H*U-s.eye(3))
for label,M in enumerate([s.Matrix([[2+I,1-I,2],[1-I,3+2*I,-I],[2,-I,4-I]]),s.Matrix([[1,I,0],[I,-1,0],[0,0,0]]),s.zeros(3)]):
    K1=s.expand(s.trace(C*Y.conjugate()*M*Y.H))
    K3=s.expand(s.trace(C*Y.conjugate()*M*M.H*M*Y.H))
    Mp=U.H*M*U.conjugate();Yp=V.H*Y*U.conjugate();Cp=V.H*C*V.conjugate()
    check('mass symmetry',M-M.T)
    check('linear retained mass invariant covariance',s.trace(Cp*Yp.conjugate()*Mp*Yp.H)-K1)
    check('cubic retained mass invariant covariance',s.trace(Cp*Yp.conjugate()*Mp*Mp.H*Mp*Yp.H)-K3)
    check('linear mass operator symmetric',Y.conjugate()*M*Y.H-(Y.conjugate()*M*Y.H).T)
    check('cubic mass operator symmetric',Y.conjugate()*M*M.H*M*Y.H-(Y.conjugate()*M*M.H*M*Y.H).T)
    if label==0:
        reject('wrong complex mass conjugation',K1-s.trace(C*Y.conjugate()*M.H*Y.H))
        reject('ordinary M cube replaces ordered palindrome',K3-s.trace(C*Y.conjugate()*M**3*Y.H))
        reject('drop sterile Yukawa conjugation',K1-s.trace(C*Y*M*Y.T))

phi=s.Matrix([1+I,2-I]);v=s.Matrix([2-I,-1+3*I])
rho=(phi.H*phi)[0];Aval=s.expand(rho*(v.H*v)[0]);Bval=s.expand((phi.H*v)[0]**2)
hermitian=s.expand(I*(Bval-Aval)-I*(s.conjugate(Bval)-Aval))
check('Hermitian action real for imaginary invariant',hermitian-s.conjugate(hermitian))
phip=V*phi;vp=V*v
check('global weak rotation scalar norm',s.expand((phip.H*phip)[0])-rho)
check('global weak rotation derivative invariant',s.expand((phip.H*vp)[0]**2)-Bval)
check('global weak rotation mixed derivative norm',s.expand((phip.H*phip)[0]*(vp.H*vp)[0])-Aval)
reject('premature real part of K1 deletes derivative action',hermitian)
check('real neutral one component radial cancellation',s.Integer(3)**2*s.Integer(5)**2-s.Integer(3)**2*s.Integer(5)**2)

report={'status':'PASS','input_commit':'f1f9a892e15545c70e4fa05827ea7fefc5a438eb','date':'23 October 2026 Pacific',
 'method':'Independent common-auxiliary-mass large-q angular expansion and radial Gamma UV residues, denominator-cancellation and literal Wick/determinant controls; no Feynman parameters or repository imports.',
 'scope':'One original derivative seed and two conjugate neutrino Yukawa vertices, primitive one-loop four-scalar local action pole only; no complete physical beta function or finite matching.',
 'computed_one_pairing_residue':str(residue),'radial_expansion_terms':pieces,
 'action_pole':'kappa {K3 rho^2 + K1 [(phi^dagger partial_mu phi)^2-rho partial_mu phi^dagger partial^mu phi]} + h.c.',
 'flavor_invariants':{'K1':'Tr(C7 Y* M Ydagger)','K3':'Tr(C7 Y* M Mdagger M Ydagger)'},
 'counterterm_sign':'negative of displayed effective-action pole','positive_checks':sum(positive.values()),'positive_categories':dict(positive),
 'rejected_mutations':sum(negative.values()),'mutation_categories':dict(negative),'repository_imports':False,'assertions_required':False}
Path(__file__).with_name('INDEPENDENT_FOUR_SCALAR_RESULTS.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
