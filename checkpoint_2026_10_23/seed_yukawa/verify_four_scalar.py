#!/usr/bin/env python3
"""Literal seed + two conjugate Yukawas: four-scalar 1PI UV pole.

Feynman-parameter/Gamma integration, explicit Weyl/Wick factors, complete
labelled local-polynomial matching, and complex-M flavor covariance.
Pinned SymPy; no project imports or removable scientific assertions.
"""
from collections import Counter
from itertools import product
from pathlib import Path
import json
import sympy as s

I=s.I; R=s.Rational
positive=Counter();negative=Counter()

def vanishes(x):
    if isinstance(x,s.MatrixBase):return all(vanishes(a) for a in x)
    x=s.expand(x)
    return x==0 or s.cancel(x)==0

def check(label,x):
    if not vanishes(x):raise RuntimeError(label+': exact identity failed: '+str(x))
    positive[label]+=1

def reject(label,x):
    if vanishes(x):raise RuntimeError(label+': mutation escaped')
    negative[label]+=1


def simplex(poly,b,c):
    """Integral over b,c>=0,b+c<=1, evaluated by exact monomial moments."""
    out=0
    for (j,k),coef in s.Poly(s.expand(poly),b,c).terms():
        out+=coef*s.factorial(j)*s.factorial(k)/s.factorial(j+k+2)
    return s.expand(out)


def main():
    # For I_n(Delta), factoring i/(16 pi^2 epsilon), Gamma recursion gives
    # I1 -> Delta, I2 -> 1, I3 -> 0. Polynomial division then fixes both
    # divergent tensor radial masters in the shifted triangle.
    def gamma_residue(integer_shift):
        if integer_shift>0:return s.Integer(0)
        return s.Integer(-1)**(-integer_shift)/s.factorial(-integer_shift)
    Delta=s.symbols('Delta')
    masters={n:s.Integer(-1)**n*gamma_residue(n-2)*Delta**(2-n)
             for n in (1,2,3,4)}
    check('Gamma tadpole residue',masters[1]-Delta)
    check('Gamma logarithmic bubble residue',masters[2]-1)
    check('Gamma UV-finite triangle residue',masters[3])
    l2_over_D3=masters[2]+Delta*masters[3]
    l4_over_D3=masters[1]+2*Delta*masters[2]+Delta**2*masters[3]
    check('shifted quadratic radial master',l2_over_D3-1)
    check('shifted quartic radial master',l4_over_D3-3*Delta)
    reject('omit radial tadpole in quartic master',l4_over_D3-2*Delta)

    b,c=s.symbols('b c');mass2=s.symbols('mass_squared')
    invariants=s.symbols('r2 s2 t2 rs rt st')
    r2,s2,t2,rs,rt,st=invariants
    G=s.Matrix([[r2,rs,rt],[rs,s2,st],[rt,st,t2]])
    dot=lambda x,y:s.expand((x.T*G*y)[0])
    r=s.Matrix([1,0,0]);v_s=s.Matrix([0,1,0]);t=s.Matrix([0,0,1])
    v=b*r+c*v_s
    rr=r-v;ss=v_s-v;tt=t-v
    delta=(1-b-c)*mass2-b*r2-c*s2+dot(v,v)
    # Full numerator [(l+R).(l+S)](l+T)^2. The cross term is
    # 2(l.(R+S))(l.T); its angular coefficient is 2/d -> 1/2.
    angular=dot(tt,tt)+dot(rr,ss)+(dot(rr,tt)+dot(ss,tt))/2
    integrated=simplex(2*(3*delta+angular),b,c)
    expected=mass2-(r2+s2-2*rs)/2+t2-(rt+st)/2
    check('complete abstract Feynman-parameter UV polynomial',integrated-expected)
    for inv in (mass2,*invariants):
        check('exhaustive invariant coefficient',s.diff(integrated,inv)-s.diff(expected,inv))
    check('simplex area',simplex(1,b,c)-R(1,2))
    check('mass-weighted simplex moment',simplex(6*(1-b-c),b,c)-1)
    bad=simplex(2*(3*delta+dot(tt,tt)+dot(rr,ss)),b,c)
    reject('drop numerator mixed angular term',integrated-bad)
    bad=simplex(2*(2*delta+angular),b,c)
    reject('wrong radial factor two instead of three',integrated-bad)
    reject('omit sterile retained-mass pole',integrated-integrated.subs(mass2,0))

    # Weyl conventions directly inherited from Cspin=diag(eps,-eps).
    eps=s.Matrix([[0,1],[-1,0]])
    pauli=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-I],[I,0]]),s.diag(1,-1)]
    sig=[s.eye(2)]+[-a for a in pauli]
    eta=(1,-1,-1,-1)
    def spin(A,B,swap=False):
        return s.expand(sum(eps[a,b_]*eps[c_,d]*eps[e,f]*eps[d,f]
                           *A[c_,b_ if swap else a]*B[e,a if swap else b_]
                           for a,b_,c_,d,e,f in product(range(2),repeat=6)))
    for mu,nu in product(range(4),repeat=2):
        target=2*eta[mu] if mu==nu else 0
        check('all closed Weyl numerator coefficients',spin(sig[mu],sig[nu])-target)
        check('exchanged active contraction sign',spin(sig[mu],sig[nu],True)+target)
    def wick_sign(pairs):
        return (-1)**sum(a<c_<b_<d or c_<a<d<b_
                        for j,(a,b_) in enumerate(pairs) for c_,d in pairs[j+1:])
    direct=[(0,2),(1,4),(3,5)];crossed=[(0,4),(1,2),(3,5)]
    check('direct six-fermion Wick sign',wick_sign(direct)-1)
    check('crossed six-fermion Wick sign',wick_sign(crossed)+1)
    # iS seed -i/2; two iB=-i and expansion1/2!, then three ordered
    # contractions -i. Active1 carries -q-k1, adding one numerator minus.
    expansion=(-I/2)*(-I)**2/2
    contractions=(-I)**3
    coefficient=s.expand(expansion*contractions*(-1)
                         *(wick_sign(direct)*2+wick_sign(crossed)*(-2)))
    check('literal closed-loop action-integrand coefficient',coefficient-1)
    reject('append an extra Majorana half',coefficient-coefficient/2)
    reject('omit one active Wick pairing',coefficient-s.expand(expansion*contractions*(-1)*2))
    reject('drop closed-chain active momentum reversal',coefficient+coefficient)

    # Local action matching after momentum conservation. Six abstract Gram
    # entries span the unrestricted degree-two polynomial of three momenta.
    gram=s.symbols('p1sq p2sq k1sq p1p2 p1k1 p2k1')
    GG=s.Matrix([[gram[0],gram[3],gram[4]],
                 [gram[3],gram[1],gram[5]],
                 [gram[4],gram[5],gram[2]]])
    dp=lambda x,y:s.expand((x.T*GG*y)[0])
    p1=s.Matrix([1,0,0]);p2=s.Matrix([0,1,0]);k1=s.Matrix([0,0,1]);k2=-p1-p2-k1
    K=k1+k2;w=k1-k2;hD=p1+k1;hX=p1+k2
    JD=s.expand(-dp(K,K)/2+dp(hD,hD)-dp(w,hD)/2)
    JX=s.expand(-dp(K,K)/2+dp(hX,hX)+dp(w,hX)/2)
    # Vertex from A=rho dphi^dagger.dphi and B=(phi^dagger dphi)^2.
    AD=-dp(k1,p1)-dp(k2,p2);AX=-dp(k1,p2)-dp(k2,p1)
    BD=BX=-2*dp(p1,p2)
    for name,left,right in [('direct',2*JD,BD-AD),('crossed',2*JX,BX-AX)]:
        check('local derivative action polynomial '+name,left-right)
        for g in gram:check('local-action exhaustive degree-two coefficient',s.diff(left-right,g))
    # Potential vertex rho² has2(deltaD+deltaX), not an unexplained weak sum.
    phi=s.symbols('phi0:2');barphi=s.symbols('barphi0:2')
    rho=sum(a*b_ for a,b_ in zip(phi,barphi))
    for a,b_,c_,d in product(range(2),repeat=4):
        weakD=int(a==c_)*int(b_==d);weakX=int(a==d)*int(b_==c_)
        deriv=s.diff(rho**2,phi[a],phi[b_],barphi[c_],barphi[d])
        check('literal four-scalar potential weak derivative',deriv-2*(weakD+weakX))
    reject('mistake fixed scalar weak indices for factor-two trace',JD-2*JD)
    badJ=expected.subs({t2:s2,rt:rs,st:s2})
    reject('P0 only on internal lepton instead of full composite',expected-badJ)

    # Strong exact Bose checks at12 generic nonexceptional points and all
    # sixteen weak assignments, retaining complex S-only coefficients.
    etaM=s.diag(*eta);fourdot=lambda a,b_:s.expand((a.T*etaM*b_)[0])
    def kernel(ps,ks,weak,c1,c3):
        aa,bb,cc,dd=weak;DD=int(aa==cc)*int(bb==dd);XX=int(aa==dd)*int(bb==cc)
        kk=ks[0]+ks[1];ww=ks[0]-ks[1];hd=ps[0]+ks[0];hx=ps[0]+ks[1]
        jd=-fourdot(kk,kk)/2+fourdot(hd,hd)-fourdot(ww,hd)/2
        jx=-fourdot(kk,kk)/2+fourdot(hx,hx)+fourdot(ww,hx)/2
        return s.expand(2*c3*(DD+XX)+2*c1*(DD*jd+XX*jx))
    def action_kernel(ps,ks,weak,c1,c3):
        aa,bb,cc,dd=weak;DD=int(aa==cc)*int(bb==dd);XX=int(aa==dd)*int(bb==cc)
        aaD=-fourdot(ks[0],ps[0])-fourdot(ks[1],ps[1])
        aaX=-fourdot(ks[0],ps[1])-fourdot(ks[1],ps[0])
        bbBoth=-2*fourdot(ps[0],ps[1])
        return s.expand(2*c3*(DD+XX)+c1*(DD*(bbBoth-aaD)+XX*(bbBoth-aaX)))
    for trial in range(12):
        ps=[s.Matrix([2+trial,1,-2,3]),s.Matrix([-3,2+trial,1,-1])]
        ks=[s.Matrix([1,-2,3+trial,2])];ks.append(-ps[0]-ps[1]-ks[0])
        for weak in product(range(2),repeat=4):
            aa,bb,cc,dd=weak;baseline=kernel(ps,ks,weak,2+I,3-2*I)
            check('phi Bose exchange',baseline-kernel(ps[::-1],ks,(bb,aa,cc,dd),2+I,3-2*I))
            check('phi-dagger Bose exchange',baseline-kernel(ps,ks[::-1],(aa,bb,dd,cc),2+I,3-2*I))
            herm=baseline+s.conjugate(kernel(ks,ps,(cc,dd,aa,bb),2+I,3-2*I))
            local=action_kernel(ps,ks,weak,2+I,3-2*I)+s.conjugate(action_kernel(ks,ps,(cc,dd,aa,bb),2+I,3-2*I))
            check('four-field Fourier action including conjugate',herm-local)

    # Complex flavor covariance in fixed, non-diagonal mass coordinates.
    C=s.Matrix([[2-I,1+2*I],[1+2*I,3+I]])
    Y=s.Matrix([[1+I,2-I],[3-2*I,-1+I]])
    M=s.Matrix([[2+I,1-2*I],[1-2*I,4+3*I]])
    U=s.diag(I,(3+4*I)/5)*s.Matrix([[R(3,5),R(4,5)],[-R(4,5),R(3,5)]])
    V=s.Matrix([[R(5,13),R(12,13)],[-R(12,13),R(5,13)]])*s.diag((5+12*I)/13,I)
    check('sterile rotation unitary',U.H*U-s.eye(2));check('active rotation unitary',V.H*V-s.eye(2))
    def sources(c7,y,m):
        return s.expand(s.trace(c7*y.conjugate()*m*y.H)),s.expand(s.trace(c7*y.conjugate()*m*m.H*m*y.H))
    b1,b3=sources(C,Y,M)
    for mm in [M,s.Matrix([[1,I],[I,-1]]),s.zeros(2)]:
        check('retained Majorana mass symmetry',mm-mm.T)
        check('anomalous cubic mass numerator symmetry',mm*mm.H*mm-(mm*mm.H*mm).T)
        check('left/right anomalous mass order identity',(mm*mm.H)*mm-mm*(mm.H*mm))
        cc=V.H*C*V.conjugate();yy=V.H*Y*U.conjugate();mp=U.H*mm*U.conjugate()
        old=sources(C,Y,mm);new=sources(cc,yy,mp)
        check('one-mass scalar source flavor covariance',new[0]-old[0])
        check('three-mass scalar source flavor covariance',new[1]-old[1])
    check('zero retained mass removes displayed sector',sum(sources(C,Y,s.zeros(2))))
    reject('wrong conjugated Majorana mass in source',b1-s.trace(C*Y.conjugate()*M.conjugate()*Y.H))
    reject('wrong noncommuting cubic mass ordering',b3-s.trace(C*Y.conjugate()*M.H*M*M*Y.H))
    reject('replace cubic mass source by linear one',b3-b1)
    reject('omit the conjugate of a complex scalar action coefficient',b1-s.conjugate(b1))

    # Constant-background real one-active/one-sterile Gaussian determinant
    # checks the absolute potential sign/factor, not the full momentum pole.
    q2,a,dd,m,cc=s.symbols('q_squared a d m C7',nonzero=True)
    mm=s.Matrix([[a,dd],[dd,m]])
    det=(q2*s.eye(2)-mm*mm).det()
    linear=s.diff(det,a).subs(a,0)
    check('constant-background determinant linear mass insertion',linear+2*dd**2*m)
    baseline=det.subs({a:0,dd:0})
    loglinear=s.cancel(linear/baseline*(-cc*q2))
    check('derivative seed mass insertion in Gaussian determinant',loglinear-2*cc*dd**2*m/(q2-m*m))
    tadpole_coefficient=s.limit((q2-m*m)*loglinear,q2,m*m)
    radial_potential=tadpole_coefficient*masters[1].subs(Delta,m*m)
    check('constant-background radial pole and displayed plus conjugate',radial_potential-2*cc*dd**2*m**3)

    report={'status':'PASS','input_commit':'f1f9a892e15545c70e4fa05827ea7fefc5a438eb',
            'dateline':'23 October 2026 Pacific',
            'method':'Actual Feynman-parameter triangle UV integration and Gamma radial masters; Weyl Wick normalization; exact off-shell local-polynomial matching.',
            'pole_core_after_factoring_i_kappa':'M_squared-(r-s)^2/2+t^2-(r+s).t/2',
            'displayed_effective_action_after_factoring_kappa':'B3 rho^2+B1[(phi_dagger partial phi)^2-rho partial_phi_dagger.partial_phi]+h.c.',
            'sources':{'B1':'Tr(C7 Y* M Ydagger)','B3':'Tr(C7 Y* M Mdagger M Ydagger)'},
            'positive_checks':sum(positive.values()),'positive_categories':dict(positive),
            'rejected_mutations':sum(negative.values()),'negative_categories':dict(negative),
            'scope':['one original symmetric half-normalized derivative seed plus exactly two conjugate neutrino Yukawa vertices',
                     'primitive one-loop four-scalar 1PI pole at zero gauge field; retained Majorana field not integrated out',
                     'unreduced local action and cancelling negative counterterm; no complete physical Wilson map, beta function, finite threshold or observable'],
            'project_imports':False,'assertions_required':False}
    Path(__file__).with_name('FOUR_SCALAR_RESULTS.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
