#!/usr/bin/env python3
"""Independent primitive seed/Yukawa four-lepton UV action.

Auxiliary-mass asymptotics and radial Gamma residues; scalar determinant
and literal Grassmann/Wick normalization; exact Clifford/Fierz/flavor checks.
No producer or repository imports and no optimization-removable assertions.
"""
from collections import Counter
from pathlib import Path
import json
import sympy as s

I=s.I; R=s.Rational
checks=Counter(); adverse=Counter()


class G:
    """Finite exterior algebra; field adjoints are specified by mate."""
    mate={}
    def __init__(self, data=None):
        self.d={k:s.expand(v) for k,v in (data or {}).items() if s.expand(v)!=0}
    @staticmethod
    def cast(x):
        return x if isinstance(x,G) else G({():s.sympify(x)})
    def __add__(self,x):
        x=G.cast(x);out=dict(self.d)
        for k,v in x.d.items():out[k]=out.get(k,0)+v
        return G(out)
    __radd__=__add__
    def __neg__(self):return G({k:-v for k,v in self.d.items()})
    def __sub__(self,x):return self+-G.cast(x)
    def __rsub__(self,x):return G.cast(x)+-self
    def __mul__(self,x):
        x=G.cast(x);out={}
        for a,v in self.d.items():
            for b,w in x.d.items():
                if set(a)&set(b):continue
                sign=(-1)**sum(i>j for i in a for j in b)
                key=tuple(sorted(a+b))
                out[key]=out.get(key,0)+sign*v*w
        return G(out)
    __rmul__=__mul__
    def __truediv__(self,x):return self*(1/s.sympify(x))
    def dagger(self):
        out=G()
        for k,v in self.d.items():
            term=G.cast(s.conjugate(v))
            for j in reversed(k):term=term*G({(self.mate[j],):1})
            out+=term
        return out


def iszero(x):
    if isinstance(x,G):return not x.d
    if isinstance(x,s.MatrixBase):return all(s.simplify(v)==0 for v in x)
    return s.simplify(x)==0


def check(name,x):
    if not iszero(x):raise RuntimeError(name+': exact comparison failed')
    checks[name]+=1


def reject(name,x):
    if iszero(x):raise RuntimeError(name+': adverse control escaped')
    adverse[name]+=1


eps=s.Matrix([[0,1],[-1,0]])
pauli=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-I],[I,0]]),s.diag(1,-1)]
z=s.zeros(2)
gamma=[s.BlockMatrix([[z,s.eye(2)],[s.eye(2),z]]).as_explicit()]
gamma += [s.BlockMatrix([[z,p],[-p,z]]).as_explicit() for p in pauli]
eta=[1,-1,-1,-1]
g5=I*gamma[0]*gamma[1]*gamma[2]*gamma[3]
PL=(s.eye(4)-g5)/2; PR=(s.eye(4)+g5)/2
Cspin=I*gamma[2]*gamma[0]


def dot(a,b):return sum(eta[mu]*a[mu]*b[mu] for mu in range(4))
def slash(a):return sum((eta[mu]*a[mu]*gamma[mu] for mu in range(4)),s.zeros(4))


def radial_uv(power,denominator,aux):
    """Coefficient of i/(16 pi^2 epsilon) in (q²)^power/(q²-aux)^den.

    Direct Gamma radial formula, with d=4-2epsilon.  The residue of
    Gamma(epsilon-m) is (-1)^m/m!, derived by its recurrence.
    """
    pole_shift=denominator-power-2
    if pole_shift>0:return s.Integer(0)
    residue=s.Integer(-1)**(-pole_shift)/s.factorial(-pole_shift)
    return s.simplify(s.Integer(-1)**(power-denominator)*aux**(2+power-denominator)
                      *s.gamma(power+2)*residue/(s.gamma(2)*s.gamma(denominator)))


def fields(n):
    L=[];bar=[]
    for f in range(n):
        lf=[];bf=[]
        for a in range(2):
            lv=[];bv=[]
            for r in range(2):
                j=((f*2+a)*2+r)*2
                G.mate[j]=j+1;G.mate[j+1]=j
                lv.append(G({(j,):1}));bv.append(G({(j+1,):1}))
            lf.append(lv);bf.append(bv)
        L.append(lf);bar.append(bf)
    return L,bar


def bil(a,b):
    return sum((eps[r,t]*a[r]*b[t] for r in range(2) for t in range(2)),G())


def scalar_pair(C,bar):
    n=C.rows
    return [[sum((C[f,g]*bil(bar[f][a],bar[g][b])
                  for f in range(n) for g in range(n)),G())
             for b in range(2)] for a in range(2)]


def P(C,S,L,bar):
    n=C.rows
    return sum((C[f,g]*S[h,j]*bil(bar[f][a],bar[g][b])*bil(L[h][a],L[j][b])
                for f in range(n) for g in range(n)
                for h in range(n) for j in range(n)
                for a in range(2) for b in range(2)),G())


def current(b,l,mu):
    # bar L gamma^mu L in the declared chiral basis.
    bs=[G(),G(),*b];ls=[*l,G(),G()]
    return sum((bs[r]*gamma[mu][r,t]*ls[t] for r in range(4) for t in range(4)),G())


def charge_left(lbar):
    return [G(),G(),*[-sum((eps[r,t]*lbar[t] for t in range(2)),G()) for r in range(2)]]


def four_scalar(b1,b2,l1,l2):
    # (bar L1 L2^c)(bar L1^c L2), both charge conjugations actual.
    bar1=[G(),G(),*b1]
    lc2=charge_left(b2)
    barlc1=[sum((l1[r]*eps[r,t] for r in range(2)),G()) for t in range(2)]+[G(),G()]
    lfour2=[*l2,G(),G()]
    return sum((bar1[r]*lc2[r] for r in range(4)),G())*sum((barlc1[r]*lfour2[r] for r in range(4)),G())


def main():
    for mu in range(4):
        for nu in range(4):
            check('Clifford coefficients',gamma[mu]*gamma[nu]+gamma[nu]*gamma[mu]
                  -2*eta[mu]*(mu==nu)*s.eye(4))
    check('charge-conjugation chiral blocks',Cspin-s.diag(eps,-eps))
    check('charge-conjugation antisymmetry',Cspin+Cspin.T)
    for mu in range(4):check('charge-conjugated gamma',Cspin*gamma[mu].T*Cspin.inv()+gamma[mu])

    # Derive the anomalous <NL NL> numerator from the large-q inverse of
    # slash q - M PR - Mdag PL, followed by the actual Majorana C^T.
    M=s.Matrix([[2+I,1-2*I],[1-2*I,4+I]])
    Y=s.Matrix([[1+I,2-I],[3-2*I,-1+I]])
    C=s.Matrix([[2-I,1+2*I],[1+2*I,3+I]])
    check('seed symmetric',C-C.T);check('Majorana mass symmetric',M-M.T)
    mom=s.Matrix([3,1,-1,0]);qnorm=dot(mom,mom)
    qslash=s.kronecker_product(s.eye(2),slash(mom))
    mass=s.kronecker_product(M,PR)+s.kronecker_product(M.H,PL)
    lc=s.kronecker_product(s.eye(2),PL)
    charge=s.kronecker_product(s.eye(2),Cspin.T)
    leading=lc*qslash*charge*lc.T/qnorm
    first_mass=I*lc*qslash*mass*qslash*charge*lc.T/qnorm**2
    eleft=s.diag(eps,s.zeros(2))
    check('same-chiral kinetic numerator vanishes',leading)
    check('complex anomalous numerator and minus epsilon',first_mass+I*s.kronecker_product(M,eleft)/qnorm)
    reject('replace anomalous M by Mdagger',first_mass+I*s.kronecker_product(M.H,eleft)/qnorm)
    # q=0 calibration against the two-variable holomorphic mass Gaussian.
    m=s.Integer(3);A=-m*eps
    check('zero-momentum Grassmann Gaussian sign',I*A.inv()-I*eps/m)
    check('zero-momentum four-component propagator agrees',PL*(-I*s.eye(4)/m)*Cspin.T*PL.T-I*eleft/m)

    # The B ordering is n_s eps_sr l_r. The internal n's straddle one
    # external l, so Wick's selected-pair sign is -1.
    order=['n1','l1','n2','l2'];left=order.index('n1');right=order.index('n2')
    wick=(-1)**(right-left-1)
    check('selected anomalous Wick crossing',wick+1)
    for r in range(2):
        for u in range(2):
            contraction=sum(eps[st,r]*eps[tt,u]*eps[st,tt] for st in range(2) for tt in range(2))
            check('actual anomalous spin contraction',contraction-eps[r,u])

    # i*seed=-i C/2 in these Weyl coordinates; two B vertices -i,
    # two scalar propagators i; anomalous propagator -i M eps.
    # The S B^2/2! term is obtained from 3 choices of S in order three.
    expansion=R(3,s.factorial(3))
    phase=(-I/2)*(-I)**2*I**2*(-I)*wick*expansion
    check('literal third-order expansion factor',expansion-R(1,2))
    check('literal ordered contraction phase',phase-R(1,4))
    reject('append closed-fermion-loop minus',phase+R(1,4))
    reject('extra Majorana propagator half',phase/2-R(1,4))

    # Independent UV integral: only the homogeneous degree-two numerator
    # and leading three denominators can be logarithmically divergent.
    aux=s.symbols('auxiliary_mass_squared',positive=True)
    residue=radial_uv(1,3,aux)
    check('actual Gamma radial UV residue',residue-1)
    check('auxiliary IR mass cancels in primitive pole',s.diff(residue,aux))
    check('external-only derivative numerator has no UV pole',radial_uv(0,3,aux))
    reject('delete derivative on internal scalar composite',residue-radial_uv(0,3,aux))
    for a in [R(1,3),R(7,2),11]:check('IR regulator witnesses',residue.subs(aux,a)-1)
    t=s.symbols('large_q_parameter')
    for trial in range(12):
        q=s.Matrix([3+trial,1,-1,0]);p=s.Matrix([1,trial-2,2,1]);k=s.Matrix([-2,1,trial,1])
        rr=s.Matrix([1,-1,trial,2]);q2=dot(q,q)
        if q2==0:raise RuntimeError('null UV witness direction')
        seed=s.Poly(s.expand(dot(t*q+rr,t*q+rr)),t)
        ds=[s.Poly(s.expand(dot(t*q+v,t*q+v)-mm),t)
            for v,mm in [(s.zeros(4,1),2),(p,2),(k,5)]]
        homogeneous=s.cancel(seed.nth(2)/s.prod(x.nth(2) for x in ds))
        check('whole-composite homogeneous loop extraction',homogeneous-1/q2**2)
        remainder=s.Poly(seed.as_expr()-seed.nth(2)*t**2,t)
        if remainder.degree()>1:raise RuntimeError('subleading numerator has divergent degree')
        checks['subleading numerator UV degree']+=1
        reject('external-only seed false UV coefficient',homogeneous)

    L,bar=fields(2)
    S=Y.conjugate()*M*Y.H
    check('mass-Yukawa pair symmetry',S-S.T)
    # Compute every flavor coefficient directly from the two B endpoints.
    endpoint=s.Matrix(2,2,lambda h,j:sum(s.conjugate(Y[h,i])*M[i,k]*s.conjugate(Y[j,k])
                                      for i in range(2) for k in range(2)))
    check('literal sterile/flavor contraction',endpoint-S)
    reject('mass conjugation omitted from left-Weyl dictionary',S-Y.conjugate()*M.conjugate()*Y.H)
    reject('Yukawa endpoint conjugation reversed',S-Y*M*Y.T)
    p=P(C,S,L,bar)
    reject('nonzero complex-flavor primitive',p)
    # Direct and crossed scalar contractions are retained before symmetry
    # is used. This is independent of the producer's seed vertex rules.
    crossed=sum((C[f,g]*S[h,j]*bil(bar[f][a],bar[g][b])*bil(L[h][b],L[j][a])
                 for f in range(2) for g in range(2) for h in range(2) for j in range(2)
                 for a in range(2) for b in range(2)),G())
    check('two scalar Wick assignments after flavor symmetry',crossed-p)
    raw=phase*(p+crossed)
    check('literal Wick action normalization',raw-p/2)
    reject('omit crossed scalar contraction',phase*p-p/2)
    reject('seed half cancelled twice',2*raw-p/2)

    # Scalar determinant cross: i/2 Tr log D gives -i/4 Tr X^2;
    # off-diagonal Nambu blocks supply exactly two traces of A B.
    # The independently integrated sterile source has B=-S(ll)/D.
    abar=scalar_pair(C,bar);bll=scalar_pair(S,L)
    for a in range(2):
        for b in range(2):
            check('seed scalar Hessian symmetry',abar[a][b]-abar[b][a])
            check('sterile-source scalar Hessian symmetry',bll[a][b]-bll[b][a])
    traceab=sum((abar[a][b]*bll[b][a] for a in range(2) for b in range(2)),G())
    check('scalar determinant weak-index trace',traceab-p)
    det_loop=(-I/4)*2*(I*residue)*traceab
    check('independent determinant versus literal Wick loop',det_loop-raw)
    reject('wrong determinant trace multiplicity',(-I/4)*(I*residue)*traceab-raw)
    loop=raw+raw.dagger();ct=-raw-raw.dagger()
    check('local action Hermiticity',loop-loop.dagger())
    check('counterterm cancels loop action',loop+ct)
    reject('counterterm has loop sign',loop-ct)
    check('Hermitian-conjugate pair exchange',p.dagger()-P(S.conjugate(),C.conjugate(),L,bar))

    jj=G();scalar4=G()
    for f in range(2):
        for g in range(2):
            for h in range(2):
                for j in range(2):
                    for a in range(2):
                        for b in range(2):
                            pair=bil(bar[f][a],bar[g][b])*bil(L[h][a],L[j][b])
                            vector=sum((eta[mu]*current(bar[f][a],L[h][a],mu)
                                        *current(bar[g][b],L[j][b],mu) for mu in range(4)),G())
                            scalar=four_scalar(bar[f][a],bar[g][b],L[h][a],L[j][b])
                            check('full Grassmann vector Fierz coefficient',vector+2*pair)
                            check('actual four-component scalar conversion',scalar+pair)
                            jj+=C[f,g]*S[h,j]*vector
                            scalar4+=C[f,g]*S[h,j]*scalar
    check('summed vector-current action',raw+jj/4)
    check('summed scalar four-component action',raw+scalar4/2)
    reject('Fierz without fermion reorder sign',jj-2*p)

    V=s.Matrix([[R(3,5),4*I/5],[4*I/5,R(3,5)]])*s.diag(I,(5+12*I)/13)
    U=s.diag((3+4*I)/5,I)*s.Matrix([[R(5,13),R(12,13)],[-R(12,13),R(5,13)]])
    check('active rotation unitary',V.H*V-s.eye(2));check('sterile rotation unitary',U.H*U-s.eye(2))
    Cp=V.H*C*V.conjugate()
    for label,mm in [('general',M),('singular',s.Matrix([[1,I],[I,-1]])),('zero',s.zeros(2))]:
        ss=Y.conjugate()*mm*Y.H
        yp=V.H*Y*U.conjugate();mp=U.H*mm*U.conjugate()
        sp=yp.conjugate()*mp*yp.H
        check('complex mass-source covariance '+label,sp-V.T*ss*V)
        lo=[[[sum((V[f,g]*L[g][a][r] for g in range(2)),G()) for r in range(2)]
             for a in range(2)] for f in range(2)]
        bo=[[[sum((s.conjugate(V[f,g])*bar[g][a][r] for g in range(2)),G()) for r in range(2)]
             for a in range(2)] for f in range(2)]
        check('actual exterior active-coordinate covariance '+label,P(C,ss,lo,bo)-P(Cp,sp,L,bar))
    check('massless primitive vanishes',P(C,s.zeros(2),L,bar))
    weak=s.Matrix([[R(3,5),4*I/5],[4*I/5,R(3,5)]])
    lw=[[[sum((weak[a,b]*L[f][b][r] for b in range(2)),G()) for r in range(2)]
          for a in range(2)] for f in range(2)]
    bw=[[[sum((s.conjugate(weak[a,b])*bar[f][b][r] for b in range(2)),G()) for r in range(2)]
          for a in range(2)] for f in range(2)]
    check('weak transformation unitary',weak.H*weak-s.eye(2))
    check('full weak-index action invariance',P(C,S,lw,bw)-p)
    oneL,onebar=fields(1)
    cone=s.Matrix([[1+I]]);yone=s.Matrix([[1+2*I,2-I]])
    sone=yone.conjugate()*M*yone.H
    pone=P(cone,sone,oneL,onebar)
    reject('one active flavor survives',pone+pone.dagger())
    check('one-flavor Hermiticity',(pone+pone.dagger())-(pone+pone.dagger()).dagger())

    result={'status':'PASS','input_commit':'f1f9a892e15545c70e4fa05827ea7fefc5a438eb',
            'dateline':'23 October 2026 Pacific',
            'method':'Independent auxiliary-mass homogeneous UV extraction and radial Gamma residue; scalar Nambu determinant cross-checked against literal Wick/exterior contractions.',
            'positive_checks':sum(checks.values()),'positive_categories':dict(checks),
            'negative_checks':sum(adverse.values()),'negative_categories':dict(adverse),
            'radial_UV_residue_over_i_kappa':str(residue),
            'anomalous_left_sterile_propagator_UV':'-i M epsilon/q^2',
            'source_matrix':'R=Y* M Ydagger, symmetric',
            'ordered_Weyl_action_loop_over_kappa':'+C_ab R_cd (barl_a^u epsilon barl_b^v)(l_c^u epsilon l_d^v)/2+h.c.',
            'ordered_Weyl_action_counterterm_over_kappa':'-C_ab R_cd (barl_a^u epsilon barl_b^v)(l_c^u epsilon l_d^v)/2+h.c.',
            'vector_current_counterterm_over_kappa':'+C_ab R_cd J_ac^mu J_bd_mu/4+h.c., with weak indices summed separately in each current',
            'repository_imports':False,'assertions_required':False,
            'scope':'Primitive one-seed two-conjugate-neutrino-Yukawa four-lepton 1PI UV action in its declared chiral representative; vector Fierz is a four-dimensional projection, not evanescent-basis closure. No total physical Wilson coefficients, complete beta function, finite matching or observable.'}
    Path(__file__).with_name('INDEPENDENT_FOUR_LEPTON_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
