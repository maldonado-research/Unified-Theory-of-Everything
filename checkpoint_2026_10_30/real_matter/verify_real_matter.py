#!/usr/bin/env python3
"""First-order real and phase matter transport; no loop integration."""
from collections import Counter
import json
from pathlib import Path
import sympy as s

I=s.I
Q=s.Rational
positive=Counter()
negative=Counter()


class Exterior:
    """Independent odd components with an explicit ordered exterior product."""
    def __init__(self, terms=None):
        self.terms={k:s.expand(v) for k,v in (terms or {}).items() if s.expand(v)!=0}
    @staticmethod
    def cast(value):
        return value if isinstance(value,Exterior) else Exterior({():s.sympify(value)})
    def __add__(self,value):
        terms=dict(self.terms)
        for k,v in Exterior.cast(value).terms.items():
            terms[k]=terms.get(k,0)+v
        return Exterior(terms)
    __radd__=__add__
    def __neg__(self):
        return Exterior({k:-v for k,v in self.terms.items()})
    def __sub__(self,value):
        return self+-Exterior.cast(value)
    def __mul__(self,value):
        terms={}
        for a,x in self.terms.items():
            for b,y in Exterior.cast(value).terms.items():
                if set(a)&set(b):
                    continue
                key=tuple(sorted(a+b))
                parity=sum(i>j for i in a for j in b)
                terms[key]=terms.get(key,0)+(-1)**parity*x*y
        return Exterior(terms)
    __rmul__=__mul__
    def __truediv__(self,value):
        return self*s.sympify(value)**-1
    def adjoint(self,conjugate_ids):
        out=Exterior()
        for key,value in self.terms.items():
            term=Exterior.cast(s.conjugate(value))
            for index in reversed(key):
                term=term*Exterior({(conjugate_ids[index],):1})
            out=out+term
        return out


def zero(value):
    if isinstance(value,Exterior):
        return all(s.simplify(v)==0 for v in value.terms.values())
    if isinstance(value,s.MatrixBase):
        return all(s.simplify(v)==0 for v in value)
    return s.simplify(value)==0


def check(name,value):
    if not zero(value):
        raise RuntimeError(name+' failed')
    positive[name]+=1


def reject(name,value):
    if zero(value):
        raise RuntimeError(name+' mutation escaped')
    negative[name]+=1


def odd(index):
    return Exterior({(index,):1})


def main():
    alpha,g,t=s.symbols('alpha g t',real=True)
    x=s.symbols('x0:2',real=True)
    y=s.symbols('y0:2',real=True)
    dx=s.symbols('dx0:2',real=True)
    dy=s.symbols('dy0:2',real=True)
    phi=s.Matrix([x[j]+I*y[j] for j in range(2)])
    dphi=s.Matrix([dx[j]+I*dy[j] for j in range(2)])
    rho=(phi.H*phi)[0]
    u=(dphi.H*phi+phi.H*dphi)[0]
    z=(phi.H*dphi)[0]
    current=I*(z-s.conjugate(z))
    kinetic=(dphi.H*dphi)[0]
    variables=list(x)+list(y)
    jets=list(dx)+list(dy)
    def derivative(value):
        return value.applyfunc(lambda f:sum(s.diff(f,q)*dq for q,dq in zip(variables,jets)))
    delta_phi=alpha*(g+I*t)*rho*phi
    ddelta_phi=derivative(delta_phi)
    kinetic_change=(ddelta_phi.H*dphi+dphi.H*ddelta_phi)[0]
    norm_change=(delta_phi.H*phi+phi.H*delta_phi)[0]
    check('scalar complex first-jet kinetic variation',kinetic_change-alpha*(g*(2*rho*kinetic+u*u)-t*u*current))
    check('scalar density changes by radial component',norm_change-2*alpha*g*rho**2)
    check('complex scalar-current decomposition',g*(z*z+s.conjugate(z)**2)-g*(u*u-current**2)/2)
    primitive=alpha*(g*((u*u-current**2)/2-2*rho*kinetic)+t*u*current)
    residual=alpha*g*(3*u*u-current**2)/2
    check('full scalar derivative packet after combined map',primitive+kinetic_change-residual)
    check('residual in independent invariant notation',residual-2*alpha*g*(u*u-z*s.conjugate(z)))
    reject('erase all scalar derivative terms',residual)
    reject('drop derivative of nonlinear radial factor',primitive+2*alpha*g*rho*kinetic-alpha*t*u*current-residual)
    reject('radial scalar factor has wrong sign',primitive-kinetic_change-residual)
    m2,lam,h3=s.symbols('m2 lam h3',real=True)
    potential_change=-(m2+2*lam*rho)*norm_change
    check('potential variation uses signed quadratic and stated quartic',potential_change+2*alpha*g*m2*rho**2+4*alpha*g*lam*rho**3)
    check('raw mass-cubic pole plus transported potential',2*alpha*h3*rho**2+potential_change-(2*alpha*(h3-g*m2)*rho**2-4*alpha*g*lam*rho**3))
    reject('drop induced scalar sextic',potential_change+2*alpha*g*m2*rho**2)
    eps=s.Matrix([[0,1],[-1,0]])
    Higgs=-eps*phi.conjugate()
    delta_Higgs=-eps*delta_phi.conjugate()
    check('epsilon inverse for conjugate Higgs',eps*Higgs.conjugate()-phi)
    check('radial and opposite phase of Higgs',delta_Higgs-alpha*(g-I*t)*rho*Higgs)
    reject('give Higgs same complex phase as phi',delta_Higgs-alpha*(g+I*t)*rho*Higgs)

    # Arbitrary independent Grassmann field jets in every weak/spin/flavor slot.
    G=s.Matrix([[2,1+I],[1-I,3]])
    A=s.Matrix([[1,2*I],[-2*I,-2]])
    F=G/2+I*A
    real_trace=s.trace(G)
    phase_trace=s.trace(A)
    b=real_trace+I*phase_trace
    D=G/2+real_trace*s.eye(2)
    B=A-phase_trace*s.eye(2)
    U=D+I*B
    check('conjugate-safe matrix generator',F.H-(G/2-I*A))
    reject('flavor order is irrelevant',G*A-A*G)
    density,dr=s.symbols('rho dr',real=True)
    lepton=[[[odd(4*a+2*w+r) for r in range(2)] for w in range(2)] for a in range(2)]
    barred=[[[odd(8+4*a+2*w+r) for r in range(2)] for w in range(2)] for a in range(2)]
    dlep=[[[odd(16+4*a+2*w+r) for r in range(2)] for w in range(2)] for a in range(2)]
    dbar=[[[odd(24+4*a+2*w+r) for r in range(2)] for w in range(2)] for a in range(2)]
    def left(matrix,field):
        return [[[sum((matrix[a,c]*field[c][w][r] for c in range(2)),Exterior()) for r in range(2)] for w in range(2)] for a in range(2)]
    def right(field,matrix):
        return [[[sum((field[c][w][r]*matrix[c,a] for c in range(2)),Exterior()) for r in range(2)] for w in range(2)] for a in range(2)]
    def scale(field,factor):
        return [[[factor*field[a][w][r] for r in range(2)] for w in range(2)] for a in range(2)]
    def plus(a,b):
        return [[[a[i][w][r]+b[i][w][r] for r in range(2)] for w in range(2)] for i in range(2)]
    delta_l=scale(left(F,lepton),alpha*density)
    delta_bar=scale(right(barred,F.H),alpha*density)
    ddelta_l=scale(plus(scale(left(F,lepton),dr),scale(left(F,dlep),density)),alpha)
    ddelta_bar=scale(plus(scale(right(barred,F.H),dr),scale(right(dbar,F.H),density)),alpha)
    pauli=[s.eye(2),s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-I],[I,0]]),s.diag(1,-1)]
    def bilinear(row,matrix,col,spin):
        return sum((row[a][w][r]*matrix[a,c]*spin[r,j]*col[c][w][j] for a in range(2) for c in range(2) for w in range(2) for r in range(2) for j in range(2)),Exterior())
    for mu,spin in enumerate(pauli):
        change=I*Q(1,2)*(bilinear(delta_bar,s.eye(2),dlep,spin)+bilinear(barred,s.eye(2),ddelta_l,spin)-bilinear(ddelta_bar,s.eye(2),lepton,spin)-bilinear(dbar,s.eye(2),delta_l,spin))
        bilateral=I*(bilinear(barred,G,dlep,spin)-bilinear(dbar,G,lepton,spin))
        expected=alpha*density*bilateral/2-alpha*dr*bilinear(barred,A,lepton,spin)
        check('literal lepton kinetic slots mu='+str(mu),change-expected)
        reject('lepton real normalization doubled mu='+str(mu),change-expected-alpha*density*bilateral/2)
    # Full weak-doublet Yukawa variations, keeping the external odd components.
    p=[1+I,2-I]
    rh=sum(s.conjugate(v)*v for v in p)
    h=list(-eps*s.Matrix(p).conjugate())
    dp=[alpha*rh*b*v for v in p]
    dh=list(-eps*s.Matrix(dp).conjugate())
    dbars=scale(right(barred,F.H),alpha*rh)
    singlets=[[[odd(32+4*k+2*a+r) for r in range(2)] for a in range(2)] for k in range(4)]
    quarkbar=[[[odd(48+4*a+2*w+r) for r in range(2)] for w in range(2)] for a in range(2)]
    Yn=s.Matrix([[1+I,2],[-1,3-I]])
    Ye=s.Matrix([[2,1-I],[I,3]])
    Yu=s.Matrix([[1,2+I],[3-I,-1]])
    Yd=s.Matrix([[2-I,-1],[1+I,4]])
    def yukawa(row,ymat,scalar,singlet):
        return sum((row[a][w][r]*ymat[a,c]*scalar[w]*singlet[c][r] for a in range(2) for c in range(2) for w in range(2) for r in range(2)),Exterior())
    dy_nu=-yukawa(dbars,Yn,p,singlets[0])-yukawa(barred,Yn,dp,singlets[0])
    dy_e=-yukawa(dbars,Ye,h,singlets[1])-yukawa(barred,Ye,dh,singlets[1])
    dy_u=-yukawa(quarkbar,Yu,dp,singlets[2])
    dy_d=-yukawa(quarkbar,Yd,dh,singlets[3])
    check('full neutrino Yukawa weak spin flavor slots',dy_nu+alpha*rh*yukawa(barred,(D-I*B)*Yn,p,singlets[0]))
    check('full charged-lepton Yukawa weak spin flavor slots',dy_e+alpha*rh*yukawa(barred,(D-I*(A+phase_trace*s.eye(2)))*Ye,h,singlets[1]))
    check('full up-quark Yukawa weak spin flavor slots',dy_u+alpha*rh*b*yukawa(quarkbar,Yu,p,singlets[2]))
    check('full down-quark Yukawa weak spin flavor slots',dy_d+alpha*rh*s.conjugate(b)*yukawa(quarkbar,Yd,h,singlets[3]))
    reject('retain neutrino phase terms but drop radial Yukawa terms',dy_nu-I*alpha*rh*yukawa(barred,B*Yn,p,singlets[0]))
    reject('drop charged-lepton spectator',dy_e)
    reject('drop quark spectators',dy_u+dy_d)

    # Probe rows and columns are independent, so cancellation between them is forbidden.
    eta=[[[odd(56+4*a+2*w+r) for r in range(2)] for w in range(2)] for a in range(2)]
    bareta=[[[odd(64+4*a+2*w+r) for r in range(2)] for w in range(2)] for a in range(2)]
    probe=bilinear(bareta,s.eye(2),delta_l,s.eye(2))+bilinear(delta_bar,s.eye(2),eta,s.eye(2))
    check('lepton independent probe transport',probe-alpha*density*(bilinear(bareta,F,lepton,s.eye(2))+bilinear(barred,F.H,eta,s.eye(2))))
    reject('drop lepton probes',probe)
    j=s.Matrix([3+I,2-2*I])
    scalar_probe=(j.H*s.Matrix(dp)+s.Matrix(dp).H*j)[0]
    check('scalar independent probe transport',scalar_probe-alpha*rh*(b*(j.H*s.Matrix(p))[0]+s.conjugate(b)*(s.Matrix(p).H*j)[0]))
    reject('drop scalar probes',scalar_probe)

    # Literal C5 action: both slots are built out of the SAME barred exterior fields.
    C5=s.Matrix([[1+I,2-I],[2-I,3+2*I]])
    xb=[[sum((barred[a][w][r]*p[w] for w in range(2)),Exterior()) for r in range(2)] for a in range(2)]
    dxb=[[sum((dbars[a][w][r]*p[w]+barred[a][w][r]*dp[w] for w in range(2)),Exterior()) for r in range(2)] for a in range(2)]
    xc=[[sum((-eps[r,j]*xb[a][j] for j in range(2)),Exterior()) for r in range(2)] for a in range(2)]
    dxc=[[sum((-eps[r,j]*dxb[a][j] for j in range(2)),Exterior()) for r in range(2)] for a in range(2)]
    def c5op(matrix,row=xb,col=xc):
        return sum((Q(1,2)*row[a][r]*matrix[a,c]*col[c][r] for a in range(2) for c in range(2) for r in range(2)),Exterior())
    change5=c5op(C5,dxb,xc)+c5op(C5,xb,dxc)
    induced=(D-I*B)*C5+C5*(D-I*B).T
    check('dimension-five two actual conjugate slots',change5-alpha*rh*c5op(induced))
    check('dimension-five transport preserves coefficient symmetry',induced-induced.T)
    check('composite generator includes scalar conjugation',U-(F+s.conjugate(b)*s.eye(2)))
    reject('dimension-five omit second slot',change5-alpha*rh*c5op((D-I*B)*C5))
    reject('dimension-five replace transpose by adjoint',change5-alpha*rh*c5op((D-I*B)*C5+C5*(D-I*B).H))
    reject('dimension-five drop all radial transport',change5+I*alpha*rh*c5op(B*C5+C5*B.T))
    # Explicit anti-involution builds h.c. from independent conjugate exterior generators.
    conjugate_ids={k:(k+100 if k<100 else k-100) for k in range(200)}
    full5=change5+change5.adjoint(conjugate_ids)
    check('full dimension-five action remains Hermitian',full5-full5.adjoint(conjugate_ids))
    check('full dimension-five transport including h.c.',full5-alpha*rh*(c5op(induced)+c5op(induced).adjoint(conjugate_ids)))
    reject('dimension-five erase conditional insertion',full5)
    # First-order truncation is algebraic: the pole packet already carries alpha.
    seed,shift=s.symbols('seed shift')
    check('shifting a one-loop contact starts at second loop order',s.diff(alpha*seed*(1+alpha*shift)**4,alpha).subs(alpha,0)-seed)
    c7,rr=s.symbols('c7 rr')
    check('seed self-insertion is beyond linear seed order',s.diff(c7*(1+alpha*c7*rr),c7).subs(c7,0)-1)
    g1,a1=s.symbols('g1 a1',real=True)
    y1,p1=s.symbols('y1 p1')
    bar1,n1=odd(80),odd(81)
    f1=g1/2+I*a1
    b1=g1+I*a1
    actual_one=-(alpha*density*s.conjugate(f1)*bar1)*y1*p1*n1-bar1*y1*(alpha*density*b1*p1)*n1
    check('one flavor retains radial Yukawa source',actual_one+Q(3,2)*alpha*density*g1*bar1*y1*p1*n1)
    reject('one-flavor phase invariance erases radial source',actual_one)
    report={
      'status':'PASS','date':'30 October 2026 UTC',
      'scope':'First-order classical transport of previous primitive matter poles, with complete declared Yukawa and probe slots and conditional independent C5; no new loop integration or full physical coefficient.',
      'positive_checks':sum(positive.values()),'rejected_mutations':sum(negative.values()),
      'positive_categories':dict(positive),'negative_categories':dict(negative),
      'no_project_imports':True,'assertions_required':False,
      'scalar_test':'Arbitrary real scalar components and first jets; exact polynomial identities.',
      'fermion_test':'Independent exterior Grassmann weak, spin and flavor slots; noncommuting Hermitian generators.',
      'hermitian_completion':'Explicit conjugate-generator anti-involution reverses Grassmann order.'}
    Path(__file__).with_name('REAL_MATTER_RESULTS.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
