#!/usr/bin/env python3
"""Selected common classical action; no new loop or physical beta function."""
from collections import Counter
import json
from pathlib import Path
import sympy as s

I=s.I
Q=s.Rational
positive=Counter()
negative=Counter()


class Exterior:
    def __init__(self,terms=None):
        self.terms={k:s.expand(v) for k,v in (terms or {}).items() if s.expand(v)!=0}
    @staticmethod
    def cast(v):
        return v if isinstance(v,Exterior) else Exterior({():s.sympify(v)})
    def __add__(self,other):
        result=dict(self.terms)
        for k,v in Exterior.cast(other).terms.items():result[k]=result.get(k,0)+v
        return Exterior(result)
    __radd__=__add__
    def __neg__(self):return Exterior({k:-v for k,v in self.terms.items()})
    def __sub__(self,other):return self+-Exterior.cast(other)
    def __mul__(self,other):
        result={}
        for a,x in self.terms.items():
            for b,y in Exterior.cast(other).terms.items():
                if set(a)&set(b):continue
                sign=(-1)**sum(i>j for i in a for j in b)
                key=tuple(sorted(a+b))
                result[key]=result.get(key,0)+sign*x*y
        return Exterior(result)
    __rmul__=__mul__
    def __truediv__(self,value):return self*s.sympify(value)**-1
    def adjoint(self):
        result=Exterior()
        for key,value in self.terms.items():
            term=Exterior.cast(s.conjugate(value))
            for k in reversed(key):term=term*odd(k+100 if k<100 else k-100)
            result=result+term
        return result


def odd(i):return Exterior({(i,):1})
def zero(value):
    if isinstance(value,Exterior):return all(s.simplify(v)==0 for v in value.terms.values())
    if isinstance(value,s.MatrixBase):return all(s.simplify(v)==0 for v in value)
    return s.simplify(value)==0
def check(name,value):
    if not zero(value):raise RuntimeError(name+' failed')
    positive[name]+=1
def reject(name,value):
    if zero(value):raise RuntimeError(name+' mutation escaped')
    negative[name]+=1
def block(a,b,c,d):return a.row_join(b).col_join(c.row_join(d))
def transform(matrix,fields):
    return [[sum((matrix[a,b]*fields[b][r] for b in range(matrix.cols)),Exterior())
             for r in range(2)] for a in range(matrix.rows)]
def kinetic(rows,matrix,fields,spin):
    return sum((rows[a][r]*matrix[a,b]*spin[r,t]*fields[b][t]
                for a in range(matrix.rows) for b in range(matrix.cols)
                for r in range(2) for t in range(2)),Exterior())
def mass(matrix,fields):
    eps=s.Matrix([[0,1],[-1,0]])
    hol=sum((-Q(1,2)*fields[a][r]*matrix[a,b]*eps[r,t]*fields[b][t]
             for a in range(matrix.rows) for b in range(matrix.cols)
             for r in range(2) for t in range(2)),Exterior())
    return hol+hol.adjoint()
def mass_variation(matrix,fields,delta):
    eps=s.Matrix([[0,1],[-1,0]])
    hol=sum((-Q(1,2)*matrix[a,b]*eps[r,t]
             *(delta[a][r]*fields[b][t]+fields[a][r]*delta[b][t])
             for a in range(matrix.rows) for b in range(matrix.cols)
             for r in range(2) for t in range(2)),Exterior())
    return hol+hol.adjoint()


def main():
    # First variation with independent P,Q and R coefficients: the original
    # counterterm and primitive loop have opposite cancelling-action signs.
    alpha,x,y=s.symbols('alpha x y',real=True)
    a,b0,c,d,p,q,r,vs=s.symbols('a b c d p q r vs',real=True)
    S=x*x*y+2*x*y*y+3*x*x
    vC=s.Matrix([a,b0]);w=s.Matrix([c,d]);variables=[x,y]
    def first(v):return sum(s.diff(S,z)*v[k] for k,z in enumerate(variables))
    original=S+alpha*(p-q+r)
    shifted=original.subs({x:x+alpha*(a-c),y:y+alpha*(b0-d)},simultaneous=True)
    actual=s.diff(shifted,alpha).subs(alpha,0)
    check('common cancelling-action first variation',actual-(p+first(vC)-q-first(w)+r))
    reject('use matter-loop map with counterterm sign unchanged',actual-(p+first(vC)-q+first(w)+r))
    vS=s.Matrix([vs,vs*x])
    sourced=p+first(vS)
    check('source-to-common bridge subtracts source map',sourced+first(vC-vS-w)-(p+first(vC)-first(w)))
    reject('apply original normalizer directly to sourced coordinates',sourced+first(vC-w)-(p+first(vC)-first(w)))
    check('one-loop inputs change only at second loop order',s.diff(alpha*(p+alpha*x)*(q+alpha*y),alpha).subs(alpha,0)-p*q)

    # Arbitrary scalar first jets, using the negative of the earlier loop map.
    xx=s.symbols('x0:2',real=True);yy=s.symbols('y0:2',real=True)
    dx=s.symbols('dx0:2',real=True);dy=s.symbols('dy0:2',real=True)
    gg,tt=s.symbols('g_scalar t_scalar',real=True)
    scalar=s.Matrix([xx[k]+I*yy[k] for k in range(2)])
    jet=s.Matrix([dx[k]+I*dy[k] for k in range(2)])
    density=(scalar.H*scalar)[0]
    gradient=(jet.H*scalar+scalar.H*jet)[0]
    current=I*(scalar.H*jet-jet.H*scalar)[0]
    norm=(jet.H*jet)[0]
    dscalar=-density*(gg+I*tt)*scalar
    djet=dscalar.applyfunc(lambda f:sum(s.diff(f,v)*j for v,j in zip(list(xx)+list(yy),list(dx)+list(dy))))
    kinetic_change=(djet.H*jet+jet.H*djet)[0]
    scalar_loop=gg*((gradient**2-current**2)/2-2*density*norm)+tt*gradient*current
    scalar_remainder=gg*(3*gradient**2-current**2)/2
    check('literal arbitrary scalar cancelling-action transport',-scalar_loop+kinetic_change+scalar_remainder)
    check('negative matter map transports scalar density',(dscalar.H*scalar+scalar.H*dscalar)[0]+2*gg*density**2)
    reject('drop nonlinear scalar gradient in cancelling action',-scalar_loop-2*gg*density*norm+tt*gradient*current+scalar_remainder)

    Y=s.Matrix([[1+I,2-I],[-1+2*I,3]])
    C7=s.Matrix([[1+2*I,2-I],[2-I,-1+I]])
    M=s.Matrix([[2+I,1-I],[1-I,4-2*I]])
    HL=(Y*Y.H).applyfunc(s.expand);HR=(Y.H*Y).applyfunc(s.expand);TY=s.trace(HR)
    Fs=(C7*Y.conjugate()).applyfunc(s.expand);Km=(Fs*M*Y.H).applyfunc(s.expand)
    G=((Km+Km.H)/2).applyfunc(s.expand);A=((Km-Km.H)/(2*I)).applyfunc(s.expand)
    g=s.simplify(s.trace(G));t=s.simplify(s.trace(A));W=G/2+I*A;b=g+I*t
    if zero(G*A-A*G) or zero(g) or zero(t):
        raise RuntimeError('Complex noncommuting flavor and nonzero scalar trace witness required')
    check('seed and retained mass symmetric',C7-C7.T+M-M.T)
    check('matter Hermitian decomposition',Km-G-I*A)
    check('matter trace convention',s.trace(Km)-b)
    u=Q(3,2);m2=Q(-2,3);momentum2=Q(19);P=momentum2+m2
    H=u*P*Fs;empty=s.zeros(2)
    Z=block(empty,H,H.H,empty)
    Ds=Z+block(u*u*G,empty,empty,empty)
    Rz=block(-HL/2,empty,empty,-HR.T)
    J=block(-u*u*W,-H/2,-H.H/2,empty)
    check('selected neutral kinetic cancellation',Ds+J+J.H)
    check('renormalizable kinetic poles remain',Ds+Rz+J+J.H-Rz)
    reject('claim entire kinetic action is canonical',Rz)
    Jwrong=block(u*u*W,-H/2,-H.H/2,empty)
    reject('forget matter cancelling-action sign',Ds+Jwrong+Jwrong.H)
    reject('discard phase part of chosen field coordinates',J+Ds/2)
    fields=[[odd(2*a+r) for r in range(2)] for a in range(4)]
    rows=[[odd(100+2*a+r) for r in range(2)] for a in range(4)]
    delta=transform(J,fields)
    # A barred row transforms with J^dagger; transform acts on a column of rows.
    drows=transform(J.conjugate(),rows)
    spin_matrices=[s.eye(2),s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-I],[I,0]]),s.diag(1,-1)]
    for mu,spin in enumerate(spin_matrices):
        actual=kinetic(drows,s.eye(4),fields,spin)+kinetic(rows,s.eye(4),delta,spin)
        check('literal exterior kinetic variation '+str(mu),actual+kinetic(rows,Ds,fields,spin))
        check('literal kinetic remainder '+str(mu),actual+kinetic(rows,Ds+Rz,fields,spin)-kinetic(rows,Rz,fields,spin))
    B=u*Y.conjugate();Mc=block(empty,B,B.T,M.conjugate())
    scalarB=-u**3*s.conjugate(b)*Y.conjugate()
    scalarM=block(empty,scalarB,scalarB.T,empty)
    dMc=J.T*Mc+Mc*J+scalarM
    dA=-(H.conjugate()*B.T+B*H.H)/2
    dB=-H.conjugate()*M.conjugate()/2-u**3*(W.T+s.conjugate(b)*s.eye(2))*Y.conjugate()
    dD=-(H.T*B+B.T*H)/2
    expected=block(dA,dB,dB.T,dD)
    check('common neutral mass block reconstruction',dMc-expected)
    check('common neutral mass remains symmetric',dMc-dMc.T)
    actual_mass=mass_variation(Mc,fields,delta)+mass(scalarM,fields)
    check('literal exterior mass and scalar background transport',actual_mass-mass(expected,fields))
    reject('omit scalar background variation',actual_mass-mass(J.T*Mc+Mc*J,fields))
    reject('replace holomorphic transpose by adjoint',actual_mass-mass(J.H*Mc+Mc*J+scalarM,fields))
    reject('halve Majorana field variation again',actual_mass-mass((J.T*Mc+Mc*J)/2+scalarM,fields))
    wrong_scalar=block(empty,-u**3*b*Y.conjugate(),(-u**3*b*Y.conjugate()).T,empty)
    reject('give holomorphic background the wrong phase',actual_mass-mass(J.T*Mc+Mc*J+wrong_scalar,fields))
    probes=[[odd(20+2*a+r) for r in range(2)] for a in range(4)]
    probe_actual=kinetic(probes,s.eye(4),delta,s.eye(2))
    check('independent probe row follows common Jacobian',probe_actual-kinetic(probes,J,fields,s.eye(2)))
    reject('omit common probe transport',probe_actual)

    # Bare-parameter versus literal-action coefficients: all slots included.
    zL=-HL/4;zR=-HR/2;zNL=-HR.T/2;zphi=-TY/2
    rY=HL*Y/4+Y*HR/2+TY*Y/2
    rM=(HR.T*M+M*HR)/2
    check('renormalizable Yukawa parameter versus action closure',rY+zL.H*Y+Y*zR+zphi*Y)
    check('renormalizable Majorana parameter versus action closure',rM+zR.T*M+M*zR)
    TM=s.trace(HR*M.H*M);rm=m2*TY-2*TM
    check('scalar mass parameter versus action coefficient',rm+2*m2*zphi+2*TM)
    deltaBren=u*(rY.conjugate()+zL.T*Y.conjugate()+Y.conjugate()*zNL+zphi*Y.conjugate())
    deltaDren=(rM+zR.T*M+M*zR).conjugate()
    check('full known neutral renormalizable mass input closure',block(empty,deltaBren,deltaBren.T,deltaDren))
    reject('mistake zero proper Yukawa pole for zero parameter pole',rY)
    reject('mistake zero proper Majorana pole for zero parameter pole',rM)
    reject('replace Majorana transpose by dagger',rM+zR.H*M+M*zR)

    # Literal two barred slots of the original seed, with its half.
    eps=s.Matrix([[0,1],[-1,0]])
    barred=rows[:2]
    xb=[[u*barred[a][r] for r in range(2)] for a in range(2)]
    dxbar=transform(zL.conjugate(),barred)
    dxb=[[u*(dxbar[a][r]+zphi*barred[a][r]) for r in range(2)] for a in range(2)]
    def conjugate_slot(row):
        return [[sum((-eps[r,k]*row[a][k] for k in range(2)),Exterior())
                 for r in range(2)] for a in range(len(row))]
    xc=conjugate_slot(xb);dxc=conjugate_slot(dxb)
    def seed(matrix,row=xb,col=xc):
        term=sum((Q(1,2)*momentum2*row[a][r]*matrix[a,c]*col[c][r]
                  for a in range(2) for c in range(2) for r in range(2)),Exterior())
        return term+term.adjoint()
    seed_change=seed(C7,dxb,xc)+seed(C7,xb,dxc)
    r7=(HL*C7+C7*HL.T)/4+TY*C7
    check('literal tree seed external-field pole',seed_change+seed(r7))
    reject('omit tree seed field transport as if already one loop',seed_change)
    reject('count field-compensation parameter part twice',seed_change+2*seed(r7))
    reject('replace seed transpose by dagger',seed_change+seed((HL*C7+C7*HL.H)/4+TY*C7))
    check('tree seed variation under inverse field normalization',-seed_change-seed(r7))

    # Conditional C5: actual active and retained-sterile slots from V=vC-w.
    C5=s.Matrix([[1+I,3-I],[3-I,2-2*I]])
    E=W.H+b*s.eye(2)
    # delta bar X = -u^3 bar nu(W^dagger+bI) -u^2 P bar N Fs^dagger/2.
    active_change=transform((-u**3*E).T,barred)
    mixed_change=transform((-u**2*P*Fs.H/2).T,rows[2:])
    dxbV=[[active_change[a][r]+mixed_change[a][r] for r in range(2)] for a in range(2)]
    dxcV=conjugate_slot(dxbV)
    def five(matrix,row=xb,col=xc):
        term=sum((Q(1,2)*row[a][r]*matrix[a,c]*col[c][r]
                  for a in range(2) for c in range(2) for r in range(2)),Exterior())
        return term+term.adjoint()
    actual5=five(C5,dxbV,xc)+five(C5,xb,dxcV)
    active5=-u*u*five(E*C5+C5*E.T)
    mixed5=five(C5,mixed_change,xc)+five(C5,xb,conjugate_slot(mixed_change))
    check('conditional C5 full active and retained-field transport',actual5-active5-mixed5)
    reject('drop conditional retained-sterile contact',actual5-active5)
    reject('keep only conditional retained-sterile contact',actual5-mixed5)
    reject('replace conditional transpose by dagger',actual5+u*u*five(E*C5+C5*E.H)-mixed5)

    # Potential sign and nonzero one-flavor radial source survive the change.
    rho,lam,greal,k3real=s.symbols('rho lam greal k3real',real=True)
    V=m2*rho+lam*rho*rho
    common_potential=-2*k3real*rho*rho+2*greal*rho*rho*s.diff(V,rho)
    check('common cancelling scalar potential',common_potential-(-2*(k3real-greal*m2)*rho*rho+4*greal*lam*rho**3))
    g1,a1=s.symbols('g1 a1',real=True)
    check('one-flavor cancelling radial Yukawa coefficient',s.conjugate(g1/2+I*a1)+(g1+I*a1)-3*g1/2)
    reject('one-flavor phase cancellation erases radial source',3*g1/2)
    report={'status':'PASS','date':'4 November 2026 Pacific',
        'input_commit':'41f9b709c0fafdefdc9e1b990caa5421ee99f7dc',
        'scope':'Selected common classical cancelling action P-Q+R, source/matter coordinate composition, neutral kinetic/mass/probe projection, known bare-parameter closure, literal tree-seed field compensation and conditional C5 partners. No new loop, full basis or physical beta function.',
        'positive_checks':sum(positive.values()),'rejected_mutations':sum(negative.values()),
        'positive_categories':dict(positive),'negative_categories':dict(negative),
        'neutral_witness':{'scalar_u':str(u),'momentum_squared':str(momentum2),
                           'real_trace_g':str(g),'imaginary_trace_t':str(t),
                           'scalar_background_Dirac_block':[[str(s.simplify(v)) for v in scalarB.row(a)] for a in range(2)]},
        'no_project_imports':True,'assertions_required':False}
    Path(__file__).with_name('COMMON_ACTION_RESULTS.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'positive_checks':report['positive_checks'],'rejected_mutations':report['rejected_mutations']}))


if __name__=='__main__':main()
