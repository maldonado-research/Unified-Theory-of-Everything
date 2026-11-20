#!/usr/bin/env python3
"""Transport the calculated primitive poles; no new loop integration.

Explicit exterior-algebra flavor/weak projection, bare-mu weights and a
source-retaining Higgs phase substitution. No project imports or asserts.
"""
from collections import Counter
from pathlib import Path
import json
import sympy as s

I=s.I; R=s.Rational
positive=Counter(); negative=Counter()

def simp(x): return s.simplify(s.expand(x))
def zero(x):
    if isinstance(x,G): return not x.d
    if isinstance(x,s.MatrixBase): return all(simp(v)==0 for v in x)
    return simp(x)==0

def check(name,x):
    if not zero(x): raise RuntimeError(name+': exact identity failed')
    positive[name]+=1

def reject(name,x):
    if zero(x): raise RuntimeError(name+': mutation escaped')
    negative[name]+=1

class G:
    def __init__(self,data=None):
        self.d={k:simp(v) for k,v in (data or {}).items() if simp(v)!=0}
    @staticmethod
    def cast(x): return x if isinstance(x,G) else G({():s.sympify(x)})
    def __add__(self,x):
        out=dict(self.d)
        for k,v in G.cast(x).d.items(): out[k]=out.get(k,0)+v
        return G(out)
    __radd__=__add__
    def __neg__(self): return G({k:-v for k,v in self.d.items()})
    def __sub__(self,x): return self+-G.cast(x)
    def __rsub__(self,x): return G.cast(x)+-self
    def __mul__(self,x):
        out={}
        for a,v in self.d.items():
            for b,w in G.cast(x).d.items():
                if set(a)&set(b): continue
                sign=(-1)**sum(i>j for i in a for j in b)
                k=tuple(sorted(a+b));out[k]=out.get(k,0)+sign*v*w
        return G(out)
    __rmul__=__mul__
    def dagger(self):
        out=G()
        for a,v in self.d.items():
            term=G.cast(s.conjugate(v))
            for j in reversed(a): term=term*G({((j+8)%16,):1})
            out+=term
        return out


def main():
    Y=s.Matrix([[1+I,2-I],[3-2*I,-1+I]])
    M=s.Matrix([[2+I,1-2*I],[1-2*I,4+I]])
    C=s.Matrix([[2-I,1+2*I],[1+2*I,3+I]])
    Q=Y.conjugate()*M*Y.H
    Q3=Y.conjugate()*M*M.H*M*Y.H
    k1=s.trace(C*Q);k3=s.trace(C*Q3)
    check('symmetric four-lepton mass-flavor source',Q-Q.T)
    check('symmetric cubic mass-flavor source',Q3-Q3.T)
    E=s.Matrix([[0,1],[-1,0]])
    def var(a,u,r,bar=False): return G({(8*int(bar)+4*a+2*u+r,):1})
    def pair(a,u,b,v,bar=False):
        return sum((E[r,t]*var(a,u,r,bar)*var(b,v,t,bar)
                    for r in range(2) for t in range(2)),G())
    for a in range(2):
        for b in range(2):
            for u in range(2):
                for v in range(2):
                    check('actual spin-singlet pair exchange',pair(a,u,b,v)-pair(b,v,a,u))
                    check('actual conjugate pair minus',pair(a,u,b,v).dagger()+pair(a,u,b,v,True))
    basis=[s.Matrix([[1,0],[0,0]]),s.Matrix([[0,1],[1,0]])/s.sqrt(2),s.Matrix([[0,0],[0,1]])]
    c=s.Matrix([s.trace(e.T*C) for e in basis]);r=s.Matrix([s.trace(e.T*Q) for e in basis]);v=r.conjugate()
    H=(c*v.H+v*c.H).applyfunc(simp)
    check('Hermitian primitive pair coefficient',H-H.H)
    for e in basis:
        for f in basis: check('orthonormal symmetric pair coordinates',s.trace(e.T*f)-int(e==f))
    qs={}
    for a,e in enumerate(basis):
        for b,f in enumerate(basis):
            qs[a,b]=sum((e[i,j]*f[u,w]*pair(i,u,j,w)
                         for i in range(2) for j in range(2) for u in range(2) for w in range(2)),G())
    raw=sum((C[a,b]*Q[c0,d]*pair(a,u,b,w,True)*pair(c0,u,d,w)
             for a in range(2) for b in range(2) for c0 in range(2) for d in range(2)
             for u in range(2) for w in range(2)),G())
    # The independently integrated pole is kappa(raw+raw^dagger)/2.
    # Actual normalized pair substitution fixes the matrix normalization.
    projected=sum((-H[a,b]*qs[a,t].dagger()*qs[b,t]
                   for a in range(3) for b in range(3) for t in range(3)),G())
    check('complete normalized flavor-weak exterior projection',raw+raw.dagger()-projected)
    reject('extra identical-pair half',raw+raw.dagger()-projected*R(1,2))
    reject('omit Hermitian-conjugate seed',raw-projected)
    check('nonzero three-flavor-pair example has rank two',H.rank()-2)
    check('three weak-triplet copies have rank six',s.kronecker_product(s.eye(3),H).rank()-6)
    # A general algebraic proof is Hx=c(vdag x)+v(cdag x): its image lies
    # in span(c,v). The following symbolic Gram identity fixes its spectrum.
    cr=s.symbols('cr0:3',real=True);ci=s.symbols('ci0:3',real=True)
    vr=s.symbols('vr0:3',real=True);vi=s.symbols('vi0:3',real=True)
    cc=s.Matrix([a+I*b for a,b in zip(cr,ci)]);vv=s.Matrix([a+I*b for a,b in zip(vr,vi)])
    hh=cc*vv.H+vv*cc.H;ov=(vv.H*cc)[0];nc=(cc.H*cc)[0];nv=(vv.H*vv)[0]
    wedge=sum((s.conjugate(cc[a]*vv[b]-cc[b]*vv[a])*(cc[a]*vv[b]-cc[b]*vv[a])
               for a in range(3) for b in range(a+1,3)),s.Integer(0))
    check('general Gram determinant is sum of squared minors',nc*nv-s.conjugate(ov)*ov-wedge)
    check('general Hermitian pair trace',s.trace(hh)-2*s.re(ov))
    check('general Hermitian pair squared trace',s.trace(hh*hh)-2*s.re(ov**2)-2*nc*nv)
    # Concrete exact eigenvalue checks in the two-dimensional support.
    a=(v.H*c)[0];g=(c.H*c)[0]*(v.H*v)[0]-s.conjugate(a)*a
    if not bool(simp(g)>0): raise RuntimeError('independent-source Gram witness failed')
    check('pair characteristic constant on support',((s.trace(H)**2-s.trace(H*H))/2)+g)
    reject('primitive pole asserted positive semidefinite',g*0-g)
    V=s.diag(I,(3+4*I)/5)*s.Matrix([[R(3,5),R(4,5)],[-R(4,5),R(3,5)]])
    U=s.Matrix([[R(5,13),R(12,13)],[-R(12,13),R(5,13)]])*s.diag((5+12*I)/13,I)
    Yp=V.H*Y*U.conjugate();Mp=U.H*M*U.conjugate();Cp=V.H*C*V.conjugate()
    Qp=(Yp.conjugate()*Mp*Yp.H).applyfunc(simp)
    Q3p=(Yp.conjugate()*Mp*Mp.H*Mp*Yp.H).applyfunc(simp)
    check('complex mass source active-sterile covariance',Qp-V.T*Q*V)
    check('cubic complex mass source active-sterile covariance',Q3p-V.T*Q3*V)
    check('scalar mass-linear flavor invariant',s.trace(Cp*Qp)-k1)
    check('scalar cubic-mass flavor invariant',s.trace(Cp*Q3p)-k3)
    D=s.Matrix(3,3,lambda a,b:s.trace(basis[a].T*V.H*basis[b]*V.conjugate())).applyfunc(simp)
    check('normalized symmetric pair representation is unitary',D.H*D-s.eye(3))
    cp=s.Matrix([s.trace(e.T*Cp) for e in basis]);vp=s.Matrix([s.trace(e.T*Qp) for e in basis]).conjugate()
    check('Hermitian primitive pole flavor-coordinate covariance',cp*vp.H+vp*cp.H-D*H*D.H)
    reject('wrong complex sterile mass conjugation',Q-Y.conjugate()*M.H*Y.H)
    # d-dimensional engineering powers for each four-field target.
    eps=s.symbols('epsilon')
    for name,opdim,dim4 in [('four-fermion',6-4*eps,-2),('four-scalar derivative',6-4*eps,-2),('four-scalar potential',4-4*eps,0)]:
        check(name+' bare mu exponent',(4-2*eps)-opdim-(dim4+2*eps))
    # Quadratic Yukawa residues linear in C have weight4; four-field
    # recipients have bare exponent2. M carries no explicit bare mu power.
    cs,cb,y,yb,m=s.symbols('C Cbar Y Ybar M');weighted=[(cs,2),(cb,2),(y,1),(yb,1)]
    def beta(residue): return s.expand(sum(w*x*s.diff(residue,x) for x,w in weighted)-2*residue)
    for name,residue in [('target seed',cs*y*yb),('four-lepton pair',R(1,2)*(cs*yb**2*m+cb*y**2*m)),('complex scalar derivative',-cs*yb**2*m),('quartic potential parameter',cs*yb**2*m**3+cb*y**2*m**3)]:
        check(name+' exact canonical-epsilon conversion',beta(residue)-2*residue)
        reject(name+' pole directly identified with beta',beta(residue)-residue)
    # Direct scalar phase change, with independent first jets. phi=x+i z.
    xs=s.symbols('x0:2',real=True);zs=s.symbols('z0:2',real=True)
    dx=s.symbols('dx0:2',real=True);dz=s.symbols('dz0:2',real=True)
    variables=list(xs)+list(zs);jets=list(dx)+list(dz)
    phi=s.Matrix([x+I*z for x,z in zip(xs,zs)]);dp=s.Matrix([x+I*z for x,z in zip(dx,dz)])
    rho=(phi.H*phi)[0];dr=sum(s.diff(rho,x)*j for x,j in zip(variables,jets))
    aa=(phi.H*dp)[0];bb=(dp.H*phi)[0];current=I*(aa-bb);A=rho*(dp.H*dp)[0]
    kr,ki=s.symbols('K_real K_imag',real=True);K=kr+I*ki
    check('complete complex scalar derivative decomposition',K*(aa**2-A)+s.conjugate(K)*(bb**2-A)-(kr*((dr**2-current**2)/2-2*A)+ki*dr*current))
    reject('replace complex derivative coefficient by its real part',K*(aa**2-A)+s.conjugate(K)*(bb**2-A)-kr*(aa**2+bb**2-2*A))
    alpha=s.symbols('alpha',real=True)
    variation=I*alpha*rho*phi
    dvariation=variation.applyfunc(lambda f:sum(s.diff(f,x)*j for x,j in zip(variables,jets)))
    check('actual Higgs phase kinetic variation',(dvariation.H*dp+dp.H*dvariation)[0]+alpha*dr*current)
    check('Higgs norm fixed by phase substitution',(variation.H*phi+phi.H*variation)[0])
    realvector=[-alpha*rho*z for z in zs]+[alpha*rho*x for x in xs]
    check('real scalar phase vector has zero local divergence',sum(s.diff(f,x) for f,x in zip(realvector,variables)))
    # Arbitrary independent interaction/probe rows retain source transport.
    w=s.Matrix([[1+I,2-I]]);probe=s.Matrix([[3-I,-2+I]])
    VY=(w*phi)[0];VYc=s.conjugate(VY)
    deltaY=-(w*variation)[0]-s.conjugate((w*variation)[0])
    check('literal pure-Yukawa source transport',deltaY-(-I*alpha*rho*VY+I*alpha*rho*VYc))
    deltaProbe=(probe*variation)[0]+s.conjugate((probe*variation)[0])
    check('literal external scalar probe transport',deltaProbe-I*alpha*rho*((probe*phi)[0]-s.conjugate((probe*phi)[0])))
    reject('drop interaction source after scalar phase change',deltaY)
    reject('drop probe source after scalar phase change',deltaProbe)
    check('imaginary loop packet canceled by its declared phase change',ki*dr*current+(dvariation.H*dp+dp.H*dvariation)[0].subs(alpha,ki))
    result={'status':'PASS','input_commit':'f1f9a892e15545c70e4fa05827ea7fefc5a438eb','date':'23 October 2026 Pacific',
        'scope':'Algebraic transport of independently integrated primitive poles and restricted unreduced Green-coordinate RG inputs; no additional loop or complete physical beta function.',
        'positive_checks':sum(positive.values()),'negative_checks':sum(negative.values()),'positive_categories':dict(positive),'negative_categories':dict(negative),
        'four_lepton_pair_rank_bound':2,'weak_triplet_multiplicity':3,'pair_source_subspace':'Sym^2(active flavor) tensor Sym^2(weak doublet)',
        'new_loop_integration':False,'repository_imports':False,'assertions_required':False}
    Path(__file__).with_name('GREEN_MIXING_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
