#!/usr/bin/env python3
"""Independent actual neutral-action/field-map check, exact Grassmann modes.

SymPy 1.14.0; no project imports. Print JSON with python -I -B -O.
All polynomials are explicitly truncated at fermion degree two. Two opposite
Fourier modes and their independent Grassmann adjoints retain Majorana mass
and kinetic terms simultaneously. Constant scalar background need not solve
its equation. This is a classical action identity, not a quantum 1PI claim.
"""
import json
import sys
import sympy as s

I,R=s.I,s.Rational
checks,mutations=[],[]
momenta,conjugates={},{}

class Poly:
    def __init__(self,data=None):
        self.d={k:s.expand(v) for k,v in (data or {}).items() if s.expand(v)!=0}
    @staticmethod
    def constant(x):return Poly({():x})
    def __add__(self,other):
        if not isinstance(other,Poly):other=Poly.constant(other)
        out=self.d.copy()
        for k,v in other.d.items():out[k]=out.get(k,0)+v
        return Poly(out)
    __radd__=__add__
    def __neg__(self):return Poly({k:-v for k,v in self.d.items()})
    def __sub__(self,other):return self+-other
    def __rsub__(self,other):return -self+other
    def __mul__(self,other):
        if not isinstance(other,Poly):other=Poly.constant(other)
        out={}
        for a,x in self.d.items():
            for b,y in other.d.items():
                if len(a)+len(b)>2 or set(a)&set(b):continue
                inv=sum(i>j for i in a for j in b)
                key=tuple(sorted(a+b))
                out[key]=out.get(key,0)+(-1)**inv*x*y
        return Poly(out)
    __rmul__=__mul__
    def __truediv__(self,x):return self*(1/s.sympify(x))
    def dag(self):
        out=Poly()
        for key,value in self.d.items():
            term=Poly.constant(s.conjugate(value))
            for a in reversed(key):term=term*Poly({(conjugates[a],):1})
            out=out+term
        return out
    def derivative(self):
        return Poly({k:-I*sum(momenta[a] for a in k)*v for k,v in self.d.items()})
    def integrated(self):
        return Poly({k:v for k,v in self.d.items() if sum(momenta[a] for a in k)==0})

def check(name,p):
    if isinstance(p,Poly):bad=p.d
    elif isinstance(p,s.MatrixBase):bad=any(s.expand(v)!=0 for v in p)
    else:bad=s.expand(p)!=0
    if bad:raise RuntimeError(name+': '+str(bad))
    checks.append(name)

def reject(name,p):
    if isinstance(p,Poly):bad=p.d
    elif isinstance(p,s.MatrixBase):bad=any(s.expand(v)!=0 for v in p)
    else:bad=s.expand(p)!=0
    if not bad:raise RuntimeError(name+': mutation was not detected')
    mutations.append(name)

z=s.zeros(2)
pauli=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-I],[I,0]]),s.diag(1,-1)]
gamma=[s.BlockMatrix([[z,s.eye(2)],[s.eye(2),z]]).as_explicit()]+[s.BlockMatrix([[z,p],[-p,z]]).as_explicit() for p in pauli]
C=I*gamma[2]*gamma[0]

def matvec(mat,v):return [sum((v[j]*mat[i,j] for j in range(len(v))),Poly()) for i in range(mat.rows)]
def addv(a,b):return [x+y for x,y in zip(a,b)]
def scale(a,c):return [x*c for x in a]
def dot(a,b):return sum((x*y for x,y in zip(a,b)),Poly())
def bar(v):return [sum((v[j].dag()*gamma[0][j,i] for j in range(4)),Poly()) for i in range(4)]
def cc(v):return matvec(C,bar(v))
def Avec(v):return scale(matvec(gamma[0],[p.derivative() for p in v]),I)
def P0vec(v):return [-p.derivative().derivative() for p in v]
def Krow(v):return [sum((p.derivative()*(-I)*gamma[0][j,i] for j,p in enumerate(v)),Poly()) for i in range(4)]
def flav(mat,fields):
    return [sumv([scale(v,mat[i,j]) for j,v in enumerate(fields)]) for i in range(mat.rows)]
def sumv(vectors):
    out=[Poly() for _ in range(4)]
    for v in vectors:out=addv(out,v)
    return out
def paired(rows,mat,cols):
    return sum((dot(row,col)*mat[a,i] for a,row in enumerate(rows) for i,col in enumerate(cols)),Poly())
def hconj(p):return p+p.dag()

def makefields(n,E):
    fields=[]
    for f in range(n):
        out=[Poly() for _ in range(4)]
        for spin in range(2):
            for mode in [-1,1]:
                a=len(momenta);b=a+1
                momenta[a]=mode*E;momenta[b]=-mode*E
                conjugates[a]=b;conjugates[b]=a
                out[spin]=out[spin]+Poly({(a,):1})
        fields.append(out)
    return fields

def run_trial(label,u,m2,lam,E,Y,Cseed,M):
    momenta.clear();conjugates.clear()
    n,ns=Y.shape
    q=makefields(n+ns,E);nu,nl=q[:n],q[n:]
    nr=[cc(v) for v in nl]
    F=Cseed*Y.conjugate()
    c=2*m2/3+lam*u*u/3
    a=m2+2*lam*u*u
    bnu=[bar(v) for v in nu]
    SL=flav(u*Y,nr);bSL=[bar(v) for v in SL]
    SNL=flav(u*Y.conjugate().T,nu)
    SR=[cc(v) for v in SNL]
    RR=[addv(x,y) for x,y in zip(flav(M,nr),SR)]
    TT=[sumv(v) for v in zip(flav(M*M,nl),flav(M,SNL),[Avec(v) for v in SR])]
    UU=[sumv(v) for v in zip(flav(M*M*M,nr),flav(M*M,SR),flav(M,[Avec(v) for v in SNL]),[P0vec(v) for v in SR])]
    pred=(paired(bnu,F,UU)*R(5,6)*u+paired(bnu,F,RR)*c*u
          +paired([Krow(v) for v in bSL],F,RR)*u/3
          -paired(bSL,F,TT)*u/6+paired(bSL,F,nl)*m2*u/2)
    red=hconj(pred).integrated()
    # The actual five nonzero placements; the scalar-Laplacian placement is zero.
    original=(paired(bnu,F,[Avec(P0vec(v)) for v in nl])*R(5,6)*u
              +paired([P0vec(v) for v in bnu],F,[Avec(v) for v in nl])*u/3
              -paired([Krow(v) for v in bnu],F,[P0vec(v) for v in nl])*u/6
              +paired(bnu,F,[Avec(v) for v in nl])*m2*u/2
              +paired([Krow(v) for v in bnu],F,nl)*m2*u/2)
    original=hconj(original).integrated()
    compact=hconj(paired(bnu,F,[Avec(addv(P0vec(v),scale(v,m2))) for v in nl])*u).integrated()
    check(label+': actual six-placement neutral packet',original-compact)
    # Evaluate the reference action itself, with both Majorana slots.
    def reference(fields):
        active,sterile=fields[:n],fields[n:]
        right=[cc(v) for v in sterile]
        major=[addv(l,r) for l,r in zip(sterile,right)]
        ba=[bar(v) for v in active];bm=[bar(v) for v in major]
        kinetic=paired(ba,s.eye(n),[Avec(v) for v in active])
        kinetic+=paired(bm,s.eye(ns),[Avec(v) for v in major])/2
        mass=-paired(bm,M,major)/2
        yuk=hconj(-paired(ba,Y,right)*u)
        return (kinetic+mass+yuk).integrated()
    def variation(shift):
        return reference([addv(v,d) for v,d in zip(q,shift)])-reference(q)-reference(shift)
    ell=flav(u*F,[addv(scale(P0vec(v),R(1,6)),scale(v,m2/2)) for v in nl])
    x=[sumv(v) for v in zip(flav(R(5,6)*u*(F*M*M).conjugate().T,nu),
          flav(u*F.conjugate().T,[addv(scale(P0vec(v),R(5,6)),scale(v,c)) for v in nu]),
          flav(u*u*F.conjugate().T*Y/6,[Avec(v) for v in nr]))]
    y=[addv(p,r) for p,r in zip(flav(R(5,6)*u*(F*M).conjugate().T,[Avec(v) for v in nu]),
                                     flav(-u*u*(F*M).conjugate().T*Y/6,nr))]
    vS=[scale(v,-1) for v in ell]+[scale(addv(xx,cc(yy)),-1) for xx,yy in zip(x,y)]
    pfield=lambda v:addv(P0vec(v),scale(v,m2))
    vC=[scale(v,-1) for v in flav(u*F/2,[pfield(v) for v in nl])]
    vC += [scale(v,-1) for v in flav(u*F.conjugate().T/2,[pfield(v) for v in nu])]
    bridge=[addv(vc,scale(vs,-1)) for vc,vs in zip(vC,vS)]
    for name,shift in [('source',vS),('normalizing',vC),('bridge',bridge)]:
        for i,v in enumerate(shift):
            check(f'{label}: {name} left chirality {i}',sum((p for p in v[2:]),Poly()))
    t=paired(bnu,F,[Avec(v) for v in nl])/6
    # Direct potential expansion through degree two, rather than invoking a
    # Hessian formula. The exterior algebra truncates only fermion degree.
    potential=lambda phi,phibar:-(phibar*phi)*m2-(phibar*phi)*(phibar*phi)*lam
    pot0=potential(Poly.constant(u),Poly.constant(u))
    tadS=(potential(Poly.constant(u)-t.dag(),Poly.constant(u)-t)-pot0).integrated()
    tadB=(potential(Poly.constant(u)+t.dag(),Poly.constant(u)+t)-pot0).integrated()
    predicted_tad=hconj(t)*a*u
    check(label+': actual scalar potential source-map Hessian',tadS-predicted_tad.integrated())
    check(label+': scalar bridge reverses tadpole Hessian',tadB+tadS)
    dS=variation(vS);dC=variation(vC);dB=variation(bridge)
    check(label+': direct reduced action vs actual full source substitution',red-original-dS-tadS)
    check(label+': composed bridge full quadratic action',red+dB+tadB-original-dC)
    check(label+': first-order map composition',dB-dC+dS)
    # Independently build canonical remaining mass terms from the holomorphic
    # transpose rule in its equivalent four-component full-field variation:
    # kinetic part of the normalizing variation must remove the actual packet.
    def kinetic(fields):
        active,sterile=fields[:n],fields[n:]
        return (paired([bar(v) for v in active],s.eye(n),[Avec(v) for v in active])
                +paired([bar(v) for v in sterile],s.eye(ns),[Avec(v) for v in sterile])).integrated()
    dk=kinetic([addv(v,d) for v,d in zip(q,vC)])-kinetic(q)-kinetic(vC)
    check(label+': actual normalization removes original kinetic packet',dk+original)
    canonical_remainder=original+dC
    reject(label+': induced mass/Yukawa remainder discarded',canonical_remainder)
    reject(label+': reverse bridge fermion sign',red-dB+tadB-original-dC)
    reject(label+': Majorana source shift half-sized',red-original-variation(vS[:n]+[scale(v,R(1,2)) for v in vS[n:]])-tadS)
    if a:
        reject(label+': omit source scalar Hessian term',red-original-dS)
        reject(label+': omit bridge scalar Hessian term',red+dB-original-dC)
        reject(label+': reverse bridge scalar sign',red+dB-tadB-original-dC)
    else:
        check(label+': stationary scalar Hessian vanishes',tadS)
        check(label+': stationary bridge needs only fermion Jacobian',red+dB-original-dC)
    reject(label+': flavor products do not commute',Y*F.conjugate().T-F*Y.conjugate().T)
    return {'label':label,'active':n,'sterile':ns,'u':str(u),'m_H_squared':str(m2),'lambda':str(lam),'energy':str(E),'tadpole_factor':str(a),'grassmann_generators':len(momenta),'original_monomials':len(original.d),'reduced_monomials':len(red.d),'scalar_hessian_monomials':len(tadS.d)}

Y=s.Matrix([[1+I,2-I],[3,I]])
Cseed=s.Matrix([[2,1+I],[1+I,3-I]])
M=s.Matrix([[2,1],[1,3]])
trials=[]
for E in [R(2,3),R(5,2)]:
    trials.append(run_trial('offshell E='+str(E),R(3,2),R(7,3),R(2,5),E,Y,Cseed,M))
trials.append(run_trial('stationary',R(3,2),-R(9,5),R(2,5),R(4,3),Y,Cseed,M))
print(json.dumps({'status':'pass','python':sys.version.split()[0],'sympy':s.__version__,
    'scope':'Actual neutral Grassmann action and source-to-normalizing composition, including offshell scalar-potential Hessian; classical first-order selected packet only.',
    'exact_checks':len(checks),'rejected_mutations':len(mutations),'trials':trials,'checks':checks,'mutations':mutations},indent=2))
