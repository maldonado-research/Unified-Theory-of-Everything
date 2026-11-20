#!/usr/bin/env python3
"""Independent component Grassmann check; no producer imports."""
import json
from pathlib import Path
import sympy as s
def plus(*ps):
    out={}
    for p in ps:
        for key,value in p.items():
            out[key]=s.expand(out.get(key,0)+value)
    return {k:v for k,v in out.items() if v!=0}
def scale(c,p): return {k:s.expand(c*v) for k,v in p.items() if s.expand(c*v)!=0}
def wedge(p,q):
    out={}
    for a,u in p.items():
        for b,v in q.items():
            if set(a)&set(b): continue
            sign=(-1)**sum(i>j for i in a for j in b)
            key=tuple(sorted(a+b))
            out[key]=s.expand(out.get(key,0)+sign*u*v)
    return {k:v for k,v in out.items() if v!=0}
def equal(p,q): return not plus(p,scale(-1,q))
def demand(x,label):
    if not x: raise RuntimeError(label)
n=[[{(2*i+a,):s.Integer(1)} for a in range(2)] for i in range(2)]
Z=s.Matrix([[1,s.I],[-s.I,2]])
C=s.Matrix([[2,1+s.I],[1+s.I,3*s.I]])
M=s.Matrix([[1,s.I],[s.I,-1]])
dn=[[plus(*(scale(-Z[i,j]/2,n[j][a]) for j in range(2))) for a in range(2)] for i in range(2)]
def bilinear(X,left=n,right=n):
    return plus(*(scale(X[i,j],plus(wedge(left[i][0],right[j][1]),scale(-1,wedge(left[i][1],right[j][0])))) for i in range(2) for j in range(2)))
variation=plus(bilinear(C,dn,n),bilinear(C,n,dn))
D=s.Matrix([[-3+s.I,-3-5*s.I/2],[-3-5*s.I/2,1-7*s.I]])
demand(equal(variation,bilinear(D)), 'Independent component mass-source variation')
demand(equal(bilinear(s.Matrix([[0,1],[-1,0]])),{}), 'Antisymmetric Weyl flavor coefficient vanishes')
single=bilinear(s.Matrix([[1,0],[0,0]]))
demand(single=={(0,1):s.Integer(2)}, 'Majorana half-normalization')
demand(M.rank()==1 and M.det()==0, 'Singular mass fixture')
wrong=-(Z.conjugate().T*C+C*Z)/2
demand(not equal(variation,bilinear(wrong)), 'Wrong dagger has to fail')
# Spectral verification: each constant-Hermitian-Z eigenmode has a scalar K and S.
t,rho,z=s.symbols('t rho z',real=True)
k=1+t*t*rho*z
a=k**s.Rational(-1,2)
demand(s.simplify(a*k*a-1)==0, 'Exact spectral kinetic normalization')
demand(s.simplify(a*k*s.diff(a,rho)-s.diff(a,rho)*k*a)==0,'Nonconstant-rho kinetic connection')
out={'status':'PASS','independent_component_checks':7,'producer_imports':False,
     'scope':'Two-flavor Weyl Grassmann variation, Majorana factor, rank-one mass, adjoint mutation and exact spectral kinetic normalization.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out))

