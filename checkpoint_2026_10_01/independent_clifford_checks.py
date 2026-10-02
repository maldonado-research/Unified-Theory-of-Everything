#!/usr/bin/env python3
"""Independent exact SymPy curvature and physical-null-direction controls."""
import json
import sympy as s
I=s.I
z=s.zeros(2); e=s.eye(2)
pauli=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-I],[I,0]]),s.diag(1,-1)]
gamma=[s.diag(1,1,-1,-1)]+[z.row_join(p).col_join((-p).row_join(z)) for p in pauli]
g5=I*gamma[0]*gamma[1]*gamma[2]*gamma[3]
pl=(s.eye(4)-g5)/2
pr=(s.eye(4)+g5)/2
metric=[1,-1,-1,-1]
def slash(p): return sum((metric[k]*p[k]*gamma[k] for k in range(4)),s.zeros(4))
def dot(p,q): return sum(metric[k]*p[k]*q[k] for k in range(4))
def adj(a): return gamma[0]*a.conjugate().T*gamma[0]
def check(cond,label):
 if not cond: raise RuntimeError(label)
checks=[]
for mu in range(4):
 for nu in range(4):
  check(gamma[mu]*gamma[nu]+gamma[nu]*gamma[mu]==2*(metric[mu] if mu==nu else 0)*s.eye(4),'Clifford')
checks.append('Clifford algebra')
# A constant noncommuting gauge connection. D_mu=partial_mu+i Gamma_mu;
# F_mu_nu=i[Gamma_mu,Gamma_nu]. Tensor axes = gauge, spin.
connections=[pauli[0],pauli[1],s.zeros(2),pauli[2]/2]
b=sum((s.kronecker_product(connections[k],gamma[k]) for k in range(4)),s.zeros(8))
gsq=sum((metric[k]*s.kronecker_product(connections[k]**2,s.eye(4)) for k in range(4)),s.zeros(8))
curv=s.zeros(8)
for mu in range(4):
 for nu in range(4):
  f=I*(connections[mu]*connections[nu]-connections[nu]*connections[mu])
  sigma=I*(gamma[mu]*gamma[nu]-gamma[nu]*gamma[mu])/2
  curv+=s.kronecker_product(f,sigma)/2
check(b*b+curv==gsq,'row curvature sign')
check(b*b-curv!=gsq,'wrong row curvature mutation rejected')
checks.append('nonabelian constant-connection row curvature')
# On-shell gauge-null-direction amplitude under paper's convention.
p=s.Matrix([2,0,0,0]); r=s.Matrix([1,0,0,1]); k=s.Matrix([s.Rational(1,2),0,0,-s.Rational(1,2)]); h=k.copy(); eps=s.Matrix([0,1,0,0])
check(p==r+h+k,'conservation')
check(dot(p,p)==4 and dot(r,r)==0 and dot(h,h)==0 and dot(k,k)==0,'shells')
check(dot(eps,k)==0,'polarization')
def kernel(pol): return slash(k)*dot(p,pol)-slash(pol)*dot(p,k)
check(kernel(k)==s.zeros(4),'Ward test')
v=kernel(eps)*pl
spinsum=s.simplify(s.trace((pl*slash(r))*v*(slash(p)+2*s.eye(4))*adj(v)))
check(spinsum==4,'spin sum 4')
# A nonzero amplitude on four free shell legs cannot be a pure free EOM or
# total derivative term. This alone says nothing about loop power counting.
checks.append('gauge-null-direction shell/ward/spin sum')
# Verify CP-complex flavor orientation with explicitly nonreal matrices.
y=s.Matrix([[1+I,2-I],[I,3],[2+3*I,-1]])
f=s.Matrix([[2-I,I],[1+I,4-2*I],[-I,3]])
m=s.diag(2,5); dy=-f*m**3
dc=f*m**2*y.T+y*m**2*f.T
schur=(dc+dy*m.inv()*y.T+y*m.inv()*dy.T).applyfunc(s.expand)
check(schur==s.zeros(3),'transpose Schur')
wrong=(dc+dy*m.inv()*y.T+y*m.inv()*dy.conjugate().T).applyfunc(s.expand)
check(wrong!=s.zeros(3),'dagger cannot substitute transpose')
checks.append('explicit complex flavor transpose orientation')
print(json.dumps({'status':'PASS','checks':checks,'spin_sum':str(spinsum),'scope':'independent exact algebra controls; no beta function or novelty claim'},indent=2))
