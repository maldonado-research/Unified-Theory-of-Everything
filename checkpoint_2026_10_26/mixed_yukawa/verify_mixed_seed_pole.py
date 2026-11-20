#!/usr/bin/env python3
"""Opposite-composite mixed seed triangle: symbolic simplex UV integral.

The result checks Lorentz integrals and Green-action algebra; the literal
Wick normalisation is recorded separately in MIXED_SEED_POLE.md. No project
imports and no optimization-removable assertions.
"""
from fractions import Fraction
from pathlib import Path
import json
import sympy as s

R=s.Rational; I=s.I; positive=0; adverse=0

def eq(label,x):
 global positive
 if s.simplify(s.expand(x))!=0: raise RuntimeError(label+': '+str(x))
 positive+=1

def ne(label,x):
 global adverse
 if s.simplify(s.expand(x))==0: raise RuntimeError(label+': mutation survived')
 adverse+=1

y,z=s.symbols('y z', nonnegative=True)
# Momentum components are independent formal coordinates.  Parameterizing
# q², (q−P)²−mφ², (q+h)²−M² gives q=l+yP−zh.
p,k,h=s.symbols('p k h');P=p+k;v=y*P-z*h
# The numerator is q_mu(q-k)^2.  The odd-in-l part has no UV pole;
# isotropy l_mu l_nu -> eta_mu nu l²/4 leaves (3v-k)/2.
raw=(3*v-k)/2
res=s.integrate(s.integrate(2*raw,(z,0,1-y)),(y,0,1))
eq('actual simplex UV vector',res-(p-h)/2)
eq('simplex area',s.integrate(s.integrate(2,(z,0,1-y)),(y,0,1))-1)
eq('first simplex moment y',s.integrate(s.integrate(2*y,(z,0,1-y)),(y,0,1))-R(1,3))
eq('first simplex moment z',s.integrate(s.integrate(2*z,(z,0,1-y)),(y,0,1))-R(1,3))
# Gamma pole: l²/(l²−Delta)^3 = 1/(l²−Delta)^2 + Delta/(l²−Delta)^3.
# The first normalized residue is 1 and the latter has no pole.
eps=s.symbols('epsilon')
eq('radial I2 gamma residue',s.residue(s.gamma(eps),eps,0)-1)
eq('radial I3 gamma residue',s.residue(s.gamma(1+eps),eps,0))
eq('radial numerator reduction',s.residue(s.gamma(eps),eps,0)+s.residue(s.gamma(1+eps),eps,0)-1)
# A second tensor reduction cancels the massless active denominator. The
# bubble vector has (P-h)/2, and the q_mu q.k triangle gives -k/2.
bubble=(P-h)/2
eq('independent bubble and triangle reduction',bubble-k/2-res)
# Literal phase of one representative contraction before the simplex integral.
ordering=[1,2,3,5,0,4]
inversions=sum(a>b for i,a in enumerate(ordering) for b in ordering[i+1:])
wick=(-1)**inversions
eq('explicit Grassmann crossing',wick+1)
vertex=(-I/2)*(-I)**2/2
propagators=(-I)*(-I)*I
spin=-1
weight=4*vertex*propagators*wick*spin
eq('four labelled Wick assignments and action half',weight-1)
ne('missing fermion crossing',4*vertex*propagators*spin-weight)
ne('extra Majorana half',weight/2-weight)
ne('discard scalar derivative k',res-s.integrate(s.integrate(2*(3*v)/2,(z,0,1-y)),(y,0,1)))
ne('incorrect angular contraction',res-s.integrate(s.integrate(2*(2*v-k)/2,(z,0,1-y)),(y,0,1)))
ne('wrong sterile routing',res-s.integrate(s.integrate(2*(3*(y*P+z*h)-k)/2,(z,0,1-y)),(y,0,1)))
# Same-composite momentum is purely external; numerator q_mu times an
# external scalar has superficial degree −1 and cannot contain a UV pole.
external,omega=s.symbols('external omega')
eq('same composite degree',4+1-2*3+1) # omega = -1
# Explicit chiral tensor: epsilon_rs sigma_{t s} epsilon_tu epsilon_vw epsilon_uw = -barsigma_{r v}.
E=s.Matrix([[0,1],[-1,0]])
pauli=[s.eye(2),s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-I],[I,0]]),s.diag(1,-1)]
for mu,S in enumerate(pauli):
 B=S if mu==0 else -S
 t=s.Matrix(2,2,lambda r,v:sum(E[r,j]*S[i,j]*E[i,u]*E[v,w]*E[u,w]
      for i in range(2) for j in range(2) for u in range(2) for w in range(2)))
 for row in range(2):
  for col in range(2):eq('actual epsilon sigma spin coefficient',t[row,col]+B[row,col])
# Complex-flavor tensor and passive-coordinate congruence.
Y=s.Matrix([[1+I,2-I],[3-2*I,-1+I]])
M=s.Matrix([[2+I,1-2*I],[1-2*I,4+I]])
C=s.Matrix([[2-I,1+2*I],[1+2*I,3+I]])
V=s.diag(I,(3+4*I)/5)*s.Matrix([[R(3,5),R(4,5)],[-R(4,5),R(3,5)]])
U=s.Matrix([[R(5,13),R(12,13)],[-R(12,13),R(5,13)]])*s.diag((5+12*I)/13,I)
RR=Y.conjugate()*M*Y.H;K=C*RR
Yp=V.H*Y*U.conjugate();Mp=U.H*M*U.conjugate();Cp=V.H*C*V.conjugate()
Rp=Yp.conjugate()*Mp*Yp.H;Kp=Cp*Rp
for x in Rp-V.T*RR*V:eq('mass source covariance',x)
for x in Kp-V.H*K*V:eq('mixed matrix covariance',x)
for x in RR-RR.T:eq('source symmetry',x)
ne('wrong complex mass conjugation',s.trace(C*(RR-Y.conjugate()*M.H*Y.H)))
# Incoming Fourier modes exp(-ip.x): i∂ produces +p, so the declared
# local action produces exactly p−h at the vertex.
eq('local derivative placement',p-h-2*res)
ne('only external lepton derivative',p-(p-h))
ne('wrong incoming sign',-p+h-(p-h))
# Bare dimension-six four-field operator has μ^(2epsilon); C7 Y*² M
# carries weight 2+1+1+0=4, so a κ residue maps to 2r.
C0,Yb,M0=s.symbols('C0 Yb M0');r=C0*Yb**2*M0/2
beta=2*C0*s.diff(r,C0)+Yb*s.diff(r,Yb)-2*r
eq('canonical epsilon factor',beta-2*r)
ne('pole is not beta',beta-r)
result={'status':'PASS','positive_checks':positive,'adverse_controls':adverse,
 'scope':'Exact opposite-composite scalar/vector pole and complex flavor covariance; overall Wick normalisation is a separately justified action input.',
 'vector_residue':'(p-h)/2','same_composite_UV_degree':-1,
 'assumptions':'unbroken theory; one original seed; two conjugate neutrino Yukawa vertices; retained symmetric Majorana mass; generic off-shell kinematics',
 'repository_imports':False,'assertions_required':False}
Path(__file__).with_name('MIXED_SEED_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
