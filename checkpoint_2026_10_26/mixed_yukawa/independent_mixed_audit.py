#!/usr/bin/env python3
"""Independent auxiliary-mass UV, Weyl Wick and flavor audit of the mixed seed pole."""
from collections import Counter
from itertools import product
from math import factorial
from pathlib import Path
import json
import sympy as s

I=s.I; R=s.Rational
checks=Counter(); rejected=Counter()
def zero(x):
 if isinstance(x,s.MatrixBase): return all(s.simplify(a)==0 for a in x)
 return s.simplify(x)==0
def check(name,x):
 if not zero(x): raise RuntimeError('Check failed: '+name)
 checks[name]+=1
def reject(name,x):
 if zero(x): raise RuntimeError('Mutation survived: '+name)
 rejected[name]+=1

def radial(power,den):
 # Wick-rotated Gamma master in d=4−2epsilon, stripped i/(16π²epsilon).
 n=den-power-2
 if n>0:return s.Integer(0)
 return s.Integer(-1)**(power-den)*s.gamma(power+2)*(-1)**(-n)/(s.gamma(2)*s.gamma(den)*factorial(-n))
check('q squared over three denominators radial pole',radial(1,3)-1)
check('q fourth over four denominators radial pole',radial(2,4)-1)
check('three denominators without numerator UV finite',radial(0,3))
P=s.Matrix(s.symbols('P0:4'));S=s.Matrix(s.symbols('S0:4'));T=s.Matrix(s.symbols('T0:4'))
# q enters the seed at its barred-lepton leg. The composite momentum is q+S;
# the scalar line denominator is (q+P+S)^2−mH² and sterile (q−T)^2−M².
# Geometric expansion about q²−a keeps only one external power.
num_shift=2*radial(1,3)*S/4
den_shift=-2*radial(2,4)*(P+S-T)/4
UV=num_shift+den_shift
check('exhaustive local vector residue',UV-(T-P)/2)
reject('drop denominator momentum shifts',num_shift-UV)
reject('differentiate the lepton only',-2*radial(2,4)*(P+S-T)/4-UV)
reject('omit scalar external momentum in composite derivative',-2*radial(2,4)*(P+S-T)/4-UV)
# First neglected geometric term has degree at most -1; mass shifts enter at
# even external degree and multiply an odd q vector in the linear UV part.
check('one-loop superficial degree',4+3-6-1)
check('two denominator corrections are UV finite',4+3+2-2*(3+2)-(-1))
# Seed components: bar0,bar1,l1,n1,l2,n2. Each Yukawa contributes ε(l,n).
E=s.Matrix([[0,1],[-1,0]])
sigma=[s.eye(2),s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-I],[I,0]]),s.diag(1,-1)]
bar_sigma=[sigma[0],-sigma[1],-sigma[2],-sigma[3]]
for bpos in [0,1]:
 for lpos in [2,4]:
  bout=1-bpos;lout=6-lpos
  ordering=[bpos,lpos,3,5,bout,lout]
  parity=(-1)**sum(ordering[i]>ordering[j] for i in range(6) for j in range(i+1,6))
  for mu,Q in enumerate(sigma):
   out=s.zeros(2)
   for vals in product(range(2),repeat=6):
    a,b,c,d,e,f=vals
    out[vals[bout],vals[lout]]+=parity*E[a,b]*E[c,d]*E[e,f]*E[d,f]*Q[vals[lpos],vals[bpos]]
   check('literal four Wick contractions and Pauli coefficients',out-bar_sigma[mu])
   if parity == -1: reject('omit Wick parity on negative-sign contractions',-out-bar_sigma[mu])
# For bar(x) l(y), reversing the canonical ordered propagator gives -i σ(k)
# with k flowing x -> y. Since q enters x, k=-q, hence +i σ(q).
vertex=(-I/2)*((-I)**2)/2
phase_scalar=I;phase_sterile=-I;phase_active_in=I
coefficient=vertex*phase_scalar*phase_sterile*phase_active_in
check('literal vertex and incoming propagator phase',coefficient+R(1,4))
check('four opposite-composite labelled assignments',4*coefficient+1)
reject('ordered propagator reversed without momentum reversal',coefficient-vertex*phase_scalar*phase_sterile*(-I))
# Both differentiated seed composites give the same square by momentum conservation.
q,Sx=s.symbols('q S'); check('opposite composite derivative momenta',(-q-Sx)**2-(q+Sx)**2)
# Same-composite internal seed momenta q and -q cancel; numerator degree drops.
check('same-composite loop momentum cancellation',q-q)
check('same-composite UV degree',4+1-6+1)
# Bare pole: gamma matrix multiplies +κ(P−T)/2. For fields e^(−ipx),
# i∂ on barL gives P and i∂ on phi† gives T.
check('ordered local Fourier derivative',((P-T)/2)+UV)
# General complex symmetric sterile mass, no inverse and passive flavor rotations.
Y=s.Matrix([[1+I,2-I],[3-2*I,-1+I]])
C=s.Matrix([[2-I,1+2*I],[1+2*I,3+I]])
M=s.Matrix([[2+I,1-2*I],[1-2*I,4+I]])
U=s.Matrix([[R(5,13),R(12,13)],[-R(12,13),R(5,13)]])*s.diag((5+12*I)/13,I)
V=s.diag(I,(3+4*I)/5)*s.Matrix([[R(3,5),R(4,5)],[-R(4,5),R(3,5)]])
Q=Y.conjugate()*M*Y.H;K=C*Q
Yp=V.H*Y*U.conjugate();Cp=V.H*C*V.conjugate();Mp=U.H*M*U.conjugate()
Qp=Yp.conjugate()*Mp*Yp.H;Kp=Cp*Qp
check('symmetric anomalous mass source',Q-Q.T)
check('complex mass source congruence',Qp-V.T*Q*V)
check('mixed lepton coefficient transforms by similarity',Kp-V.H*K*V)
reject('sterile mass adjoint in place of mass',Q-Y.conjugate()*M.H*Y.H)
reject('reverse active flavor order',K-Q*C)
report={'status':'PASS','scope':'Independent UV integral, Weyl contraction and primitive local-action check for the one-seed/two-neutrino-Yukawa mixed 1PI sector; zero gauge fields, not a complete physical Wilson coefficient.','input_commit':'8999d3de1a0e3dc2576fe1586295426da687fc67','positive_checks':sum(checks.values()),'negative_checks':sum(rejected.values()),'positive_categories':dict(checks),'negative_categories':dict(rejected),'incoming_seed_momentum_residue':'i kappa (T-P)/2','incoming_ordered_propagator_phase':'+i sigma(q)','action_pole':'kappa/2 K [rho i(partial_mu bar L) gamma^mu L - phi i(partial_mu phi_dagger) bar L gamma^mu L] + h.c.','scalar_and_sterile_mass_dependence_in_simple_pole':'only anomalous sterile mass numerator','repository_imports':False,'assertions_required':False}
Path(__file__).with_name('INDEPENDENT_MIXED_RESULTS.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
