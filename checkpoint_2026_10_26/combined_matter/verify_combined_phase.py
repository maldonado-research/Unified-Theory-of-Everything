#!/usr/bin/env python3
"""Correlated scalar/lepton phase transport; no new loop integration."""
from collections import Counter
from pathlib import Path
import json
import sympy as s
I=s.I; R=s.Rational
pos=Counter();neg=Counter()
def z(x):
 if isinstance(x,G):return not x.d
 if isinstance(x,s.MatrixBase):return all(s.simplify(v)==0 for v in x)
 return s.simplify(x)==0
def check(n,x):
 if not z(x):raise RuntimeError(n+' failed')
 pos[n]+=1
def reject(n,x):
 if z(x):raise RuntimeError(n+' mutation escaped')
 neg[n]+=1
class G:
 def __init__(self,d=None):self.d={k:s.expand(v) for k,v in (d or {}).items() if s.expand(v)!=0}
 @staticmethod
 def cast(v):return v if isinstance(v,G) else G({():s.sympify(v)})
 def __add__(self,v):
  d=dict(self.d)
  for k,x in G.cast(v).d.items():d[k]=d.get(k,0)+x
  return G(d)
 __radd__=__add__
 def __neg__(self):return G({k:-v for k,v in self.d.items()})
 def __sub__(self,v):return self+-G.cast(v)
 def __mul__(self,v):
  out={}
  for a,x in self.d.items():
   for b,y in G.cast(v).d.items():
    if set(a)&set(b):continue
    k=tuple(sorted(a+b));out[k]=out.get(k,0)+(-1)**sum(i>j for i in a for j in b)*x*y
  return G(out)
 __rmul__=__mul__
def main():
 alpha,rho,dr=s.symbols('alpha rho dr',real=True)
 A=s.Matrix([[2,1+2*I],[1-2*I,-1]]);t=s.trace(A);B=A-t*s.eye(2)
 Yn=s.Matrix([[1+I,2],[3-I,-1]]);Ye=s.Matrix([[2,1-I],[-1+2*I,3]])
 C5=s.Matrix([[1+I,2-I],[2-I,3+2*I]])
 check('Hermitian phase generator',A-A.H)
 # Arbitrary jets, including unrelated complex entries, certify product rule.
 l=s.Matrix(s.symbols('l0:2'));dl=s.Matrix(s.symbols('dl0:2'));bar=s.Matrix(1,2,s.symbols('b0:2'));dbar=s.Matrix(1,2,s.symbols('db0:2'))
 delta_l=I*alpha*rho*A*l; delta_bar=-I*alpha*rho*bar*A
 ddelta_l=I*alpha*(dr*A*l+rho*A*dl)
 kinetic=(delta_bar*I*dl+bar*I*ddelta_l)[0]
 check('lepton kinetic variation cancels anti-Hermitian current',kinetic+alpha*dr*(bar*A*l)[0])
 reject('omit derivative of density',kinetic)
 x=s.symbols('x0:2',real=True);y=s.symbols('y0:2',real=True);dx=s.symbols('dx0:2',real=True);dy=s.symbols('dy0:2',real=True)
 p=s.Matrix([xj+I*yj for xj,yj in zip(x,y)]);dp=s.Matrix([xj+I*yj for xj,yj in zip(dx,dy)])
 rh=(p.H*p)[0];drh=(dp.H*p+p.H*dp)[0];J=I*((p.H*dp)[0]-(dp.H*p)[0])
 delta=I*alpha*t*rh*p; dd=delta.applyfunc(lambda f:sum(s.diff(f,q)*v for q,v in zip(list(x)+list(y),list(dx)+list(dy))))
 check('scalar kinetic variation cancels trace current',(dd.H*dp+dp.H*dd)[0]+alpha*t*drh*J)
 check('scalar norm invariant',(delta.H*p+p.H*delta)[0])
 weak_eps=s.Matrix([[0,1],[-1,0]])
 Higgs=-weak_eps*p.conjugate(); dHiggs=-weak_eps*delta.conjugate()
 check('Higgs doublet inversion',weak_eps*Higgs.conjugate()-p)
 check('conjugate Higgs phase from actual inversion',dHiggs+I*alpha*t*rh*Higgs)
 j=s.Matrix(s.symbols('j0:2'))
 scalar_probe=(j.H*delta+delta.H*j)[0]
 check('scalar probe source transport',scalar_probe-I*alpha*t*rh*((j.H*p)[0]-(p.H*j)[0]))
 reject('drop scalar probe transport',scalar_probe)
 # Flavor-row source transport with independent external currents.
 n=s.Matrix(s.symbols('n0:2'));e=s.Matrix(s.symbols('e0:2'));u=s.Matrix(s.symbols('u0:2'));d=s.Matrix(s.symbols('d0:2'));q=s.Matrix(1,2,s.symbols('q0:2'))
 Yu=s.Matrix([[1,2+I],[I,-1]]);Yd=s.Matrix([[3-I,1],[2,1+I]])
 Yn0=(bar*Yn*n)[0];Ye0=(bar*Ye*e)[0];Yu0=(q*Yu*u)[0];Yd0=(q*Yd*d)[0]
 # Treat the scalar and conjugate scalar as independent factors; phi -> +it phi, H -> -it H.
 dn=-(delta_bar*Yn*n)[0]-I*alpha*rho*t*Yn0
 de=-(delta_bar*Ye*e)[0]+I*alpha*rho*t*Ye0
 du=-I*alpha*rho*t*Yu0;ddown=I*alpha*rho*t*Yd0
 check('neutrino Yukawa source variation',dn-I*alpha*rho*(bar*B*Yn*n)[0])
 check('charged lepton Yukawa source variation',de-I*alpha*rho*(bar*(A+t*s.eye(2))*Ye*e)[0])
 check('up quark Yukawa source variation',du+I*alpha*rho*t*Yu0)
 check('down quark Yukawa source variation',ddown-I*alpha*rho*t*Yd0)
 reject('drop charged-lepton spectator',de)
 reject('drop quark spectator',du+ddown)
 # Independent Grassmann spinor and weak components reproduce the full dimension-five transport.
 eps=s.Matrix([[0,1],[-1,0]]); phi=[1+I,2-I]; density=sum(s.conjugate(v)*v for v in phi)
 bl=[[[G({(4*a+2*w+r,):1}) for r in range(2)] for w in range(2)] for a in range(2)]
 dbl=[[[-I*alpha*density*sum((bl[b][w][r]*A[b,a] for b in range(2)),G()) for r in range(2)] for w in range(2)] for a in range(2)]
 dphi=[I*alpha*density*t*v for v in phi]
 xb=[[sum((bl[a][w][r]*phi[w] for w in range(2)),G()) for r in range(2)] for a in range(2)]
 dxb=[[sum((dbl[a][w][r]*phi[w]+bl[a][w][r]*dphi[w] for w in range(2)),G()) for r in range(2)] for a in range(2)]
 xc=[[sum((-eps[r,j]*xb[a][j] for j in range(2)),G()) for r in range(2)] for a in range(2)]
 dxc=[[sum((-eps[r,j]*dxb[a][j] for j in range(2)),G()) for r in range(2)] for a in range(2)]
 def op(c,rows=xb,cols=xc):return sum((R(1,2)*c[a,b]*rows[a][r]*cols[b][r] for a in range(2) for b in range(2) for r in range(2)),G())
 direct=op(C5,dxb,xc)+op(C5,xb,dxc)
 induced=B*C5+C5*B.T
 check('literal dimension-five action variation',direct+I*alpha*density*op(induced))
 reject('omit dimension-five times dimension-six field-change contact',direct)
 reject('replace transpose by adjoint on charge-conjugate slot',direct+I*alpha*density*op(B*C5+C5*B.H))
 # The cancellation criterion holds precisely on the image of Y; test valid seeds.
 Y1=s.Matrix([[1],[0]]); M=s.Matrix([[1]]);C=s.diag(I,2);K=C*Y1.conjugate()*M*Y1.H;A1=(K-K.H)/(2*I);t1=s.trace(A1)
 check('rank-one compatible seed has nonzero generator',A1-s.diag(1,0))
 check('rank-one Yukawa support satisfies phase cancellation',(A1-t1*s.eye(2))*Y1)
 reject('rank-one cancellation implies zero generator',A1)
 Yfull=s.eye(2); Cf=I*s.eye(2);Kf=Cf*Yfull.conjugate()*s.eye(2)*Yfull.H;Af=(Kf-Kf.H)/(2*I);tf=s.trace(Af)
 reject('multi-flavor full-rank Yukawa phase cancels automatically',(Af-tf*s.eye(2))*Yfull)
 one=s.Matrix([[1]]);C1=s.Matrix([[I]]);K1=C1*one*one*one;Aone=(K1-K1.H)/(2*I)
 check('one-flavor neutrino and neutral composite are invariant',(Aone-s.trace(Aone)*s.eye(1))*one)
 # Probe action delta = i alpha rho (bar eta A L - bar L A eta).
 eta=s.Matrix(s.symbols('eta0:2'));bareta=s.Matrix(1,2,s.symbols('be0:2'))
 probe=(bareta*delta_l+delta_bar*eta)[0]
 check('lepton probe source transport',probe-I*alpha*rho*((bareta*A*l)[0]-(bar*A*eta)[0]))
 reject('drop lepton probe transport',probe)
 report={'status':'PASS','date':'26 October 2026 Pacific','scope':'Classical first-order correlated phase transport using previous pole inputs; no loop integration, complete physical basis or observable.','positive_checks':sum(pos.values()),'rejected_mutations':sum(neg.values()),'positive_categories':dict(pos),'negative_categories':dict(neg),'no_project_imports':True,'assertions_required':False}
 Path(__file__).with_name('COMBINED_PHASE_RESULTS.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps(report,indent=2))
if __name__=='__main__':main()
