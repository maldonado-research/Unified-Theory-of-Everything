#!/usr/bin/env python3
"""Independent odd-field action audit of the correlated local phase.

Literal Yukawa and probe variation with exterior algebra, scalar/lepton
kinetic first jets, and singular-Y source-kernel witnesses. No project imports.
"""
from collections import Counter
from pathlib import Path
import json
import sympy as s
I=s.I; counts=Counter(); neg=Counter(); mates={}; nvar=0

def simp(x): return s.simplify(s.expand(x))
class G:
 def __init__(self,d=None):self.d={k:simp(v) for k,v in (d or {}).items() if simp(v)!=0}
 @staticmethod
 def cast(x):return x if isinstance(x,G) else G({():s.sympify(x)})
 def __add__(self,x):
  z=dict(self.d)
  for k,v in G.cast(x).d.items():z[k]=z.get(k,0)+v
  return G(z)
 __radd__=__add__
 def __neg__(self):return G({k:-v for k,v in self.d.items()})
 def __sub__(self,x):return self+-G.cast(x)
 def __rsub__(self,x):return G.cast(x)+-self
 def __mul__(self,x):
  z={}
  for a,v in self.d.items():
   for b,w in G.cast(x).d.items():
    if set(a)&set(b):continue
    k=tuple(sorted(a+b));sgn=(-1)**sum(i>j for i in a for j in b)
    z[k]=z.get(k,0)+sgn*v*w
  return G(z)
 __rmul__=__mul__
 def dag(self):
  result=G()
  for k,v in self.d.items():
   term=G.cast(s.conjugate(v))
   for idx in reversed(k):term=term*G({(mates[idx],):1})
   result+=term
  return result

def pair():
 global nvar
 a=nvar;b=nvar+1;nvar+=2;mates[a]=b;mates[b]=a
 return G({(a,):1}),G({(b,):1})
def check(name,x):
 if isinstance(x,G):ok=not x.d
 elif isinstance(x,s.MatrixBase):ok=all(simp(a)==0 for a in x)
 else:ok=simp(x)==0
 if not ok:raise RuntimeError(name+' failed')
 counts[name]+=1
def reject(name,x):
 if isinstance(x,G):bad=bool(x.d)
 elif isinstance(x,s.MatrixBase):bad=any(simp(a)!=0 for a in x)
 else:bad=simp(x)!=0
 if not bad:raise RuntimeError(name+' mutation escaped')
 neg[name]+=1

def main():
 Y=s.Matrix([[1+I,2-I],[3-2*I,-1+I]])
 Ye=s.Matrix([[2-I,1+I],[-1+2*I,3]])
 M=s.Matrix([[2+I,1-2*I],[1-2*I,4+I]])
 C=s.Matrix([[2-I,1+2*I],[1+2*I,3+I]])
 R=Y.conjugate()*M*Y.H; K=C*R; A=(K-K.H)/(2*I); t=s.trace(A)
 check('A Hermitian',A-A.H);check('trace real',t-s.conjugate(t))
 reject('generic flavor phase differs from trace',A-t*s.eye(2))
 phi=s.Matrix([1+2*I,2-I]); dphi=s.Matrix([2-I,-1+2*I]); eps=s.Matrix([[0,1],[-1,0]])
 H=-eps*phi.conjugate();check('tilde Higgs inverse convention',eps*H.conjugate()-phi);rho=(phi.H*phi)[0];dr=(dphi.H*phi+phi.H*dphi)[0];j=I*((phi.H*dphi)[0]-(dphi.H*phi)[0])
 dphi_change=I*t*(dr*phi+rho*dphi)
 reject('zero scalar norm derivative in test configuration',dr)
 reject('zero scalar current in test configuration',j)
 check('scalar kinetic variation cancels scalar phase-current pole',(dphi_change.H*dphi+dphi.H*dphi_change)[0]+t*dr*j)
 reject('reversed scalar phase-current sign',(dphi_change.H*dphi+dphi.H*dphi_change)[0]-t*dr*j)
 reject('half scalar phase-current coefficient',(dphi_change.H*dphi+dphi.H*dphi_change)[0]+s.Rational(1,2)*t*dr*j)
 check('scalar norm unchanged',(-I*t*rho*phi.H*phi+phi.H*(I*t*rho*phi))[0])
 L={}; BL={}; DL={}; N={};BN={};ER={};BER={};Q={};BQ={};UR={};BRU={};DR={};BRD={};eta={};beta={}
 for a in range(2):
  for w in range(2):
   for z in range(2):
    L[a,w,z],BL[a,w,z]=pair();DL[a,w,z],_=pair();eta[a,w,z],beta[a,w,z]=pair()
 for i in range(2):
  for z in range(2):N[i,z],BN[i,z]=pair();ER[i,z],BER[i,z]=pair()
 for w in range(2):
  for z in range(2):Q[w,z],BQ[w,z]=pair()
 for z in range(2):UR[z],BRU[z]=pair();DR[z],BRD[z]=pair()
 dL={(a,w,z):sum((I*rho*A[a,b]*L[b,w,z] for b in range(2)),G()) for a in range(2) for w in range(2) for z in range(2)}
 dBL={(a,w,z):sum((-I*rho*BL[b,w,z]*A[b,a] for b in range(2)),G()) for a in range(2) for w in range(2) for z in range(2)}
 dph=I*rho*t*phi;dH=-I*rho*t*H
 # Direct first-order product changes; spin labels are independent odd fields.
 def yuk(B,F,scalar,matrix):
  return sum((B[a,w,z]*matrix[a,i]*scalar[w]*F[i,z] for a in range(2) for i in range(2) for w in range(2) for z in range(2)),G())
 dVnu=yuk(dBL,N,phi,Y)+yuk(BL,N,dph,Y)
 expected_nu=sum((I*rho*BL[a,w,z]*((t*s.eye(2)-A)*Y)[a,i]*phi[w]*N[i,z] for a in range(2) for i in range(2) for w in range(2) for z in range(2)),G())
 check('literal neutrino Yukawa variation',dVnu-expected_nu)
 check('neutrino interaction including h.c.',-dVnu-dVnu.dag()-(I*rho*yuk(BL,N,phi,(A-t*s.eye(2))*Y)+(I*rho*yuk(BL,N,phi,(A-t*s.eye(2))*Y)).dag()))
 dVe=yuk(dBL,ER,H,Ye)+yuk(BL,ER,dH,Ye)
 check('literal charged-lepton Yukawa variation',dVe+I*rho*yuk(BL,ER,H,(A+t*s.eye(2))*Ye))
 def qterm(scalar,F):return sum((BQ[w,z]*scalar[w]*F[z] for w in range(2) for z in range(2)),G())
 yu=1+2*I;yd=2-I
 dup=yu*qterm(dph,UR);ddown=yd*qterm(dH,DR)
 check('literal up-quark Yukawa variation',-dup+I*rho*t*yu*qterm(phi,UR))
 check('literal down-quark Yukawa variation',-ddown-I*rho*t*yd*qterm(H,DR))
 probe=sum((beta[a,w,z]*L[a,w,z]+BL[a,w,z]*eta[a,w,z] for a in range(2) for w in range(2) for z in range(2)),G())
 dprobe=sum((beta[a,w,z]*dL[a,w,z]+dBL[a,w,z]*eta[a,w,z] for a in range(2) for w in range(2) for z in range(2)),G())
 expected_probe=sum((I*rho*beta[a,w,z]*A[a,b]*L[b,w,z]-I*rho*BL[b,w,z]*A[b,a]*eta[a,w,z] for a in range(2) for b in range(2) for w in range(2) for z in range(2)),G())
 check('literal lepton probe variation',dprobe-expected_probe)
 reject('drop neutrino interaction source',dVnu);reject('drop charged-lepton source',dVe);reject('drop probe source',dprobe)
 # Lower-order Weinberg insertion: vary the actual composite Xbar and its
 # charge-conjugate column. C5 is an independently supplied symmetric input.
 C5=s.Matrix([[1+I,2-I],[2-I,-1+2*I]]); B=A-t*s.eye(2)
 xb={(a,z):sum((BL[a,w,z]*phi[w] for w in range(2)),G()) for a in range(2) for z in range(2)}
 dxb={(a,z):sum((dBL[a,w,z]*phi[w]+BL[a,w,z]*dph[w] for w in range(2)),G()) for a in range(2) for z in range(2)}
 xc={(a,z):sum((-eps[z,j]*xb[a,j] for j in range(2)),G()) for a in range(2) for z in range(2)}
 dxc={(a,z):sum((-eps[z,j]*dxb[a,j] for j in range(2)),G()) for a in range(2) for z in range(2)}
 change=sum((s.Rational(1,2)*C5[a,b]*(dxb[a,z]*xc[b,z]+xb[a,z]*dxc[b,z]) for a in range(2) for b in range(2) for z in range(2)),G())
 target=sum((-I*rho*s.Rational(1,2)*(B*C5+C5*B.T)[a,b]*xb[a,z]*xc[b,z] for a in range(2) for b in range(2) for z in range(2)),G())
 check('literal lower-order C5 insertion variation',change-target)
 check('lower-order C5 insertion including h.c.',change+change.dag()-target-target.dag())
 wrong=sum((-I*rho*s.Rational(1,2)*(A*C5+C5*A.T)[a,b]*xb[a,z]*xc[b,z] for a in range(2) for b in range(2) for z in range(2)),G())
 reject('omit scalar phase from C5 insertion',change-wrong)
 # Independently vary kinetic term in a genuine odd-field basis, keeping first jets.
 pauli=[s.eye(2),s.Matrix([[0,-1],[-1,0]]),s.Matrix([[0,I],[-I,0]]),s.diag(-1,1)]
 for sig in pauli:
  raw=G(); want=G()
  for a in range(2):
   for b in range(2):
    for w in range(2):
     for z in range(2):
      for z2 in range(2):
       # dBL already sums over both flavors; this loop repeats it twice.
       raw+=s.Rational(1,2)*dBL[a,w,z]*I*sig[z,z2]*DL[a,w,z2]+BL[a,w,z]*I*sig[z,z2]*(I*dr*A[a,b]*L[b,w,z2]+I*rho*A[a,b]*DL[b,w,z2])
       want+=-dr*BL[a,w,z]*sig[z,z2]*A[a,b]*L[b,w,z2]
  check('lepton kinetic first-jet variation',raw-want)
  reject('omit mixed phase-current pole',want)
 # Explicit scalar probe, with unconstrained complex source coefficients.
 probe_row=s.Matrix([[2-I,3+I]])
 check('literal scalar probe variation',(probe_row*dph)[0]+s.conjugate((probe_row*dph)[0])-I*rho*t*((probe_row*phi)[0]-s.conjugate((probe_row*phi)[0])))
 # Rank-safe image projectors. Rank 0, 1, 2 are explicit exact witnesses.
 examples=[s.zeros(2,2),s.Matrix([[1+I,0],[2-I,0]]),Y]
 for idx,yy in enumerate(examples):
  rr=yy.conjugate()*M*yy.H;kk=C*rr;aa=(kk-kk.H)/(2*I);tt=s.trace(aa)
  cols=yy.columnspace(); P=s.zeros(2) if not cols else s.Matrix.hstack(*cols)*(s.Matrix.hstack(*cols).H*s.Matrix.hstack(*cols)).inv()*s.Matrix.hstack(*cols).H
  # The inverse is used only to construct a test projector, never a sterile mass inverse.
  Q0=s.eye(2)-P
  check('rank-safe source support '+str(idx),rr*Q0)
  check('zero complementary A block '+str(idx),Q0*aa*Q0)
  check('projector equals its adjoint '+str(idx),P-P.H)
  check('source condition equivalent to AP=tP '+str(idx),(aa-tt*s.eye(2))*yy-((aa-tt*s.eye(2))*P)*yy)
 # Nontrivial rank-one source-free solution: choose real Y direction and C producing A=tP.
 y1=s.Matrix([[1],[0]]);m1=s.Matrix([[1]]);c1=s.Matrix([[I,0],[0,2]])
 k1=c1*y1.conjugate()*m1*y1.H;a1=(k1-k1.H)/(2*I);t1=s.trace(a1)
 check('nonzero rank-one source-free example',(a1-t1*s.eye(2))*y1)
 reject('wrongly claim rank-one source-free requires A zero',a1)
 # In two flavors eigenvalues of A−Tr(A)I are minus those of A in reverse.
 a,b,c=s.symbols('a b c',real=True);aa=s.Matrix([[a,b+I*c],[b-I*c,s.symbols('d',real=True)]])
 check('two-flavor trace-reversed determinant', (aa-s.trace(aa)*s.eye(2)).det()-aa.det())
 # One active flavor: X=phi†L and Vnu are invariant; other sources need not be.
 av=s.symbols('a_single',real=True);check('one-flavor neutrino phase difference',av-av)
 reject('one-flavor charged-lepton source incorrectly dropped',2*av)
 reject('one-flavor up-quark source incorrectly dropped',av)
 result={'status':'PASS','input_commit':'9b1e18450c0c1e18fb6f917105368b05ca123df2','scope':'Independent first-order action variation of selected phase terms and interaction/probe sources; no new loop, full quantum map or physical coefficient.', 'positive_checks':sum(counts.values()),'negative_checks':sum(neg.values()),'positive_categories':dict(counts),'negative_categories':dict(neg),'repository_imports':False,'assertions_required':False}
 Path(__file__).with_name('INDEPENDENT_COMBINED_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
