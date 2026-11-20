#!/usr/bin/env python3
"""Off-shell Fourier/Grassmann certificate for a mixed-source identity.

This is an action-coordinate calculation. It does not compute a loop pole.
"""
from fractions import Fraction as F
from collections import Counter
import json
from pathlib import Path

pos=Counter(); neg=Counter(); I=1j
class P:
    def __init__(self,d=None): self.d={k:v for k,v in (d or {}).items() if v != 0}
    @staticmethod
    def cast(v): return v if isinstance(v,P) else P({(0,()):v})
    def __add__(self,v):
        out=dict(self.d)
        for k,x in P.cast(v).d.items(): out[k]=out.get(k,0)+x
        return P(out)
    __radd__=__add__
    def __neg__(self): return P({k:-x for k,x in self.d.items()})
    def __sub__(self,v): return self + -P.cast(v)
    def __rsub__(self,v): return P.cast(v)+-self
    def __mul__(self,v):
        out={}
        for (p,a),x in self.d.items():
            for (q,b),y in P.cast(v).d.items():
                if set(a)&set(b):continue
                sign=(-1)**sum(i>j for i in a for j in b)
                key=(p+q,tuple(sorted(a+b)));out[key]=out.get(key,0)+sign*x*y
        return P(out)
    __rmul__=__mul__
    def D(self):return P({(p,g):-I*p*x for (p,g),x in self.d.items()})
    def integral(self):return P({(p,g):x for (p,g),x in self.d.items() if p==0})
    def dagger(self):
        out=P()
        for (p,g),x in self.d.items():
            term=P({(-p,()):x.conjugate()})
            for j in reversed(g):
                jj=(j+4 if j<4 else j-4 if j<8 else j+4 if j<12 else j-4)
                term=term*P({(0,(jj,)):1})
            out+=term
        return out

def check(label,v):
    if P.cast(v).d:raise RuntimeError(label+': identity failed')
    pos[label]+=1

def reject(label,v):
    if not P.cast(v).d:raise RuntimeError(label+': incorrect formula survived')
    neg[label]+=1

def main():
    k=[[1+I,2-I],[-1+2*I,3-I]]
    fields=[[P({(1+2*a+u,(2*a+u,)):1,(-2-a-u,(2*a+u,)):2-I}) for u in range(2)] for a in range(2)]
    bar=[[v.dagger() for v in row] for row in fields]
    src=[[P({(2+3*a+u,(8+2*a+u,)):2+I,(-1-a-2*u,(8+2*a+u,)):1-I}) for u in range(2)] for a in range(2)]
    barsrc=[[v.dagger() for v in row] for row in src]
    phi=[P({(1+u,()):1+I,(-3-u,()):2-I}) for u in range(2)]
    phibar=[v.dagger() for v in phi]
    rho=sum((phibar[u]*phi[u] for u in range(2)),P())
    A=sum((phibar[u].D()*phi[u] for u in range(2)),P())
    B=sum((phibar[u]*phi[u].D() for u in range(2)),P())
    W=A-B
    check('scalar product rule',rho.D()-A-B)
    total=P();reduced=P();residual=P();boundary=P()
    for a in range(2):
        for b in range(2):
            J=sum((bar[a][u]*fields[b][u] for u in range(2)),P())
            rr=sum((bar[a][u].D()*fields[b][u] for u in range(2)),P())
            ss=sum((bar[a][u]*fields[b][u].D() for u in range(2)),P())
            bls=sum((barsrc[a][u]*fields[b][u] for u in range(2)),P())
            lbs=sum((bar[a][u]*src[b][u] for u in range(2)),P())
            ebarL=I*rr+bls;bare=I*ss-lbs
            check('lepton-current product rule',J.D()-rr-ss)
            raw=I*(A*J-rho*rr)
            red=I*F(1,2)*W*J+rho*(F(3,2)*bls-F(1,2)*lbs)
            rem=-rho*F(1,2)*(3*ebarL+bare)
            edge=I*F(1,2)*(rho*J).D()
            check('local off-shell source and boundary identity',raw-red-rem-edge)
            check('integrated off-shell source identity',(raw-red-rem).integral())
            total+=k[a][b]*raw;reduced+=k[a][b]*red;residual+=k[a][b]*rem;boundary+=k[a][b]*edge
    check('complex-flavor off-shell identity',total-reduced-residual-boundary)
    check('complete Hermitian-conjugate identity',total.dagger()-reduced.dagger()-residual.dagger()-boundary.dagger())
    check('barred equation is minus adjoint',sum(((I*bar[a][u].D()+barsrc[a][u])+(I*fields[a][u].D()-src[a][u]).dagger() for a in range(2) for u in range(2)),P()))
    # Variation of a first-order action: delta S=delta bar L e - ebar delta L.
    deltal=[[sum((rho*F(1,2)*(3*k[a][b]-k[b][a].conjugate())*fields[b][u] for b in range(2)),P()) for u in range(2)] for a in range(2)]
    deltabar=[[v.dagger() for v in row] for row in deltal]
    variation=P()
    for a in range(2):
        for u in range(2):
            e=I*fields[a][u].D()-src[a][u]
            ebar=I*bar[a][u].D()+barsrc[a][u]
            variation+=deltabar[a][u]*e-ebar*deltal[a][u]
    check('conjugate-safe field variation gives both residuals',variation-residual-residual.dagger())
    reject('omit all interaction sources',(total-(I*F(1,2)*sum((k[a][b]*W*sum((bar[a][u]*fields[b][u] for u in range(2)),P()) for a in range(2) for b in range(2)),P()))-residual-boundary).integral())
    reject('wrong row-derivative coefficient',(total-reduced-residual*F(2,3)-boundary).integral())
    reject('discard scalar-current term',(total-reduced+I*F(1,2)*sum((k[a][b]*W*sum((bar[a][u]*fields[b][u] for u in range(2)),P()) for a in range(2) for b in range(2)),P())-residual-boundary).integral())
    reject('use independent nonconjugate fermion changes',(variation-residual).integral())
    out={'status':'PASS','date':'26 October 2026 Pacific','positive_checks':sum(pos.values()),'rejected_mutations':sum(neg.values()),'positive_categories':dict(pos),'negative_categories':dict(neg),'scope':'Formal coefficient-independent off-shell source and boundary identity with a conjugate-safe first-order fermion field change; not a loop integration or complete EFT reduction.','no_project_imports':True,'assertions_required':False}
    Path(__file__).with_name('MIXED_SOURCE_IDENTITY_RESULTS.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
