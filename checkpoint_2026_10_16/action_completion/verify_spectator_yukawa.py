#!/usr/bin/env python3
"""Exact spectator-Yukawa transport under the two declared field maps.

Standalone exterior-polynomial substitution, with actual scalar/lepton shifts.
No project imports and no optimization-removable assertions.
"""
from collections import Counter
from pathlib import Path
import json
import sympy as s

I=s.I
positive=Counter()
negative=Counter()
conjugates={}


class E:
    def __init__(self, data=None):
        self.d={k:s.expand(v) for k,v in (data or {}).items() if v!=0}
        self.d={k:v for k,v in self.d.items() if v!=0}
    @staticmethod
    def scalar(x):
        return x if isinstance(x,E) else E({():s.sympify(x)})
    def __add__(self,other):
        other=E.scalar(other);out=dict(self.d)
        for k,v in other.d.items():out[k]=out.get(k,0)+v
        return E(out)
    __radd__=__add__
    def __neg__(self):return E({k:-v for k,v in self.d.items()})
    def __sub__(self,other):return self+-E.scalar(other)
    def __rsub__(self,other):return E.scalar(other)+-self
    def __mul__(self,other):
        other=E.scalar(other);out={}
        for a,x in self.d.items():
            for b,y in other.d.items():
                if set(a)&set(b):continue
                sign=(-1)**sum(i>j for i in a for j in b)
                key=tuple(sorted(a+b))
                out[key]=out.get(key,0)+sign*x*y
        return E(out)
    __rmul__=__mul__
    def __truediv__(self,other):return self*s.Rational(1,other)
    def star(self):
        out=E()
        for key,c in self.d.items():
            term=E.scalar(s.conjugate(c))
            for k in reversed(key):term=term*E({(conjugates[k],):1})
            out=out+term
        return out


def field(name):
    a=len(conjugates);conjugates[a]=a+1;conjugates[a+1]=a
    return E({(a,):1}),E({(a+1,):1})


def check(name,value):
    if E.scalar(value).d:raise RuntimeError(name+': exact identity failed')
    positive[name]+=1


def reject(name,value):
    if not E.scalar(value).d:raise RuntimeError(name+': bad formula escaped')
    negative[name]+=1


def closed(x):return x+x.star()


def mul3(a,b,c):
    # Coefficient of kappa from direct polynomial substitution.
    return a[1]*b[0]*c[0]+a[0]*b[1]*c[0]+a[0]*b[0]*c[1]


def main():
    eps=s.Matrix([[0,1],[-1,0]])
    sig=[s.eye(2),s.Matrix([[0,1],[1,0]]),
         s.Matrix([[0,-I],[I,0]]),s.diag(1,-1)]
    F=s.Matrix([[1+I,2-I],[3-2*I,-1+3*I]])
    Ye=s.Matrix([[1+2*I,3-I],[2-I,-1+I]])
    Yd=s.Matrix([[2-I,1+I],[-2+I,3+I]])
    Yu=s.Matrix([[1-I,2+I],[4-I,-2+3*I]])
    modes=[]
    for trial in range(3):
        phi=s.Matrix([1+I*(trial+1),2-trial+I])
        H=-eps*phi.conjugate()
        dp=[s.Matrix([2+mu+I*(trial+1),1+trial-I*(mu+1)])
            for mu in range(4)]
        p=[3+trial,1,-1,2]
        p2=p[0]**2-sum(x*x for x in p[1:]);m2=trial-2
        # AN_L has right chirality: lower Weyl block of slash(p) N_L.
        B=p[0]*sig[0]-sum((p[k]*sig[k] for k in range(1,4)),s.zeros(2))
        N=[];Nb=[]
        for i in range(2):
            pair=[field(('N',trial,i,r)) for r in range(2)]
            N.append([x[0] for x in pair]);Nb.append([x[1] for x in pair])
        bl=[[[field(('L',trial,a,w,r))[1] for r in range(2)]
              for w in range(2)] for a in range(2)]
        bq=[[[field(('Q',trial,a,w,r))[1] for r in range(2)]
              for w in range(2)] for a in range(2)]
        rights={name:[[field((name,trial,b,r))[0] for r in range(2)]
                      for b in range(2)] for name in ['e','d','u']}
        AN=[[sum((B[r,c]*N[i][c] for c in range(2)),E())
             for r in range(2)] for i in range(2)]
        barAN=[[AN[i][r].star() for r in range(2)] for i in range(2)]
        t=[sum((F[a,i]*bl[a][w][r]*AN[i][r]/s.Integer(6)
                for a in range(2) for i in range(2) for r in range(2)),E())
           for w in range(2)]
        ellbar=[[[E() for r in range(2)] for w in range(2)] for a in range(2)]
        cbar=[[[E() for r in range(2)] for w in range(2)] for a in range(2)]
        for a in range(2):
            for w in range(2):
                for r in range(2):
                    for i in range(2):
                        ellbar[a][w][r]+=s.conjugate(F[a,i])*(
                            -I*s.Rational(1,3)*sum((s.conjugate(dp[mu][w])*
                              barAN[i][c]*sig[mu][c,r]
                              for mu in range(4) for c in range(2)),E())
                            +s.conjugate(phi[w])*(s.Rational(1,6)*p2+
                              s.Rational(1,2)*m2)*Nb[i][r])
                        cbar[a][w][r]+=-s.Rational(1,2)*(p2+m2)*s.conjugate(F[a,i])*s.conjugate(phi[w])*Nb[i][r]
        dphi=[-x.star() for x in t]
        dH=[sum((eps[w,v]*t[v] for v in range(2)),E()) for w in range(2)]
        scalar_source=E()
        for w in range(2):
            for a in range(2):
                for b in range(2):
                    for r in range(2):
                        scalar_source+=sum((-t[w]*bl[a][v][r]*Ye[a,b]*eps[v,w]*rights['e'][b][r]
                                            -t[w]*bq[a][v][r]*Yd[a,b]*eps[v,w]*rights['d'][b][r]
                                            for v in range(2)),E())
                        # Up-type S_phi=-bar u_R Y_u^dagger Q, preserving order.
                        scalar_source+=t[w]*rights['u'][b][r].star()*s.conjugate(Yu[a,b])*bq[a][w][r].star()
        by_species={}
        norm_species={}
        bridge_species={}
        fermion_source=E()
        for name,Y,row,scalar,dscalar in [
            ('e',Ye,bl,H,dH),('d',Yd,bq,H,dH),('u',Yu,bq,phi,dphi)]:
            direct=E();expected=E();normal=E();bridge=E()
            for a in range(2):
                for b in range(2):
                    for w in range(2):
                        for r in range(2):
                            v=-ellbar[a][w][r] if name=='e' else E()
                            vc=cbar[a][w][r] if name=='e' else E()
                            right=rights[name][b][r]
                            direct+=-Y[a,b]*mul3((row[a][w][r],v),(E.scalar(scalar[w]),dscalar[w]),(right,E()))
                            normal+=-Y[a,b]*mul3((row[a][w][r],vc),(E.scalar(scalar[w]),E()),(right,E()))
                            bridge+=-Y[a,b]*mul3((row[a][w][r],vc-v),(E.scalar(scalar[w]),-dscalar[w]),(right,E()))
                            expected+=-Y[a,b]*(v*scalar[w]*right+row[a][w][r]*dscalar[w]*right)
                            if name=='e':fermion_source+=ellbar[a][w][r]*Y[a,b]*H[w]*right
            check(name+' direct source-map transport',direct-expected)
            check(name+' local-normalizing transport zero',normal)
            check(name+' coordinate-bridge spectator cancellation',bridge+direct)
            check(name+' conjugate transport',(direct-expected).star())
            check(name+' adjoint involution',direct.star().star()-direct)
            reject(name+' omitted source-map spectator',closed(direct))
            reject(name+' reversed scalar-shift sign',closed(sum((row[a][w][r]*Y[a,b]*dscalar[w]*rights[name][b][r]
                          for a in range(2) for b in range(2) for w in range(2) for r in range(2)),E())))
            by_species[name]=direct
            norm_species[name]=normal
            bridge_species[name]=bridge
        total=sum(by_species.values(),E())
        derivative_contact=E()
        FY=F.conjugate().T*Ye
        for i in range(2):
            for b in range(2):
                for mu in range(4):
                    weak=sum(s.conjugate(dp[mu][w])*H[w] for w in range(2))
                    for c in range(2):
                        for r in range(2):
                            derivative_contact+=-I*s.Rational(1,3)*FY[i,b]*weak*barAN[i][c]*sig[mu][c,r]*rights['e'][b][r]
        check('differentiated scalar charged-lepton source',fermion_source-derivative_contact)
        reject('charged-lepton derivative adjoint sign',fermion_source+derivative_contact)
        check('literal variational scalar and lepton sources',closed(total-fermion_source-scalar_source))
        reject('missing scalar source',closed(scalar_source))
        reject('missing derivative lepton source',closed(fermion_source))
        reject('Majorana half inserted into t',closed(scalar_source*s.Rational(1,2)))
        wrong_adjoint=E()
        for key,c in total.d.items():
            term=E.scalar(s.conjugate(c))
            for k in key:term=term*E({(conjugates[k],):1})
            wrong_adjoint+=term
        reject('conjugate fermion order unreversed',total.star()-wrong_adjoint)
        check('scalar adjoint involution',sum((x.star().star()-x for x in t),E()))
        modes.append({'trial':trial,'source_map_monomials':{k:len(v.d) for k,v in by_species.items()},
                      'total_with_conjugate_monomials':len(closed(total).d)})
    report={'status':'PASS','input_commit':'5f9889b76661faa155c681a5ecfac5d036398bc8',
            'scope':'Classical first-order spectator charged-lepton and quark Yukawa transport under the declared source and local-normalizing maps. Not total physical coefficient matching.',
            'positive_checks':sum(positive.values()),'positive_categories':dict(positive),
            'negative_checks':sum(negative.values()),'negative_categories':dict(negative),
            'trials':modes,'arithmetic':'exact Gaussian rational SymPy coefficients and exterior Grassmann generators',
            'repository_imports':False,'assertions_required':False}
    Path(__file__).with_name('SPECTATOR_YUKAWA_RESULTS.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
