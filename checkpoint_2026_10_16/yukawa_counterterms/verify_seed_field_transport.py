#!/usr/bin/env python3
"""Seed counterterm from the independently computed pure-Yukawa field poles.

Direct four-field, varying-scalar Fourier/Grassmann substitution; formal
canonical-epsilon RG conversion. No project imports or removable assertions.
"""
from collections import Counter
from pathlib import Path
import json
import sympy as s

I=s.I;R=s.Rational
positive=Counter();negative=Counter()


class P:
    def __init__(self,data=None):
        self.d={k:s.expand(v) for k,v in (data or {}).items() if s.expand(v)!=0}
    @staticmethod
    def scalar(x):return x if isinstance(x,P) else P({((0,0),()):s.sympify(x)})
    def __add__(self,x):
        x=P.scalar(x);out=dict(self.d)
        for k,v in x.d.items():out[k]=out.get(k,0)+v
        return P(out)
    __radd__=__add__
    def __neg__(self):return P({k:-v for k,v in self.d.items()})
    def __sub__(self,x):return self+-P.scalar(x)
    def __rsub__(self,x):return P.scalar(x)+-self
    def __mul__(self,x):
        x=P.scalar(x);out={}
        for (p,a),v in self.d.items():
            for (q,b),w in x.d.items():
                if set(a)&set(b):continue
                sign=(-1)**sum(i>j for i in a for j in b)
                k=(tuple(i+j for i,j in zip(p,q)),tuple(sorted(a+b)))
                out[k]=out.get(k,0)+sign*v*w
        return P(out)
    __rmul__=__mul__
    def p0(self):return P({k:(k[0][0]**2-k[0][1]**2)*v for k,v in self.d.items()})
    def integral(self):return P({k:v for k,v in self.d.items() if k[0]==(0,0)})


def zero(x):
    if isinstance(x,P):return not x.d
    if isinstance(x,s.MatrixBase):return all(s.expand(v)==0 for v in x)
    return s.expand(x)==0


def check(name,x):
    if not zero(x):raise RuntimeError(name+': exact identity failed')
    positive[name]+=1


def reject(name,x):
    if zero(x):raise RuntimeError(name+': mutation escaped')
    negative[name]+=1


def main():
    eps=s.Matrix([[0,1],[-1,0]])
    Y=s.Matrix([[1+I,2-I],[3-2*I,-1+I]])
    C=s.Matrix([[2-I,1+2*I],[1+2*I,3+I]])
    HL=Y*Y.H;HR=Y.H*Y;T=s.trace(HR)
    zL=-HL/2;zphi=-T  # pole coefficients, common kappa stripped
    rc=(HL*C+C*HL.T)/4+T*C
    check('symmetric seed field-pole coefficient',rc-rc.T)
    check('constant field coefficient formula',(zL*C+C*zL.T)/2+zphi*C+rc)
    for trial in range(2):
        phi=[P({((1+trial,1),()):1+I,((-1-trial,-1),()):2-I}),
             P({((2,0),()):2+I,((-2,0),()):-1+2*I})]
        bl=[[[P() for r in range(2)] for a in range(2)] for alpha in range(2)]
        index=0
        for alpha in range(2):
            for a in range(2):
                for r in range(2):
                    bl[alpha][a][r]=P({((3+alpha+a,1+r),(index,)):1,
                                      ((-3-alpha-a,-1-r),(index+1,)):1+I})
                    index+=2
        dbl=[[[sum((zL[g,alpha]*bl[g][a][r]/s.Integer(2)
                      for g in range(2)),P()) for r in range(2)]
                for a in range(2)] for alpha in range(2)]
        dphi=[zphi*x/2 for x in phi]
        xbar=[[sum((bl[alpha][a][r]*phi[a] for a in range(2)),P())
               for r in range(2)] for alpha in range(2)]
        dxbar=[[sum((dbl[alpha][a][r]*phi[a]+bl[alpha][a][r]*dphi[a]
                     for a in range(2)),P()) for r in range(2)] for alpha in range(2)]
        xc=[[sum((-eps[r,t]*xbar[alpha][t] for t in range(2)),P())
             for r in range(2)] for alpha in range(2)]
        dxc=[[sum((-eps[r,t]*dxbar[alpha][t] for t in range(2)),P())
              for r in range(2)] for alpha in range(2)]
        def op(coeff,rows=xbar,cols=xc):
            return sum((R(1,2)*coeff[a,b]*rows[a][r]*cols[b][r].p0()
                        for a in range(2) for b in range(2) for r in range(2)),P())
        direct=op(C,dxbar,xc)+op(C,xbar,dxc)
        check('literal four-field bare seed substitution',direct+op(rc))
        check('integrated bare seed substitution',direct.integral()+op(rc).integral())
        check('field-specific parameter pole cancels field substitution',direct+op(rc))
        check('scalar/Box symmetric coefficient',op(rc-rc.T).integral())
        reject('omit the field-induced seed pole',direct.integral())
        reject('miss one Higgs half-field contribution',(direct+op(rc-T*C/2)).integral())
        reject('count each Higgs field twice',(direct+op(rc+T*C)).integral())
        reject('adjoint instead of transpose on right seed slot',
               (direct+op((HL*C+C*HL.H)/4+T*C)).integral())
        reject('extra Majorana half on seed parameter pole',(direct+op(rc/2)).integral())
        # Applying P0 only to the lepton, then multiplying phi, loses scalar
        # and mixed derivatives in the actual composite operator.
        bad_cols=[[sum((-eps[r,t]*phi[a]*bl[alpha][a][t].p0()
                        for a in range(2) for t in range(2)),P()) for r in range(2)]
                  for alpha in range(2)]
        bad=sum((R(1,2)*rc[a,b]*xbar[a][r]*bad_cols[b][r]
                 for a in range(2) for b in range(2) for r in range(2)),P())
        reject('delete derivatives on Higgs inside Xc',(op(rc)-bad).integral())
    # Formal d=4-2epsilon bookkeeping: Y0=mu^epsilon(...),
    # C70=mu^(2epsilon)(...). Y and Ybar are independent formal variables.
    y=s.Matrix(2,2,s.symbols('y0:4'));yb=s.Matrix(2,2,s.symbols('yb0:4'))
    cs=s.symbols('c0:3');c=s.Matrix([[cs[0],cs[1]],[cs[1],cs[2]]])
    hl=y*yb.T;hr=yb.T*y;t=s.trace(hr)
    q=(hl*c+c*hl.T)/4+t*c
    weighted=[(x,1) for x in list(y)+list(yb)]+[(x,2) for x in cs]
    def canonical_beta(residue,bare_power):
        return residue.applyfunc(lambda f:s.expand(sum(w*x*s.diff(f,x) for x,w in weighted)-bare_power*f))
    beta=canonical_beta(q,2)
    check('formal epsilon residue-to-RG conversion',beta-2*q)
    check('field anomalous-dimension form',beta-(hl*c+c*hl.T)/2-2*t*c)
    reject('UV pole identified directly with beta',beta-q)
    reject('scalar anomalous dimension counted once',beta-(hl*c+c*hl.T)/2-t*c)
    # Dimensional analysis in d fixes the bare mu exponent independently:
    # op=2*(3/2-eps)+2*(1-eps)+2=7-4eps.
    e=s.symbols('epsilon')
    check('seed bare mu exponent',(4-2*e)-(2*(R(3,2)-e)+2*(1-e)+2)-(-3+2*e))
    V=s.diag(I,(3+4*I)/5)*s.Matrix([[R(3,5),R(4,5)],[-R(4,5),R(3,5)]])
    U=s.Matrix([[R(5,13),R(12,13)],[-R(12,13),R(5,13)]])*s.diag((5+12*I)/13,I)
    Yp=V.H*Y*U.conjugate();Cp=V.H*C*V.conjugate()
    hp=Yp*Yp.H;tp=s.trace(Yp.H*Yp)
    rcp=(hp*Cp+Cp*hp.T)/4+tp*Cp
    check('unitary active basis',V.H*V-s.eye(2))
    check('unitary sterile basis',U.H*U-s.eye(2))
    check('seed field-pole flavor covariance',rcp-V.H*rc*V.conjugate())
    check('sterile-basis trace covariance',tp-T)
    reject('transpose discarded in complex basis',rcp-(hp*Cp+Cp*hp)/4-tp*Cp)
    result={'status':'PASS','input_commit':'0770c304c7e05f213dd796d3fd5a194107f38858',
            'scope':'External-field contribution to seed parameter pole and its canonical-epsilon RG input; not the full seed anomalous dimension.',
            'positive_checks':sum(positive.values()),'positive_categories':dict(positive),
            'negative_checks':sum(negative.values()),'negative_categories':dict(negative),
            'source_inputs':{'delta_ZL_over_kappa':'-YYdagger/2','delta_Zphi_over_kappa':'-Tr(YdaggerY)'},
            'new_loop_integration':False,'repository_imports':False,'assertions_required':False}
    Path(__file__).with_name('SEED_FIELD_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


P.__truediv__=lambda self,x:self*(1/s.sympify(x))
if __name__=='__main__':main()
