#!/usr/bin/env python3
"""Exact checks for the selected source-preserving joint phase obstruction.

This is finite-dimensional flavor algebra applied to existing loop inputs;
it performs no loop integration and imports no project files.
"""
from collections import Counter
from pathlib import Path
import json
import sympy as s
I=s.I; Q=s.Rational
ok=Counter(); bad=Counter()

def zero(x):
    if isinstance(x,s.MatrixBase): return all(s.simplify(s.expand(z))==0 for z in x)
    return s.simplify(s.expand(x))==0

def check(name,x):
    if not zero(x): raise RuntimeError(name+': identity failed')
    ok[name]+=1

def reject(name,x):
    if zero(x): raise RuntimeError(name+': mutation survived')
    bad[name]+=1

def projector(Y):
    cols=Y.columnspace()
    if not cols:return s.zeros(Y.rows)
    B=s.Matrix.hstack(*cols)
    return (B*(B.H*B).inv()*B.H).applyfunc(s.simplify)

def quantities(Y,M,C):
    R=Y.conjugate()*M*Y.H
    K=C*R
    A=(K-K.H)/(2*I)
    a=s.trace(A)
    P=projector(Y)
    return R,K,A,a,P

def main():
    cases=[
      ('zero',s.zeros(3,2),s.eye(2),s.eye(3)),
      ('one flavor',s.Matrix([[1+I]]),s.Matrix([[2]]),s.Matrix([[1+I]])),
      ('rank1 source-free',s.Matrix([[1,0],[0,0],[0,0]]),s.diag(2,1),s.diag(I,0,0)),
      ('rank1 source survives',s.Matrix([[1,0],[0,0],[0,0]]),s.diag(2,1),s.Matrix([[I,1,0],[1,0,0],[0,0,0]])),
      ('rank2 source-free',s.Matrix([[1,0],[0,1],[0,0]]),s.eye(2),s.eye(3)),
      ('rank2 source survives',s.Matrix([[1,0],[0,1],[0,0]]),s.Matrix([[2,1],[1,3]]),I*s.eye(3)),
      ('rank3 singular mass',s.eye(3),s.diag(1,0,0),s.diag(I,0,0)),
      ('complex rectangular',s.Matrix([[1+I,2-I],[3,-I],[1-2*I,-2]]),s.Matrix([[2+I,1-2*I],[1-2*I,4+I]]),s.Matrix([[2-I,1+2*I,I],[1+2*I,3+I,2], [I,2,1-I]]))]
    for label,Y,M,C in cases:
        n=Y.rows; E=s.eye(n); R,K,A,a,P=quantities(Y,M,C); r=Y.rank()
        check('symmetric Majorana and seed inputs',M-M.T);check('symmetric Majorana and seed inputs',C-C.T)
        check('orthogonal Yukawa support projector',P*P-P);check('orthogonal Yukawa support projector',P-P.H)
        check('Yukawa support',P*Y-Y);check('R annihilates complement',R*(E-P));check('K annihilates complement',K*(E-P))
        check('anti-Hermitian flavor part is Hermitian',A-A.H);check('complement block of A vanishes',(E-P)*A*(E-P))
        check('trace relation',a-(s.trace(K)-s.conjugate(s.trace(K)))/(2*I))
        check('source criterion equivalence residual',(A-a*E)*Y-(A*P-a*P)*Y)
        sourcefree=zero((A-a*E)*Y)
        if sourcefree:
            check('source-free implies full A equals trace times support',A-a*P)
            check('trace rank constraint',(r-1)*a)
            if r>=2:check('rank >= 2 source-free implies A=0',A)
        if r==0:check('zero Yukawa implies zero A',A)
        if r>=2:
            if sourcefree != zero(A): raise RuntimeError('rank >= 2 source-free equivalence failed: '+label)
            check('rank >= 2 equivalence in both directions',s.Integer(sourcefree)-s.Integer(zero(A)))
        elif r==1:
            if sourcefree != zero(A-a*P): raise RuntimeError('rank-one source-free equivalence failed: '+label)
            check('rank-one equivalence in both directions',s.Integer(sourcefree)-s.Integer(zero(A-a*P)))
        Vcase=s.eye(n)
        Vcase[0,0]=I
        if n>=2:
            mix=s.eye(n); mix[0,0]=mix[1,1]=Q(3,5); mix[0,1]=mix[1,0]=4*I/5
            Vcase=Vcase*mix
        ns=Y.cols; Ucase=s.eye(ns); Ucase[0,0]=I
        if ns>=2:
            mix=s.eye(ns); mix[0,0]=mix[1,1]=Q(5,13); mix[0,1]=mix[1,0]=12*I/13
            Ucase=Ucase*mix
        Ycase=Vcase.H*Y*Ucase.conjugate(); Mcase=Ucase.H*M*Ucase.conjugate(); Ccase=Vcase.H*C*Vcase.conjugate()
        _,_,Acase,acase,Pcase=quantities(Ycase,Mcase,Ccase)
        check('all-rank rotated support covariance',Pcase-Vcase.H*P*Vcase)
        check('all-rank rotated phase covariance',Acase-Vcase.H*A*Vcase)
        check('all-rank rotated source covariance',(Acase-acase*s.eye(n))*Ycase-Vcase.H*(A-a*E)*Y*Ucase.conjugate())
        if sourcefree != zero((Acase-acase*s.eye(n))*Ycase): raise RuntimeError('source-free property changed by rotation: '+label)
        check('all-rank rotated equivalence',s.Integer(sourcefree)-s.Integer(zero((Acase-acase*s.eye(n))*Ycase)))
        if label=='one flavor':
            check('one-flavor neutrino source cancels',(A-a*E)*Y)
            reject('one-flavor charged-lepton spectator survives',(A+a*E)*s.ones(1,1))
            reject('omitting scalar phase produces false neutrino source',A*Y)
        if label=='rank1 source-free':
            check('rank-one nonzero phase can be source-free',(A-a*E)*Y)
            reject('rank-one source-free example has nonzero A',A)
        if label in ('rank1 source survives','rank2 source survives','rank3 singular mass'):
            reject('generic selected neutrino source survives',(A-a*E)*Y)
    Y=cases[-1][1];M=cases[-1][2];C=cases[-1][3];R,K,A,a,P=quantities(Y,M,C)
    V=s.eye(3);V[0,0]=Q(3,5);V[0,1]=4*I/5;V[1,0]=4*I/5;V[1,1]=Q(3,5)
    U=s.Matrix([[Q(5,13),12*I/13],[12*I/13,Q(5,13)]])
    Yp=V.H*Y*U.conjugate();Mp=U.H*M*U.conjugate();Cp=V.H*C*V.conjugate()
    Rp,Kp,Ap,ap,Pp=quantities(Yp,Mp,Cp)
    check('active rotation unitary',V.H*V-s.eye(3));check('sterile rotation unitary',U.H*U-s.eye(2))
    check('R covariance',Rp-V.T*R*V);check('K covariance',Kp-V.H*K*V)
    check('phase matrix covariance',Ap-V.H*A*V);check('trace scalar invariant',ap-a)
    check('Yukawa support covariance',Pp-V.H*P*V)
    check('neutrino source covariance',(Ap-ap*s.eye(3))*Yp-V.H*(A-a*s.eye(3))*Y*U.conjugate())
    Rwrong=Y.conjugate()*M*Y.T
    reject('wrong Yukawa adjoint ordering',Rwrong*(s.eye(3)-P))
    # General Hermitian block proof for r=2, n=3: AP=aP forces both
    # off-block entries to zero, and trace A=a forces (r-1)a=0.
    aa,xr,xi,yr,yi=s.symbols('a xr xi yr yi',real=True)
    Ab=s.Matrix([[aa,0,xr+I*xi],[0,aa,yr+I*yi],[xr-I*xi,yr-I*yi,0]])
    Pb=s.diag(1,1,0)
    check('rank-two support trace identity',s.trace(Ab)-2*aa)
    check('rank-two remaining trace condition',s.trace(Ab)-aa-aa)
    reject('rank-two off-block source omitted',(Ab-aa*s.eye(3))*Pb)
    result={'status':'PASS','scope':'Exact algebraic obstruction for the selected first-order joint scalar/lepton phase, using existing seed/Yukawa input; not an independent loop or full physical coefficient.',
      'positive_checks':sum(ok.values()),'negative_checks':sum(bad.values()),'positive_categories':dict(ok),'negative_categories':dict(bad),
      'source':'(A-tr(A) I)Y','spectator':'(A+tr(A) I)Y_e',
      'rank_2_or_more_source_free_iff':'A=0','rank_1_source_free_iff':'A=tr(A) P',
      'projector_note':'The implementation uses an inverse only for a chosen independent column Gram matrix; no sterile mass inverse is used.',
      'repository_imports':False,'assertions_required':False}
    Path(__file__).with_name('YUKAWA_PHASE_OBSTRUCTION_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
