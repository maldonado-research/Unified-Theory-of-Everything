#!/usr/bin/env python3
"""Exact complex-Majorana mass recursion and source-coordinate covariance.

Standalone SymPy coefficient checks, with no project imports or assertions.
The fields/sources are arbitrary chiral Fourier jets, not on-shell waves.
"""
from collections import Counter
from pathlib import Path
import json
import sympy as s

I=s.I
positive=Counter();negative=Counter()


def clean(x):
    if isinstance(x,s.MatrixBase):return x.applyfunc(s.expand)
    return s.expand(x)


def iszero(x):
    x=clean(x)
    return all(v==0 for v in x) if isinstance(x,s.MatrixBase) else x==0


def check(name,x):
    if not iszero(x):raise RuntimeError(name+': exact identity failed')
    positive[name]+=1


def reject(name,x):
    if iszero(x):raise RuntimeError(name+': mutation escaped')
    negative[name]+=1


def jets(n,offset,left):
    return s.Matrix(n,4,lambda a,k:
        (offset+2*a+k+I*(offset-a+2*k)) if (k<2)==left else 0)


def main():
    z=s.zeros(2);one=s.eye(2)
    pauli=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-I],[I,0]]),s.diag(1,-1)]
    gam=[s.BlockMatrix([[z,one],[one,z]]).as_explicit()]
    gam += [s.BlockMatrix([[z,t],[-t,z]]).as_explicit() for t in pauli]
    eta=[1,-1,-1,-1]
    for mu in range(4):
        for nu in range(4):
            check('Clifford convention',gam[mu]*gam[nu]+gam[nu]*gam[mu]-2*(eta[mu] if mu==nu else 0)*s.eye(4))
    Y=s.Matrix([[1+I,2-I],[3-2*I,-1+I],[2+3*I,4-I]])
    C=s.Matrix([[2,1+I,3-I],[1+I,4-I,2*I],[3-I,2*I,-1+2*I]])
    F=C*Y.conjugate()
    check('symmetric input seed',C-C.T)
    qL=jets(2,1,True);qR=jets(2,3,False)
    SR=jets(2,5,False);SLN=jets(2,7,True)
    light=jets(3,2,True);lightc=jets(3,4,False)
    masses=[('generic complex symmetric',s.Matrix([[2+I,1+2*I],[1+2*I,3-I]])),
            ('singular complex symmetric',s.Matrix([1+I,2-I])*s.Matrix([1+I,2-I]).T),
            ('zero',s.zeros(2))]
    rot=s.Matrix([[s.Rational(3,5),s.Rational(4,5)],[-s.Rational(4,5),s.Rational(3,5)]])
    Us=[s.diag(I,(3+4*I)/5)*rot,rot*s.diag((5+12*I)/13,I)]
    V=s.diag((3+4*I)/5,I,1)*s.Matrix([[1,0,0],[0,s.Rational(5,13),s.Rational(12,13)],[0,-s.Rational(12,13),s.Rational(5,13)]])
    check('active unitary coordinate change',V.H*V-s.eye(3))
    momenta=[(2,1,-1,1),(3,-2,1,1),(s.Rational(5,3),s.Rational(2,3),1,0)]
    cases=[]
    for label,M in masses:
        check('symmetric Majorana mass',M-M.T)
        for p in momenta:
            slash=sum((eta[k]*p[k]*gam[k] for k in range(4)),s.zeros(4))
            p2=sum(eta[k]*p[k]**2 for k in range(4))
            A=lambda v:clean(v*slash.T)
            check('Fourier A squared',slash*slash-p2*s.eye(4))
            eR=clean(A(qL)-M*qR-SR)
            eL=clean(A(qR)-M.H*qL-SLN)
            R=clean(M*qR+SR)
            T=clean(M*M.H*qL+M*SLN+A(SR))
            UU=clean(M*M.H*M*qR+M*M.H*SR+M*A(SLN)+A(A(SR)))
            check('first sourced recurrence',A(qL)-R-eR)
            check('second sourced recurrence',A(A(qL))-T-M*eL-A(eR))
            check('third sourced recurrence',A(A(A(qL)))-UU-M*M.H*eR-M*A(eL)-A(A(eR)))
            if M.rank():
                wrongT=M*M*qL+M*SLN+A(SR)
                reject('M squared used for complex mass',A(A(qL))-wrongT-M*eL-A(eR))
                wrongU=M*M*M*qR+M*M*SR+M*A(SLN)+A(A(SR))
                reject('real-mass cubic rule reused',A(A(A(qL)))-wrongU-M*M.H*eR-M*A(eL)-A(A(eR)))
                reject('adjoint missing from opposite chiral equation',eL-(A(qR)-M*qL-SLN))
            for U in Us:
                check('sterile unitary coordinate change',U.H*U-s.eye(2))
                Mp=clean(U.H*M*U.conjugate())
                Yp=clean(V.H*Y*U.conjugate())
                Cp=clean(V.H*C*V.conjugate())
                Fp=clean(Cp*Yp.conjugate())
                check('mass congruence remains symmetric',Mp-Mp.T)
                check('selected F coordinate change',Fp-V.H*F*U)
                lp=clean(U.H*qL);rp=clean(U.T*qR)
                srp=clean(U.H*SR);slnp=clean(U.T*SLN)
                erp=clean(A(lp)-Mp*rp-srp)
                elp=clean(A(rp)-Mp.H*lp-slnp)
                Rp=clean(Mp*rp+srp)
                Tp=clean(Mp*Mp.H*lp+Mp*slnp+A(srp))
                Up=clean(Mp*Mp.H*Mp*rp+Mp*Mp.H*srp+Mp*A(slnp)+A(A(srp)))
                for name,old,new in [('eR',eR,erp),('R',R,Rp),('T',T,Tp),('U',UU,Up)]:
                    check(name+' sterile covariance',new-U.H*old)
                    check(name+' contracted active covariance',Fp*new-V.H*F*old)
                check('opposite-chirality source covariance',slnp-U.T*SLN)
                check('opposite-chirality residual covariance',elp-U.T*eL)
                # Restore the actual light composite sources independently.
                actual_sr=Y.T*lightc;actual_sln=Y.H*light
                actual_srp=Yp.T*(V.T*lightc)
                actual_slnp=Yp.H*(V.H*light)
                check('literal active source S_R covariance',actual_srp-U.H*actual_sr)
                check('literal active source S_LN covariance',actual_slnp-U.T*actual_sln)
                for name,old,new in [
                    ('active neutral Box',F*Y.T,Fp*Yp.T),
                    ('active mass-squared',F*M*M.H*Y.T,Fp*Mp*Mp.H*Yp.T)]:
                    check(name+' congruence',new-V.H*old*V.conjugate())
                    check(name+' symmetric coefficient congruence',(new+new.T)/2-V.H*((old+old.T)/2)*V.conjugate())
                G=F*M*Y.H;Gp=Fp*Mp*Yp.H
                check('active LNC derivative similarity',Gp-V.H*G*V)
                check('active LNC Hermitian projection',(Gp+Gp.H)-V.H*(G+G.H)*V)
                check('retained-mass tensor transport',Fp*Mp-V.H*F*M*U.conjugate())
                HN=Y.H*F;HNp=Yp.H*Fp
                check('sterile scalar/tensor source congruence',HNp-U.T*HN*U)
                check('selected sterile source symmetry',HNp-HNp.T)
                # Integrated field-map rows: bar x transforms with U and
                # bar y with U*, hence x+y^c has the correct left family.
                W=s.Matrix([[1+I,2-I]]);Vrow=s.Matrix([[3-I,-2+I]])
                K0W=s.Matrix([[2+I,1-3*I]]);K02W=s.Matrix([[-1+I,4+2*I]])
                Jrow=s.Matrix([[2-I,3+I]])
                bx=s.Rational(5,6)*(W*M*M.H+K02W)+s.Rational(7,9)*W+Jrow
                by=s.Rational(5,6)*K0W*M-s.Rational(1,6)*Vrow*M
                bxp=s.Rational(5,6)*(W*U*Mp*Mp.H+K02W*U)+s.Rational(7,9)*W*U+Jrow*U
                byp=s.Rational(5,6)*K0W*U*Mp-s.Rational(1,6)*Vrow*U*Mp
                check('integrated x field-map row covariance',bxp-bx*U)
                check('integrated y field-map row covariance',byp-by*U.conjugate())
                if M.rank():
                    reject('unitary similarity substituted for Majorana congruence',Mp-U.H*M*U)
                    wrong=Fp*Mp*Mp*Yp.T
                    reject('active mass-squared real rule not covariant',wrong-V.H*(F*M*M*Y.T)*V.conjugate())
            cases.append({'mass':label,'rank':M.rank(),'p_squared':str(p2)})
    result={'status':'PASS','input_commit':'4ed8560e019dfe74ed6c5178b0f990fe3ccb3008',
            'scope':'Complex symmetric retained mass: chiral sourced recurrences, ordered residuals, integrated field-map rows, and active coefficient covariance. No physical matching or running.',
            'positive_checks':sum(positive.values()),'positive_categories':dict(positive),
            'negative_checks':sum(negative.values()),'negative_categories':dict(negative),
            'mass_momentum_cases':cases,'sterile_unitaries':2,'active_flavors':3,'sterile_flavors':2,
            'arithmetic':'exact Gaussian rational SymPy matrices','repository_imports':False,'assertions_required':False}
    Path(__file__).with_name('COMPLEX_MASS_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
