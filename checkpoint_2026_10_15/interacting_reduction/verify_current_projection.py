#!/usr/bin/env python3
"""Exact weak-index current and broken-basis field-shift controls.

This standalone program checks action-coordinate identities, not loop integrals
or physical observables. It imports no project code. All numbers are exact.
"""
from pathlib import Path
import json
import sympy as s

I=s.I
T=[s.Matrix([[0,1],[1,0]])/2,
   s.Matrix([[0,-I],[I,0]])/2,s.diag(1,-1)/2]
eye=s.eye(2)
counts={}
mutations={}


def check(name,left,right):
    difference=left-right
    if isinstance(difference,s.MatrixBase):
        ok=all(s.simplify(v)==0 for v in difference)
    else:ok=s.simplify(difference)==0
    if not ok:raise RuntimeError(name+': '+str(difference))
    counts[name]=counts.get(name,0)+1


def reject(name,wrong,right):
    difference=wrong-right
    if isinstance(difference,s.MatrixBase):
        detected=any(s.simplify(v)!=0 for v in difference)
    else:detected=s.simplify(difference)!=0
    if not detected:raise RuntimeError('Undetected mutation: '+name)
    mutations[name]=mutations.get(name,0)+1


def bilinear(row,mat,col):return (row*mat*col)[0]


def main():
    # Prove the fundamental completeness tensor component by component.
    for a in range(2):
        for b in range(2):
            for c in range(2):
                for d in range(2):
                    tensor=sum(t[a,b]*t[c,d] for t in T)
                    rhs=s.Rational(1,2)*int(a==d and b==c)-s.Rational(1,4)*int(a==b and c==d)
                    check('Pauli completeness components',tensor,rhs)
    for n in range(1,9):
        phi=s.Matrix([n+I,2*n-1-I*n])
        dphi=s.Matrix([1-2*I*n,n*n+I])
        b=s.Matrix([[n-2*I,1+I*n]])
        # M_cd stands for the already ordered doublet current bar f_c gamma f_d.
        # No Grassmann field is exchanged in the completeness contraction.
        M=s.Matrix([[n+I,2-I*n],[-1+I*n,n*n-3*I]])
        j=(b*phi)[0]
        jI=[bilinear(b,t,phi) for t in T]
        fcurrent=[sum(t[c,d]*M[c,d] for c in range(2) for d in range(2)) for t in T]
        direct=sum(jI[k]*fcurrent[k] for k in range(3))
        exchange=sum(b[a]*phi[bb]*M[bb,a] for a in range(2) for bb in range(2))/2
        singlet=j*s.trace(M)/4
        check('ordered doublet-current contraction',direct,exchange-singlet)
        reject('missing weak singlet subtraction',exchange,direct)
        rho=(phi.conjugate().T*phi)[0]
        d_rho=(phi.conjugate().T*dphi+dphi.conjugate().T*phi)[0]
        scalarI=[I*(bilinear(phi.conjugate().T,t,dphi)-bilinear(dphi.conjugate().T,t,phi)) for t in T]
        direct_scalar=sum(jI[k]*scalarI[k] for k in range(3))
        reduced_scalar=I*rho*(b*dphi)[0]/2-I*j*d_rho/4
        check('scalar weak-current contraction',direct_scalar,reduced_scalar)
        reject('missing scalar norm derivative',I*rho*(b*dphi)[0]/2,direct_scalar)
        # A genuinely non-diagonal rational SU(2) rotation of both slots.
        U=((1-n*n)*eye+2*I*n*(2*T[1]))/(1+n*n)
        check('finite weak rotation unitarity',U.conjugate().T*U,eye)
        phip=U*phi;bp=b*U.conjugate().T;dp=U*dphi
        jIp=[bilinear(bp,t,phip) for t in T]
        scalarIp=[I*(bilinear(phip.conjugate().T,t,dp)-bilinear(dp.conjugate().T,t,phip)) for t in T]
        check('weak scalar contraction covariance',sum(jIp[k]*scalarIp[k] for k in range(3)),direct_scalar)
        check('singlet current covariance',(bp*phip)[0],j)
    # Treat the current and its Hermitian conjugate as independent symbols,
    # preserving the complex charged-current orientation.
    j0,j0c,jc,jcc,u=s.symbols('j0 j0c jc jcc u',commutative=True)
    g,gp=s.symbols('g gp',positive=True)
    gz=s.sqrt(g*g+gp*gp);sw=gp/gz;cw=g/gz
    dB=-gp*(j0+j0c)/12
    dW=[g*(jc+jcc)/12,I*g*(jc-jcc)/12,g*(j0+j0c)/12]
    dA=cw*dB+sw*dW[2]
    dZ=-sw*dB+cw*dW[2]
    dWp=(dW[0]-I*dW[1])/s.sqrt(2)
    dWm=(dW[0]+I*dW[1])/s.sqrt(2)
    check('neutral photon field shift',dA,0)
    check('neutral Z field shift',dZ,gz*(j0+j0c)/12)
    check('positive charged field shift',dWp,g*jc/(6*s.sqrt(2)))
    check('negative charged field shift',dWm,g*jcc/(6*s.sqrt(2)))
    reject('hypercharge-only photon shift',cw*dB,dA)
    reject('wrong weak generator normalization',cw*dB+2*sw*dW[2],dA)
    reject('swapped charged-current orientation',g*jcc/(6*s.sqrt(2)),dWp)
    # Check the full unitary rotation of the source/current contraction, not
    # just its neutral projection. Gauge currents are independently retained.
    JY,J1,J2,J3=s.symbols('JY J1 J2 J3',real=True)
    jy=j0;ji=[jc/2,I*jc/2,j0/2]
    raw=-(gp*gp*(-s.Rational(1,2))*jy*JY+g*g*sum(ji[k]*[J1,J2,J3][k] for k in range(3)))/6
    JZ=cw*g*J3-sw*gp*JY
    reduced=-gz*j0*JZ/12-g*g*jc*(J1+I*J2)/12
    check('full neutral and charged source-current projection',raw,reduced)
    reject('free weak equation drops charged source',-gz*j0*JZ/12,raw)
    report={'status':'PASS','input_commit':'d822e9f4d70ad5b455436f40abc3c2b45909bbef',
            'scope':'Weak-index completeness, scalar and ordered fermion current descendants, finite weak covariance, real gauge-shift and source-current rotations. No loops or physical running.',
            'positive_checks':sum(counts.values()),'positive_categories':counts,
            'mutation_rejections':sum(mutations.values()),'mutation_categories':mutations,
            'exact_arithmetic':True,'repository_imports':False,'assertions_required':False}
    Path('CURRENT_PROJECTION_RESULTS.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
