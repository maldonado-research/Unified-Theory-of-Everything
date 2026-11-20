#!/usr/bin/env python3
"""Scalar first-jet metric and linearized curvature of the retained pole packet."""
from collections import Counter
from pathlib import Path
import json
import sympy as s

positive=Counter(); adverse=Counter()
def zero(value):
    if isinstance(value,s.MatrixBase):return all(s.expand(v)==0 for v in value)
    return s.expand(value)==0

def check(name,value):
    if not zero(value):raise RuntimeError(name+' failed: '+str(value))
    positive[name]+=1

def reject(name,value):
    if zero(value):raise RuntimeError(name+' mutation escaped')
    adverse[name]+=1

def first_order(value,h):
    return value.applyfunc(lambda v:s.expand(v).coeff(h,0)+h*s.expand(v).coeff(h,1))

def curvature(metric_correction,q):
    n=len(q)
    return {(i,j,k,l):s.expand((s.diff(metric_correction[i,l],q[j],q[k])
              +s.diff(metric_correction[j,k],q[i],q[l])
              -s.diff(metric_correction[i,k],q[j],q[l])
              -s.diff(metric_correction[j,l],q[i],q[k]))/2)
            for i in range(n) for j in range(n) for k in range(n) for l in range(n)}

def main():
    h=s.symbols('h',real=True)
    q=s.Matrix(s.symbols('q0:4',real=True));dq=s.Matrix(s.symbols('u0:4',real=True))
    J=s.diag(s.Matrix([[0,-1],[1,0]]),s.Matrix([[0,-1],[1,0]]))
    eye=s.eye(4);v=J*q;radius2=q.dot(q)
    phi=s.Matrix([(q[0]+s.I*q[1])/s.sqrt(2),(q[2]+s.I*q[3])/s.sqrt(2)])
    dphi=s.Matrix([(dq[0]+s.I*dq[1])/s.sqrt(2),(dq[2]+s.I*dq[3])/s.sqrt(2)])
    rho=(phi.H*phi)[0];drho=(dphi.H*phi+phi.H*dphi)[0]
    current=s.I*((phi.H*dphi)[0]-(dphi.H*phi)[0]);kinetic=(dphi.H*dphi)[0]
    check('complex structure antisymmetry',J.T+J)
    check('complex structure square',J*J+eye)
    check('literal scalar density',rho-radius2/2)
    check('literal density derivative',drho-q.dot(dq))
    check('literal scalar current',current-(q.T*J*dq)[0])
    check('literal scalar kinetic term',kinetic-dq.dot(dq)/2)
    H_old=q*q.T-v*v.T-radius2*eye
    H_new=3*q*q.T-v*v.T
    old_packet=drho**2/2-current**2/2-2*rho*kinetic
    new_packet=(3*drho**2-current**2)/2
    check('old scalar packet gives original metric',old_packet-(dq.T*H_old*dq)[0]/2)
    check('remaining scalar packet gives transformed metric',new_packet-(dq.T*H_new*dq)[0]/2)
    # Actual radial coordinate substitution q_old=q+h*rho*q, including its full Jacobian.
    xi=radius2*q/2
    old_coordinates=q+h*xi
    jacobian=old_coordinates.jacobian(q)
    old_metric=eye+h*H_old
    substituted=old_metric.subs(dict(zip(q,old_coordinates)),simultaneous=True)
    pulled_back=first_order(jacobian.T*substituted*jacobian,h)
    check('literal radial pullback',pulled_back-eye-h*H_new)
    check('metric change is radial Lie derivative',H_new-H_old-xi.jacobian(q)-xi.jacobian(q).T)
    reject('radial map with omitted density derivative',2*(radius2/2)*eye-(H_new-H_old))
    reject('half radial map',s.Rational(1,2)*(xi.jacobian(q)+xi.jacobian(q).T)-(H_new-H_old))
    # Both metrics have vanishing first derivatives at q=0. The exact curvature
    # at that point equals the displayed linearized tensor (Christoffel squares vanish).
    R_old=curvature(H_old,q);R_new=curvature(H_new,q)
    origin={x:0 for x in q}
    check('old metric correction vanishes at origin',H_old.subs(origin))
    check('new metric correction vanishes at origin',H_new.subs(origin))
    for n,H in [('old',H_old),('new',H_new)]:
        for a in range(4):check(n+' metric first jets vanish',H.diff(q[a]).subs(origin))
    for i in range(4):
        for j in range(4):
            for k in range(4):
                for l in range(4):
                    key=(i,j,k,l)
                    delta=lambda a,b:int(a==b)
                    expression=3*(delta(i,k)*delta(j,l)-delta(i,l)*delta(j,k)) \
                        +2*J[i,j]*J[k,l]+J[i,k]*J[j,l]-J[i,l]*J[j,k]
                    check('all curvature entries match closed tensor',R_new[key]-expression)
                    check('all curvature entries survive radial map',R_new[key]-R_old[key])
                    check('curvature antisymmetry first pair',R_new[key]+R_new[j,i,k,l])
                    check('curvature antisymmetry last pair',R_new[key]+R_new[i,j,l,k])
                    check('curvature pair exchange',R_new[key]-R_new[k,l,i,j])
                    check('first Bianchi identity',R_new[key]+R_new[i,k,l,j]+R_new[i,l,j,k])
    ricci=s.Matrix(4,4,lambda j,l:sum(R_new[i,j,i,l] for i in range(4)))
    check('Ricci tensor',ricci-12*eye)
    check('scalar curvature',s.trace(ricci)-48)
    check('complex plane sectional curvature',R_new[0,1,0,1]-6)
    check('totally real plane sectional curvature',R_new[0,2,0,2]-3)
    # Symbolic arbitrary plane: no orthogonality assumption required here.
    X=s.Matrix(s.symbols('x0:4',real=True));Y=s.Matrix(s.symbols('y0:4',real=True))
    contracted=sum(R_new[i,j,k,l]*X[i]*Y[j]*X[k]*Y[l]
                   for i in range(4) for j in range(4) for k in range(4) for l in range(4))
    area2=X.dot(X)*Y.dot(Y)-X.dot(Y)**2;omega=(X.T*J*Y)[0]
    check('arbitrary plane curvature polynomial',contracted-3*(area2+omega**2))
    theta=s.symbols('theta',real=True)
    # A unit plane interpolating from totally real to a complex line.
    plane_y=s.Matrix([0,s.cos(theta),s.sin(theta),0])
    plane_x=s.Matrix([1,0,0,0])
    check('interpolating plane norm',s.simplify(plane_y.dot(plane_y)-1))
    check('interpolating plane orthogonality',plane_x.dot(plane_y))
    # trigonometric identities require simplify, while all other checks are polynomial.
    plane_value=s.simplify(contracted.subs(dict(zip(list(X)+list(Y),list(plane_x)+list(plane_y)))))
    check('interpolating sectional curvature',s.expand_trig(plane_value-3*(1+s.cos(theta)**2)))
    R_drop_current=curvature(3*q*q.T,q)
    R_bad_radial=curvature(2*q*q.T-v*v.T,q)
    R_drop_density=curvature(q*q.T-v*v.T,q)
    reject('erase residual curvature',R_new[0,2,0,2])
    reject('omit scalar current metric term',R_drop_current[0,1,0,1]-R_new[0,1,0,1])
    reject('wrong coefficient of density gradient',R_bad_radial[0,2,0,2]-R_new[0,2,0,2])
    reject('omit original radial kinetic term',R_drop_density[0,2,0,2]-R_old[0,2,0,2])
    reject('reverse curvature sign',-R_new[0,1,0,1]-6)
    reject('identify complex and totally real planes',R_new[0,1,0,1]-R_new[0,2,0,2])
    check('zero selected coefficient is flat',(h*R_new[0,1,0,1]).subs(h,0))
    report={'status':'PASS','calculation':'First-order scalar field-space metric and curvature of the selected real matter pole; no new loop integration.',
        'conventions':{'phi':'(q0+i q1,q2+i q3)/sqrt(2)','h':'kappa Re Tr(C7 Y_nu* M Y_nu^dagger)',
        'R_ijkl':'(h_il,jk+h_jk,il-h_ik,jl-h_jl,ik)/2 at origin',
        'sectional':'R(X,Y,X,Y)/(X^2 Y^2-(X dot Y)^2)'},
        'results':{'old_metric_correction':'h[q q^T-(Jq)(Jq)^T-q^2 I]',
          'new_metric_correction':'h[3q q^T-(Jq)(Jq)^T]',
          'unit_orthogonal_sectional':'3h[1+(X^T J Y)^2]',
          'complex_line_sectional':'6h','totally_real_sectional':'3h','Ricci_at_origin':'12h I','scalar_at_origin':'48h'},
        'positive_checks':len(positive),'component_equalities_evaluated':sum(positive.values()),'rejected_adverse_controls':sum(adverse.values()),
        'positive_groups':dict(positive),'adverse_groups':dict(adverse),
        'scope':'Local derivative-free scalar field reparametrization in the selected zero-external-gauge, first-order packet only; not a full EFT EOM quotient, physical amplitude, spacetime curvature or novelty result.'}
    (Path(__file__).resolve().parent/'SCALAR_GEOMETRY_RESULTS.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'positive_checks':report['positive_checks'],'rejected_adverse_controls':report['rejected_adverse_controls']}))

if __name__=='__main__':main()
