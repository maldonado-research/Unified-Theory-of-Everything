"""Exact rational checks for ONE declared derivative-seed mixed bubble.

This is not a complete dim-7 matching/RGE calculation. The analytic DR pole
and rotational-average identities are stated in ROUTED_BUBBLE.md.
No external packages; ordinary and optimized Python must both enforce checks.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import random

OUT=Path(__file__).resolve().parent
checks=[]


def check(name, ok):
    if not ok:
        raise RuntimeError(name)
    checks.append(name)


def integral_polynomial(cs):
    return sum((F(c)/F(j+1) for j,c in enumerate(cs)), F(0))


def dot(a,b):
    return a[0]*b[0]-sum(x*y for x,y in zip(a[1:],b[1:]))


def add(a,b):
    return [x+y for x,y in zip(a,b)]


def mul(c,a):
    return [c*x for x in a]


def pole_tensor_reduction(p,h,m2):
    p2,ph,h2=dot(p,p),dot(p,h),dot(h,h)
    # UV ∫l_mu/D=p_mu/2 and ∫l_mu l_nu/D=
    # (m2/4-p2/12)g_mu_nu+p_mu p_nu/3.
    I1=mul(F(1,2),p)
    I2h=add(mul(m2/F(4)-p2/F(12),h),mul(ph/F(3),p))
    # Cancel l² before integration; the shifted massive tadpole is m2*p.
    return add(add(mul(m2,p),mul(-2,I2h)),mul(2*p2+h2,I1))


def pole_uv_angular(p,h,m2):
    p2,ph,h2=dot(p,p),dot(p,h),dot(h,h)
    # Independently isolate degree -4 of the large-loop integrand.
    # <n_mu n_nu>=g/4; <n_mu n_nu n_rho n_sigma>=sym(g*g)/24.
    cubic_3=mul(F(8*3,24)*p2,p)
    cubic_2=mul(F(4,4)*(m2-p2),p)
    quadratic_2=mul(F(-8,24),add(mul(p2,h),mul(2*ph,p)))
    quadratic_1=mul(F(-2,4)*(m2-p2),h)
    linear_1=mul(F(2,4)*(2*p2+h2),p)
    total=[F(0)]*4
    for term in [cubic_3,cubic_2,quadratic_2,quadratic_1,linear_1]:
        total=add(total,term)
    return total


def pole_closed(p,h,m2):
    p2,ph,h2=dot(p,p),dot(p,h),dot(h,h)
    return add(mul(p2+m2-F(2,3)*ph+h2/2,p),mul(p2/6-m2/2,h))


def main():
    check('Feynman x integral',integral_polynomial([0,1])==F(1,2))
    check('Feynman x squared integral',integral_polynomial([0,0,1])==F(1,3))
    check('Feynman Delta metric mass term',integral_polynomial([0,F(1,2)])==F(1,4))
    check('Feynman Delta metric momentum term',integral_polynomial([0,F(-1,2),F(1,2)])==F(-1,12))
    eps=[[0,1],[-1,0]]
    for a in range(2):
        for c in range(2):
            crossed=sum(eps[a][d]*eps[b][c]*eps[b][d] for b in range(2) for d in range(2))
            direct=sum(eps[a][c]*eps[b][d]*eps[b][d] for b in range(2) for d in range(2))
            check('SU2 crossed %d %d'%(a,c),crossed==eps[a][c])
            check('SU2 direct %d %d'%(a,c),direct==2*eps[a][c])
    rng=random.Random(20260906)
    for j in range(400):
        p=[F(rng.randint(-9,9),rng.randint(1,5)) for _ in range(4)]
        h=[F(rng.randint(-9,9),rng.randint(1,5)) for _ in range(4)]
        m2=F(rng.randint(0,7),rng.randint(1,5))
        a=pole_tensor_reduction(p,h,m2)
        b=pole_uv_angular(p,h,m2)
        c=pole_closed(p,h,m2)
        check('exact independent residue routes %d'%j,a==b==c)
        # Simultaneously reversing the two physical vectors is a routing parity control.
        check('momentum reversal %d'%j,pole_closed(mul(-1,p),mul(-1,h),m2)==mul(-1,c))
    p=[F(3),F(0),F(0),F(0)]
    h=[F(1),F(1),F(0),F(0)]
    actual=pole_closed(p,h,F(0))
    naive=mul(F(3,2)*dot(p,p),p)
    check('off-shell nonparallel counterexample',actual==[F(45,2),F(3,2),F(0),F(0)] and actual!=naive)
    # In the free massless light-leg projection: h²=r²=0,p=r+h,
    # bar u_r slash r=0.  The p-slash coefficient is 5*p²/6.
    projected=F(1)-F(1,3)+F(1,6)
    check('massless projected kinematic coefficient',projected==F(5,6))
    check('projection differs from naive p2 multiplication',projected/F(3,2)==F(5,9))
    # Do not infer any beta-function replacement from this scalar projection.
    result={'status':'EXACT_RATIONAL_CHECKS_PASS','check_count':len(checks),
        'scope':'UV residue of declared one-W1 derivative seed mixed bubble, before complete sourced EOM/gauge completion',
        'pole_removed_common_factor':'i/(16*pi^2*epsilon), d=4-2epsilon',
        'basis':['p2*p_mu','(p.h)*p_mu','h2*p_mu','p2*h_mu','mH2*p_mu','mH2*h_mu'],
        'coefficients':['1','-2/3','1/2','1/6','1','-1/2'],
        'off_shell_counterexample':{'p':[str(x) for x in p],'h':[str(x) for x in h],
            'mH2':'0','actual':[str(x) for x in actual],'naive':[str(x) for x in naive]},
        'light_massless_free_EOM_projection':{'coefficient_of_p2_slash_p':'5/6','naive_coefficient':'3/2','ratio':'5/9','not_a_beta_function':True},
        'not_completed':['full dimension-seven Green-basis matching','sourced EOM and all induced operators','full physical beta functions','finite thresholds','phenomenological prediction']}
    (OUT/'ROUTED_POLE_CHECKS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
