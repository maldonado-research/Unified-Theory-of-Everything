#!/usr/bin/env python3
"""Exact seven-diagram vertex/routing audit; standard library, no local imports.

This checks rational unintegrated kernels, not loop integration. Complex Dirac
matrices use an independent 8 by 8 real representation over Fraction.
"""
from fractions import Fraction as Q
from random import Random
import json

DIM = 8
ETA = (1, -1, -1, -1)
ZERO = tuple(tuple(Q(0) for _ in range(DIM)) for _ in range(DIM))
UNIT = tuple(tuple(Q(i == j) for j in range(DIM)) for i in range(DIM))


def realify(rows):
    """Entries are integers or pairs (real integer, imaginary integer)."""
    out = [[Q(0) for _ in range(DIM)] for _ in range(DIM)]
    for i in range(4):
        for j in range(4):
            x = rows[i][j]
            a, b = x if isinstance(x, tuple) else (x, 0)
            out[2*i][2*j] = out[2*i+1][2*j+1] = Q(a)
            out[2*i][2*j+1] = Q(-b)
            out[2*i+1][2*j] = Q(b)
    return tuple(map(tuple, out))


G = (
    realify([[1,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,-1]]),
    realify([[0,0,0,1],[0,0,1,0],[0,-1,0,0],[-1,0,0,0]]),
    realify([[0,0,0,(0,-1)],[0,0,(0,1),0],
             [0,(0,1),0,0],[(0,-1),0,0,0]]),
    realify([[0,0,1,0],[0,0,0,-1],[-1,0,0,0],[0,1,0,0]]),
)
IMAG = realify([[(0,1) if i == j else 0 for j in range(4)] for i in range(4)])


def add(*ms):
    return tuple(tuple(sum(m[i][j] for m in ms) for j in range(DIM))
                 for i in range(DIM))


def scale(a, m):
    return tuple(tuple(a*x for x in row) for row in m)


def mul(*ms):
    out = UNIT
    for m in ms:
        # Skip exact zeros; this changes no arithmetic or test criterion.
        out = tuple(tuple(sum(out[i][k]*m[k][j] for k in range(DIM)
                              if out[i][k] and m[k][j])
                          for j in range(DIM)) for i in range(DIM))
    return out


def transpose_complex(m):
    return tuple(tuple(m[2*(j//2)+i%2][2*(i//2)+j%2]
                       for j in range(DIM)) for i in range(DIM))


def vadd(*vs):
    return tuple(sum(v[i] for v in vs) for i in range(4))


def vsub(a, b):
    return tuple(x-y for x,y in zip(a,b))


def vs(a, v):
    return tuple(a*x for x in v)


def dot(a, b):
    return sum(s*x*y for s,x,y in zip(ETA,a,b))


def slash(a):
    return add(*(scale(s*x, g) for s,x,g in zip(ETA,a,G)))


def comm(a,b):
    return add(mul(a,b),scale(-1,mul(b,a)))


ZV = (Q(0),)*4
C = mul(IMAG,G[2],G[0])
CINV = scale(-1,C)
G5 = mul(IMAG,*G)
PL = scale(Q(1,2),add(UNIT,scale(-1,G5)))


def one_b(p,h,k,e,q,m2):
    t = vsub(p,k)
    qk = vsub(q,k)
    dp = dot(vsub(q,p),vsub(q,p))-m2
    dt = dot(vsub(q,t),vsub(q,t))-m2
    sh = 2*dot(t,t)+dot(vsub(q,h),vsub(q,h))
    sl = 2*dot(t,t)+dot(vsub(qk,h),vsub(qk,h))
    scalar = scale(sh*dot(vsub(vs(2,vsub(p,q)),k),e)/(dot(q,q)*dp*dt),slash(q))
    fermion = scale(sl/(dot(q,q)*dot(qk,qk)*dp),mul(slash(qk),slash(e),slash(q)))
    return add(scalar,scale(-1,fermion))


def diagrams(p,h,k1,k2,e1,e2,q,m2,wrong_seed=False,reverse_ff=False):
    k = vadd(k1,k2)
    t = vsub(p,k)
    s = vsub(p,q)
    def seed(a):
        b = vsub(vsub(q,h),a)
        return 2*dot(t,t)+dot(b,b)
    def fd(a):
        b = vsub(q,a)
        return dot(b,b)
    def sd(a):
        b = vsub(s,a)
        return dot(b,b)-m2
    terms = []
    for ki,kj,ei,ej in [(k1,k2,e1,e2),(k2,k1,e2,e1)]:
        chain = [slash(vsub(q,k)),slash(ej),slash(vsub(q,ki)),slash(ei),slash(q)]
        if reverse_ff:
            chain.reverse()
        terms.append(scale(seed(ZV if wrong_seed else k)/(fd(ZV)*fd(ki)*fd(k)*sd(ZV)),mul(*chain)))
    for ki,kj,ei,ej in [(k1,k2,e1,e2),(k2,k1,e2,e1)]:
        a = dot(vsub(vs(2,s),ki),ei)
        b = dot(vsub(vsub(vs(2,s),vs(2,ki)),kj),ej)
        terms.append(scale(seed(ZV)*a*b/(fd(ZV)*sd(ZV)*sd(ki)*sd(k)),slash(q)))
    for ki,kj,ei,ej in [(k1,k2,e1,e2),(k2,k1,e2,e1)]:
        a = dot(vsub(vs(2,s),kj),ej)
        terms.append(scale(seed(ki)*a/(fd(ZV)*fd(ki)*sd(ZV)*sd(kj)),
                           mul(slash(vsub(q,ki)),slash(ei),slash(q))))
    terms.append(scale(seed(ZV)*dot(e1,e2)/(fd(ZV)*sd(ZV)*sd(k)),slash(q)))
    return terms


def combine(terms,seagull=-2,mixed=-1,drop=None):
    weights = [1,1,1,1,mixed,mixed,seagull]
    if drop is not None:
        weights[drop]=0
    return add(*(scale(w,t) for w,t in zip(weights,terms)))


def ward_rhs(p,h,k1,k2,e2,q,m2):
    """Four exact integrands; the last pair is a loop-shift difference."""
    pm = vsub(p,k1)
    external = add(one_b(p,vadd(h,k1),k2,e2,q,m2),
                   scale(-1,one_b(p,h,k2,e2,q,m2)))
    shift = add(one_b(pm,h,k2,e2,q,m2),
                scale(-1,one_b(pm,h,k2,e2,vsub(q,k1),m2)))
    return external,shift


def residue_one(p,h,k,e,m2):
    r = vsub(vsub(p,h),k)
    cov = add(scale(-dot(vadd(vs(2,r),k),e)/3+dot(vadd(vs(2,h),k),e)/6,slash(p)),
              scale(dot(p,p)/6-m2/2,slash(e)))
    transverse = add(scale(Q(1,4),mul(comm(slash(k),slash(e)),slash(r))),
                     scale(Q(1,3),add(scale(dot(k,e),slash(k)),scale(-dot(k,k),slash(e)))))
    return add(cov,transverse)


def residue_two(p,k1,k2,e1,e2):
    return add(scale(dot(e1,e2),slash(p)),
               scale(Q(-1,4),add(mul(comm(slash(k1),slash(e1)),slash(e2)),
                                  mul(comm(slash(k2),slash(e2)),slash(e1)))))


def main():
    checks = {}
    def check(name,actual,expected):
        if actual != expected:
            raise RuntimeError('FAILED: '+name)
        checks[name] = checks.get(name,0)+1
    for a in range(4):
        for b in range(4):
            check('Clifford algebra',add(mul(G[a],G[b]),mul(G[b],G[a])),
                  scale(2*ETA[a] if a == b else 0,UNIT))
        check('charge conjugation gamma',mul(C,transpose_complex(G[a]),CINV),scale(-1,G[a]))
    check('chirality projector',mul(PL,PL),PL)
    check('baseline loop phase',(1j)*(1j)*(-1j)*(1j),-1)
    check('sequential insertion relative phase',1j*(-1j),1)
    check('seagull extra propagator relative phase',1j*(2j),-2)
    check('physical FF charge product',(-1)*(-1),1)
    check('physical mixed charge product',(-1)*(1),-1)
    killed = {name:0 for name in ['missing seagull','half seagull','wrong seagull sign',
                                  'wrong mixed charges','wrong FF seed','reversed FF chain',
                                  'discarding shift before integration','one missing FF ordering']}
    rng = Random(710152026)
    accepted = 0
    attempts = 0
    while accepted < 24:
        attempts += 1
        def rv():
            return tuple(Q(rng.randint(-7,7),rng.randint(1,3)) for _ in range(4))
        p,h,k1,k2,e1,e2,q = [rv() for _ in range(7)]
        m2 = Q(rng.choice([-5,-2,0,3,7]),3)
        try:
            ordinary = diagrams(p,h,k1,k2,e1,e2,q,m2)
            contracted = diagrams(p,h,k1,k2,k1,e2,q,m2)
            external,shift = ward_rhs(p,h,k1,k2,e2,q,m2)
            rhs = add(external,shift)
            swapped = diagrams(p,h,k2,k1,e2,e1,q,m2)
            contracted2 = diagrams(p,h,k2,k1,k2,e1,q,m2)
            rhs2 = add(*ward_rhs(p,h,k2,k1,e1,q,m2))
        except ZeroDivisionError:
            continue
        accepted += 1
        check('Bose identity exact integrand',combine(ordinary),combine(swapped))
        for a,b in [(0,1),(2,3),(4,5),(6,6)]:
            check('Bose diagram pairing',ordinary[a],swapped[b])
        check('first Ward identity with shift',combine(contracted),rhs)
        check('second Ward identity with shift',combine(contracted2),rhs2)
        check('first Ward identity projected',mul(combine(contracted),PL),mul(rhs,PL))
        check('scalar seagull functional factor',dot(vadd(e1,e2),vadd(e1,e2))-dot(e1,e1)-dot(e2,e2),2*dot(e1,e2))
        # With A_i=2s-2(previous scalar emissions)-k_i, scalar contraction
        # of the two ordered graphs differs from its endpoint terms by +2SG.
        sg_contracted = contracted[6]
        s = vsub(p,q)
        k = vadd(k1,k2)
        t = vsub(p,k)
        seed0 = 2*dot(t,t)+dot(vsub(q,h),vsub(q,h))
        d0 = dot(s,s)-m2
        d1 = dot(vsub(s,k1),vsub(s,k1))-m2
        d2 = dot(vsub(s,k2),vsub(s,k2))-m2
        dk = dot(vsub(s,k),vsub(s,k))-m2
        endpoint = scale(seed0/dot(q,q)*(dot(vsub(vsub(vs(2,s),vs(2,k1)),k2),e2)/(dk*d1)
                         -dot(vsub(vs(2,s),k2),e2)/(d2*d0)),slash(q))
        check('scalar ordered contraction needs factor two',
              add(contracted[2],contracted[3],scale(-1,endpoint)),scale(2,sg_contracted))
        # Ordered two-insertion charge conjugation, before loop integration.
        ordinary_chain = mul(slash(vs(-1,q)),slash(e1),slash(vsub(k1,q)),slash(e2),slash(vsub(k,q)))
        charged_chain = mul(slash(vsub(q,k)),scale(-1,slash(e2)),slash(vsub(q,k1)),scale(-1,slash(e1)),slash(q))
        check('ordered two-insertion charge conjugation',
              mul(C,transpose_complex(ordinary_chain),CINV),charged_chain)
        check('local action Ward identity',residue_two(p,k1,k2,k1,e2),
              add(residue_one(p,vadd(h,k1),k2,e2,m2),scale(-1,residue_one(p,h,k2,e2,m2))))
        check('local action Bose identity',residue_two(p,k1,k2,e1,e2),residue_two(p,k2,k1,e2,e1))
        mutations = {
            'missing seagull':combine(contracted,seagull=0),
            'half seagull':combine(contracted,seagull=-1),
            'wrong seagull sign':combine(contracted,seagull=2),
            'wrong mixed charges':combine(contracted,mixed=1),
            'wrong FF seed':combine(diagrams(p,h,k1,k2,k1,e2,q,m2,wrong_seed=True)),
            'reversed FF chain':combine(diagrams(p,h,k1,k2,k1,e2,q,m2,reverse_ff=True)),
            'discarding shift before integration':external,
            'one missing FF ordering':combine(contracted,drop=1),
        }
        for name,bad in mutations.items():
            if bad != rhs:
                killed[name] += 1
    for name,count in killed.items():
        if count == 0:
            raise RuntimeError('Negative control survived: '+name)
    print(json.dumps({'status':'passed','arithmetic':'exact Fraction; independent 8x8 realified Dirac matrices',
                      'accepted_generic_points':accepted,'attempted_points':attempts,
                      'checks':checks,'total_positive_checks':sum(checks.values()),
                      'mutations_detected_at_points':killed,
                      'scope':'seven unintegrated kernels, action phases, charge ordering, Bose and Ward identities; no loop integration'},
                     indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
