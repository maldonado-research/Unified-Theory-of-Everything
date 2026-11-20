#!/usr/bin/env python3
"""Exact controls of two declared 1PI gauge integrals, using only stdlib.

The Gaussian-rational matrix representation follows the public October 1
Clifford checker. Loop poles are reduced here by denominator cancellation
and rank-two/three moments. The optional gauge_critic.py uses Feynman-simplex moments;
verify_gauge_asymptotic.py separately uses a large-momentum expansion.
This script is import-safe and never writes files.
"""
from dataclasses import dataclass
from fractions import Fraction as F
import json
import random


@dataclass(frozen=True)
class G:
    re: F = F(0)
    im: F = F(0)

    def __post_init__(self):
        if isinstance(self.re, float) or isinstance(self.im, float):
            raise TypeError("floating point is forbidden in exact controls")
        object.__setattr__(self, "re", F(self.re))
        object.__setattr__(self, "im", F(self.im))

    @staticmethod
    def coerce(x):
        return x if isinstance(x, G) else G(x)

    def __add__(self, x):
        x = self.coerce(x)
        return G(self.re+x.re, self.im+x.im)

    __radd__ = __add__

    def __neg__(self):
        return G(-self.re, -self.im)

    def __sub__(self, x):
        return self + -self.coerce(x)

    def __rsub__(self, x):
        return self.coerce(x) + -self

    def __mul__(self, x):
        x = self.coerce(x)
        return G(self.re*x.re-self.im*x.im,
                 self.re*x.im+self.im*x.re)

    __rmul__ = __mul__

    def conj(self):
        return G(self.re, -self.im)


I = G(0, 1)
Z = G()


def mat(rows):
    return [[G.coerce(x) for x in row] for row in rows]


ONE = mat([[int(a == b) for b in range(4)] for a in range(4)])
ZERO = mat([[0]*4 for _ in range(4)])


def add(*args):
    return [[sum((a[i][j] for a in args), Z) for j in range(4)]
            for i in range(4)]


def scale(c, a):
    return [[c*x for x in row] for row in a]


def mul(a, b):
    return [[sum((a[i][k]*b[k][j] for k in range(4)), Z)
             for j in range(4)] for i in range(4)]


def trace(a):
    return sum((a[i][i] for i in range(4)), Z)


def adjoint(a):
    return [[a[j][i].conj() for j in range(4)] for i in range(4)]


gamma = [
    mat([[1,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,-1]]),
    mat([[0,0,0,1],[0,0,1,0],[0,-1,0,0],[-1,0,0,0]]),
    mat([[0,0,0,-I],[0,0,I,0],[0,I,0,0],[-I,0,0,0]]),
    mat([[0,0,1,0],[0,0,0,-1],[-1,0,0,0],[0,1,0,0]]),
]
metric = (1, -1, -1, -1)
g5 = scale(I, mul(mul(mul(gamma[0], gamma[1]), gamma[2]), gamma[3]))
PL = scale(F(1, 2), add(ONE, scale(-1, g5)))


def vadd(a, b):
    return tuple(x+y for x, y in zip(a, b))


def vsub(a, b):
    return tuple(x-y for x, y in zip(a, b))


def dot(a, b):
    return sum(metric[j]*a[j]*b[j] for j in range(4))


def slash(v):
    return add(*(scale(metric[j]*v[j], gamma[j]) for j in range(4)))


def triple(a, nu, b):
    return mul(mul(slash(a), gamma[nu]), slash(b))


def antisym(a, nu, b):
    return add(triple(a, nu, b), scale(-a[nu], slash(b)),
               scale(-b[nu], slash(a)), scale(dot(a,b), gamma[nu]))


def contract(v, tensor):
    return add(*(scale(metric[j]*v[j], tensor[j]) for j in range(4)))


def bubble(p, h, m2):
    return add(scale(dot(p,p)+m2-F(2,3)*dot(p,h)+dot(h,h)/2, slash(p)),
               scale(dot(p,p)/6-m2/2, slash(h)))


def fermion(p, h, k, m2):
    """Closed denominator-cancellation expression in the derivation."""
    t, v, s = vsub(p,k), vadd(h,k), vadd(p,k)
    c = -F(5,6)*dot(t,t)-dot(v,v)/2+dot(v,s)/3-m2/2
    return [add(scale(F(1,3), triple(t,nu,t)),
                scale(F(1,2), triple(t,nu,k)),
                scale(-F(1,6), triple(v,nu,s)),
                scale(-F(1,6), triple(s,nu,v)),
                scale(F(1,2), triple(k,nu,v)), scale(c,gamma[nu]))
            for nu in range(4)]


def scalar(p, h, k, m2):
    t = vsub(p,k)
    u = vsub(vadd(p,p),k)
    c = dot(k,k)/6-m2+dot(h,u)/3-dot(h,h)/2-dot(t,t)
    return [add(scale(-k[nu]/6, slash(k)), scale(-u[nu]/6, slash(h)),
                scale(h[nu]/3, slash(u)), scale(c,gamma[nu]))
            for nu in range(4)]


def local_representative(p, h, k, m2):
    """Public O_J one-B kernel in the common Ward normalization."""
    r = vsub(vsub(p,h),k)
    return [add(scale(-(2*r[nu]+k[nu])/3+(2*h[nu]+k[nu])/6,slash(p)),
                scale(dot(p,p)/6-m2/2,gamma[nu])) for nu in range(4)]


def transverse(p, h, k):
    r = vsub(vsub(p,h),k)
    return [add(scale(r[nu]/2,slash(k)), scale(-dot(r,k)/2,gamma[nu]),
                scale(F(1,2),antisym(k,nu,r)),
                scale(k[nu]/3,slash(k)), scale(-dot(k,k)/3,gamma[nu]))
            for nu in range(4)]


def compact_transverse(p, h, k):
    r = vsub(vsub(p,h),k)
    return [add(scale(F(1,4),mul(add(mul(slash(k),gamma[nu]),
                                    scale(-1,mul(gamma[nu],slash(k)))),slash(r))),
                scale(k[nu]/3,slash(k)), scale(-dot(k,k)/3,gamma[nu]))
            for nu in range(4)]


def t2(a, b):
    return F(metric[a],4) if a == b else F(0)


def t3(a, b, c, shiftsum):
    return F(1,12)*((metric[a]*shiftsum[c] if a == b else 0)
            +(metric[a]*shiftsum[b] if a == c else 0)
            +(metric[b]*shiftsum[a] if b == c else 0))


def fermion_moments(p, h, k, m2):
    """Direct moment contraction after q² denominator cancellation."""
    t, v, s = vsub(p,k), vadd(h,k), vadd(p,k)
    kval = dot(v,v)+2*dot(t,t)
    answer = []
    for nu in range(4):
        terms = [scale(F(1,2),triple(t,nu,k))]
        for a in range(4):
            for b in range(4):
                b2 = t[a]*t[b]/3
                if a == b:
                    b2 += metric[a]*(m2/4-dot(t,t)/12)
                coefficient = b2+kval*t2(a,b)
                coefficient -= 2*sum(metric[c]*v[c]*t3(a,b,c,s)
                                     for c in range(4))
                terms.append(scale(metric[a]*metric[b]*coefficient,
                                   mul(mul(gamma[a],gamma[nu]),gamma[b])))
                # +2 kslash gamma^nu qslash (q.v); triangle rank two.
                terms.append(scale(2*metric[a]*metric[b]*v[b]*t2(a,b),
                                   mul(mul(slash(k),gamma[nu]),gamma[a])))
        answer.append(add(*terms))
    return answer


def scalar_moments(p, h, k, m2):
    """Equal-mass bubble plus direct rank-two/three triangle moments."""
    t, u = vsub(p,k), vsub(vadd(p,p),k)
    kval = dot(h,h)+2*dot(t,t)
    answer = []
    for nu in range(4):
        terms = [scale(-k[nu]/6,slash(k)),
                 scale(dot(k,k)/6-m2,gamma[nu])]
        for a in range(4):
            coefficient = -2*kval*t2(a,nu)
            coefficient += sum(metric[c]*h[c]*
                               (-2*u[nu]*t2(a,c)+4*t3(a,nu,c,u))
                               for c in range(4))
            terms.append(scale(metric[a]*coefficient,gamma[a]))
        answer.append(add(*terms))
    return answer


def main():
    counts = {"matrix_equalities":0, "scalar_equalities":0,
              "negative_controls":0}

    def check(condition, label, kind="matrix_equalities"):
        counts[kind] += 1
        if not condition:
            raise RuntimeError(label)

    for a in range(4):
        for b in range(4):
            rhs = scale(2*metric[a] if a == b else 0,ONE)
            check(add(mul(gamma[a],gamma[b]),mul(gamma[b],gamma[a])) == rhs,
                  "Clifford relation")

    rng = random.Random(20261015)
    cases = []
    for case in range(12):
        vectors = [tuple(F(rng.randrange(-5,6),rng.randrange(1,5))
                         for _ in range(4)) for _ in range(3)]
        m2 = F(rng.randrange(-4,5),rng.randrange(1,5))
        if case == 0:
            vectors[2] = (F(0),)*4
        if case == 1:
            m2 = F(0)
        cases.append((*vectors,m2))

    for case, (p,h,k,m2) in enumerate(cases):
        rf, rs = fermion(p,h,k,m2), scalar(p,h,k,m2)
        gf, gs = fermion_moments(p,h,k,m2), scalar_moments(p,h,k,m2)
        gj, dt = local_representative(p,h,k,m2), transverse(p,h,k)
        dc = compact_transverse(p,h,k)
        for nu in range(4):
            check(rf[nu] == gf[nu], "fermion rank moments")
            check(rs[nu] == gs[nu], "scalar rank moments")
            check(add(rs[nu],scale(-1,rf[nu])) == add(gj[nu],dt[nu]),
                  "new transverse matching")
            check(dt[nu] == dc[nu], "compact transverse identity")
            massdiff = add(scalar(p,h,k,m2+1)[nu],
                           scale(-1,fermion(p,h,k,m2+1)[nu]),
                           scale(-1,local_representative(p,h,k,m2+1)[nu]))
            check(massdiff == dt[nu], "no extra transverse scalar-mass term")
        t = vsub(p,k)
        delta = 2*(dot(t,t)-dot(p,p))
        common = add(bubble(t,h,m2),scale(-delta/2,slash(p)))
        check(contract(k,rf) == add(common,scale(-1,bubble(p,vadd(h,k),m2))),
              "fermion Ward contraction including auxiliary bubble")
        check(contract(k,rs) == add(common,scale(-1,bubble(p,h,m2))),
              "scalar Ward contraction including auxiliary bubble")
        ward = add(bubble(p,vadd(h,k),m2),scale(-1,bubble(p,h,m2)))
        combined = [add(rs[n],scale(-1,rf[n])) for n in range(4)]
        check(contract(k,combined) == ward, "opposite-charge Ward contraction")
        check(contract(k,gj) == ward, "covariant representative Ward contraction")
        check(contract(k,dt) == ZERO, "transverse remainder Ward contraction")

    p,h,k,m2 = cases[-1]
    rf,rs = fermion(p,h,k,m2),scalar(p,h,k,m2)
    ward = add(bubble(p,vadd(h,k),m2),scale(-1,bubble(p,h,m2)))
    for label, bad in [
        ("omit scalar attachment",[scale(-1,x) for x in rf]),
        ("omit fermion attachment",rs),
        ("give same signed charges",[add(rs[n],rf[n]) for n in range(4)]),
    ]:
        check(contract(k,bad) != ward,label,"negative_controls")
    actual = [add(rs[n],scale(-1,rf[n])) for n in range(4)]
    gj,dt = local_representative(p,h,k,m2),transverse(p,h,k)
    check(actual != gj,"drop all transverse terms","negative_controls")
    check(contract(k,gj) == ward,
          "dropping the whole transverse remainder is invisible to Ward")
    r = vsub(vsub(p,h),k)
    check(actual != [add(gj[n],dt[n],scale(-F(1,2),antisym(k,n,r)))
                     for n in range(4)],"drop antisymmetric gamma term",
          "negative_controls")
    check(actual != [add(gj[n],dt[n],scale(-k[n]/3,slash(k)),
                        scale(dot(k,k)/3,gamma[n])) for n in range(4)],
          "drop derivative of field strength","negative_controls")

    p = tuple(map(F,(2,0,0,0)))
    r = tuple(map(F,(1,0,0,1)))
    k = (F(1,2),F(0),F(0),F(-1,2))
    h = k
    eps = tuple(map(F,(0,1,0,0)))
    check(p == vadd(vadd(r,h),k),"on-shell momentum conservation",
          "scalar_equalities")
    check(dot(r,r)==dot(h,h)==dot(k,k)==dot(k,eps)==0,
          "on-shell light kinematics","scalar_equalities")
    dep = contract(eps,transverse(p,h,k))
    expected = add(scale(dot(r,eps),slash(k)),scale(-dot(r,k),slash(eps)))
    check(mul(slash(r),dep) == mul(slash(r),expected),
          "outgoing free-lepton transverse projection")
    vertex = mul(dep,PL)
    diracadj = mul(mul(gamma[0],adjoint(vertex)),gamma[0])
    spin = trace(mul(mul(mul(mul(PL,slash(r)),vertex),
                         add(slash(p),scale(2,ONE))),diracadj))
    check(spin == G(4),"nonzero on-shell spin sum","scalar_equalities")
    print(json.dumps({"status":"PASS", "trials":len(cases),
                      "checks":counts, "on_shell_spin_sum":str(spin.re),
                      "ward_invisible_mutation":"omitting the entire transverse remainder passes Ward but fails the independently reduced loop residue",
                      "scope":"declared two 1PI integrals and local transverse residue; no physical beta function"},
                     indent=2))


if __name__ == "__main__":
    main()
