#!/usr/bin/env python3
"""Exact Gaussian-rational controls for a sourced sterile EOM identity.

Standard library only. This checks algebra, not a loop anomalous dimension.
Explicit failures remain active under python -O. No input files or network.
"""
from dataclasses import dataclass
from fractions import Fraction
import json
import random


@dataclass(frozen=True)
class Q:
    re: Fraction = Fraction(0)
    im: Fraction = Fraction(0)

    def __post_init__(self):
        object.__setattr__(self, "re", Fraction(self.re))
        object.__setattr__(self, "im", Fraction(self.im))

    def __add__(self, other):
        other = asq(other)
        return Q(self.re + other.re, self.im + other.im)

    __radd__ = __add__

    def __neg__(self):
        return Q(-self.re, -self.im)

    def __sub__(self, other):
        return self + -asq(other)

    def __rsub__(self, other):
        return asq(other) + -self

    def __mul__(self, other):
        other = asq(other)
        return Q(self.re * other.re - self.im * other.im,
                 self.re * other.im + self.im * other.re)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = asq(other)
        den = other.re * other.re + other.im * other.im
        if not den:
            raise ZeroDivisionError
        return self * Q(other.re / den, -other.im / den)

    def conj(self):
        return Q(self.re, -self.im)


def asq(value):
    return value if isinstance(value, Q) else Q(value)


def mat(rows):
    return [[asq(x) for x in row] for row in rows]


def tr(a):
    return [list(row) for row in zip(*a)]


def dag(a):
    return [[x.conj() for x in row] for row in tr(a)]


def add(*args):
    return [[sum((a[i][j] for a in args), Q())
             for j in range(len(args[0][0]))] for i in range(len(args[0]))]


def scale(c, a):
    return [[asq(c) * x for x in row] for row in a]


def mul(a, *args):
    for b in args:
        if len(a[0]) != len(b):
            raise ValueError("incompatible matrix dimensions")
        a = [[sum((a[i][k] * b[k][j] for k in range(len(b))), Q())
              for j in range(len(b[0]))] for i in range(len(a))]
    return a


def diag(values):
    return [[asq(v) if i == j else Q() for j in range(len(values))]
            for i, v in enumerate(values)]


def zero(a):
    return all(x == Q() for row in a for x in row)


counts = {}


def check(label, condition):
    if not condition:
        raise RuntimeError("check failed: " + label)
    counts[label] = counts.get(label, 0) + 1


def same(label, a, b):
    check(label, a == b)


def randmat(rng, m, n):
    return [[Q(Fraction(rng.randint(-5, 5), rng.randint(1, 4)),
               Fraction(rng.randint(-5, 5), rng.randint(1, 4)))
             for _ in range(n)] for _ in range(m)]


def main():
    rng = random.Random(20261002)
    unit2 = diag([1, 1])
    sig = [unit2, mat([[0, 1], [1, 0]]),
           mat([[0, Q(0, -1)], [Q(0, 1), 0]]), diag([1, -1])]
    sigbar = [sig[0]] + [scale(-1, s) for s in sig[1:]]
    # Contract p_mu = (p0,-p1,-p2,-p3), metric (+---).
    sigma_controls = 0
    for mu in range(4):
        for nu in range(4):
            rhs = scale(2 * ([1, -1, -1, -1][mu] if mu == nu else 0), unit2)
            same("sigma anticommutator", add(mul(sig[mu], sigbar[nu]),
                 mul(sig[nu], sigbar[mu])), rhs)
            sigma_controls += 1
    masses = [2, 5]
    m = diag(masses)
    m2 = mul(m, m)
    m3 = mul(m2, m)
    mi = diag([Fraction(1, x) for x in masses])
    for trial in range(48):
        y, f = randmat(rng, 3, 2), randmat(rng, 3, 2)
        j, jb = randmat(rng, 3, 2), randmat(rng, 3, 2)
        n, nb = randmat(rng, 2, 2), randmat(rng, 2, 2)
        p = [3 + trial % 3, 1, 1, 0]
        p2 = p[0] ** 2 - sum(x * x for x in p[1:])
        d = add(scale(p[0], sig[0]), *[scale(-p[k], sig[k]) for k in range(1, 4)])
        db = add(scale(p[0], sigbar[0]), *[scale(-p[k], sigbar[k]) for k in range(1, 4)])
        same("Dirac square", mul(d, db), scale(p2, unit2))
        # Spinors are rows within each flavor: D nbar is nbar D^T.
        e = add(mul(nb, tr(d)), scale(-1, mul(m, n)), scale(-1, mul(tr(y), j)))
        eb = add(mul(n, tr(db)), scale(-1, mul(m, nb)), scale(-1, mul(dag(y), jb)))
        source1 = mul(m2, tr(y), j)
        source2 = mul(m, dag(y), jb, tr(d))
        source3 = scale(p2, mul(tr(y), j))
        lhs = add(scale(p2, mul(nb, tr(d))), scale(-1, mul(m3, n)),
                  scale(-1, source1), scale(-1, source2), scale(-1, source3))
        rhs = add(scale(p2, e), mul(m2, e), mul(m, eb, tr(d)))
        same("off-shell sourced identity", lhs, rhs)
        # Solve both sourced EOM at nonexceptional p^2, with arbitrary sources.
        inverse = diag([Fraction(1, p2 - x * x) for x in masses])
        ns = mul(inverse, add(mul(m, tr(y), j), mul(dag(y), jb, tr(d))))
        nbs = mul(mi, add(mul(ns, tr(db)), scale(-1, mul(dag(y), jb))))
        es = add(mul(nbs, tr(d)), scale(-1, mul(m, ns)), scale(-1, mul(tr(y), j)))
        ebs = add(mul(ns, tr(db)), scale(-1, mul(m, nbs)), scale(-1, mul(dag(y), jb)))
        check("sourced solution satisfies both EOM", zero(es) and zero(ebs))
        exact = scale(p2, mul(nbs, tr(d)))
        terms = [mul(m3, ns), source1, source2, source3]
        same("identity on sourced solution", exact, add(*terms))
        check("free cubic substitution rejected", exact != terms[0])
        for index in range(1, 4):
            check("each missing source rejected", exact != add(*[t for k, t in enumerate(terms) if k != index]))
        dy = scale(-1, mul(f, m3))
        dc = add(mul(f, m2, tr(y)), mul(y, m2, tr(f)))
        ds = add(mul(dy, mi, tr(y)), mul(y, mi, tr(dy)))
        check("correlated Schur variation cancels", zero(add(dc, ds)))
        check("dropping Weinberg contact rejected", not zero(ds))
        check("conjugate transpose mutation rejected",
              not zero(add(dc, mul(dy, mi, tr(y)), mul(y, mi, dag(dy)))))
        # Exact finite triangular redefinition, including its quadratic contact.
        a = randmat(rng, 2, 3)
        c0 = randmat(rng, 3, 3)
        c = add(c0, tr(c0))
        yp = add(y, mul(tr(a), m))
        cp = add(c, scale(-1, mul(y, a)), scale(-1, mul(tr(a), tr(y))),
                 scale(-1, mul(tr(a), m, a)))
        old = add(c, mul(y, mi, tr(y)))
        new = add(cp, mul(yp, mi, tr(yp)))
        same("exact finite Schur invariance", old, new)
        check("missing quadratic contact rejected",
              old != add(new, mul(tr(a), m, a)))
    result = {
        "status": "PASS",
        "scope": "sourced EOM and basis-change algebra; no loop coefficient or physical beta function",
        "arithmetic": "exact Gaussian rationals, Python standard library",
        "flavors": {"light": 3, "sterile": 2},
        "independent_trials": 48,
        "checks": counts,
        "total_checks": sum(counts.values()),
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
