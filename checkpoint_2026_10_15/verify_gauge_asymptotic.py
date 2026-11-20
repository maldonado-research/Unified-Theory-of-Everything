#!/usr/bin/env python3
"""Independent homogeneous-UV extraction for the two declared gauge triangles.

Standard library only; Q(i) matrices and exact diagonal-metric angular moments.
No producer module is imported. Scalars named m2 are signed squared masses.
The UV angular operation uses an auxiliary IR separation: it is not evaluation
of a scaleless integral as zero. Its radial logarithm is normalized to C_UV.
"""
from fractions import Fraction as F
from math import prod
import json
import random


class Q:
    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        if isinstance(a, Q):
            self.a, self.b = a.a, a.b
        else:
            self.a, self.b = F(a), F(b)

    def __add__(self, other):
        other = Q(other)
        return Q(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self):
        return Q(-self.a, -self.b)

    def __sub__(self, other):
        return self + (-Q(other))

    def __mul__(self, other):
        other = Q(other)
        return Q(self.a * other.a - self.b * other.b,
                 self.a * other.b + self.b * other.a)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = Q(other)
        norm = other.a * other.a + other.b * other.b
        return self * Q(other.a / norm, -other.b / norm)

    def __eq__(self, other):
        other = Q(other)
        return self.a == other.a and self.b == other.b

    def conj(self):
        return Q(self.a, -self.b)

    def __bool__(self):
        return bool(self.a or self.b)


def matrix(rows):
    return tuple(tuple(Q(x) for x in row) for row in rows)


ZERO = matrix([[0] * 4 for _ in range(4)])
ID = matrix([[int(i == j) for j in range(4)] for i in range(4)])


def add(a, b):
    return tuple(tuple(a[i][j] + b[i][j] for j in range(4)) for i in range(4))


def scale(c, a):
    return tuple(tuple(c * x for x in row) for row in a)


def mul(a, b):
    return tuple(tuple(sum((a[i][k] * b[k][j] for k in range(4)), Q())
                       for j in range(4)) for i in range(4))


def dagger(a):
    return tuple(tuple(a[j][i].conj() for j in range(4)) for i in range(4))


def trace(a):
    return sum((a[i][i] for i in range(4)), Q())


I = Q(0, 1)
G = (
    matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]]),
    matrix([[0, 0, 0, 1], [0, 0, 1, 0], [0, -1, 0, 0], [-1, 0, 0, 0]]),
    matrix([[0, 0, 0, -I], [0, 0, I, 0], [0, I, 0, 0], [-I, 0, 0, 0]]),
    matrix([[0, 0, 1, 0], [0, 0, 0, -1], [-1, 0, 0, 0], [0, 1, 0, 0]]),
)
METRIC = (1, -1, -1, -1)
E0 = (0, 0, 0, 0)
UNIT = tuple(tuple(int(i == j) for i in range(4)) for j in range(4))


def dot(a, b):
    return sum(METRIC[i] * a[i] * b[i] for i in range(4))


def vadd(a, b):
    return tuple(x + y for x, y in zip(a, b))


def vsub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def slash(a):
    return sum_matrices(scale(METRIC[i] * a[i], G[i]) for i in range(4))


def sum_matrices(items):
    out = ZERO
    for item in items:
        out = add(out, item)
    return out


def padd(a, b):
    out = dict(a)
    for exp, value in b.items():
        out[exp] = add(out.get(exp, ZERO), value)
    return {exp: value for exp, value in out.items() if value != ZERO}


def pmul(a, b):
    out = {}
    for ea, va in a.items():
        for eb, vb in b.items():
            exp = tuple(x + y for x, y in zip(ea, eb))
            out[exp] = add(out.get(exp, ZERO), mul(va, vb))
    return {exp: value for exp, value in out.items() if value != ZERO}


def pscale(c, a):
    return {exp: scale(c, value) for exp, value in a.items()}


def constant(c):
    return {E0: scale(c, ID)}


def qdot(a):
    return {UNIT[i]: scale(METRIC[i] * a[i], ID) for i in range(4)}


QS = {UNIT[i]: scale(METRIC[i], G[i]) for i in range(4)}
Q2 = {tuple(2 * x for x in UNIT[i]): scale(METRIC[i], ID) for i in range(4)}


def degree(a, n):
    return {exp: value for exp, value in a.items() if sum(exp) == n}


def moment(exp):
    """<q^exp>/(q²)^(degree/2), at d=4, using full Wick pairings."""
    if any(x % 2 for x in exp):
        return F(0)
    n = sum(exp) // 2
    pairings = prod(prod(range(1, x, 2)) for x in exp)
    signs = prod(METRIC[i] ** (exp[i] // 2) for i in range(4))
    return F(signs * pairings, prod(4 + 2 * j for j in range(n)))


def angular(a):
    return sum_matrices(scale(moment(exp), value) for exp, value in a.items())


def residue(numerator, a, b, ma2, mb2):
    """Log term of N/[q² ((q-a)²-ma²) ((q-b)²-mb²)]."""
    qa, qb = qdot(a), qdot(b)
    linear = pscale(2, padd(qa, qb))
    constant2 = constant(-dot(a, a) - dot(b, b) + ma2 + mb2)
    quadratic = pscale(4, padd(padd(pmul(qa, qa), pmul(qb, qb)), pmul(qa, qb)))
    return sum_matrices((angular(degree(numerator, 2)),
                         angular(pmul(degree(numerator, 3), linear)),
                         angular(pmul(degree(numerator, 4), constant2)),
                         angular(pmul(degree(numerator, 4), quadratic))))


def triangles(p, h, k, m2, reverse=False):
    t = vsub(p, k)
    seed_f = padd(padd(Q2, pscale(-2, qdot(vadd(h, k)))),
                  constant(dot(vadd(h, k), vadd(h, k)) + 2 * dot(t, t)))
    seed_s = padd(padd(Q2, pscale(-2, qdot(h))),
                  constant(dot(h, h) + 2 * dot(t, t)))
    qmk = padd(QS, {E0: scale(-1, slash(k))})
    fs, ss = [], []
    for nu in range(4):
        left, right = (QS, qmk) if reverse else (qmk, QS)
        nf = pmul(seed_f, pmul(pmul(left, {E0: G[nu]}), right))
        current = padd(constant(2 * p[nu] - k[nu]), {UNIT[nu]: scale(-2, ID)})
        ns = pmul(seed_s, pmul(QS, current))
        fs.append(residue(nf, k, p, 0, m2))
        ss.append(residue(ns, p, t, m2, m2))
    return fs, ss


def bubble(p, h, m2):
    return add(scale(dot(p, p) + m2 - F(2, 3) * dot(p, h) + F(1, 2) * dot(h, h), slash(p)),
               scale(F(1, 6) * dot(p, p) - F(1, 2) * m2, slash(h)))


def representative(p, h, k, m2):
    return [add(scale(-F(2, 3) * p[nu] + h[nu] + F(1, 2) * k[nu], slash(p)),
                scale(F(1, 6) * dot(p, p) - F(1, 2) * m2, G[nu])) for nu in range(4)]


def transverse(p, h, k):
    r = vsub(vsub(p, h), k)
    return [add(scale(F(1, 4), mul(add(mul(slash(k), G[nu]), scale(-1, mul(G[nu], slash(k)))), slash(r))),
                scale(F(1, 3), add(scale(k[nu], slash(k)), scale(-dot(k, k), G[nu])))) for nu in range(4)]


def contract(k, mats):
    return sum_matrices(scale(METRIC[i] * k[i], mats[i]) for i in range(4))


def check(condition, label):
    if not condition:
        raise RuntimeError(label)


def main():
    check(all(add(mul(G[i], G[j]), mul(G[j], G[i])) == scale(2 * METRIC[i] * int(i == j), ID)
              for i in range(4) for j in range(4)), "Clifford algebra")
    check(angular(Q2) == ID, "quadratic angular trace")
    check(angular(pmul(Q2, Q2)) == ID, "quartic angular trace")
    check(angular(pmul(pmul(Q2, Q2), Q2)) == ID, "sextic angular trace")
    rng = random.Random(20261015)
    samples, reverse_distinct, ward_blind = 0, 0, 0
    omitted_attachment_failures, dropped_transverse_failures = 0, 0
    while samples < 24:
        p, h, k = [tuple(F(rng.randint(-4, 4)) for _ in range(4)) for _ in range(3)]
        r = vsub(vsub(p, h), k)
        if any(dot(v, v) == 0 for v in (p, h, k, r, vsub(p, k))):
            continue
        m2 = F(rng.randint(-5, 5))
        fs, ss = triangles(p, h, k, m2)
        total = [add(ss[i], scale(-1, fs[i])) for i in range(4)]
        gj, dt = representative(p, h, k, m2), transverse(p, h, k)
        ward = add(bubble(p, vadd(h, k), m2), scale(-1, bubble(p, h, m2)))
        check(total == [add(gj[i], dt[i]) for i in range(4)], "independent UV tensor")
        check(contract(k, total) == ward, "Ward identity")
        check(contract(k, dt) == ZERO, "transversality")
        check(contract(k, gj) == ward, "representative Ward identity")
        omitted_attachment_failures += contract(k, ss) != ward
        dropped_transverse_failures += total != gj
        if samples < 4:
            flipped, _ = triangles(p, h, k, m2, reverse=True)
            changed = [add(ss[i], scale(-1, flipped[i])) for i in range(4)]
            reverse_distinct += changed != total
            ward_blind += contract(k, changed) == ward
        samples += 1
    check(omitted_attachment_failures == samples, "omitted attachment mutation")
    check(dropped_transverse_failures == samples, "dropped transverse mutation")
    check(reverse_distinct == 4 and ward_blind == 4, "Ward-blind spinor reversal control")
    p, r = tuple(map(F, (2, 0, 0, 0))), tuple(map(F, (1, 0, 0, 1)))
    h = k = (F(1, 2), F(0), F(0), F(-1, 2))
    eps = tuple(map(F, (0, 1, 0, 0)))
    check(p == vadd(vadd(r, h), k), "on-shell witness routing")
    check(dot(r, r) == dot(h, h) == dot(k, k) == dot(k, eps) == 0, "on-shell witness kinematics")
    gamma5 = scale(I, mul(mul(mul(G[0], G[1]), G[2]), G[3]))
    pl, pr = scale(F(1, 2), add(ID, scale(-1, gamma5))), scale(F(1, 2), add(ID, gamma5))
    de = contract(eps, transverse(p, h, k))
    projected = add(scale(dot(r, eps), slash(k)), scale(-dot(r, k), slash(eps)))
    check(mul(slash(r), add(de, scale(-1, projected))) == ZERO, "free-lepton projection")
    vertex = mul(de, pl)
    adjoint = mul(mul(G[0], dagger(vertex)), G[0])
    spin_sum = trace(mul(mul(mul(slash(r), vertex), add(slash(p), scale(2, ID))), adjoint))
    check(spin_sum == 4, "nonzero on-shell trace")
    print(json.dumps({"status": "passed", "arithmetic": "exact Q(i), stdlib Fraction",
                      "independent_generic_kinematics": samples,
                      "matrix_components_compared": 16 * 4 * samples,
                      "omitted_attachment_mutations_detected": omitted_attachment_failures,
                      "dropped_transverse_mutations_detected": dropped_transverse_failures,
                      "spinor_reversal_changes_transverse_but_passes_Ward": ward_blind,
                      "on_shell_spin_sum": str(spin_sum.a),
                      "scope": "declared reduced integral UV residues, not physical beta"}, indent=2))


if __name__ == "__main__":
    main()
