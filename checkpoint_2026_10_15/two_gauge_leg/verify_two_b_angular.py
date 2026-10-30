#!/usr/bin/env python3
"""Independent exact large-q extraction of the seven two-B one-loop graphs.

Only Python's standard library is used. No producer or repository module is
imported. Clifford blades implement Cl(1,3) over Q; UV extraction uses the
homogeneous radial logarithm and exact d=4 angular moments, with an implicit
infrared separator. This never identifies an unseparated scaleless integral
with a UV pole. All returned poles strip i/(16*pi**2*epsilon).
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import prod
import argparse
import json
import random

ETA = (1, -1, -1, -1)
Z = (0, 0, 0, 0)
E = tuple(tuple(int(i == j) for i in range(4)) for j in range(4))


def va(*vectors):
    return tuple(sum(v[i] for v in vectors) for i in range(4))


def vs(a, b):
    return tuple(a[i] - b[i] for i in range(4))


def vm(c, a):
    return tuple(c * x for x in a)


def dot(a, b):
    return sum(ETA[i] * a[i] * b[i] for i in range(4))


@lru_cache(None)
def blade(a, b):
    """Ordered gamma^0...gamma^3 blades; return product sign and mask."""
    inversions = sum((b & ((1 << i) - 1)).bit_count()
                     for i in range(4) if a & (1 << i))
    sign = (-1) ** inversions
    sign *= prod(ETA[i] for i in range(4) if a & b & (1 << i))
    return sign, a ^ b


def clean(a):
    return {k: v for k, v in a.items() if v}


def add(*items):
    out = {}
    for item in items:
        for k, v in item.items():
            out[k] = out.get(k, 0) + v
    return clean(out)


def scale(c, a):
    return clean({k: c * v for k, v in a.items()})


def mul(a, b):
    out = {}
    for ia, ca in a.items():
        for ib, cb in b.items():
            sign, k = blade(ia, ib)
            out[k] = out.get(k, 0) + sign * ca * cb
    return clean(out)


def chain(*items):
    out = {0: 1}
    for item in items:
        out = mul(out, item)
    return out


def slash(v):
    return clean({1 << i: ETA[i] * v[i] for i in range(4)})


def comm(a, b):
    return add(mul(a, b), scale(-1, mul(b, a)))


def pc(a):
    return {(Z, mask): value for mask, value in a.items()}


def qdot(v):
    return clean({(E[i], 0): ETA[i] * v[i] for i in range(4)})


QS = {(E[i], 1 << i): ETA[i] for i in range(4)}
Q2 = {(vm(2, E[i]), 0): ETA[i] for i in range(4)}


def pmul(a, b):
    out = {}
    for (ea, ia), ca in a.items():
        for (eb, ib), cb in b.items():
            sign, mask = blade(ia, ib)
            key = (va(ea, eb), mask)
            out[key] = out.get(key, 0) + sign * ca * cb
    return clean(out)


def pchain(*items):
    out = {(Z, 0): 1}
    for item in items:
        out = pmul(out, item)
    return out


@lru_cache(None)
def moment(exp):
    """Normalized Minkowski angular moment from Wick pairings at d=4."""
    if any(x % 2 for x in exp):
        return F(0)
    n = sum(exp) // 2
    pairings = prod(prod(range(1, x, 2)) for x in exp)
    signs = prod(ETA[i] ** (exp[i] // 2) for i in range(4))
    return F(signs * pairings, prod(4 + 2 * j for j in range(n)))


def angular(poly):
    out = {}
    for (exp, mask), value in poly.items():
        out[mask] = out.get(mask, 0) + value * moment(exp)
    return clean(out)


def jmul(a, b):
    """Top two homogeneous numerator degrees, descending in q."""
    return pmul(a[0], b[0]), add(pmul(a[0], b[1]), pmul(a[1], b[0]))


def jchain(*items):
    out = ({(Z, 0): 1}, {})
    for item in items:
        out = jmul(out, item)
    return out


def seed_jet(v):
    # Full seed is q^2 - 2 q.v + v^2 + 2 t^2. The constant is two
    # degrees below its top, and hence below the logarithmic UV threshold.
    return Q2, scale(-2, qdot(v))


def fprop_jet(shift):
    return QS, pc(scale(-1, slash(shift)))


def scalar_vertex_jet(u, eps):
    return scale(-2, qdot(eps)), pc({0: dot(u, eps)})


def vertex_jet(eps):
    return pc(slash(eps)), {}


def pole(jet, shifts, denominators):
    # Box: N4/(q^2)^4 + 2(q.sum shifts)N5/(q^2)^5.
    # Triangle: N2/(q^2)^3 + 2(q.sum shifts)N3/(q^2)^4.
    expected = 2 * denominators - 3
    if any(sum(exp) != expected for exp, _ in jet[0]):
        raise RuntimeError("unexpected leading numerator degree")
    if any(sum(exp) != expected - 1 for exp, _ in jet[1]):
        raise RuntimeError("unexpected subleading numerator degree")
    return add(angular(jet[1]),
               angular(pmul(jet[0], scale(2, qdot(va(*shifts))))))


NAMES = ("FF12", "FF21", "SS12", "SS21", "FS12", "FS21", "SG")
WEIGHTS = (1, 1, 1, 1, -1, -1, -2)


def diagrams(p, h, k1, k2, e1, e2):
    K = va(k1, k2)
    t = vs(p, K)
    out = []
    # Index j labels the B emitted first, adjacent to the Yukawa vertex.
    for kj, ki, ej, ei in ((k1, k2, e1, e2), (k2, k1, e2, e1)):
        nj = jchain(seed_jet(va(h, K)), fprop_jet(K), vertex_jet(ei),
                    fprop_jet(kj), vertex_jet(ej), fprop_jet(Z))
        out.append(pole(nj, (Z, kj, K, p), 4))
    for kj, ki, ej, ei in ((k1, k2, e1, e2), (k2, k1, e2, e1)):
        u, w = vs(vm(2, p), kj), vs(vs(vm(2, p), vm(2, kj)), ki)
        nj = jchain(seed_jet(h), fprop_jet(Z),
                    scalar_vertex_jet(u, ej), scalar_vertex_jet(w, ei))
        out.append(pole(nj, (Z, p, vs(p, kj), t), 4))
    for ki, kj, ei, ej in ((k1, k2, e1, e2), (k2, k1, e2, e1)):
        nj = jchain(seed_jet(va(h, ki)), fprop_jet(ki), vertex_jet(ei),
                    fprop_jet(Z), scalar_vertex_jet(vs(vm(2, p), kj), ej))
        out.append(pole(nj, (Z, ki, p, vs(p, kj)), 4))
    # SG includes eps1.eps2 but not its separate phase/identical-field -2.
    nj = jchain(seed_jet(h), fprop_jet(Z))
    out.append(scale(dot(e1, e2), pole(nj, (Z, p, t), 3)))
    return out


def total(parts, weights=WEIGHTS):
    return add(*(scale(w, part) for w, part in zip(weights, parts)))


def closed_diagrams(p, h, k1, k2, e1, e2):
    """Individual analytic angular results, separately from target action."""
    out = []
    K = va(k1, k2)
    for b, a, eb, ea in ((k1, k2, e1, e2), (k2, k1, e2, e1)):
        d = vs(va(b, p), h)
        out.append(add(scale(F(1, 2), chain(slash(K), slash(ea), slash(eb))),
                       scale(F(1, 2), chain(slash(eb), slash(b), slash(ea))),
                       scale(-F(1, 6), add(chain(slash(d), slash(ea), slash(eb)),
                                           chain(slash(eb), slash(d), slash(ea)),
                                           chain(slash(ea), slash(eb), slash(d))))))
    for b, a, eb, ea in ((k1, k2, e1, e2), (k2, k1, e2, e1)):
        u = vs(vm(2, p), b)
        w = vs(vs(vm(2, p), vm(2, b)), a)
        d = vs(vs(vs(vm(3, p), vm(2, b)), a), h)
        out.append(add(scale(-F(1, 2) * dot(u, eb), slash(ea)),
                       scale(-F(1, 2) * dot(w, ea), slash(eb)),
                       scale(F(1, 3), add(scale(dot(ea, eb), slash(d)),
                                           scale(dot(d, ea), slash(eb)),
                                           scale(dot(d, eb), slash(ea))))))
    for a, b, ea, eb in ((k1, k2, e1, e2), (k2, k1, e2, e1)):
        u, d = vs(vm(2, p), b), va(h, vm(-2, p), b)
        out.append(add(scale(-F(1, 2) * dot(u, eb), slash(ea)),
                       scale(F(1, 2), chain(slash(a), slash(ea), slash(eb))),
                       scale(F(1, 3), add(scale(dot(ea, eb), slash(d)),
                                           scale(dot(ea, d), slash(eb)),
                                           scale(-2 * dot(eb, d), slash(ea))))))
    out.append(scale(F(1, 2) * dot(e1, e2), slash(va(vm(2, p), vm(-1, K), vm(-1, h)))))
    return out


def action(p, k1, k2, e1, e2):
    return add(scale(dot(e1, e2), slash(p)),
               scale(-F(1, 4), add(mul(comm(slash(k1), slash(e1)), slash(e2)),
                                   mul(comm(slash(k2), slash(e2)), slash(e1)))))


def one_b(p, h, k, eps, m2):
    """Previously audited one-B pole; used only for longitudinal controls."""
    r = vs(vs(p, h), k)
    gj = add(scale(-F(1, 3) * dot(va(vm(2, r), k), eps)
                   + F(1, 6) * dot(va(vm(2, h), k), eps), slash(p)),
             scale((dot(p, p) / 6 - m2 / 2), slash(eps)))
    dt = add(scale(F(1, 4), mul(comm(slash(k), slash(eps)), slash(r))),
             scale(F(1, 3), add(scale(dot(k, eps), slash(k)),
                                 scale(-dot(k, k), slash(eps)))))
    return add(gj, dt)


CHECKS = 0


def check(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


def encode(a):
    return {str(k): str(v) for k, v in sorted(a.items())}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=24)
    parser.add_argument("--skip-basis", action="store_true")
    args = parser.parse_args()
    for a, b in product(range(4), repeat=2):
        check(add(mul(slash(E[a]), slash(E[b])), mul(slash(E[b]), slash(E[a])))
              == clean({0: 2 * ETA[a] * int(a == b)}), "Clifford generators")
    for n in (1, 2, 3):
        check(angular(pchain(*(Q2 for _ in range(n)))) == {0: 1}, "angular trace")
    # Scalar seagull: +i times quadratic kinetic action, two field derivatives.
    # Relative replacement of one propagator is (i)*(2i)=-2; ordinary i*(-i)=+1.
    check(1j * (2j) == -2 and 1j * (-1j) == 1, "relative vertex phases")
    check((1j) * (1j) * (-1j) * (1j) == -1, "zero-B source phase")
    rng = random.Random(0x2B20261015)
    mutations = {name: 0 for name in ("wrong_seagull_sign", "drop_FF21", "drop_SG",
                                      "same_internal_charges", "omit_sigma_contact",
                                      "opposite_counterterm_sign")}
    witness = None
    for sample in range(args.samples):
        p, h, k1, k2, e1, e2 = [tuple(F(rng.randint(-7, 7), rng.randint(1, 4))
                                      for _ in range(4)) for _ in range(6)]
        parts = diagrams(p, h, k1, k2, e1, e2)
        analytic = closed_diagrams(p, h, k1, k2, e1, e2)
        for name, raw, reduced in zip(NAMES, parts, analytic):
            check(raw == reduced, name + " individual analytic residue")
        raw, target = total(parts), action(p, k1, k2, e1, e2)
        check(raw == target, "full unprojected two-B action match")
        check(total(diagrams(p, h, k2, k1, e2, e1)) == raw, "Bose symmetry")
        for ki, kj, ej in ((k1, k2, e2), (k2, k1, e1)):
            contracted = total(diagrams(p, h, ki, kj, ki, ej))
            for m2 in (F(-5, 3), F(0), F(11, 7)):
                ward = add(one_b(p, va(h, ki), kj, ej, m2),
                           scale(-1, one_b(p, h, kj, ej, m2)))
                check(contracted == ward, "one-leg Ward descent")
        double = total(diagrams(p, h, k1, k2, k1, k2))
        check(double == scale(dot(k1, k2), slash(p)), "double Ward contraction")
        check(add(scale(-1, raw), target) == {}, "literal loop plus counterterm")
        mutations["wrong_seagull_sign"] += total(parts, (1, 1, 1, 1, -1, -1, 2)) != target
        mutations["drop_FF21"] += total(parts, (1, 0, 1, 1, -1, -1, -2)) != target
        mutations["drop_SG"] += total(parts, (1, 1, 1, 1, -1, -1, 0)) != target
        mutations["same_internal_charges"] += total(parts, (1, 1, 1, 1, 1, 1, -2)) != target
        mutations["omit_sigma_contact"] += raw != scale(dot(e1, e2), slash(p))
        mutations["opposite_counterterm_sign"] += add(scale(-1, raw), scale(-1, target)) != {}
        if witness is None:
            witness = {"inputs": {n: [str(x) for x in v] for n, v in
                                    zip(("p", "h", "k1", "k2", "eps1", "eps2"),
                                        (p, h, k1, k2, e1, e2))},
                       "diagrams": {n: encode(v) for n, v in zip(NAMES, parts)},
                       "sum": encode(raw), "action": encode(target)}
    for name, detections in mutations.items():
        check(detections == args.samples, name + " detected at every generic point")
    # The residue is degree one jointly in p,h,k1,k2 and bilinear in eps1,eps2.
    # Checking all 4*4*4*4=256 basis coefficients proves the local polynomial
    # identity, independent of the random rational checks above.
    basis = 0
    if not args.skip_basis:
        for slot, component, a, b in product(range(4), repeat=4):
            momenta = [Z, Z, Z, Z]
            momenta[slot] = E[component]
            p, h, k1, k2 = momenta
            parts = diagrams(p, h, k1, k2, E[a], E[b])
            check(total(parts) == action(p, k1, k2, E[a], E[b]), "complete coefficient basis")
            for raw, reduced in zip(parts, closed_diagrams(p, h, k1, k2, E[a], E[b])):
                check(raw == reduced, "individual complete coefficient basis")
            basis += 1
    print(json.dumps({"status": "PASS", "route": "homogeneous-large-q angular extraction",
                      "arithmetic": "exact rational Clifford blades Cl(1,3)",
                      "repository_imports": False, "checks": CHECKS,
                      "generic_rational_samples": args.samples,
                      "complete_polynomial_basis_points": basis,
                      "diagrams_per_evaluation": 7, "mutations_detected": mutations,
                      "result": "R_BB = eps1.eps2 slash(p) - ([slash(k1),slash(eps1)]slash(eps2) + [slash(k2),slash(eps2)]slash(eps1))/4",
                      "full_loop_phase": "-i F (gprime*yL)^2/(16*pi^2*epsilon)",
                      "scope": "specified single-seed one-loop 1PI sector only",
                      "first_exact_witness": witness}, indent=2))


if __name__ == "__main__":
    main()
