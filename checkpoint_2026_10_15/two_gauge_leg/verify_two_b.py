#!/usr/bin/env python3
"""Exact two-external-B UV residues from seven literal routed diagrams.

Standard library only; no repository modules are imported. The engine computes
Feynman-simplex rank-2/3 triangle and rank-4/5 box moments in (+---). Those
moments act on the unreduced loop polynomial, independently of the compact
diagram formulas checked below. All checks remain active under python -O.
"""
from fractions import Fraction as F
from itertools import product
from random import Random
from math import factorial
import argparse
import json

ETA = (1, -1, -1, -1)
ZERO = (F(0),)*4
E0 = (0, 0, 0, 0)


def vs(*vectors):
    return tuple(sum(v[i] for v in vectors) for i in range(4))


def scale(c, v):
    return tuple(c*x for x in v)


def dot(a, b):
    return sum(ETA[i]*a[i]*b[i] for i in range(4))


def blade_product(a, b):
    """Ordered Clifford basis gamma^0,...,gamma^3; {gamma^a,gamma^b}=2g^ab."""
    sign = 1
    for i in range(4):
        if a >> i & 1:
            if (b & ((1 << i)-1)).bit_count() % 2:
                sign = -sign
            if b >> i & 1:
                sign *= ETA[i]
    return a ^ b, sign


MUL = {(a, b): blade_product(a, b) for a in range(16) for b in range(16)}


def plus(*polys):
    out = {}
    for poly in polys:
        for key, val in poly.items():
            out[key] = out.get(key, F(0))+val
    return {key: val for key, val in out.items() if val}


def times(c, poly):
    return {key: c*val for key, val in poly.items() if c*val}


def mul(*polys):
    out = {(E0, 0): F(1)}
    for poly in polys:
        nxt = {}
        for (a, am), av in out.items():
            for (b, bm), bv in poly.items():
                mask, sign = MUL[am, bm]
                key = (tuple(x+y for x, y in zip(a, b)), mask)
                nxt[key] = nxt.get(key, F(0))+sign*av*bv
        out = {key: val for key, val in nxt.items() if val}
    return out


def const(c):
    return {(E0, 0): F(c)} if c else {}


def slash(v):
    return {(E0, 1 << i): F(ETA[i]*v[i]) for i in range(4) if v[i]}


QCOMP = [{(tuple(int(i == j) for j in range(4)), 0): F(1)} for i in range(4)]
QSLASH = plus(*(mul(QCOMP[i], slash(tuple(int(i == j) for j in range(4))))
                for i in range(4)))


def qdot(v):
    return plus(*(times(ETA[i]*v[i], QCOMP[i]) for i in range(4)))


def square_shift(v):
    return plus(*(times(ETA[i], mul(plus(QCOMP[i], const(-v[i])),
                                    plus(QCOMP[i], const(-v[i]))))
                  for i in range(4)))


def qslash_shift(v):
    return plus(QSLASH, times(-1, slash(v)))


def pair_metric(exponents):
    if any(x % 2 for x in exponents):
        return F(0)
    out = F(1)
    for i, x in enumerate(exponents):
        out *= ETA[i]**(x//2)
        for n in range(1, x, 2):
            out *= n
    return out


def simplex_moment(exponents, shifts):
    """Only potentially divergent ranks in these diagrams; no angular method.

    For N=3, rank 2 = g/4 and rank 3 = sym(g S)/12.
    For N=4, rank 4 = sym(gg)/24 and rank 5 = sym(gg S)/96.
    S=sum(shifts). Coefficients follow from simplex volume 1/(N-1)!
    and first moment 1/N!, with the isotropic 1/2^j factor.
    """
    n, rank = len(shifts), sum(exponents)
    minimum = 2*n-4
    if rank < minimum:
        return F(0)
    if rank not in (minimum, minimum+1) or n not in (3, 4):
        raise ValueError((n, rank))
    pairs = n-2
    if rank == minimum:
        return pair_metric(exponents)/F(2**pairs*factorial(n-1))
    total = vs(*shifts)
    out = F(0)
    for i, count in enumerate(exponents):
        if count:
            reduced = list(exponents)
            reduced[i] -= 1
            out += count*total[i]*pair_metric(tuple(reduced))
    return out/F(2**pairs*factorial(n))


def residue(poly, shifts):
    out = {}
    for (powers, mask), coefficient in poly.items():
        val = coefficient*simplex_moment(powers, shifts)
        key = (E0, mask)
        out[key] = out.get(key, F(0))+val
    return {key: val for key, val in out.items() if val}


def diagrams(p, h, k1, k2, e1, e2):
    """Seven independently routed integrands, before physical charge weights."""
    K = vs(k1, k2)
    t = vs(p, scale(-1, K))
    out = {}
    for i, ki, kj, ei, ej in [(1, k1, k2, e1, e2), (2, k2, k1, e2, e1)]:
        seed_ff = plus(const(2*dot(t, t)), square_shift(vs(h, K)))
        chain_ff = mul(qslash_shift(K), slash(ej), qslash_shift(ki), slash(ei), QSLASH)
        out[f'FF{i}'] = residue(mul(seed_ff, chain_ff), [ZERO, ki, K, p])
        seed_ss = plus(const(2*dot(t, t)), square_shift(h))
        u = vs(scale(2, p), scale(-1, ki))
        v = vs(scale(2, p), scale(-2, ki), scale(-1, kj))
        v1 = plus(const(dot(u, ei)), times(-2, qdot(ei)))
        v2 = plus(const(dot(v, ej)), times(-2, qdot(ej)))
        out[f'SS{i}'] = residue(mul(seed_ss, QSLASH, v1, v2),
                                [ZERO, p, vs(p, scale(-1, ki)), t])
        seed_fs = plus(const(2*dot(t, t)), square_shift(vs(h, ki)))
        uj = vs(scale(2, p), scale(-1, kj))
        vj = plus(const(dot(uj, ej)), times(-2, qdot(ej)))
        out[f'FS{i}'] = residue(mul(seed_fs, qslash_shift(ki), slash(ei), QSLASH, vj),
                                [ZERO, ki, p, vs(p, scale(-1, kj))])
    seed_sg = plus(const(2*dot(t, t)), square_shift(h))
    out['SG'] = residue(mul(seed_sg, QSLASH), [ZERO, p, t])
    return out


def compact_diagrams(p, h, k1, k2, e1, e2):
    """Analytic denominator-cancellation reduction, independent of raw moments."""
    K = vs(k1, k2)
    out = {}
    for i, ki, kj, ei, ej in [(1, k1, k2, e1, e2), (2, k2, k1, e2, e1)]:
        Ei, Ej, Ki = slash(ei), slash(ej), slash(ki)
        w = vs(p, ki, scale(-1, h))
        W = slash(w)
        out[f'FF{i}'] = plus(
            times(F(1, 2), plus(mul(slash(K), Ej, Ei), mul(Ei, Ki, Ej))),
            times(F(-1, 6), plus(mul(W, Ej, Ei), mul(Ei, W, Ej), mul(Ej, Ei, W))))
        w = vs(scale(3, p), scale(-2, ki), scale(-1, kj), scale(-1, h))
        u = vs(scale(2, p), scale(-1, ki))
        v = vs(scale(2, p), scale(-2, ki), scale(-1, kj))
        out[f'SS{i}'] = plus(
            times(F(1, 3), plus(times(dot(ei, ej), slash(w)), times(dot(w, ei), Ej),
                               times(dot(w, ej), Ei))),
            times(F(-1, 2), plus(times(dot(u, ei), Ej), times(dot(v, ej), Ei))))
        u = vs(scale(2, p), scale(-1, kj))
        w = vs(u, scale(-1, h))
        W = slash(w)
        out[f'FS{i}'] = plus(times(F(1, 2), mul(Ki, Ei, Ej)),
                             times(F(-1, 6), plus(mul(Ej, Ei, W), mul(W, Ei, Ej))),
                             times(dot(w, ej)/3-dot(u, ej)/2, Ei))
    out['SG'] = times(F(1, 2), slash(vs(scale(2, p), scale(-1, K), scale(-1, h))))
    return out


def combine(residues, e1, e2):
    return plus(residues['FF1'], residues['FF2'], residues['SS1'], residues['SS2'],
                times(-1, residues['FS1']), times(-1, residues['FS2']),
                times(-2*dot(e1, e2), residues['SG']))


def prediction(p, h, k1, k2, e1, e2):
    E1, E2, K1, K2 = slash(e1), slash(e2), slash(k1), slash(k2)
    return plus(times(dot(e1, e2), slash(p)),
                times(F(-1, 4), plus(mul(plus(mul(K1, E1), times(-1, mul(E1, K1))), E2),
                                      mul(plus(mul(K2, E2), times(-1, mul(E2, K2))), E1))))


def serialize(poly):
    return {str(mask): str(value) for (_, mask), value in sorted(poly.items())}


def run_checks(generic_count=12):
    counts = {}
    def equal(a, b, group):
        if a != b:
            raise RuntimeError(f'{group}: {serialize(plus(a, times(-1, b)))}')
        counts[group] = counts.get(group, 0)+1
    def reject(a, b, group):
        if a == b:
            raise RuntimeError(f'Mutation survived: {group}')
        counts[group] = counts.get(group, 0)+1
    basis = [tuple(F(i == j) for j in range(4)) for i in range(4)]
    for i, a in enumerate(basis):
        for j, b in enumerate(basis):
            equal(plus(mul(slash(a), slash(b)), mul(slash(b), slash(a))),
                  const(2*dot(a, b)), 'clifford_metric')
    # The residue is linear in p,h,k1,k2 and bilinear in e1,e2 by UV power
    # counting. Exhausting its 4*4*4*4=256 input basis coefficients proves the
    # polynomial identity, rather than merely sampling generic momenta.
    basis_cases = 0
    for vector_slot, mu, a, b in product(range(4), repeat=4):
        momenta = [ZERO]*4
        momenta[vector_slot] = basis[mu]
        args = (*momenta, basis[a], basis[b])
        raw = diagrams(*args)
        compact = compact_diagrams(*args)
        for name in raw:
            equal(raw[name], compact[name], 'diagram_basis_coefficients')
        equal(combine(raw, *args[-2:]), prediction(*args), 'summed_basis_coefficients')
        basis_cases += 1
    rng = Random(15327)
    generic = []
    for _ in range(generic_count):
        args = tuple(tuple(F(rng.randint(-7, 7), rng.randint(1, 5)) for _ in range(4))
                     for _ in range(6))
        p, h, k1, k2, e1, e2 = args
        raw = diagrams(*args)
        compact = compact_diagrams(*args)
        for name in raw:
            equal(raw[name], compact[name], 'generic_diagram_residues')
        total = combine(raw, e1, e2)
        equal(total, prediction(*args), 'generic_sum')
        swapped = diagrams(p, h, k2, k1, e2, e1)
        equal(total, combine(swapped, e2, e1), 'bose_symmetry')
        # Each deliberate corruption must fail the complete off-shell match.
        reject(plus(total, times(2*dot(e1, e2), raw['SG'])), total, 'omit_seagull')
        reject(plus(total, times(4*dot(e1, e2), raw['SG'])), total, 'reverse_seagull')
        reject(plus(total, times(2, plus(raw['FS1'], raw['FS2']))), total, 'equal_line_charges')
        reject(plus(total, times(-1, raw['FF2'])), total, 'omit_fermion_ordering')
        reject(times(dot(e1, e2), slash(p)), total, 'omit_sigma_connection')
        generic.append({'momenta': [[str(x) for x in v] for v in args],
                        'diagrams': {name: serialize(poly) for name, poly in raw.items()},
                        'sum': serialize(total)})
    return {'status': 'PASS', 'method': 'exact Feynman-simplex moments plus independent denominator cancellation',
            'normalization': 'loop = -i F (gprime*yL)^2 R2/(16*pi^2*epsilon)',
            'basis_cases': basis_cases, 'generic_cases': generic_count, 'checks': counts,
            'total_checks': sum(counts.values()), 'samples': generic[:2]}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json', action='store_true')
    parser.add_argument('--generic-count', type=int, default=12)
    args = parser.parse_args()
    result = run_checks(args.generic_count)
    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print(f"PASS: {result['total_checks']} exact checks; {result['basis_cases']} complete polynomial-basis cases, "
              f"{result['generic_cases']} generic rational cases; seven diagrams integrated independently.")
