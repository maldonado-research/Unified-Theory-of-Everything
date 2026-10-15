#!/usr/bin/env python3
"""Two weak fields and mixed BW: exact ordered-color UV residue checker.

PROVENANCE: the Lorentz polynomial, Clifford, moment, routed-integrand and
compact individual-residue helpers below are adapted from
checkpoint_2026_10_15/two_gauge_leg/verify_two_b.py at source commit
5c77bc5243b32da7b97a0ab86ea8571cb1cff5c0. They are reused methodology, not
a new independent derivation of the primitive moments or seven Lorentz
residues. The routed integrands retain only the crossed seed; the direct
seed is separately checked UV finite. Other new work begins at the
Gaussian/color helpers after prediction().
The new calculation evaluates their ordered weak-index combinations and
checks the commutator residue and full nonabelian action polynomial.

Original engine description:

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
    """Seven crossed integrands, before weak-index contractions."""
    K = vs(k1, k2)
    t = vs(p, scale(-1, K))
    out = {}
    for i, ki, kj, ei, ej in [(1, k1, k2, e1, e2), (2, k2, k1, e2, e1)]:
        seed_ff = square_shift(vs(h, K))
        chain_ff = mul(qslash_shift(K), slash(ej), qslash_shift(ki), slash(ei), QSLASH)
        out[f'FF{i}'] = residue(mul(seed_ff, chain_ff), [ZERO, ki, K, p])
        seed_ss = square_shift(h)
        u = vs(scale(2, p), scale(-1, ki))
        v = vs(scale(2, p), scale(-2, ki), scale(-1, kj))
        v1 = plus(const(dot(u, ei)), times(-2, qdot(ei)))
        v2 = plus(const(dot(v, ej)), times(-2, qdot(ej)))
        out[f'SS{i}'] = residue(mul(seed_ss, QSLASH, v1, v2),
                                [ZERO, p, vs(p, scale(-1, ki)), t])
        seed_fs = square_shift(vs(h, ki))
        uj = vs(scale(2, p), scale(-1, kj))
        vj = plus(const(dot(uj, ej)), times(-2, qdot(ej)))
        out[f'FS{i}'] = residue(mul(seed_fs, qslash_shift(ki), slash(ei), QSLASH, vj),
                                [ZERO, ki, p, vs(p, scale(-1, kj))])
    seed_sg = square_shift(h)
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


# New ordered weak-field calculation. All Gaussian arithmetic is rational.
def comm(a, b):
    return plus(mul(a, b), times(-1, mul(b, a)))


def commutator_residue(residues):
    return plus(residues['FF1'], times(-1, residues['FF2']),
                times(-1, residues['SS1']), residues['SS2'],
                residues['FS1'], times(-1, residues['FS2']))


def commutator_formula(p, h, k1, k2, e1, e2):
    r = vs(p, scale(-1, h), scale(-1, k1), scale(-1, k2))
    E1, E2 = slash(e1), slash(e2)
    connection = times(F(1, 4), plus(mul(comm(slash(k1), E1), E2),
                                   times(-1, mul(comm(slash(k2), E2), E1))))
    curvature = times(F(1, 2), mul(comm(E1, E2), slash(r)))
    divergence = times(F(2, 3), plus(
        times(dot(vs(scale(2, k1), k2), e2), E1),
        times(-dot(vs(k1, scale(2, k2)), e1), E2),
        times(dot(e1, e2), slash(vs(k2, scale(-1, k1))))))
    return plus(connection, curvature, divergence)


def z(x=0, y=0):
    return (F(x), F(y))


def za(*zs):
    return tuple(sum(v[i] for v in zs) for i in range(2))


def zm(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])


Z = z()
WID = ((z(1), Z), (Z, z(1)))
WZERO = ((Z, Z), (Z, Z))
GENERATORS = (
    ((Z, z(F(1, 2))), (z(F(1, 2)), Z)),
    ((Z, z(0, F(-1, 2))), (z(0, F(1, 2)), Z)),
    ((z(F(1, 2)), Z), (Z, z(F(-1, 2)))),
    ((z(F(-1, 2)), Z), (Z, z(F(-1, 2)))),
)


def wa(*ms):
    return tuple(tuple(za(*(m[i][j] for m in ms)) for j in range(2)) for i in range(2))


def wt(m):
    return tuple(tuple(m[j][i] for j in range(2)) for i in range(2))


def ws(c, m):
    return tuple(tuple(zm(z(c), x) for x in row) for row in m)


def wm(*ms):
    out = WID
    for m in ms:
        out = tuple(tuple(za(*(zm(out[i][k], m[k][j]) for k in range(2)))
                          for j in range(2)) for i in range(2))
    return out


def weak_contraction(C, S, crossed):
    """Literal seed indices, no pre-contracted matrix shortcut."""
    out = [[Z, Z], [Z, Z]]
    for a, c, b, d, e, f in product(range(2), repeat=6):
        seed = (a == d and b == c) if crossed else (a == c and b == d)
        if seed and e == f:
            out[a][c] = za(out[a][c], zm(C[b][e], S[d][f]))
    return tuple(map(tuple, out))


def diagram_colors(X1, X2, crossed=True, missing_transpose=False):
    trans = (lambda x: x) if missing_transpose else wt
    C1, C2 = ws(-1, trans(X1)), ws(-1, trans(X2))
    return {
        'FF1': weak_contraction(wm(C2, C1), WID, crossed),
        'FF2': weak_contraction(wm(C1, C2), WID, crossed),
        'SS1': weak_contraction(WID, wm(X2, X1), crossed),
        'SS2': weak_contraction(WID, wm(X1, X2), crossed),
        'FS1': weak_contraction(C1, X2, crossed),
        'FS2': weak_contraction(C2, X1, crossed),
        # The -1 is the extra propagator times seagull relative phase.
        'SG': ws(-1, weak_contraction(WID, wa(wm(X1, X2), wm(X2, X1)), crossed)),
    }


def tensor(color, poly):
    return {(a, c, mask): zm(color[a][c], z(value))
            for a, c in product(range(2), repeat=2)
            for (_, mask), value in poly.items() if color[a][c] != Z and value}


def ta(*ts):
    out = {}
    for t in ts:
        for key, val in t.items():
            out[key] = za(out.get(key, Z), val)
    return {key: val for key, val in out.items() if val != Z}


def weighted(residues, colors, e1, e2):
    return ta(*(tensor(colors[name], times(dot(e1, e2) if name == 'SG' else 1, value))
                for name, value in residues.items()))


def full_prediction(args, X1, X2):
    return ta(tensor(ws(F(1, 2), wa(wm(X1, X2), wm(X2, X1))), prediction(*args)),
              tensor(ws(F(1, 2), wa(wm(X1, X2), ws(-1, wm(X2, X1)))),
                     commutator_formula(*args)))


def one_residue(p, h, k, e, m2=F(0)):
    r = vs(p, scale(-1, h), scale(-1, k))
    return plus(
        times(-dot(vs(scale(2, r), k), e)/3 + dot(vs(scale(2, h), k), e)/6, slash(p)),
        times(dot(p, p)/6-m2/2, slash(e)),
        times(F(1, 4), mul(comm(slash(k), slash(e)), slash(r))),
        times(F(1, 3), plus(times(dot(k, e), slash(k)), times(-dot(k, k), slash(e)))))


def ward_rhs(args, X1, X2, m2):
    p, h, k1, k2, e1, e2 = args
    return ta(
        tensor(wm(X2, X1), one_residue(p, vs(h, k1), k2, e2, m2)),
        tensor(ws(-1, wm(X1, X2)), one_residue(p, h, k2, e2, m2)),
        tensor(wa(wm(X1, X2), ws(-1, wm(X2, X1))),
               one_residue(p, h, vs(k1, k2), e2, m2)))


def serialize(poly):
    return {str(mask): str(value) for (_, mask), value in sorted(poly.items())}


def run_checks(generic_count=12):
    counts = {}
    def equal(a, b, group):
        if a != b:
            raise RuntimeError(f'{group}: actual={a}; expected={b}')
        counts[group] = counts.get(group, 0)+1
    def reject(a, b, group):
        if a == b:
            raise RuntimeError(f'Mutation survived: {group}')
        counts[group] = counts.get(group, 0)+1
    # All three weak generators and hypercharge, including imaginary tau2.
    for X1, X2 in product(GENERATORS, repeat=2):
        c = diagram_colors(X1, X2)
        expected = {'FF1': wm(X1, X2), 'FF2': wm(X2, X1),
                    'SS1': wm(X2, X1), 'SS2': wm(X1, X2),
                    'FS1': ws(-1, wm(X2, X1)), 'FS2': ws(-1, wm(X1, X2)),
                    'SG': ws(-1, wa(wm(X1, X2), wm(X2, X1)))}
        direct = diagram_colors(X1, X2, crossed=False)
        for name in c:
            equal(c[name], expected[name], 'literal_crossed_weak_indices')
            trace = za(expected[name][0][0], expected[name][1][1])
            equal(direct[name], ((trace, Z), (Z, trace)), 'literal_direct_weak_indices')
    # UV-finite direct numerator: box degree <=3, SG triangle degree 1.
    for n, degree in ((4, 3), (3, 1)):
        for powers in product(range(degree+1), repeat=4):
            if sum(powers) <= degree:
                equal(simplex_moment(powers, [ZERO]*n), 0, 'direct_seed_uv_power_count')
    basis = [tuple(F(i == j) for j in range(4)) for i in range(4)]
    basis_cases = 0
    for vector_slot, mu, a, b in product(range(4), repeat=4):
        momenta = [ZERO]*4
        momenta[vector_slot] = basis[mu]
        args = (*momenta, basis[a], basis[b])
        raw = diagrams(*args)
        # Reuse old Lorentz integration as an input method, explicitly counted
        # separately from the new noncommutative/color conclusions.
        compact = compact_diagrams(*args)
        for name in raw:
            equal(raw[name], compact[name], 'inherited_lorentz_diagram_reproduction')
        equal(commutator_residue(raw), commutator_formula(*args), 'new_commutator_basis_coefficients')
        equal(combine(raw, *args[-2:]), prediction(*args), 'inherited_abelian_reproduction')
        for i, j in product(range(4), repeat=2):
            X1, X2 = GENERATORS[i], GENERATORS[j]
            group = 'new_mixed_BW_basis_coefficients' if (i == 3) != (j == 3) else (
                'inherited_BB_embedding' if i == j == 3 else 'new_WW_basis_coefficients')
            equal(weighted(raw, diagram_colors(X1, X2), *args[-2:]),
                  full_prediction(args, X1, X2), group)
        basis_cases += 1
    rng = Random(151020269)
    samples = []
    for _ in range(generic_count):
        args = tuple(tuple(F(rng.randint(-7, 7), rng.randint(1, 5)) for _ in range(4))
                     for _ in range(6))
        p, h, k1, k2, e1, e2 = args
        raw = diagrams(*args)
        equal(commutator_residue(raw), commutator_formula(*args), 'new_generic_commutator')
        swapped = (p, h, k2, k1, e2, e1)
        equal(commutator_formula(*args), times(-1, commutator_formula(*swapped)), 'commutator_bose_antisymmetry')
        for i, j in product(range(4), repeat=2):
            X1, X2 = GENERATORS[i], GENERATORS[j]
            equal(full_prediction(args, X1, X2), full_prediction(swapped, X2, X1), 'full_bose_symmetry')
            contracted = (p, h, k1, k2, k1, e2)
            equal(full_prediction(contracted, X1, X2), ward_rhs(args, X1, X2, F(-7, 3)),
                  'nonabelian_ward_with_lower_gauge_term')
        X1, X2 = GENERATORS[0], GENERATORS[1]
        total = weighted(raw, diagram_colors(X1, X2), e1, e2)
        reject({}, total, 'mutation_omit_commutator')
        reject(weighted(raw, diagram_colors(X1, X2, missing_transpose=True), e1, e2),
               total, 'mutation_omit_conjugate_transpose')
        bad = dict(raw)
        bad['FF1'], bad['FF2'] = bad['FF2'], bad['FF1']
        reject(weighted(bad, diagram_colors(X1, X2), e1, e2), total, 'mutation_reverse_FF_color_order')
        E1, E2 = slash(e1), slash(e2)
        r = vs(p, scale(-1, h), scale(-1, k1), scale(-1, k2))
        curvature = times(F(1, 2), mul(comm(E1, E2), slash(r)))
        reject(plus(commutator_formula(*args), times(-2, curvature)),
               commutator_formula(*args), 'mutation_wrong_nonlinear_curvature_sign')
        reject(plus(commutator_formula(*args), times(-1, curvature)),
               commutator_formula(*args), 'mutation_omit_nonlinear_sigma_curvature')
        divergence = times(F(2, 3), plus(
            times(dot(vs(scale(2, k1), k2), e2), E1),
            times(-dot(vs(k1, scale(2, k2)), e1), E2),
            times(dot(e1, e2), slash(vs(k2, scale(-1, k1))))))
        reject(plus(commutator_formula(*args), times(-1, divergence)),
               commutator_formula(*args), 'mutation_omit_curvature_divergence')
        # A distinct weak pair has zero anticommutator: use equal generators
        # to test the physical scalar seagull and mixed BW to test its factor.
        X1, X2 = GENERATORS[0], GENERATORS[0]
        total = weighted(raw, diagram_colors(X1, X2), e1, e2)
        colors = diagram_colors(X1, X2)
        colors['SG'] = WZERO
        reject(weighted(raw, colors, e1, e2), total, 'mutation_omit_scalar_seagull')
        X1, X2 = GENERATORS[3], GENERATORS[1]
        total = weighted(raw, diagram_colors(X1, X2), e1, e2)
        colors = diagram_colors(X1, X2)
        colors['SG'] = ws(F(1, 2), colors['SG'])
        reject(weighted(raw, colors, e1, e2), total, 'mutation_half_mixed_seagull')
        samples.append({'inputs_p_h_k1_k2_e1_e2': [[str(x) for x in v] for v in args],
                        'commutator_residue': serialize(commutator_residue(raw)),
                        'abelian_residue': serialize(combine(raw, e1, e2)),
                        'individual_residues': {name: serialize(value) for name, value in raw.items()}})
    return {'status': 'PASS', 'source_commit': '5c77bc5243b32da7b97a0ab86ea8571cb1cff5c0',
            'method': 'attributed Feynman-simplex Lorentz integration; new explicit ordered weak indices and commutator sum',
            'not_independent_of': 'two_gauge_leg/verify_two_b.py primitive moments and Lorentz helper implementation',
            'normalization': 'loop=-i F [g1 g2 kernel]/(16*pi^2*epsilon); generators include yL=-1/2 in BW',
            'basis_cases': basis_cases, 'generic_cases': generic_count, 'checks': counts,
            'total_checks': sum(counts.values()), 'samples': samples[:2]}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json', action='store_true')
    parser.add_argument('--generic-count', type=int, default=12)
    args = parser.parse_args()
    result = run_checks(args.generic_count)
    print(json.dumps(result, indent=2, sort_keys=True) if args.json else
          f"PASS: {result['total_checks']} exact checks; {result['basis_cases']} exhaustive Lorentz coefficients.")


