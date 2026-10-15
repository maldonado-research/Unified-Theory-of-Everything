#!/usr/bin/env python3
"""Independent simplex-moment check of the two routed hypercharge residues.

Metric (+---); R = (16*pi**2/i) Res_{epsilon=0}.  The original numerator
polynomials are shifted and contracted as rank <= 4 Lorentz tensors.  No
denominator cancellation or large-q expansion is used.  Requires SymPy.
Run: python gauge_critic.py
"""

from functools import lru_cache
from math import factorial
import sys

import sympy as sp
from sympy.physics.matrices import mgamma


if sys.flags.optimize:
    raise RuntimeError("Scientific checks require assertions; run without Python -O or -OO.")


R = sp.Rational
METRIC = (1, -1, -1, -1)
GAMMA = tuple(mgamma(mu) for mu in range(4))
ZERO = sp.zeros(4)
x, y = sp.symbols("x y")
ell = sp.symbols("ell0:4")
p = sp.Matrix(sp.symbols("p0:4"))
h = sp.Matrix(sp.symbols("h0:4"))
k = sp.Matrix(sp.symbols("k0:4"))
m2 = sp.symbols("m2")


def dot(a, b):
    return sum(METRIC[mu] * a[mu] * b[mu] for mu in range(4))


def slash(a):
    return sum((METRIC[mu] * a[mu] * GAMMA[mu] for mu in range(4)), ZERO)


def clean(matrix):
    return matrix.applyfunc(sp.expand)


def assert_equal(actual, expected, label):
    difference = clean(actual - expected)
    assert difference == ZERO, (label, difference)


def nonzero(matrix):
    return clean(matrix) != ZERO


def contract(vector, tensor):
    return sum((METRIC[mu] * vector[mu] * tensor[mu] for mu in range(4)), ZERO)


@lru_cache(maxsize=None)
def simplex_moment(a, b):
    # The Feynman denominator factor 2! is included: <1> = 1.
    return R(2 * factorial(a) * factorial(b), factorial(a + b + 2))


def simplex_average(expression):
    polynomial = sp.Poly(sp.expand(expression), x, y)
    return sp.expand(sum(coef * simplex_moment(*powers)
                         for powers, coef in polynomial.terms()))


def interval_average(expression):
    return sp.expand(sum(coef / (power[0] + 1)
                         for power, coef in sp.Poly(sp.expand(expression), x).terms()))


def loop_moment(powers, delta, denominator_power):
    """UV pole moments after the shift, before parameter averaging.

    For denominator power 3, rank 0 is UV finite, rank 2 is g/4,
    rank 4 is delta*(gg+gg+gg)/8.  For power 2, rank 0 is 1 and
    rank 2 is delta*g/2.  Odd monomials vanish.  The diagonal metric
    makes all remaining mixed odd-index monomials vanish as well.
    """
    rank = sum(powers)
    if any(power % 2 for power in powers):
        return sp.S.Zero
    if rank == 0:
        return sp.S.One if denominator_power == 2 else sp.S.Zero
    if rank == 2:
        mu = powers.index(2)
        return METRIC[mu] * (delta / 2 if denominator_power == 2 else R(1, 4))
    if rank == 4 and denominator_power == 3:
        if 4 in powers:
            return 3 * delta / 8
        indices = [mu for mu, power in enumerate(powers) if power == 2]
        return METRIC[indices[0]] * METRIC[indices[1]] * delta / 8
    raise ValueError((powers, denominator_power))


def shifted_uv(numerator, delta, denominator_power=3, scalar_factor=1):
    # Keep the scalar and Clifford factors separate until the Lorentz-moment
    # selection: expanding the discarded rank-zero quartic coefficients is
    # expensive and contributes no UV pole for denominator power three.
    scalar_terms = sp.Poly(scalar_factor, *ell).terms()

    def reduce_entry(expression):
        contracted = sp.S.Zero
        for matrix_powers, matrix_coefficient in sp.Poly(expression, *ell).terms():
            for scalar_powers, scalar_coefficient in scalar_terms:
                powers = tuple(a + b for a, b in zip(matrix_powers, scalar_powers))
                moment = loop_moment(powers, delta, denominator_power)
                if moment != 0:
                    contracted += matrix_coefficient * scalar_coefficient * moment
        average = simplex_average if denominator_power == 3 else interval_average
        return average(contracted)
    return numerator.applyfunc(reduce_entry)


def triangles():
    """Start with each exact numerator and its own Feynman shift."""
    result_f, result_s = [], []
    t = p - k
    a_f = x * k + y * p
    delta_f = y * m2 - x * dot(k, k) - y * dot(p, p) + dot(a_f, a_f)
    q_f = sp.Matrix(ell) + a_f
    v_f = 2 * dot(t, t) + dot(q_f - h - k, q_f - h - k)

    a_s = x * p + y * t
    delta_s = (x + y) * m2 - x * dot(p, p) - y * dot(t, t) + dot(a_s, a_s)
    q_s = sp.Matrix(ell) + a_s
    v_s = 2 * dot(t, t) + dot(q_s - h, q_s - h)
    for nu in range(4):
        result_f.append(shifted_uv(slash(q_f - k) * GAMMA[nu] * slash(q_f), delta_f, scalar_factor=v_f))
        result_s.append(shifted_uv(slash(q_s) * (2 * p - 2 * q_s - k)[nu], delta_s, scalar_factor=v_s))
    return result_f, result_s


def bubble(momentum, spectator):
    """Compute T directly from the shifted two-propagator numerator."""
    q = sp.Matrix(ell) + x * momentum
    delta = x * m2 - x * (1 - x) * dot(momentum, momentum)
    scalar = 2 * dot(momentum, momentum) + dot(q - spectator, q - spectator)
    return shifted_uv(slash(q), delta, denominator_power=2, scalar_factor=scalar)


def gamma3(a, b, nu):
    """a_mu b_rho gamma^[mu nu rho], with unit-weight antisymmetrization."""
    return (slash(a) * GAMMA[nu] * slash(b) - a[nu] * slash(b)
            - b[nu] * slash(a) + dot(a, b) * GAMMA[nu])


def closed_form():
    pp, kk, hh = dot(p, p), dot(k, k), dot(h, h)
    pk, ph, kh = dot(p, k), dot(p, h), dot(k, h)
    f0 = -m2 / 2 - kk - R(7, 6) * pp + R(5, 2) * pk - R(5, 6) * kh + R(2, 3) * ph - hh / 2
    s0 = -m2 - pp + 2 * pk - R(5, 6) * kk + R(2, 3) * ph - kh / 3 - hh / 2
    d = 2 * p - k
    result_f, result_s = [], []
    for nu in range(4):
        result_f.append(f0 * GAMMA[nu] + R(2, 3) * p[nu] * slash(p)
                        - (k[nu] * slash(p) + p[nu] * slash(k)) / 2
                        + (k[nu] * slash(h) + h[nu] * slash(k)) / 6
                        - (p[nu] * slash(h) + h[nu] * slash(p)) / 3
                        + (gamma3(k, h, nu) - gamma3(k, p, nu)) / 2)
        result_s.append(s0 * GAMMA[nu] - k[nu] * slash(k) / 6
                        - d[nu] * slash(h) / 6 + h[nu] * slash(d) / 3)
    return result_f, result_s


def seed_gauge():
    coefficient = -R(2, 3) * p + h + k / 2
    return [coefficient[nu] * slash(p) + (dot(p, p) / 6 - m2 / 2) * GAMMA[nu]
            for nu in range(4)]


def remainder_form():
    r = p - h - k
    return [(r[nu] * slash(k) - dot(r, k) * GAMMA[nu] + gamma3(k, r, nu)) / 2
            + (k[nu] * slash(k) - dot(k, k) * GAMMA[nu]) / 3
            for nu in range(4)]


def on_shell_checks(remainder):
    cases = [
        ((5, 3, 4, 0), (5, 0, 3, 4), (0, 4, 0, 0), (2, 1, -2, 3)),
        ((1, 0, 0, 1), (1, 0, 0, -1), (0, 1, 0, 0), (3, -2, 4, 1)),
    ]
    rows_tested = 0
    for r_value, k_value, epsilon_value, h_value in cases:
        rr, kk, epsilon, hh = map(sp.Matrix, (r_value, k_value, epsilon_value, h_value))
        pp = rr + hh + kk
        assert dot(rr, rr) == dot(kk, kk) == dot(kk, epsilon) == 0
        substitutions = dict(zip(list(p) + list(h) + list(k), list(pp) + list(hh) + list(kk)))
        actual = contract(epsilon, [entry.subs(substitutions) for entry in remainder])
        target = dot(rr, epsilon) * slash(kk) - dot(rr, kk) * slash(epsilon)
        for left_kernel_vector in slash(rr).T.nullspace():
            row = left_kernel_vector.T
            assert clean(row * (actual - target)) == sp.zeros(1, 4)
            rows_tested += 1
    return rows_tested


def run():
    f, s = triangles()
    f_expected, s_expected = closed_form()
    for nu in range(4):
        assert_equal(f[nu], f_expected[nu], f"F[{nu}] simplex vs closed form")
        assert_equal(s[nu], s_expected[nu], f"S[{nu}] simplex vs closed form")
    print("PASS: all 8 triangle components from exact rank<=4 shifted numerators")

    t_ph = bubble(p, h)
    bubble_expected = ((m2 + dot(p, p) - R(2, 3) * dot(p, h) + dot(h, h) / 2) * slash(p)
                       + (dot(p, p) / 6 - m2 / 2) * slash(h))
    assert_equal(t_ph, bubble_expected, "independent three-point bubble")
    t = p - k
    shifted_t = bubble(t, h)
    t_phk = bubble(p, h + k)
    i_p = slash(p) / 2
    contact = 2 * (dot(p, p) - dot(t, t)) * i_p
    assert_equal(contract(k, f), shifted_t - t_phk + contact, "individual fermion Ward identity")
    assert_equal(contract(k, s), shifted_t - t_ph + contact, "individual scalar Ward identity")
    internal = [s[nu] - f[nu] for nu in range(4)]
    ward = t_phk - t_ph
    assert_equal(contract(k, internal), ward, "relative charge QL=-QH")
    print("PASS: independent bubble and both individual Ward identities; QL=-QH sum")

    gauge = seed_gauge()
    remainder = [clean(internal[nu] - gauge[nu]) for nu in range(4)]
    remainder_expected = remainder_form()
    for nu in range(4):
        assert_equal(remainder[nu], remainder_expected[nu], f"transverse remainder[{nu}]")
        assert not remainder[nu].has(m2)
    assert_equal(contract(k, gauge), ward, "local representative Ward identity")
    assert_equal(contract(k, remainder), ZERO, "off-shell transversality")
    print("PASS: exact transverse remainder and cancellation of all m2 terms")

    count = on_shell_checks(remainder)
    print(f"PASS: real transverse gauge leg and massless left lepton ({count} exact kernel rows)")

    # Mutation checks intentionally expose which checks do or do not constrain
    # a transverse term; merely checking the Ward identity is insufficient.
    assert nonzero(contract(k, [s[nu] + f[nu] for nu in range(4)]) - ward)
    assert nonzero(contract(k, [-f[nu] for nu in range(4)]) - ward)
    assert nonzero(contract(k, [-entry for entry in internal]) - ward)
    no_mass = [internal[nu] + m2 * GAMMA[nu] / 2 for nu in range(4)]
    assert nonzero(contract(k, no_mass) - ward)
    print("PASS mutations: wrong relative charge, missing scalar, common sign, omitted mass rejected by Ward anchor")

    r = p - h - k
    no_gamma3 = [remainder[nu] - gamma3(k, r, nu) / 2 for nu in range(4)]
    no_divergence = [remainder[nu] - (k[nu] * slash(k) - dot(k, k) * GAMMA[nu]) / 3
                     for nu in range(4)]
    for mutation in (no_gamma3, no_divergence):
        assert_equal(contract(k, mutation), ZERO, "transverse mutation still passes Ward")
        assert any(nonzero(mutation[nu] - remainder[nu]) for nu in range(4))
    try:
        on_shell_checks(no_gamma3)
    except AssertionError:
        pass
    else:
        raise AssertionError("omitting gamma3 must fail the on-shell dipole coefficient")
    on_shell_checks(no_divergence)
    print("PASS mutations: gamma3 omission evades Ward but fails residue/on-shell coefficient")
    print("PASS mutations: divergence omission evades Ward and real-B checks but fails off-shell residue")
    print("All checks passed. No physical beta function or action-level counterterm sign inferred.")


if __name__ == "__main__":
    run()
