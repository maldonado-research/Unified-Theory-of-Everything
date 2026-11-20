"""Exact one-loop tensor UV residues in Euclidean dimensional regularization.

For denominators D_i=(q-p_i)^2+m_i^2 and d=4-2*epsilon, this module
returns R=16*pi**2*Res_{epsilon=0} of the integral with measure
d**d q/(2*pi)**d.  ``masses2`` contains squared masses, not masses.

Only UV poles are represented.  Exceptional massless kinematics can have IR
poles as well; this polynomial UV functional does not infer or subtract those
IR poles.  External vectors and open components are four-dimensional, which
is sufficient for the simple UV residue; finite evanescent terms are omitted.
UVIntegral includes no Minkowski continuation; MinkowskiUVIntegral provides
the explicitly documented (+---), overall-i-stripped continuation.  Neither
class supplies couplings, vertex factors, or diagram symmetry factors.  The
calculation is a reusable standard-method helper, not a claim of a new
physical result.

Example::

    import sympy as s
    q = s.symbols('q0:4')
    loop = UVIntegral([(0,0,0,0), (1,2,0,0)], [2,3])
    loop.component((0, 1))
    loop.polynomial_residue(q[0]*q[1], q)

Run ``python tensor_engine.py --self-test`` for independent and algebraic
checks, including a large-momentum angular-average implementation.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
from itertools import combinations_with_replacement, product
from math import factorial

import sympy as sp


def _exact(value):
    result = sp.sympify(value)
    if result.has(sp.Float):
        raise ValueError('Use exact integers, Rational values, or symbols; no floats.')
    return result


def _dot(a, b):
    return sum((x*y for x, y in zip(a, b)), sp.S.Zero)


def simplex_integral(polynomial, variables):
    """Integrate a polynomial on x_i>=0, sum_i x_i=1, without a prefactor.

    Each monomial prod(x_i**a_i) integrates to
    prod(a_i!)/(sum(a_i)+N-1)!.  In particular, volume=1/(N-1)!.
    Coefficients may be exact symbolic expressions independent of x_i.
    """
    variables = tuple(variables)
    if not variables:
        raise ValueError('A simplex needs at least one barycentric variable.')
    poly = sp.Poly(_exact(polynomial), *variables)
    result = sp.S.Zero
    for exponents, coefficient in poly.terms():
        numerator = sp.prod(factorial(a) for a in exponents)
        result += coefficient*sp.Rational(numerator,
                                         factorial(sum(exponents)+len(variables)-1))
    return sp.expand(result)


class UVIntegral:
    """Cached exact UV functional for an arbitrary positive denominator count.

    ``component((mu, nu, ...))`` gives the residue for q_mu*q_nu*...
    and ``polynomial_residue(f, q)`` handles any polynomial numerator.
    For speed, exact numeric external vectors are preferable; symbolic
    masses and momenta are supported, but high-rank symbolic expansions grow.
    """

    def __init__(self, shifts, masses2):
        self.shifts = tuple(tuple(_exact(x) for x in p) for p in shifts)
        self.masses2 = tuple(_exact(m) for m in masses2)
        self.N = len(self.shifts)
        if self.N < 1 or len(self.masses2) != self.N:
            raise ValueError('Provide at least one shift and one squared mass per shift.')
        if any(len(p) != 4 for p in self.shifts):
            raise ValueError('Each Euclidean external vector needs four components.')
        # Dummy objects prevent collisions with symbols supplied by callers.
        self.x = tuple(sp.Dummy(f'x{i}') for i in range(self.N))
        self.P = tuple(sum((x*p[mu] for x, p in zip(self.x, self.shifts)),
                           sp.S.Zero) for mu in range(4))
        self.Delta = sp.expand(sum((x*(m+_dot(p, p)) for x, p, m in
                                   zip(self.x, self.shifts, self.masses2)), sp.S.Zero)
                               - _dot(self.P, self.P))
        self._moments = {}
        self._components = {}

    def _moment(self, delta_power, powers):
        key = (delta_power, powers)
        if key not in self._moments:
            polynomial = self.Delta**delta_power
            for P, exponent in zip(self.P, powers):
                if exponent:
                    polynomial *= P**exponent
            self._moments[key] = simplex_integral(polynomial, self.x)
        return self._moments[key]

    def component(self, indices=()):
        """Return an exact symmetric tensor component; indices range from 0 to 3."""
        indices = tuple(sorted(indices))
        if any(not isinstance(mu, int) or not 0 <= mu < 4 for mu in indices):
            raise ValueError('Tensor indices must be integers from 0 through 3.')
        if indices in self._components:
            return self._components[indices]
        counts = tuple(indices.count(mu) for mu in range(4))
        answer = sp.S.Zero
        # Pairings between unlike component labels vanish.  For each label,
        # n! / ((n-2k)! 2**k k!) counts distinct unnormalized metric pairings.
        for pair_counts in product(*(range(n//2+1) for n in counts)):
            k = sum(pair_counts)
            j = k+2-self.N
            if j < 0:
                continue
            powers = tuple(n-2*a for n, a in zip(counts, pair_counts))
            multiplicity = 1
            for n, a in zip(counts, pair_counts):
                multiplicity *= factorial(n)//(factorial(n-2*a)*2**a*factorial(a))
            coefficient = sp.Rational((-1)**j*multiplicity, 2**k*factorial(j))
            answer += coefficient*self._moment(j, powers)
        answer = sp.expand(answer)
        self._components[indices] = answer
        return answer

    def polynomial_residue(self, numerator, q):
        """Apply this UV functional linearly to an exact polynomial in four q's."""
        q = tuple(q)
        if len(q) != 4 or len(set(q)) != 4 or not all(x.is_Symbol for x in q):
            raise ValueError('q must contain four distinct symbols.')
        if any(expr.has(*q) for row in self.shifts for expr in row) or any(
                expr.has(*q) for expr in self.masses2):
            raise ValueError('External shifts and squared masses must not depend on q.')
        poly = sp.Poly(_exact(numerator), *q)
        answer = sp.S.Zero
        for exponents, coefficient in poly.terms():
            indices = tuple(mu for mu, degree in enumerate(exponents) for _ in range(degree))
            answer += coefficient*self.component(indices)
        return sp.expand(answer)

    def tensor(self, rank):
        """Return only independent components, indexed by sorted index tuples."""
        if not isinstance(rank, int) or rank < 0:
            raise ValueError('rank must be a nonnegative integer.')
        return {indices: self.component(indices) for indices in
                combinations_with_replacement(range(4), rank)}

    def without(self, denominator):
        """Return the integral after cancelling one denominator (N must exceed 1)."""
        if self.N <= 1:
            raise ValueError('A denominator-free polynomial is outside this API.')
        if not 0 <= denominator < self.N:
            raise ValueError('Denominator index is out of range.')
        return UVIntegral(self.shifts[:denominator]+self.shifts[denominator+1:],
                          self.masses2[:denominator]+self.masses2[denominator+1:])


class MinkowskiUVIntegral:
    """Time-first (+---) wrapper, with the overall Wick-rotation i stripped.

    Denominators are (q-p_i)^2-m_i^2+i0.  The returned normalization is
    (16*pi**2/i)*Res_epsilon integral, with no vertex/propagator i factors.
    Wick rotation gives p_E=(-i*p_0,p_1,p_2,p_3), q_0=i*q_E0,
    and one factor (-1) per denominator.  This is an algebraic continuation
    of UV polynomials and does not evaluate threshold/branch-cut finite parts.
    """

    def __init__(self, shifts, masses2):
        self.shifts = tuple(tuple(_exact(a) for a in p) for p in shifts)
        self.masses2 = tuple(_exact(m) for m in masses2)
        if any(len(p) != 4 for p in self.shifts):
            raise ValueError('Each Minkowski external vector needs four components.')
        self.euclidean = UVIntegral([(-sp.I*p[0], *p[1:]) for p in self.shifts],
                                    self.masses2)
        self.N = self.euclidean.N

    def component(self, indices=()):
        indices = tuple(indices)
        return sp.expand((-1)**self.N*sp.I**indices.count(0)
                         *self.euclidean.component(indices))

    def polynomial_residue(self, numerator, q):
        q = tuple(q)
        if len(q) != 4:
            raise ValueError('q must contain four distinct symbols.')
        rotated = _exact(numerator).subs(q[0], sp.I*q[0], simultaneous=True)
        return sp.expand((-1)**self.N*self.euclidean.polynomial_residue(rotated, q))

    def tensor(self, rank):
        if not isinstance(rank, int) or rank < 0:
            raise ValueError('rank must be a nonnegative integer.')
        return {indices: self.component(indices) for indices in
                combinations_with_replacement(range(4), rank)}

    def without(self, denominator):
        if self.N <= 1:
            raise ValueError('A denominator-free polynomial is outside this API.')
        if not 0 <= denominator < self.N:
            raise ValueError('Denominator index is out of range.')
        return MinkowskiUVIntegral(self.shifts[:denominator]+self.shifts[denominator+1:],
                                   self.masses2[:denominator]+self.masses2[denominator+1:])


@lru_cache(maxsize=None)
def _sphere_monomial(exponents):
    """Normalized angular mean on S^3, used only by independent checks."""
    if any(a % 2 for a in exponents):
        return sp.S.Zero
    degree = sum(exponents)//2
    odd_factorials = 1
    for a in exponents:
        for r in range(1, a, 2):
            odd_factorials *= r
    return sp.Rational(odd_factorials, 2**degree*factorial(degree+1))


def asymptotic_residue(shifts, masses2, indices):
    """Independent large-|q| expansion and S^3 angular average.

    This reference uses no Feynman parameters and no simplex moments.
    Expand prod_i(1-2*t*u.p_i+t**2*(p_i**2+m_i**2))**(-1).
    The coefficient of t**(rank+4-2*N), multiplied by prod(u_mu),
    has precisely the radial logarithmic UV degree.  Its angular mean is R.
    """
    shifts = tuple(tuple(_exact(a) for a in p) for p in shifts)
    masses2 = tuple(_exact(m) for m in masses2)
    degree = len(indices)+4-2*len(shifts)
    if degree < 0:
        return sp.S.Zero
    u = sp.symbols('_u0:4')
    series = [sp.S.One]+[sp.S.Zero]*degree
    for p, mass2 in zip(shifts, masses2):
        linear = 2*_dot(u, p)
        constant = _dot(p, p)+mass2
        inverse = [sp.S.One]
        for j in range(1, degree+1):
            inverse.append(sp.expand(linear*inverse[j-1]
                                     -(constant*inverse[j-2] if j >= 2 else 0)))
        series = [sp.expand(sum((series[j]*inverse[n-j] for j in range(n+1)),
                                sp.S.Zero)) for n in range(degree+1)]
    angular_polynomial = sp.Poly(series[degree]*sp.prod(u[mu] for mu in indices), *u)
    return sp.expand(sum((coefficient*_sphere_monomial(exponents)
                          for exponents, coefficient in angular_polynomial.terms()),
                         sp.S.Zero))


def self_test():
    """Meaningful producer checks plus deliberately incorrect negative controls."""
    checks = 0

    def equal(actual, expected, label):
        nonlocal checks
        difference = sp.expand(actual-expected)
        if difference != 0:
            raise AssertionError(f'{label}: difference {difference}')
        checks += 1

    def unequal(actual, wrong, label):
        nonlocal checks
        if sp.expand(actual-wrong) == 0:
            raise AssertionError(f'Negative control unexpectedly passed: {label}')
        checks += 1

    q = sp.symbols('q0:4')
    a, b = sp.symbols('a b', positive=True)
    zero = (0, 0, 0, 0)
    p = (sp.Rational(1, 2), -1, sp.Rational(2, 3), 2)
    r = (-2, sp.Rational(1, 3), 1, -1)
    t = (1, -2, 3, sp.Rational(1, 5))

    x = sp.symbols('x0:3')
    equal(simplex_integral(1, x), sp.Rational(1, 2), 'simplex volume')
    equal(simplex_integral(x[0]**2*x[1], x), sp.Rational(1, 60), 'simplex monomial')
    tadpole = UVIntegral([zero], [a])
    equal(tadpole.component(), -a, 'massive scalar tadpole')
    equal(tadpole.component((0, 0)), a**2/4, 'massive rank-two tadpole')
    equal(tadpole.component((0, 0, 1, 1)), -a**3/24, 'massive rank-four tadpole')
    bubble = UVIntegral([zero, zero], [a, b])
    equal(bubble.component(), 1, 'massive bubble scalar')
    equal(bubble.component((0, 0)), -(a+b)/4, 'massive bubble rank two')
    equal(bubble.component((0, 0, 1, 1)), (a*a+a*b+b*b)/24,
          'massive bubble rank four')
    equal(bubble.component((0, 0, 1, 1, 2, 2)),
          -(a**3+a*a*b+a*b*b+b**3)/192, 'massive bubble rank six')

    triangle = UVIntegral([zero, p, r], [2, 3, 5])
    equal(triangle.component(), 0, 'scalar triangle UV convergence')
    equal(triangle.component((0, 0)), sp.Rational(1, 4), 'triangle rank-two pole')
    equal(triangle.component((0, 1)), 0, 'triangle rank-two off-diagonal')
    box = UVIntegral([zero, p, r, t], [2, 3, 5, 7])
    equal(box.component((0, 0, 1, 1)), sp.Rational(1, 24), 'box metric pair coefficient')
    equal(box.component((0, 0, 0, 0)), sp.Rational(1, 8), 'box unnormalized pairings')
    equal(box.component((0, 0, 1)), 0, 'box subleading rank UV convergence')
    unequal(box.component((0, 0, 0, 0)), sp.Rational(1, 24),
            'wrong normalized rank-four symmetrization')

    # A denominator insertion cancels exactly, including masses and shifted q.
    # These compare two different denominator counts, not two rewrites of one
    # tensor formula.  Bubble rank six and triangle rank four are exercised.
    for loop, f in [(triangle, 1), (triangle, q[0]),
                    (triangle, q[0]*q[2]),
                    (UVIntegral([p, r], [a, b]), q[0]**2*q[1]**2)]:
        for j, (shift, mass2) in enumerate(zip(loop.shifts, loop.masses2)):
            D = sum((q[mu]-shift[mu])**2 for mu in range(4))+mass2
            equal(loop.polynomial_residue(D*f, q),
                  loop.without(j).polynomial_residue(f, q),
                  f'denominator cancellation N={loop.N}, j={j}, f={f}')

    # Explicit trace and contracted momentum differences reproduce reduced
    # denominators.  This catches signs of P, masses, and metric factors.
    f = q[1]*q[2]
    pj = triangle.shifts[1]
    trace = sum(triangle.polynomial_residue(q[mu]**2*f, q) for mu in range(4))
    trace_expected = (triangle.without(1).polynomial_residue(f, q)
                      +2*triangle.polynomial_residue(_dot(q, pj)*f, q)
                      -(_dot(pj, pj)+triangle.masses2[1])
                      *triangle.polynomial_residue(f, q))
    equal(trace, trace_expected, 'rank-four massive triangle trace')
    for f in [1, q[0], q[0]*q[1], q[2]**3]:
        i, j = 0, 2
        pi, pj = triangle.shifts[i], triangle.shifts[j]
        c = (_dot(pi, pi)+triangle.masses2[i]
             -_dot(pj, pj)-triangle.masses2[j])
        expected = (c*triangle.polynomial_residue(f, q)
                    -triangle.without(i).polynomial_residue(f, q)
                    +triangle.without(j).polynomial_residue(f, q))/2
        equal(triangle.polynomial_residue(_dot(q, tuple(x-y for x, y in zip(pi, pj)))*f, q),
              expected, f'momentum contraction f={f}')

    # Covariance under shifts of loop routing, plus denominator permutation.
    translation = (sp.Rational(2, 5), -1, 0, sp.Rational(3, 7))
    moved = UVIntegral([tuple(x+y for x, y in zip(shift, translation))
                        for shift in triangle.shifts], triangle.masses2)
    permuted = UVIntegral(triangle.shifts[::-1], triangle.masses2[::-1])
    for indices in [(), (0,), (0, 0), (1, 2), (0, 1, 1), (0, 0, 1, 2), (2, 2, 2, 2)]:
        shifted_numerator = sp.prod(q[mu]+translation[mu] for mu in indices)
        equal(moved.component(indices), triangle.polynomial_residue(shifted_numerator, q),
              f'loop-routing covariance {indices}')
        equal(permuted.component(indices), triangle.component(indices),
              f'denominator permutation {indices}')

    # Independent UV asymptotics, evaluated for several topologies/ranks.
    reference_cases = [([p], [2], [(0,), (0, 1), (0, 0, 1, 1)]),
                       ([p, r], [2, 3], [(), (0,), (0, 1), (0, 0, 1),
                                          (0, 0, 1, 2), (0, 0, 1, 1, 2, 2)]),
                       ([zero, p, r], [2, 3, 5], [(0,), (0, 0), (0, 1, 1),
                                                    (0, 0, 1, 1), (0, 1, 2, 3)]),
                       ([zero, p, r, t], [2, 3, 5, 7], [(0, 0, 1), (0, 0, 1, 1)]),
                       ([zero, p, r, t, zero], [2, 3, 5, 7, 11],
                        [(0, 0, 1, 1), (0, 0, 1, 1, 2, 2)])]
    for shifts, masses2, components in reference_cases:
        loop = UVIntegral(shifts, masses2)
        for indices in components:
            equal(loop.component(indices), asymptotic_residue(shifts, masses2, indices),
                  f'independent angular UV expansion N={loop.N}, indices={indices}')
    independent = asymptotic_residue(triangle.shifts, triangle.masses2, (0, 0))
    unequal(-triangle.component((0, 0)), independent, 'reversed UV pole sign')
    unequal(triangle.component((0, 0))/2, independent, 'spurious tensor normalization')

    # Independent analytic signature checks anchor the Wick-rotation factors.
    mink = MinkowskiUVIntegral([zero, p, r], [2, 3, 5])
    equal(mink.component((0, 0)), sp.Rational(1, 4), 'Minkowski temporal triangle')
    equal(mink.component((1, 1)), -sp.Rational(1, 4), 'Minkowski spatial triangle')
    equal(MinkowskiUVIntegral([zero], [a]).component(), a,
          'Minkowski scalar tadpole sign')
    for j, (shift, mass2) in enumerate(zip(mink.shifts, mink.masses2)):
        D = (q[0]-shift[0])**2-sum((q[mu]-shift[mu])**2 for mu in range(1, 4))-mass2
        equal(mink.polynomial_residue(D*q[0]*q[1], q),
              mink.without(j).polynomial_residue(q[0]*q[1], q),
              f'Minkowski denominator cancellation j={j}')

    try:
        UVIntegral([zero], [1.5])
    except ValueError:
        checks += 1
    else:
        raise AssertionError('Inexact floating-point input was accepted.')
    print(f'PASS: {checks} exact checks, including 18 independent angular UV cases '
          'and 3 negative controls.')
    return checks


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test', action='store_true', help='Run exact validation checks.')
    arguments = parser.parse_args()
    if arguments.self_test:
        self_test()
    else:
        parser.print_help()
