"""Verify the fully contracted same-weak-component on-shell primitive triangle."""
from fractions import Fraction
from itertools import permutations
import json
from pathlib import Path

import mpmath as mp


def require(condition, label):
    if not condition:
        raise ArithmeticError(label)


def close(left, right, label, tol='1e-32'):
    require(abs(left-right) <= mp.mpf(tol) * (1+abs(right)), label)


def spectral(s, a, b):
    if s <= 4*a:
        return mp.mpf(0)
    beta = mp.sqrt(1-4*a/s)
    rapidity_ratio = s*(1+beta)**2/(4*a)
    return (mp.mpf('0.5')-a/s) * mp.log((b+a*rapidity_ratio)/(b+a/rapidity_ratio))


def cut_from_angles(s, a, b, cos_theta):
    """Average the two seed assignments after computing each with actual momenta."""
    beta = mp.sqrt(1-4*a/s)
    h, v = s/2-a, s*beta/2
    sin_theta = mp.sqrt(1-cos_theta**2)

    def at_polar(u, direction):
        transverse = sin_theta*mp.sqrt(1-u*u)
        # Four equidistant azimuths integrate a constant plus one first harmonic exactly.
        values = []
        for cosine in [1, 0, -1, 0]:
            numerator = -h + direction*v*(cos_theta*u+transverse*cosine)
            denominator = -b-h+v*u
            values.append(numerator/denominator)
        return sum(values)/4

    return [beta*mp.quad(lambda u: at_polar(u, sign), [-1, 1])/2 for sign in (1, -1)]


def parameter_subtraction(s, a, b):
    """Upper lip of the integral, evaluating its y logarithm before integrating x."""
    endpoints = [mp.mpf(0), mp.mpf(1)]
    if s >= 4*a:
        beta = mp.sqrt(1-4*a/s)
        roots = [(1-beta)/2, (1+beta)/2]
        endpoints = sorted(set([*endpoints, *roots]))
    else:
        roots = []

    def across_y(x):
        first = a*x+b*(1-x)
        slope = b-a+s*x
        last = s*(x-roots[0])*(x-roots[1]) if roots else a-s*x*(1-x)
        if last == 0:
            return 0  # Endpoint value has zero integration measure.
        if slope == 0:
            return (1-x)/first
        if last > 0:
            # Avoid cancellation near the removable zero of the slope.
            if abs((1-x)*slope) < first/2:
                return -mp.log1p(-(1-x)*slope/first)/slope
            return mp.log(first/last)/slope
        return (mp.log(first/abs(last))+mp.j*mp.pi)/slope

    triangle = mp.quad(across_y, endpoints)
    vacuum = mp.quad(lambda u: u/(a*u+b*(1-u)), [0, 1])
    return (s/2-a)*triangle+a*vacuum


def dispersive_subtraction(s, a, b):
    """Once-subtracted cut integral in t=4a/s', including its upper-lip pole."""
    edge = 4*a
    def density(t):
        return spectral(edge/t, a, b)

    if s < edge:
        return mp.quad(lambda t: s*density(t)/(edge-s*t), [0, 1])
    if s == edge:
        # t=1-u². The endpoint singularity from 1/(1-t) is integrable.
        def at_threshold(u):
            if not u:
                return 2*a/(a+b)
            if u == 1:
                return 0  # Integrable logarithm at the infinite-energy endpoint.
            denom = (1-u)*(b*(1+u)+a*(1-u))
            return (1+u*u)*mp.log1p(4*a*u/denom)/(2*u)
        return mp.quad(at_threshold, [0, 1])
    pivot = edge/s
    at_pole = density(pivot)
    def regular(t):
        if abs(t-pivot) < mp.mpf('1e-24'):
            return -mp.diff(density, pivot)
        return (density(t)-at_pole)/(pivot-t)
    principal = mp.quad(regular, [0, pivot, 1])+at_pole*mp.log(pivot/(1-pivot))
    return principal+mp.j*mp.pi*at_pole


def dot(p, q):
    return p[0]*q[0]-sum(x*y for x,y in zip(p[1:], q[1:]))


def difference(p, q):
    return tuple(x-y for x,y in zip(p,q))


def eps(p,q):
    return p[0]*q[1]-p[1]*q[0]


def ordered_matrix_element(out, incoming):
    """Apply c†_o1 c†_o2 c_i1 c_i2 to |1,2> and project on |3,4>."""
    state = [1,2]
    sign = 1
    for mode in reversed(incoming):
        index = state.index(mode)
        sign *= (-1)**index
        state.pop(index)
    for mode in reversed(out):
        index = sum(i < mode for i in state)
        sign *= (-1)**index
        state.insert(index, mode)
    require(state == [3,4], 'outgoing state order')
    return sign


def check_literal():
    # Seed spin sign and half, two signed Yukawas, two scalar propagators,
    # sterile propagator, open-fermion crossing, identical-vertex expansion.
    wick = (-1j/2)*(-1j)**2*1j**2*(-1j)*(-1)/2
    require(wick == 0.25, 'literal single Wick coefficient')
    out_slots = list(permutations((3, 4)))
    in_slots = list(permutations((1, 2)))
    pairings = list(permutations((0, 1)))
    u = {1:(0,1), 2:(1,0), 3:(-1,3), 4:(3,1)}
    ordered_spin = eps(u[3],u[4])*eps(u[2],u[1])
    count = {3: 0, 4: 0}
    for out in out_slots:
        for incoming in in_slots:
            for pairing in pairings:
                # k belongs to the scalar ending on p1.
                count[out[pairing[incoming.index(1)]]] += 1
                spin = eps(u[out[0]],u[out[1]])*eps(u[incoming[0]],u[incoming[1]])
                require(ordered_matrix_element(out,incoming)*spin == ordered_spin,
                        'external sign in ordered S convention')
    require(count == {3: 4, 4: 4} and sum(count.values())*wick == 2,
            'eight labeled contractions with two different integrals')

    P = (10,0,0,0)
    p3, p4 = (5,3,0,4), (5,-3,0,-4)
    k, a = (2,1,2,3), Fraction(3, 2)
    s = dot(P, P)
    ns = [dot(difference(k, p), difference(k, p)) for p in (p3, p4)]
    d1, d2 = dot(k,k)-a, dot(difference(P,k),difference(P,k))-a
    require(sum(ns) == d1+d2+2*a-s and ns[0] != ns[1], 'off-shell loop numerator')
    require(dot(p3,p3) == 0 and dot(p4,p4) == 0, 'massless external states')
    # For p3 at theta with cos=4/5, choose normalized negative-helicity u3,u4.
    h, c = mp.sqrt(mp.mpf('0.1')), mp.sqrt(mp.mpf('0.9'))
    close(10**2*eps((-h,c),(c,h))*eps((1,0),(0,1)), -s, 'explicit spin S')
    return count


def display(value):
    return mp.nstr(value, 22)


def main():
    mp.mp.dps = 55
    counts = check_literal()
    checks = []
    for av, bv, sv in [(1,1,-2), (1,16,3), (1,1,4), (1,16,4),
                       (1,1,8), (1,16,8), (2,1,12), (2,5,9)]:
        a,b,s = map(mp.mpf, (av,bv,sv))
        direct = parameter_subtraction(s,a,b)
        dispersion = dispersive_subtraction(s,a,b)
        close(direct, dispersion, 'real and imaginary dispersive replay')
        row = dict(a=av, b=bv, s=sv, subtracted_real=display(mp.re(direct)),
                   subtracted_imag=display(mp.im(direct)))
        if s == 4*a:
            close(mp.im(direct), 0, 'threshold real boundary value')
            if b == a:
                close(direct, mp.mpf('0.5')+mp.pi**2/8, 'equal mass exact threshold')
        if s > 4*a:
            rho = spectral(s,a,b)
            close(mp.im(direct)/mp.pi, rho, 'Feynman delta Jacobian and cut sign')
            routed = cut_from_angles(s,a,b,mp.mpf('0.6'))
            require(routed[0] != routed[1], 'individual Wick routes have different angular kernels')
            close(sum(routed)/2, rho, 'two-particle angular cut')
            row['routed_cut_densities_at_cos_theta_0.6'] = list(map(display,routed))
            row['cut_density'] = display(rho)
        checks.append(row)
    unity = mp.mpf(1)
    beta = mp.mpf('1e-15')
    s = 4/(1-beta*beta)
    for b in (unity, 16*unity):
        close(spectral(s,unity,b)/beta, 1/(1+b), 'threshold onset', '1e-23')
    pole_weights = [1*(1+0j)**2, 4*(.5j)**2]
    require(pole_weights == [1,-1] and sum(pole_weights)==0, 'complex coherent mass weights')
    density_difference = spectral(8*unity,unity,unity)-spectral(8*unity,unity,16*unity)
    reduced = parameter_subtraction(8*unity,unity,unity)-parameter_subtraction(8*unity,unity,16*unity)
    close(mp.im(reduced), mp.pi*density_difference, 'physical unequal-mass projection')
    close((spectral(s,unity,unity)-spectral(s,unity,16*unity))/beta,
          mp.mpf(15)/34,'two-mass threshold cancellation test','1e-23')
    require(density_difference>0 and sum(abs(x) for x in pole_weights) != 0, 'cut survives pole zero')
    results = {
        'scope': 'same charged weak component, massless LL to LL, one primitive triangle in the unbroken gaugeless theory',
        'single_scalar_assignment_factor_before_integration': '1/4',
        'external_contractions_by_seed_leg_matching_p1':counts,
        'direct_flavor_selected_coefficient_of_ordered_S_K': '2 C_ee m_i (Y*_{mu i})²/(16 pi²)',
        'cut_definition': 'rho=(K(s+i0)-K(s-i0))/(2 pi i)',
        'subtraction': 'K(s+i0)-K(0)',
        'test_points': checks,
        'two_mass_example': {
            'a':1, 'sterile_masses':[1,4], 'pole_weights':[1,-1],
            'limit_rho_difference_over_beta': '15/34',
            's':8, 'rho_difference':display(density_difference),
            'subtracted_real_difference':display(mp.re(reduced)),
            'subtracted_imag_difference':display(mp.im(reduced)),
        },
    }
    Path(__file__).with_name('ON_SHELL_TRIANGLE_RESULTS.json').write_text(
        json.dumps(results, indent=2)+'\n', encoding='utf-8')
    print('PASS: literal Wick coefficients, whole-composite routing, on-shell spin, angular cuts and subtracted dispersion')


if __name__ == '__main__':
    main()
