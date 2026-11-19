"""Check mass-resolved flavor cancellation in the ordered subtracted triangle."""
from pathlib import Path
import json

import mpmath as mp
import sympy as sp


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def shape(s, a, m):
    """Euclidean kernel and its real continuation above the first threshold."""
    b = m * m
    center = m / (mp.sqrt(a) + m)
    cuts = [mp.mpf(0), center, mp.mpf(1)]
    if s < 0:
        distance = s + (mp.sqrt(a) + m)**2
        width = mp.sqrt(distance * center * (1 - center)) / (mp.sqrt(a) + m)
        cuts += [center + i * width for i in (-100, -1, 1, 100)
                 if 0 < center + i * width < 1]

    def integrand(x):
        d = a * x + b * (1 - x)
        q = x * (1 - x) * s / d
        return -mp.log1p(q) + a * x / d * q / (1 + q)

    return mp.quad(integrand, sorted(cuts))


def slope(a, m):
    b = m * m
    return -mp.quad(lambda x: b * x * (1-x)**2 / (a*x + b*(1-x))**2, [0, 1])


def main():
    mp.mp.dps = 55
    a, m = sp.symbols("a m", positive=True)
    x = sp.symbols("x", real=True)
    # Use A for the scalar mass, so a=A².
    A = sp.symbols("A", positive=True)
    threshold = (A+m)**2
    center = m/(A+m)
    d = A**2*x + m**2*(1-x)
    c = x*(1-x)
    require(sp.expand(d-c*threshold-(A+m)**2*(x-center)**2).simplify() == 0,
            "threshold minimum or Feynman-parameter routing")
    require(sp.simplify(d.subs(x, center) / c.subs(x, center) - threshold) == 0,
            "threshold is the interior minimum")
    h = a*x + m**2*(1-x)
    derivative = -m**2*x*(1-x)**2/h**2
    require(sp.simplify(derivative.subs(a, m**2).integrate((x, 0, 1))
                        + 1/(12*m**2)) == 0, "equal-mass slope")
    # Numerical approach from the real analytic interval, on the physical sheet.
    branch_controls = []
    for scalar_mass, sterile_mass in [(1, 1), (2, 3), (3, 1)]:
        A0, m0 = map(mp.mpf, (scalar_mass, sterile_mass))
        x0 = m0 / (A0 + m0)
        T = (A0 + m0)**2
        lead = -mp.pi * A0**2 * x0 / ((A0 + m0) * mp.sqrt(x0*(1-x0)))
        delta = T * mp.mpf('1e-12')
        estimated = mp.sqrt(delta) * shape(-T+delta, A0**2, m0)
        error = abs(estimated/lead-1)
        require(error < mp.mpf('1e-5'), "inverse-square-root coefficient")
        require(abs(estimated/(-lead)-1) > 1, "wrong-sign control")
        branch_controls.append({'scalar_mass':scalar_mass, 'sterile_mass':sterile_mass,
                                'relative_error_at_delta_over_T_1e-12':mp.nstr(error, 8)})

    # One active flavor, two sterile Takagi modes and a real seed.
    ybar = sp.Matrix([[1, sp.I/2]])
    mass = sp.diag(1, 4)
    residues = [ybar[i]**2 * mass[i, i] for i in range(2)]
    require(residues == [1, -1] and sum(residues) == 0, "coherent pole cancellation")
    wrong = (ybar * mass * ybar.H)[0]
    require(wrong == 2, "incoherent modulus-square mutation must lift pole zero")
    d1, d4 = (slope(mp.mpf(1), mp.mpf(j)) for j in (1, 4))
    diff_slope = d1-d4
    require(diff_slope < -mp.mpf(5)/96 and d4 < 0, "mass-resolved slope bound")
    unit = mp.mpf(1)
    response = shape(unit, unit, unit)-shape(unit, unit, 4*unit)
    require(response < 0, "nonzero finite witness at equal Euclidean and scalar squared scales")
    step = mp.mpf('1e-8')
    estimated = (shape(step, unit, unit)-shape(step, unit, 4*unit))/step
    require(abs(estimated-diff_slope) < mp.mpf('1e-8'), "finite difference slope")
    # The conjugate routing at the same external momenta gives a linear function.
    # Derive the strict mass ordering of the sum of both slopes.
    t, b = sp.symbols("t b", positive=True)
    radial = 1/(t+b) - b*t/(t+b)**3
    require(sp.cancel(sp.diff(radial, b) + (2*t**2+b**2)/(t+b)**4) == 0,
            "mass ordering of the paired slope")
    def mirror_derivative(bvalue):
        return mp.quad(lambda x: x/(x+bvalue*(1-x)), [0, 1])

    mirror_diff = mirror_derivative(1) - mirror_derivative(16)
    exact_mirror = unit/2 - (16*mp.log(16)-15)/225
    require(abs(mirror_diff-exact_mirror) < mp.mpf('1e-50'), "second routing integral")
    paired_slope = diff_slope + mirror_diff
    paired_response = response + mirror_diff  # at s=1
    require(paired_slope > 0 and paired_response > 0, "two-orientation example")
    require(diff_slope-mirror_diff < 0, "wrong conjugate-pair sign control")
    # A group of equal, nonzero masses can cancel at every momentum.
    degenerate = sp.Matrix([[1, sp.I]])
    require((degenerate * degenerate.T)[0] == 0, "degenerate group cancels")
    rotation = sp.Matrix([[3,4],[-4,3]])/5
    rotated = degenerate*rotation
    require(sp.expand((rotated*rotated.T)[0]) == 0, "real orthogonal Takagi freedom")
    # Replacing complex coherent coefficients with squared magnitudes fails.
    require(sp.simplify((rotated*rotated.H)[0]) == 2, "incoherent degenerate-group control")
    results = {
        'scope':'ordered off-shell triangle and its conjugate orientation at the same labeled momenta in the one-flavor example',
        'normalization':'a=1, Takagi masses 1 and 4, conjugated Yukawa row (1,i/2), C=1; 1/(32 pi²) stripped',
        'pole_weights':['1','-1'], 'summed_pole_weight':'0',
        'slope_at_zero':mp.nstr(diff_slope, 21),
        'rigorous_slope_upper_bound':'-5/96 (strict)',
        'finite_response_at_Q_squared_1':mp.nstr(response, 21),
        'conjugate_routing_slope_difference':mp.nstr(mirror_diff, 21),
        'paired_slope_at_zero':mp.nstr(paired_slope, 21),
        'paired_finite_response_at_Q_squared_1':mp.nstr(paired_response, 21),
        'threshold_limits':branch_controls,
    }
    Path(__file__).with_name('MASS_RESOLUTION_RESULTS.json').write_text(json.dumps(results, indent=2)+'\n')
    print('PASS: threshold geometry, three branch-coefficient limits, coherent flavor and finite witnesses')


if __name__ == '__main__':
    main()
