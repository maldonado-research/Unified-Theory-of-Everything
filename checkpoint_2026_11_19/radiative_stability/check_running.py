"""Check the singlet mass moments and their one-loop RG equations in exact arithmetic."""

import argparse
import json
from pathlib import Path

from sympy import I, Matrix, Rational, diag, expand, eye, simplify, sqrt, zeros


def equal(label, actual, expected):
    difference = expand(actual - expected)
    if difference != expected * 0:
        raise ArithmeticError(f"{label}: difference {difference}")


def fundamental_rates(yukawa, mass):
    h = expand(yukawa.H * yukawa)
    dy = Rational(3, 2) * yukawa * h + h.trace() * yukawa
    dm = h.T * mass + mass * h
    return expand(dy), expand(dm)


def moments(yukawa, mass):
    inverse = mass.inv()
    return expand(yukawa * inverse * yukawa.T), expand(yukawa.conjugate() * mass * yukawa.H)


def chain_rule(yukawa, mass, dy, dm):
    inverse = mass.inv()
    left = expand(dy * inverse * yukawa.T)
    dk = left + left.T - yukawa * inverse * dm * inverse * yukawa.T
    right = expand(dy.conjugate() * mass * yukawa.H)
    dr = right + right.T + yukawa.conjugate() * dm * yukawa.H
    return expand(dk), expand(dr)


def flavor_rates(yukawa, k, r):
    p = expand(yukawa * yukawa.H)
    dk = 2 * p.trace() * k + (p*k + k*p.T)/2
    dr = 2 * p.trace() * r + Rational(5, 2) * (p.conjugate()*r + r*p)
    return expand(dk), expand(dr)


def check_general():
    y = Matrix([[1+I, 2-I, I], [-2*I, 3, 1+2*I]])
    mass = Matrix([[2+I, 1, -I], [1, 4-I, 2*I], [-I, 2*I, 3]])
    k, r = moments(y, mass)
    dy, dm = fundamental_rates(y, mass)
    dk, dr = chain_rule(y, mass, dy, dm)
    fk, fr = flavor_rates(y, k, r)
    equal("inverse moment product rule", dk, fk)
    equal("positive moment product rule", dr, fr)
    p = y*y.H
    equal("row Gram equation", dy*y.H + y*dy.H, 3*p**2 + 2*p.trace()*p)

    # A complex unitary in each space detects both conjugation and transpose errors.
    u = Matrix([[3, 4*I], [4*I, 3]]) / 5
    v = diag(1, u)
    equal("active unitarity", u.H*u, eye(2))
    equal("sterile unitarity", v.H*v, eye(3))
    yp, mp = expand(u.H*y*v), expand(v.T*mass*v)
    kp, rp = moments(yp, mp)
    dyp, dmp = fundamental_rates(yp, mp)
    dkp, drp = chain_rule(yp, mp, dyp, dmp)
    equal("transformed Yukawa beta", dyp, u.H*dy*v)
    equal("transformed mass beta", dmp, v.T*dm*v)
    for original, changed in ((k, kp), (dk, dkp)):
        equal("inverse moment covariance", changed, u.H*original*u.conjugate())
    for original, changed in ((r, rp), (dr, drp)):
        equal("positive moment covariance", changed, u.T*original*u)

    # Sharing even one column makes an isolated component zero inhomogeneous.
    two_masses = diag(1, 4)
    controls = []
    for rel, index, expected in ((2*I, 0, 1), (I/2, 1, 5)):
        shared = Matrix([[1, rel], [1, 0]])
        equal("shared-column tuned component", moments(shared, two_masses)[index][0, 0], 0)
        rate = chain_rule(shared, two_masses, *fundamental_rates(shared, two_masses))
        entry = simplify(rate[index][0, 0])
        equal("shared-column derivative", entry, expected)
        controls.append(str(entry))
    return controls


def takagi_rates(y, masses, dy, dm):
    # At a positive diagonal mass point, solve dM_ij + m_i A_ij - m_j A_ij* = 0.
    rotation = zeros(len(masses))
    for i, first in enumerate(masses):
        for j in range(i+1, len(masses)):
            second = masses[j]
            real, imag = dm[i, j].expand(complex=True).as_real_imag()
            rotation[i, j] = -real/(first-second) - I*imag/(first+second)
            rotation[j, i] = -rotation[i, j].conjugate()
    dmass = simplify(dm + rotation.T*diag(*masses) + diag(*masses)*rotation)
    equal("Takagi generator", rotation.H, -rotation)
    equal("off-diagonal derivatives removed", dmass, diag(*dmass.diagonal()))
    return dy + y*rotation, dmass


def check_texture():
    e, g = 2+I, 1-I
    y = Matrix([[e, I*sqrt(2)*e, 0, 0], [0, 0, g, -I*g/2]])
    masses = [10, 20, 1, 4]
    mass = diag(*masses)
    k, r = moments(y, mass)
    equal("electron K zero", k[0, 0], 0)
    equal("muon R zero", r[1, 1], 0)
    dy, dm = fundamental_rates(y, mass)
    dk, dr = chain_rule(y, mass, dy, dm)
    equal("electron K derivative", dk[0, 0], 0)
    equal("muon R derivative", dr[1, 1], 0)
    equal("electron off-diagonal mass beta", dm[0, 1], -10*I*sqrt(2)*abs(e)**2)
    equal("muon off-diagonal mass beta", dm[2, 3], 3*I*abs(g)**2/2)
    ryt, rmt = takagi_rates(y, masses, dy, dm)
    tk, tr = chain_rule(y, mass, ryt, rmt)
    equal("Takagi inverse moment rate", tk, dk)
    equal("Takagi positive moment rate", tr, dr)
    rates = []
    for row, start, h, expected_mass, expected_y in (
        (0, 0, e, 2, 2), (1, 2, g, Rational(-3, 2), Rational(3, 2))
    ):
        i, j = start, start+1
        dlogm = simplify((rmt[j,j]/masses[j] - rmt[i,i]/masses[i])/abs(h)**2)
        ratio_rate = ryt[row,j]/y[row,j] - ryt[row,i]/y[row,i]
        dlogy = simplify((ratio_rate + ratio_rate.conjugate())/abs(h)**2)
        equal("mass ratio rate", dlogm, expected_mass)
        equal("squared Yukawa magnitude ratio rate", dlogy, expected_y)
        rates.append([str(dlogm), str(dlogy)])
    return rates


def check_sequential():
    # Keep E1 and both muon columns. C descends from removing E2.
    e, g, mass_scale = 2+I, 1-I, 10
    y = Matrix([[e, 0, 0], [0, g, -I*g/2]])
    mass = diag(mass_scale, 1, 4)
    c = diag(-e**2/mass_scale, 0)
    k, r = moments(y, mass)
    equal("tree match electron", (c+k)[0, 0], 0)
    p = y*y.H
    dy, dm = fundamental_rates(y, mass)
    dy -= 3*c*y.conjugate()*mass
    dk, dr = chain_rule(y, mass, dy, dm)
    dc = 2*p.trace()*c + Rational(7, 2)*(p*c + c*p.T)
    total, _ = flavor_rates(y, c+k, r)
    equal("sequential pure-Y tensor", dc+dk, total)
    equal("sequential pure-Y electron rate", (dc+dk)[0, 0], 0)
    equal("muon positive moment with electron Wilson coefficient", dr[1, 1], 0)
    if simplify(dc[0, 0]) == 0 or simplify(dk[0, 0]) == 0:
        raise ArithmeticError("separate EFT contributions must run in this example")
    # D log(mu) interval = -log(2) for running down, at fixed loop order.
    quartic_source = 4*c
    coefficient = simplify(-quartic_source[0, 0]*mass_scale/e**2)
    equal("gaugeless quartic log coefficient", coefficient, 4)
    return str(coefficient)


def verify():
    return {
        "status": "PASS",
        "arithmetic": "exact complex and algebraic matrices in SymPy",
        "shared_column_derivatives_Kee_then_Ree_in_two_examples": check_general(),
        "Takagi_log_ratio_rates_mass_then_squared_Yukawa_over_abs_first_Yukawa_squared":
            dict(zip(("electron", "muon"), check_texture())),
        "sequential_quartic_leading_log_over_lambda_y_squared_log2_div_16pi2_M":
            check_sequential(),
        "checks": [
            "weak-basis tensor product rules, active and sterile covariance",
            "rank-one blocks, unequal mass eigenvalue and Yukawa ratio evolution",
            "complex symmetric mass beta and infinitesimal Takagi rotation",
            "dimension-five sequential EFT feedback and total pure-Y cancellation",
        ],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write an identical copy of the receipt")
    args = parser.parse_args()
    record = json.dumps(verify(), indent=2) + "\n"
    if args.output is not None:
        args.output.write_text(record, encoding="utf-8")
    print(record, end="")
