#!/usr/bin/env python3
"""Sew the two scalar trees in CM four-vectors, then project the box cut."""
import argparse
import json
import math
from pathlib import Path


def gauss_legendre(n):
    """Generate the symmetric Gauss-Legendre rule on [-1,1]."""
    rule = []
    for i in range(n // 2):
        z = math.cos(math.pi * (i + 0.75) / (n + 0.5))
        for _ in range(24):
            prev, cur = 1.0, z
            for k in range(2, n + 1):
                prev, cur = cur, ((2*k - 1)*z*cur - (k - 1)*prev)/k
            derivative = n*(prev - z*cur)/(1 - z*z)
            correction = cur/derivative
            if abs(correction) < 3e-16:
                break
            z -= correction
        w = 2/((1 - z*z)*derivative**2)
        rule.extend([(-z, w), (z, w)])
    return sorted(rule)


def transfer(p, k):
    return (p[0] - k[0])**2 - sum((p[i] - k[i])**2 for i in range(1, 4))


def tree(p, pair, species):
    """Amputated LL -> scalar pair with both identical-scalar assignments."""
    return math.fsum(
        residue/(transfer(p, k) - mass**2)
        for mass, residue in species for k in pair
    )


def scalar_sphere(s, rule, azimuths):
    energy = math.sqrt(s)/2
    momentum = math.sqrt(energy**2 - 1)
    p = (energy, 0., 0., energy)
    points = []
    for x, weight in rule:
        transverse = momentum*math.sqrt(1 - x*x)
        for j in range(azimuths):
            phi = 2*math.pi*(j + 0.5)/azimuths
            k1 = (energy, transverse*math.cos(phi), transverse*math.sin(phi), momentum*x)
            k2 = (energy, -k1[1], -k1[2], -k1[3])
            pair = (k1, k2)
            points.append((pair, weight/(2*azimuths)*tree(p, pair, [(1., 1.), (4., -1.)])))
    return points


def sewn_cut(s, cosine, points):
    energy = math.sqrt(s)/2
    p3 = (energy, energy*math.sqrt(1 - cosine**2), 0., energy*cosine)
    electron = [(10., 10.), (20., -40.)]
    angular = math.fsum(w*tree(p3, pair, electron) for pair, w in points)
    phase_space = math.sqrt(1 - 4/s)/(8*math.pi)
    # B34* B12 = s by choice of external phases.
    return s * (0.5/2) * phase_space * angular


def analytic(s):
    beta = math.sqrt(1 - 4/s)
    h, v = s/2 - 1, s*beta/2
    def j(b):
        return math.log1p(2*v/(b + 1/(h+v)))/(2*v)

    mu = j(1) - j(16)
    electron = 10*(j(100) - 4*j(400))
    absorptive = s*beta*mu*electron/(8*math.pi)
    return absorptive, mu, electron


def cut_by_combined_denominators(s, c, parameter_rule):
    """Combine two spacelike transfers before integrating the scalar sphere."""
    beta = math.sqrt(1 - 4/s)
    h, v = s/2 - 1, s*beta/2

    def kernel(b, d, z):
        terms = []
        for x, weight in parameter_rule:
            u = (1 + x)/2
            shift = h + u*b + (1-u)*d
            length2 = 1 - 2*u*(1-u)*(1-z)
            terms.append(weight/2/(shift**2 - v**2*length2))
        return math.fsum(terms)

    box_sum = math.fsum(
        w_mu*w_e*2*(kernel(b, d, c) + kernel(b, d, -c))
        for b, w_mu in [(1., 1.), (16., -1.)]
        for d, w_e in [(100., 10.), (400., -40.)]
    )
    return s*beta*box_sum/(32*math.pi)


def rounded(value):
    """Keep reproducible digits in the portable receipt."""
    return float(format(value, ".10g"))


def verify():
    rule = gauss_legendre(28)
    parameter_rule = gauss_legendre(160)
    rows = []
    for s in [5., 6., 8., 9., 10., 11.]:
        points = scalar_sphere(s, rule, 32)
        cuts = [(c, w/2, sewn_cut(s, c, points)) for c, w in rule]
        numeric = math.fsum(w*cut for _, w, cut in cuts)
        exact, mu, electron = analytic(s)
        if not math.isclose(numeric, exact, rel_tol=4.e-12):
            raise ArithmeticError(f"Projected sewing failed at {s}: {numeric} vs {exact}")
        if not (mu > 0 and electron < 0):
            raise ArithmeticError("Projected tree signs failed")
        # Verify the fixed-angle paired boxes carry higher waves.
        quadrupole = math.fsum(w*(3*c*c - 1)/2*cut for c, w, cut in cuts)
        mean_square = math.fsum(w*cut**2 for _, w, cut in cuts)
        if mean_square <= numeric**2 or quadrupole == 0:
            raise ArithmeticError("Expected a higher-wave contribution")
        floor = exact**2/(32*math.pi*s)
        at_right_angle = sewn_cut(s, 0., points)
        combined = cut_by_combined_denominators(s, 0., parameter_rule)
        if not math.isclose(at_right_angle, combined, rel_tol=5.e-11):
            raise ArithmeticError(f"Fixed-angle box sewing failed at {s}")
        row = dict(s_over_m2=s, cut_at_c0=at_right_angle, m_F_mu_mean=mu, m_F_e_mean=electron,
                         J0_cut_analytic=exact, J0_cut_from_boxes=numeric,
                         J2_moment_from_boxes=quadrupole,
                         m2_sigma_floor_J0=floor,
                         m2_sigma_floor_full_cut=mean_square/(32*math.pi*s))
        rows.append({name: rounded(number) for name, number in row.items()})
    # Check dimensional factors and the coefficient of beta^2 at opening.
    near, _, _ = analytic(4.00000004)
    coefficient = 4*(15/34 * 30/(101*401))**2/(2048*math.pi**3)
    threshold_check = near**2/(32*math.pi*4.00000004*(1-4/4.00000004))
    if not math.isclose(threshold_check, coefficient, rel_tol=1.e-7):
        raise ArithmeticError("Threshold normalization failed")
    return dict(m=1, M_over_m=10, scalar_mass_over_m=1,
                squared_yukawa_strengths=[1, 1], B34_conjugate_times_B12_over_s=1,
                polar_order=len(rule), azimuthal_order=32,
                m2_floor_J0_over_beta2_at_threshold=rounded(coefficient), rows=rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write the JSON receipt to this file")
    args = parser.parse_args()
    payload = json.dumps(verify(), indent=2, allow_nan=False) + "\n"
    if args.output is None:
        print(payload, end="")
    else:
        args.output.write_text(payload, encoding="utf-8")


if __name__ == "__main__":
    main()
