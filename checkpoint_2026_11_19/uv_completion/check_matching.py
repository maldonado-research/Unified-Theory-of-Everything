"""Tree matching of the two Majorana pairs in original action coordinates."""

import argparse
import json
import math
from fractions import Fraction as F
from pathlib import Path


def equal(label, actual, expected):
    if actual != expected:
        raise ArithmeticError(f"{label}: {actual} versus {expected}")


def near(label, actual, expected):
    if abs(actual - expected) > 2e-13 * max(1, abs(expected)):
        raise ArithmeticError(f"{label}: {actual} versus {expected}")


def electron_poles(z, mass):
    """Squared Yukawas divided by y**2 are +1 and -2."""
    return mass / (mass**2 - z) - 4 * mass / (4 * mass**2 - z)


def moment(n, mass):
    return sum(coupling / pole_mass**(2*n + 1)
               for coupling, pole_mass in [(1, mass), (-2, 2*mass)])


def check_transfer(mass):
    c7 = moment(1, mass)
    c9 = moment(2, mass)
    for q in (F(0), F(1, 3), F(2), F(3), 10*mass):
        z = -q**2
        # For a spatial transfer along the third axis, the spin-up sector
        # of slash(q)-M is [[-M,-q],[q,-M]]. Extract minus its inverse's
        # upper diagonal entry by the two-by-two cofactor formula.
        from_inverse = F(0)
        for coeff, pole_mass in [(1, mass), (-2, 2*mass)]:
            k11, k12, k21, k22 = -pole_mass, -q, q, -pole_mass
            from_inverse -= coeff * k22 / (k11*k22 - k12*k21)
        full = electron_poles(z, mass)
        equal("free inverse versus pole sum", from_inverse, full)
        x = -z / mass**2
        denominator = 4 + 5*x + x**2
        factor = 4 / denominator
        equal("exact transfer", full, c7*z*factor)
        equal("dimension-seven remainder", c7*z - full,
              c7*z*x*(5+x)/denominator)
        equal("dimension-nine remainder", full - c7*z - c9*z**2,
              c7*z*x**2*(21 + 5*x)/(4*denominator))
        if not (0 < factor <= 1 and 1-factor <= 5*x/4
                and 0 <= factor - (1 - 5*x/4) <= 21*x**2/16):
            raise ArithmeticError(f"spacelike error bounds at x={x}")


def check_projection(mass):
    a = F(1)
    c7, c9 = moment(1, mass), moment(2, mass)
    traces = []
    for beta in [F(1, 3), F(3, 5), F(4, 5)]:
        s = 4*a / (1-beta**2)
        h, v = s/2-a, s*beta/2
        for cosine in [F(-1), F(0), F(2, 3), F(1)]:
            t, u = -h+v*cosine, -h-v*cosine
            equal("two assignments", t + u, 2*a-s)
            equal("normalized d7 tree", c7*(t+u)/2, c7*(a-s/2))
            equal("d9 projection", c9*(t*t+u*u)/2, c9*(h*h+v*v*cosine**2))
            projected = (electron_poles(t, mass) + electron_poles(u, mass))/2
            multiplier = projected / (-h*c7)
            xmax = (h+v)/mass**2
            if not (4/(4+5*xmax+xmax**2) <= multiplier <= 1):
                raise ArithmeticError("two-assignment transfer bound")
        traces.append((s, -h*c7))
    for (s1, f1), (s2, f2) in zip(traces, traces[1:]):
        equal("physical slope", (f2-f1)/(s2-s1), -c7/2)


def check_complex_trees():
    y, g = 2+1j, 1-1j
    e_columns, e_masses = (y, 1j*math.sqrt(2)*y), (10, 20)
    mu_columns, mu_masses = (g, -1j*g/2), (1, 4)
    mu_weights = [m*entry.conjugate()**2 for m, entry in zip(mu_masses, mu_columns)]
    equal("muon ultraviolet weight sum", sum(mu_weights), 0j)
    equal("muon phase", mu_weights[0], g.conjugate()**2)
    equal("retained muon inverse mass sum",
          sum(entry.conjugate()**2/m for m, entry in zip(mu_masses, mu_columns)),
          15*g.conjugate()**2/16)
    b12, b34 = 1-3j, 2+4j
    mu_routes, e_routes = (-1.0, -5.0), (-2.0, -4.0)
    # Sew both physical trees, including the two different scalar assignments.
    def incoming(columns, masses, routes, spin):
        return spin*sum(m*entry.conjugate()**2/(z-m*m)
                        for z in routes for m, entry in zip(masses, columns))
    direct = (incoming(e_columns, e_masses, e_routes, b34).conjugate()
              * incoming(mu_columns, mu_masses, mu_routes, b12))
    reduced_e = y**2 * sum(electron_poles(z, 10) for z in e_routes)
    reduced_mu = sum(w/(z-m*m) for z in mu_routes for m, w in zip(mu_masses, mu_weights))
    near("ordered W and literal positive action", direct,
         (-b34.conjugate()*b12)*reduced_e*reduced_mu)
    inverse = (incoming(mu_columns, mu_masses, mu_routes, b12).conjugate()
               * incoming(e_columns, e_masses, e_routes, b34))
    near("reciprocal phase", inverse, direct.conjugate())
    # Preserve sensitivity to conjugating just the electron source.
    for coefficient in [reduced_e, reduced_mu, reduced_e*reduced_mu]:
        if abs(coefficient - coefficient.conjugate()) < 1e-7:
            raise ArithmeticError("conjugation control became real")


def verify():
    for mass in [F(1), F(7, 3), F(10)]:
        for n in range(8):
            equal("mass moment", moment(n, mass), (1-F(1, 4)**n)/mass**(2*n+1))
        check_transfer(mass)
        check_projection(mass)
        delta = F(1, 7)
        c5 = 1/mass - 2*(1+delta)/(2*mass)
        c7 = 1/mass**3 - 2*(1+delta)/(2*mass)**3
        equal("detuned d5", c5, -delta/mass)
        equal("detuned d7", c7, (3-delta)/(4*mass**3))
    check_complex_trees()
    return {
        "status": "PASS",
        "electron_C5_C7_C9_in_units_y2_over_mass_power":
            [str(moment(n, F(1))) for n in range(3)],
        "slope_in_units_y2_over_M3": str(-moment(1, F(1))/2),
        "muon_residues_over_m_times_conjugate_g_squared": [1, -1],
        "spacelike_ratio_at_minus_z_equal_M_squared":
            str(electron_poles(F(-1), F(1)) / (-moment(1, F(1)))),
        "checks": [
            "positive-mass moments n=0..7, independent spacelike Dirac inverse",
            "exact exchange and rational remainders through dimensions seven and nine",
            "two scalar routings, energy slope, angular term and uniform transfer bound",
            "complex Yukawas, ordered spin factors and reciprocal tree product",
            "fractional detuning of second electron Yukawa square",
        ],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="also write the JSON receipt here")
    args = parser.parse_args()
    receipt = json.dumps(verify(), indent=2) + "\n"
    if args.output is not None:
        args.output.write_text(receipt, encoding="utf-8")
    print(receipt, end="")
