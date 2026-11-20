"""Verify the internal sterile self-energy term in the J=0 scalar-pair cut."""

import argparse
import json
from pathlib import Path

from mpmath import mp


def require_close(result, expected, tolerance, label):
    if not abs(result - expected) <= tolerance * abs(expected):
        raise ArithmeticError(f"{label}: {result} vs {expected}")


def bubble(q, log_scale=0):
    """Twice the finite vector-bubble coefficient; m=1 and q=-p²."""
    return 2 * log_scale + 2 + 1 / q - ((1 + q) / q)**2 * mp.log1p(q)


def channels():
    return {
        "electron": ((10, 20), (1, mp.j * mp.sqrt(2))),
        "muon": ((1, 4), (1, -mp.j / 2)),
    }


def transfer_terms(q, masses, row):
    denominators = [m**2 + q for m in masses]
    tree = sum(m * y**2 / d for m, y, d in zip(masses, row, denominators))
    kinetic = sum(abs(y)**2 / d for y, d in zip(row, denominators))
    coefficient = -2 * q * kinetic * tree / (16 * mp.pi**2)
    return tree.real, (coefficient * bubble(q)).real, (2 * coefficient).real


def cut(s, rule):
    beta = mp.sqrt(1 - mp.mpf(4) / s)
    h = mp.mpf(s) / 2 - 1
    v = s * beta / 2
    means = {}
    for flavor, (masses, row) in channels().items():
        mean = [0, 0, 0]
        for x, weight in rule:
            terms = transfer_terms(h - v * x, masses, row)
            for k, value in enumerate(terms):
                mean[k] += weight * value / 2
        log_tree = sum(m * y**2 * mp.log1p(2 * v / (m**2 + h - v)) / (2 * v)
                       for m, y in zip(masses, row))
        require_close(mean[0], log_tree, mp.mpf("1e-23"), "tree angular mean")
        means[flavor] = mean
    e, mu = means["electron"], means["muon"]
    norm = s * beta / (8 * mp.pi)
    return {"leading_over_gamma": norm * e[0] * mu[0],
            "electron": norm * e[1] * mu[0], "muon": norm * e[0] * mu[1],
            "d_ln_mu_electron": norm * e[2] * mu[0],
            "d_ln_mu_muon": norm * e[0] * mu[2]}


def propagator_tree(q, mass, row, kinetic_shift):
    """Direct chiral inverse at spacelike p; row is outgoing in these coordinates."""
    n = mass.rows
    sqrt_q = mp.sqrt(q)
    H = row.H * row
    K = mp.zeros(2 * n)
    K[:n, :n] = -mass.H
    K[n:, n:] = -mass
    K[:n, n:] = sqrt_q * (mp.eye(n) + kinetic_shift * H)
    K[n:, :n] = -sqrt_q * (mp.eye(n) + kinetic_shift * H.T)
    inverse = K**-1
    return -(row * inverse[n:, n:] * row.T)[0]


def verify_propagator():
    """Check the insertion sign, flavor phases, finite bubble and mass/Y RG."""
    angle = mp.mpf(".47")
    phase = mp.exp(mp.j * mp.mpf(".63"))
    U = mp.matrix([[mp.cos(angle), phase * mp.sin(angle)],
                   [-phase.conjugate() * mp.sin(angle), mp.cos(angle)]])
    U *= mp.diag([mp.exp(mp.j / 5), mp.exp(-mp.j / 7)])
    active_phase = mp.exp(mp.mpf(".43") * mp.j)
    step = mp.mpf("1e-13")
    loop = 16 * mp.pi**2

    for masses, couplings in channels().values():
        mass = U.T * mp.diag(masses) * U
        row = active_phase * mp.matrix([couplings]) * U
        H = row.H * row
        dy = row * H
        dm = H.T * mass + mass * H
        for q in map(mp.mpf, (".2", "2", "7")):
            tree, finite, slope = transfer_terms(q, masses, couplings)
            # Integrate the original Feynman parameter with the doublet trace.
            integrated = -2 * mp.quad(lambda x: (1-x) * mp.log((1-x)*(1+x*q)), [0, 1])
            require_close(bubble(q), integrated, mp.mpf("1e-37"), "bubble integral")
            core = propagator_tree(q, mass, row, 0)
            require_close(core, active_phase**2 * tree, mp.mpf("1e-37"), "basis change")
            def kinetic_derivative(f):
                plus = propagator_tree(q, mass, row, step * f / loop)
                minus = propagator_tree(q, mass, row, -step * f / loop)
                return (plus - minus) / (2 * step)
            require_close(kinetic_derivative(integrated), active_phase**2 * finite,
                          mp.mpf("1e-21"), "inverse-propagator insertion")
            shifted = lambda t: propagator_tree(q, mass + t * dm, row + t * dy, 0)
            implicit = (shifted(step) - shifted(-step)) / (2 * step * loop)
            require_close(implicit, -active_phase**2 * slope,
                          mp.mpf("1e-21"), "sterile part of running")
            # Changing ln(mu) changes only the subtracted scalar bubble at this order.
            L = mp.mpf(".27")
            actual_slope = (kinetic_derivative(bubble(q, L)) -
                            kinetic_derivative(bubble(q, -L))) / (2 * L)
            require_close(actual_slope, -implicit, mp.mpf("1e-21"), "scale cancellation")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    mp.dps = 48
    verify_propagator()
    rules = [list(zip(*mp.gauss_quadrature(order, "legendre"))) for order in (40, 64)]
    rows = []
    for s in (6, 8, 10):
        coarse, result = [cut(s, rule) for rule in rules]
        for key, value in result.items():
            require_close(coarse[key], value, mp.mpf("1e-20"), "quadrature refinement")
        rows.append({"s_over_m2": s, **{k: mp.nstr(v, 16) for k, v in result.items()}})
    receipt = {
        "sector": "two-scalar cut with one renormalized internal sterile self-energy insertion",
        "scheme": "MSbar Yukawas and four sterile masses at mu=m, fixed scalar pole mass m=1",
        "normalization": "A0_insert/gamma=|y_e|^2 electron + |g|^2 muon; B34* B12=s",
        "scale": "explicit derivative with respect to ln(mu) at fixed renormalized inputs",
        "checks": {"precision_decimal_digits": mp.dps, "angular_orders": [40, 64],
                   "angular_refinement_relative_tolerance": "1e-20",
                   "propagator_and_RG_relative_tolerance": "1e-21",
                   "spacelike_transfers_over_m2": ["0.2", "2", "7"],
                   "coordinate_test": "complex sterile unitary and external row rephasing"},
        "rows": rows}
    data = json.dumps(receipt, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(data, encoding="utf-8")
    else:
        print(data, end="")


if __name__ == "__main__":
    main()
