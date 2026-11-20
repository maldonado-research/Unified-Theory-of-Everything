"""Charged-scalar Yukawa continuum at a stable, unit-residue pole."""

import argparse
import json
from itertools import permutations
from pathlib import Path

import mpmath as mp


def near(label, value, target, tolerance=mp.mpf("1e-34")):
    if abs(value - target) > tolerance * max(abs(target), mp.mpf("1e-6")):
        raise ArithmeticError(f"{label}: {value} vs {target}")


def rho(u, a, b):
    """Continuum with |Y|^2/(16 pi^2) removed, per unit mass squared."""
    if u < b:
        return 0
    if b == a:
        return 1/u
    return ((u-b)/(u-a))**2/u


def stieltjes(z, a, b):
    # v=b/u keeps the infinite spectral interval compact.
    def integrand(v):
        ratio = 1 if a == b else (1-v)/(1-(a/b)*v)
        return ratio**2/(b-z*v)
    return mp.quad(integrand, [0, 1])


def loop_primitive(z, b):
    if z == b:
        return mp.mpf(0)
    if z == 0:
        return b
    return -(z-b)**2*mp.log1p(-z/b)/z


def subtracted_loop(z, a, b):
    slope = 0 if a == b else mp.diff(lambda x: loop_primitive(x, b), a)
    return loop_primitive(z, b)-loop_primitive(a, b)-(z-a)*slope


def parameter_loop(z, a, b):
    # Tensor reduction leaves (z-b)*B0 and an affine subtraction.
    def bubble(v, w):
        return -mp.log1p(-v*w/b)

    def integrand(v):
        at_pole = bubble(v, a)
        slope_term = 0 if a == b else (a-b)*v/(b-v*a)
        return (z-b)*bubble(v, z)-(a-b)*at_pole-(z-a)*(at_pole+slope_term)

    return mp.quad(integrand, [0, 1])


def spin_trace(u, b):
    """Evaluate tr[p.sigma k.barsigma] for a two-body CM configuration."""
    direction = mp.matrix([[2, -1-2j], [-1+2j, -2]])/3
    momentum = (u-b)/(2*mp.sqrt(u))
    sterile_energy = (u+b)/(2*mp.sqrt(u))
    p_sigma = momentum*(mp.eye(2)+direction)
    k_barsigma = sterile_energy*mp.eye(2)+momentum*direction
    result = p_sigma*k_barsigma
    return result[0, 0]+result[1, 1]


def spectator_sites():
    # The left slots belong to (electron 0, electron 1, spectator).
    # Permuting right slots assigns their scalar partners (muon 0, muon 1,
    # spectator). Fermion contractions already sew the two spectator vertices.
    counts = {pairing: 0 for pairing in permutations((0, 1))}
    vacuum_count = 0
    for partners in permutations((0, 1, 2)):
        if partners[2] == 2:
            vacuum_count += 1
            continue
        reduced = tuple(partners[2] if partner == 2 else partner for partner in partners[:2])
        counts[reduced] += 1
    if sorted(counts.values()) != [2, 2] or vacuum_count != 2:
        raise ArithmeticError("spectator-connected scalar Wick assignments")
    return list(counts.values()), vacuum_count


def verify():
    a = mp.mpf(1)
    masses = list(map(mp.mpf, [1, 16, 100, 400]))
    max_relative_error = mp.mpf(0)
    examples = []
    for b in masses:
        for scale in ["1.0001", "1.6", "11"]:
            u = b*mp.mpf(scale)
            spin = spin_trace(u, b)
            near("chiral spin trace", spin, u-b)
            phase_space = (u-b)/(8*mp.pi*u)
            # Two-point optical theorem and -Im D/pi fix the absolute weight.
            optical_density = spin*phase_space/(2*mp.pi*(u-a)**2)
            near("absolute cut", optical_density, rho(u, a, b)/(16*mp.pi**2))
            upper = loop_primitive(mp.mpc(u, "1e-45"), b)
            near("upper rim from tensor-reduced loop",
                 mp.im(upper)/(mp.pi*(u-a)**2), rho(u, a, b), mp.mpf("1e-29"))
        for z in map(mp.mpf, ["-35", "-2", "-0.125", "0"]):
            spectral = stieltjes(z, a, b)
            closed = subtracted_loop(z, a, b)/(z-a)**2
            parameters = parameter_loop(z, a, b)/(z-a)**2
            near("subtracted analytic self-energy", spectral, closed)
            near("subtracted Feynman parameter", spectral, parameters)
            if b == a:
                expected = 1/a if z == 0 else -mp.log1p(-z/a)/z
                near("threshold pole integral", spectral, expected)
            if spectral <= 0:
                raise ArithmeticError("Euclidean continuum must be positive")
            max_relative_error = max(max_relative_error, abs(closed/spectral-1),
                                     abs(parameters/spectral-1))
            if z == -2:
                examples.append({"b_over_a": int(b), "a_times_stieltjes_at_z_minus_2a": float(spectral)})

    # The charged channel excludes the other weak component. In the H
    # variables the contractions include both antisymmetric weak tensors.
    epsilon = mp.matrix([[0, 1], [-1, 0]])
    weak = epsilon*epsilon.T
    near("charged weak component", weak[0, 0], 1)
    near("off-diagonal weak component", weak[0, 1], 0)
    y, g = mp.mpc(2, -3), mp.mpc(3, 1)
    rows = [(y, 1j*mp.sqrt(2)*y), (g, -1j*g/2)]
    weights = [[float(abs(k/row[0])**2) for k in row] for row in rows]
    for row, reference in zip(weights, [[1, 2], [1, .25]]):
        for value, target in zip(row, reference):
            near("modulus weights", value, target)

    # The positive continuum at b=a gives a threshold branch while the pole
    # residue is the coefficient of the simple pole.
    for distance in [mp.mpf("1e-4"), mp.mpf("1e-12")]:
        z = a-distance
        near("residue approach at cut endpoint", distance*stieltjes(z, a, a),
             -distance*mp.log(distance/a)/z)
    counts, disconnected = spectator_sites()
    return {
        "status": "PASS",
        "pole_mass_squared_a_in_mass_units": 1,
        "density_unit": "|Y|^2/(16*pi^2)",
        "spectral_integral": "integral_b^infty du rho_reduced(u)/(u-z)",
        "samples": examples,
        "absolute_weight_factors_electron_then_muon": weights,
        "spectator_connected_insertions_per_outer_pairing": counts,
        "spectator_disconnected_assignments": disconnected,
        "spacelike_and_zero_max_relative_difference": float(max_relative_error),
        "checks": [
            "explicit two-component spin trace and integrated distinguishable phase space",
            "absorptive upper rim and per-component scalar spectral normalization",
            "twice-subtracted tensor-reduced closed loop and Feynman parameter against continuum",
            "stable pole at a branch endpoint, weak-index normalization, complex texture moduli",
            "spectator scalar Wick contractions after imposing flavor and charge",
        ],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    options = parser.parse_args()
    with mp.workdps(65):
        report = json.dumps(verify(), indent=2) + "\n"
    if options.output:
        options.output.write_text(report, encoding="utf-8")
    print(report, end="")
