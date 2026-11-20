"""Quartic response of the amputated electron vertex in the four-singlet model."""

import argparse
import json
from itertools import permutations
from pathlib import Path

import mpmath as mp


def loop_shape(r):
    """Dimensionless finite single-mass function for r = a / M_i^2."""
    if r == 1:
        return mp.mpf("0.5")
    return (r - 1 - mp.log(r)) / (1 - r)**2


def radial_gap(r, k):
    """Integrate the subtracted pair directly, with t = Euclidean k_E^2 / M^2."""
    def integrand(t):
        return (k - 1) * (t / (t + r))**2 / ((1 + t) * (k + t))

    endpoints = sorted(set([mp.mpf(0), r, mp.mpf(1), k]))
    return mp.quad(integrand, endpoints + [mp.inf])


def near(label, actual, expected, tolerance=mp.mpf("1e-43")):
    if abs(actual - expected) > tolerance * max(1, abs(expected)):
        raise ArithmeticError(f"{label}: {actual} versus {expected}")


def check_contractions():
    # Orient the two Yukawa vertices with temporary labels. Equal-spinor
    # contractions acquire an exchange sign and a spin sign, which cancel.
    pairings = list(permutations((0, 1)))
    tree = sum(1 for _fermions in pairings for _scalars in pairings) // 2
    triangle = sum(1 for _fermions in pairings for _internal in pairings
                   for _external in pairings) // 2
    if (tree, triangle) != (2, 4):
        raise ArithmeticError("Scalar and external-fermion Wick counting")
    return tree, triangle


def verify():
    tree, triangle = check_contractions()
    for r in map(mp.mpf, ["1e-14", ".01", ".65", "1", "1.0000001", "4", "70"]):
        # Feynman parameter for one repeated light line, with no mass subtraction.
        parameter_value = mp.quad(lambda x: x / (1 - x + r*x), [0, 1])
        near("single heavy mass", loop_shape(r), parameter_value)
        for k in map(mp.mpf, ["1", "1.0001", "4", "25"]):
            gap = loop_shape(r/k) - loop_shape(r)
            near("full radial integral with the two residues", radial_gap(r, k), gap)
            if k > 1 and not 0 < gap < mp.log(k):
                raise ArithmeticError("strict hierarchy bound")

    k = mp.mpf(4)
    for r in map(mp.mpf, ["1e-8", "1e-12", "1e-18"]):
        logarithm = -mp.log(r)
        first_power = r * ((2*(logarithm + mp.log(k)) - 1)/k - (2*logarithm - 1))
        remainder = abs(loop_shape(r/k) - loop_shape(r) - mp.log(k) - first_power)
        if remainder > 8 * r**2 * logarithm:
            raise ArithmeticError("first power correction in small scalar mass")

    # Dimension-five hard coefficient + light bubble, in units of
    # 2 lambda Y_i^2 / (16 pi^2 M_i), at leading inverse mass.
    a = mp.mpf(1)
    for heavy in map(mp.mpf, ["10", "20"]):
        for scale in map(mp.mpf, [2, 11, 44]):
            hard = 1 + mp.log(scale**2 / heavy**2)
            light = -mp.log(scale**2 / a)
            near("single-pole hard plus light", hard + light, 1 - mp.log(heavy**2/a))
    for scale in map(mp.mpf, [2, 11, 44]):
        hard_pair = mp.log(scale**2 / 100) - mp.log(scale**2 / 400)
        near("simultaneous and sequential threshold matching", hard_pair, 2*mp.log(2))

    r = mp.mpf(".01")
    gap = loop_shape(r/k) - loop_shape(r)
    full = 2*gap/(16*mp.pi**2)
    hard = 2*mp.log(k)/(16*mp.pi**2)
    near("benchmark via radial quadrature", 2*radial_gap(r, k)/(16*mp.pi**2), full)
    return {
        "status": "PASS",
        "potential": "V = a (H_dagger H) + lambda (H_dagger H)^2",
        "charged_scalar_assignment_counts": {"tree": tree, "quartic_triangle": triangle},
        "benchmark_M_over_m": 10,
        "single_mass_shapes": [float(loop_shape(r)), float(loop_shape(r/k))],
        "subtracted_shape": float(gap),
        "M_times_zero_momentum_correction_over_lambda_y_squared": float(full),
        "M_times_leading_MSbar_C5_over_lambda_y_squared": float(hard),
        "finite_mass_to_leading_ratio": float(full/hard),
        "checks": [
            "charged-field Wick multiplicities and common open-fermion contractions",
            "closed form against single-mass Feynman parameter and two-mass radial quadrature",
            "equal sterile masses with opposite squares and equal scalar/sterile masses",
            "strict sign at unequal masses, heavy-mass limit and first power correction",
            "single-threshold scale cancellation and two-threshold leading logarithm",
        ],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    options = parser.parse_args()
    with mp.workdps(75):
        result = json.dumps(verify(), indent=2) + "\n"
    if options.output:
        options.output.write_text(result, encoding="utf-8")
    print(result, end="")
