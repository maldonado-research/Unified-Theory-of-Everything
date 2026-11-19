"""Checks on the electron threshold and on the single-quartic graph count."""

import argparse
import itertools
import json
import math
from collections import Counter
from pathlib import Path


def close(label, actual, expected, tolerance=2e-12):
    if abs(actual - expected) > tolerance * max(1, abs(expected)):
        raise ArithmeticError(f"{label}: {actual} versus {expected}")


def connected(edges):
    neighbors = [set() for _ in range(4)]
    for u, v in edges:
        neighbors[u].add(v)
        neighbors[v].add(u)
    seen = {0}
    frontier = [0]
    while frontier:
        new = neighbors[frontier.pop()] - seen
        seen.update(new)
        frontier.extend(new)
    return len(seen) == 4


def yukawa_contractions():
    # The conjugate vertex has label 3. Its three internal legs H,L,N
    # attach to h,l,n; the two remaining N legs have equal chirality.
    counts = Counter()
    for h, l, n in itertools.product(range(3), repeat=3):
        a, b = sorted({0, 1, 2} - {n})
        edges = [(3, h), (3, l), (3, n), (a, b)]
        if not connected(edges):
            counts["disconnected"] += 1
            continue
        if any(not connected(edges[:i] + edges[i + 1:]) for i in (0, 1)):
            counts["light_reducible"] += 1
            continue
        if n == h or n == l:
            raise ArithmeticError("mixed sterile line has no external pair at endpoint")
        counts["light_irreducible_with_mixed_sterile_bridge"] += 1
    expected = {"disconnected": 3, "light_reducible": 12,
                "light_irreducible_with_mixed_sterile_bridge": 12}
    if counts != expected:
        raise ArithmeticError(counts)
    return dict(sorted(counts.items()))


def charged_scalar_contractions():
    # Both lists contain two external-Yukawa scalar slots, then two quartic
    # slots. Match the positive to negative slots in all 4! possible orders.
    counts = Counter()
    for negative in itertools.permutations(range(4)):
        quartic_self_pairs = sum(negative[p] >= 2 for p in (2, 3))
        counts[quartic_self_pairs] += 1
    if counts != {0: 4, 1: 16, 2: 4}:
        raise ArithmeticError(counts)
    # Amplitude convention: iM, four Yukawa vertices identical to the box;
    # stripping the quartic coupling, return the product of propagator,
    # quartic-vertex, and Wick-rotated loop phases, divided by i.
    box_factor = (1j**4) * (1j * (-1)**4) / 1j
    double_triangle_factor = (1j**6) * (-1j * counts[0]) * (1j * (-1)**3)**2 / 1j
    close("sewn quartic relative to the two-box normalization",
          double_triangle_factor / box_factor, -4)
    return {str(i): counts[i] for i in sorted(counts)}


def matching_contributions(masses, columns, mu):
    """Equation (88) in Zhang and Zhou v3, restricted to electron entry.

    Values are coefficients of 16 pi^2 times the correction; the last
    entry is also divided by lambda in V=lambda (H^dag H)^2.
    """
    tree_row = [sum(v[0] * v[f] / m for m, v in zip(masses, columns))
                for f in (0, 1)]
    higgs_field = 0.
    lepton_row = [0j, 0j]
    quartic = 0j
    for m, v in zip(masses, columns):
        log = 2 * math.log(mu / m)
        higgs_field += 0.5 * (1 + 2*log) * sum(abs(z)**2 for z in v)
        for f in (0, 1):
            lepton_row[f] += (3 + 2*log) * v[0] * v[f].conjugate() / 4
        quartic += 2 * (1 + log) * v[0]**2 / m
    return (-higgs_field * tree_row[0],
            -sum(z * c for z, c in zip(lepton_row, tree_row)), quartic)


def check_texture():
    y, g = 0.75 - 0.5j, -0.25 + 1.25j
    mass = 10.
    masses = (mass, 2*mass, 1., 4.)
    columns = [(y, 0), (1j*math.sqrt(2)*y, 0), (0, g), (0, -1j*g/2)]
    expected = 2 * math.log(4) * y**2 / mass
    for mu in (0.7, 10., math.sqrt(2)*10, 20., 113.):
        for end in (2, 4):
            higgs, lepton, quartic = matching_contributions(
                masses[:end], columns[:end], mu)
            close("Higgs field electron entry", higgs, 0.)
            close("lepton field electron entry", lepton, 0.)
            close("common-scale quartic threshold", quartic, expected)
    # If Y_e2 is detuned both canonical factors genuinely respond.
    shifted = list(columns)
    shifted[1] = (1.1 * shifted[1][0], 0)
    h, l, _ = matching_contributions(masses, shifted, 30.)
    if min(abs(h), abs(l)) < 1e-3:
        raise ArithmeticError("detuned control insensitive to field normalization")
    return 2 * math.log(4)


def adaptive_simpson(function, tol=2e-13):
    def panel(a, b, f_a, f_mid, f_b, coarse, target, depth):
        mid = (a+b)/2
        left_mid, right_mid = function((a+mid)/2), function((mid+b)/2)
        left = (mid-a)*(f_a + 4*left_mid + f_mid)/6
        right = (b-mid)*(f_mid + 4*right_mid + f_b)/6
        error = left+right-coarse
        if abs(error) <= 15*target:
            return left+right+error/15
        if depth == 0:
            raise ArithmeticError("integration did not converge")
        return (panel(a, mid, f_a, left_mid, f_mid, left, target/2, depth-1)
                + panel(mid, b, f_mid, right_mid, f_b, right, target/2, depth-1))
    f0, fm, f1 = function(0.), function(0.5), function(1.)
    return panel(0., 1., f0, fm, f1, (f0+4*fm+f1)/6, tol, 28)


def triangle_zero(b, a):
    if a == b:
        return 1 / (2*a)
    return (b*math.log(b/a)-b+a)/(b-a)**2


def zero_transfer_integral():
    samples = []
    for r in (0, .01, .25, 1., 4.):
        # x = k_E^2/(k_E^2+M^2); zero-transfer graph uses positive a
        # and M. r=0 is the limit of the mass-weighted sum taken together.
        def integrand(x):
            ratio = x/(x+r*(1-x)) if r else 1.
            return -3*ratio**2/(4-3*x)
        numeric = adaptive_simpson(integrand)
        analytic = (triangle_zero(1., r) - 4*triangle_zero(4., r)
                    if r else -math.log(4))
        close("Euclidean integration vs closed form", numeric, analytic)
        if numeric >= 0.:
            raise ArithmeticError("electron moment sign")
        samples.append({"a_over_M2": r, "M_times_S_e_zero": round(numeric, 13)})
    return samples


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {
        "status": "PASS",
        "yukawa_LLHH_graph_assignment_counts": yukawa_contractions(),
        "charged_four_lepton_contractions_by_quartic_self_pairs":
            charged_scalar_contractions(),
        "double_triangle_prefactor_relative_to_positive_simplex_box": -4,
        "dimension_five_quartic_shift_in_units_lambda_y2_over_16pi2_M":
            round(check_texture(), 13),
        "finite_scalar_mass_electron_triangle": zero_transfer_integral(),
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
