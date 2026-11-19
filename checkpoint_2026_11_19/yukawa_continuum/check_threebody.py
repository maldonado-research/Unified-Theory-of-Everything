#!/usr/bin/env python3
"""Compare invariant-mass sewing and direct three-body phase space, in m=1 units."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss


def unit_rule(order):
    x, w = leggauss(order)
    return (x + 1) / 2, w / 2


def f_mu(q):
    return 15 / ((1 + q) * (16 + q))


def f_e(q):
    return -30 * q / ((100 + q) * (400 + q))


def log_kernel(b, h, v):
    if v == 0:
        return 1 / (b + h)
    return math.log1p(2 * v / (b + h - v)) / (2 * v)


def spectral_cut(s, b, order=40):
    """Reduced cuts at c=0 and J=0 for a unit-strength scalar continuum."""
    root_s, r = math.sqrt(s), math.sqrt(b)
    threshold = (1 + r) ** 2
    if s <= threshold:
        return np.zeros(2)
    span = (s - threshold) * (root_s - 1 + r) / (root_s + 1 + r)
    ts, weights = unit_rule(order)
    xs, x_weights = unit_rule(order)  # integrand even
    values = np.zeros(2)
    for t, weight in zip(ts, weights):
        gap = span * t**2
        offset = span * (1 - t**2)
        u = b + offset
        h = (s - 1 - u) / 2
        v = math.sqrt(gap * (4 * root_s + gap)) / 2
        # At B=a the numerator and the squared pole denominator cancel exactly.
        pole_ratio = 1.0 if b == 1 else (offset / (b - 1 + offset)) ** 2
        rho = pole_ratio / (16 * math.pi**2 * u)
        mu_tree = f_mu(h - v * xs) + f_mu(h + v * xs)
        transverse_sq = v**2 * (1 - xs**2)
        electron_azimuth_mean = 2 * (10 / np.sqrt((100 + h)**2 - transverse_sq)
                                     - 40 / np.sqrt((400 + h)**2 - transverse_sq))
        right_angle = np.dot(x_weights, mu_tree * electron_azimuth_mean)
        mu_bar = log_kernel(1, h, v) - log_kernel(16, h, v)
        e_bar = 10 * log_kernel(100, h, v) - 40 * log_kernel(400, h, v)
        values += weight * (2 * span * t) * rho * (2 * v) / (16 * math.pi) * np.array(
            [right_angle, 4 * mu_bar * e_bar])
    return values


def recoil_cut(s, b, order=20):
    """Integrate the pole scalar's 3-momentum with both tree transfers in CM coordinates."""
    mass = math.sqrt(b)
    if math.sqrt(s) <= 1 + mass:
        return np.zeros(2)
    pmax = math.sqrt((s - (1 + mass)**2) * (s - (1 - mass)**2)) / (2 * math.sqrt(s))
    ps, p_weights = unit_rule(order)
    xs, x_weights = leggauss(order)
    phis = (np.arange(2 * order) + 0.5) * math.pi / order
    # First is the physical right angle. The other axes project the amplitude onto J=0.
    axes = np.r_[0.0, xs]
    scalar_cos_e = (np.sqrt(1 - axes[:, None, None]**2)
                    * np.sqrt(1 - xs[None, :, None]**2) * np.cos(phis)
                    + axes[:, None, None] * xs[None, :, None])
    accum = np.zeros(len(axes))
    for p_unit, weight in zip(ps, p_weights):
        p = p_unit * pmax
        e_pole = math.sqrt(1 + p**2)
        e_pair = math.sqrt(s) - e_pole
        u = e_pair**2 - p**2

        def transfers(cosine):
            return (1 - math.sqrt(s) * (e_pole - p * cosine),
                    u - math.sqrt(s) * (e_pair + p * cosine))

        def tree(cosine, masses_and_residues):
            return sum(w / (z - d) for z in transfers(cosine) for d, w in masses_and_residues)

        muon = tree(xs, [(1, 1), (16, -1)])
        electron = tree(scalar_cos_e, [(100, 10), (400, -40)])
        angular = np.sum(electron.mean(axis=2) * muon * x_weights / 2, axis=1)
        # d^3k/[(2π)^3 2E], angle-integrated decay, spin trace, propagator and 1/2 unitarity.
        recoil = p**2 / (4 * math.pi**2 * e_pole)
        decay_over_prop_sq = (1 if b == 1 else ((u - b) / (u - 1))**2) / (8 * math.pi * u)
        accum += weight * pmax * recoil * decay_over_prop_sq * s / 2 * angular
    return np.array([accum[0], np.dot(x_weights / 2, accum[1:])])


def threshold_coefficient(b):
    r = math.sqrt(b)
    trees = f_mu(r) * f_e(r)
    if b == 1:
        return trees / (96 * math.pi**3)
    kappa = r / (1 + r)
    return (trees * math.sqrt(r) * kappa**3 * (16 / 105)
            / (32 * math.pi**3 * b * (b - 1)**2))


def check_close(observed, expected, rtol, label):
    if not np.allclose(observed, expected, rtol=rtol, atol=0):
        raise ArithmeticError(f"{label}: {observed} vs {expected}")


def tree_routings(incoming_flavor):
    """Enum single external scalar and N1 slots for three Yukawa vertices."""
    slots = {0, 1, 2}  # incoming lepton 1, incoming lepton 2, outgoing muon
    flavors = [incoming_flavor, incoming_flavor, "mu"]
    counts = {"decay": 0, "exchange": 0}
    for scalar_vertex in (0, 1):  # charge fixes the scalar to an incoming vertex
        for n_vertex in sorted(slots):
            if n_vertex == scalar_vertex or flavors[n_vertex] != "mu":
                continue  # disconnected or forbidden by the N1 coupling
            left, right = sorted(slots - {n_vertex})
            if flavors[left] == flavors[right]:
                counts["decay" if n_vertex == 2 else "exchange"] += 1
    return counts


def rounded(value):
    return float(f"{value:.10g}")


def verify():
    counts = {flavor: tree_routings(flavor) for flavor in ("mu", "e")}
    if counts != {"mu": {"decay": 2, "exchange": 2}, "e": {"decay": 2, "exchange": 0}}:
        raise ArithmeticError("three-body tree routing count")
    rows = []
    for s in (6, 8, 10):
        cuts = spectral_cut(s, 1)
        direct = recoil_cut(s, 1)
        check_close(direct, cuts, 8e-11, "spectral and three-body sewing")
        check_close(spectral_cut(s, 1, 60), cuts, 5e-12, "spectral refinement")
        if not np.all(cuts < 0) or cuts[0] == cuts[1]:
            raise ArithmeticError("cut signs or angular dependence")
        for b in (16, 100, 400):
            check_close(spectral_cut(s, b), np.zeros(2), 0, "closed channel")
        rows.append(dict(s_over_m2=s, cut_c0=rounded(cuts[0]), cut_J0=rounded(cuts[1]),
                         direct_phase_space_c0=rounded(direct[0]),
                         direct_phase_space_J0=rounded(direct[1])))
    thresholds = []
    for label, b, weight in (("N1", 1, 1), ("N2", 16, 0.25), ("E1", 100, 1), ("E2", 400, 2)):
        threshold = (1 + math.sqrt(b))**2
        power = 1.5 if b == 1 else 3.5
        offset = 1e-5
        leading = threshold_coefficient(b)
        check_close(spectral_cut(threshold + offset, b) / offset**power,
                    leading, 2e-5, "threshold coefficient")
        thresholds.append(dict(column=label, B_over_m2=b, squared_coupling_weight=weight,
                               s_threshold_over_m2=int(threshold), onset_power=power,
                               coefficient_per_column_unit_strength=rounded(leading)))
    # A nondegenerate decay threshold tests the uncancelled self-energy numerator.
    check_close(recoil_cut(30, 16), spectral_cut(30, 16), 1e-10, "above N2 threshold")
    return dict(scope="three-particle cut of two scalar-line self-energy insertions",
                phase="B34_conjugate_times_B12 = s", mass_unit="m=1", M_over_m=10,
                coefficient="amplitude divided by gamma*abs(g)^2 at the three benchmark energies",
                angular_normalization="J0 is the external-cosine mean",
                charged_phi_mu_N1_tree_counts=counts, columns=thresholds, rows=rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = json.dumps(verify(), indent=2, allow_nan=False) + "\n"
    if args.output is None:
        print(payload, end="")
    else:
        args.output.write_text(payload, encoding="utf-8")


if __name__ == "__main__":
    main()
