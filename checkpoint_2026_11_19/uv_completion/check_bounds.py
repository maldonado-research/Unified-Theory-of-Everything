"""Fixed-transfer dispersion for the finite box and uniform error on its cut."""
import argparse
import json
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import mpmath as mp

mp.mp.dps = 42


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


@lru_cache(None)
def unit_rule(order):
    nodes, weights = mp.gauss_quadrature(order, "legendre")
    return tuple(((x + 1) / 2, w / 2) for x, w in zip(nodes, weights))


def angular_primitive(ratio):
    """Integral of (1+ratio*v^2)^(-1) on [0,1]."""
    if ratio == 0:
        return mp.mpf(1)
    if ratio < 0:
        root = mp.sqrt(-ratio)
        return mp.atanh(root) / root
    root = mp.sqrt(ratio)
    return mp.atan(root) / root


def reduced_amplitude(energy, cosine, order):
    """m=1, upper rim, spinor B34*B12=s, unit base Yukawas."""
    beta = mp.sqrt(1 - 4 / energy)
    pole = beta * (-mp.atanh(beta) + mp.j * mp.pi / 2)
    total = 0
    for transfer in (-energy * (1 - cosine) / 2, -energy * (1 + cosine) / 2):
        for b, mu_weight in ((1, 1), (16, -1)):
            for B, e_weight in ((100, 10), (400, -40)):
                terms = []
                # b_lambda spans quite unequal masses. Integrate uniformly in
                # log b_lambda to resolve the light endpoint at each order.
                log_range = mp.log(mp.mpf(B) / b)
                for u, weight in unit_rule(order):
                    mean_mass2 = b * mp.exp(u * log_range)
                    lam = (B - mean_mass2) / (B - b)
                    jacobian = mean_mass2 * log_range / (B - b)
                    L = (1 + mean_mass2)**2
                    H = 4 * (mean_mass2 - transfer * lam * (1 - lam)) - L
                    terms.append(weight * jacobian * (pole + angular_primitive(H / L))
                                 / (L + H * beta**2))
                total += mu_weight * e_weight * mp.fsum(terms)
    return total / (2 * mp.pi**2)  # s/(16*pi^2) times 8/s


def observables(energy, order, angular_order):
    energy = mp.mpf(energy)
    # The two identical assignments make the reduced amplitude even in cos θ.
    samples = [(w, reduced_amplitude(energy, x, order)) for x, w in unit_rule(angular_order)]
    center = reduced_amplitude(energy, 0, order)
    mean_im = mp.fsum(w * amp.imag for w, amp in samples)
    phase_space = 1 / (32 * mp.pi * energy)
    return {
        "central_amplitude_real": center.real,
        "central_amplitude_imaginary": center.imag,
        "scaled_full_rate": phase_space * mp.fsum(w * abs(amp)**2 for w, amp in samples),
        "scaled_absorptive_rate": phase_space * mp.fsum(w * amp.imag**2 for w, amp in samples),
        "scaled_absorptive_J0_floor": phase_space * mean_im**2,
    }


def cut_j0(energy):
    # An independent average of the *tree* propagators, with both t/u routings.
    s = mp.mpf(energy)
    beta = mp.sqrt(1 - 4 / s)
    n, h = s / 2 - 1, s * beta / 2
    def average(b):
        return mp.atanh(h / (b + n)) / h
    mu = average(1) - average(16)
    e = 10 * (average(100) - 4 * average(400))
    return s * beta * mu * e / (8 * mp.pi)


def cut_error():
    # z/M² < (xi-2)/100 <= 9/100 for 4 < xi <= 11.
    x = Fraction(9, 100)
    inverse = (1 + x) * (1 + x/4)
    linear_error = 1 - inverse * (1 - 5*x/4)
    require(linear_error == (21*x**2 + 5*x**3)/16, "dimension-nine remainder")
    return {
        "range_s_over_mass_squared": "(4,11]",
        "max_minus_transfer_over_M_squared": str(x),
        "lower_exact_over_dimension_seven_cut": str(1/inverse),
        "upper_dimension_seven_error_over_exact_cut": str(inverse-1),
        "upper_dimension_seven_plus_nine_error_over_exact_cut": str(linear_error),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    results = {}
    for s in (6, 10):
        coarse, fine = observables(s, 24, 10), observables(s, 38, 16)
        for key in fine:
            require(abs(coarse[key]/fine[key]-1) < mp.mpf("3e-16"), f"quadrature at {s}: {key}")
        analytic_floor = cut_j0(s)**2 / (32 * mp.pi * s)
        require(abs(fine["scaled_absorptive_J0_floor"]/analytic_floor-1) < mp.mpf("1e-19"),
                "dispersion versus projected tree cut")
        require(fine["scaled_full_rate"] > fine["scaled_absorptive_rate"]
                > fine["scaled_absorptive_J0_floor"] > 0, "positive partial rates")
        results[str(s)] = {key: mp.nstr(value, 16) for key, value in fine.items()}
    direct = json.loads(Path(__file__).with_name("BOX_DIRECT.json").read_text())
    for row in direct["rows"]:
        result = results[str(row["s_over_m2"])]
        last = row["refinement"][-1]
        comparisons = (("central_amplitude_real", last["central_re"]),
                       ("central_amplitude_imaginary", last["central_im"]),
                       ("scaled_full_rate", last["m2_sigma"]),
                       ("scaled_absorptive_rate", last["m2_sigma_from_imaginary_part"]),
                       ("scaled_absorptive_J0_floor", row["m2_sigma_floor_J0"]))
        for key, reference in comparisons:
            require(abs(mp.mpf(result[key]) / reference - 1) < mp.mpf("6e-9"),
                    f"independent direct simplex calculation: {key}")
    require({row["s_over_m2"] for row in direct["rows"]} == {6, 10},
            "both independent energy benchmarks required")
    external = json.loads(Path(__file__).with_name("QCDLOOP.json").read_text())
    require([row["s_over_m2"] for row in external["rows"]] == [6, 10],
            "missing external energies")
    max_external_difference = mp.mpf(0)
    for row in external["rows"]:
        result = results[str(row["s_over_m2"])]
        rates = row["angle_order_20"]
        comparisons = (("central_amplitude_real", row["central_real"]),
                       ("central_amplitude_imaginary", row["central_imaginary"]),
                       ("scaled_full_rate", rates["full"]),
                       ("scaled_absorptive_rate", rates["absorptive"]),
                       ("scaled_absorptive_J0_floor", rates["J0_absorptive"]))
        for key, reference in comparisons:
            difference = abs(mp.mpf(result[key]) / mp.mpf(str(reference)) - 1)
            max_external_difference = max(max_external_difference, difference)
    require(max_external_difference < mp.mpf("1e-13"), "QCDLoop external reference")
    report = {"QCDLoop_max_relative_difference": mp.nstr(max_external_difference, 4),
              "model": "m=1, sterile masses (1,4;10,20), base Yukawas divided out",
              "method": "fixed spacelike transfer, unsubtracted spectral integral",
              "cross_section_units": "m^2 sigma / |y_e y_mu|^4",
              "upper_rim": results, "cut_error": cut_error()}
    data = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(data)
    else:
        print(data, end="")


if __name__ == "__main__":
    main()
