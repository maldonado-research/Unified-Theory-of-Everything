"""Check the profiled reciprocal-rate bound with independent angular/energy quadratures."""

import json
from pathlib import Path

from mpmath import mp


def close(label, value, expected, relative="1e-26"):
    tolerance = mp.mpf(relative) * max(abs(expected), mp.mpf("1e-18"))
    if abs(value - expected) > tolerance:
        raise ArithmeticError(f"{label}: {value} versus {expected}")


def n(xi):
    return xi / 2 - 1


def log_gap(xi):
    beta = mp.sqrt(1 - 4 / xi)
    h = xi * beta / 2
    return 2 / xi * (mp.atanh(h / (1 + n(xi))) - mp.atanh(h / (16 + n(xi))))


def rate_kernel(xi, gap):
    return xi * gap**2 / (2048 * mp.pi**3)


def scan_moments(exposure):
    return [exposure(lambda xi: rate_kernel(xi, log_gap(xi)) * n(xi)**i) for i in range(3)]


def profile(moments):
    s0, s1, s2 = moments
    return s1 / s0, s2 - s1**2 / s0


def verify():
    mp.dps = 58
    points = [mp.mpf(6), mp.mpf(10)]
    bins = [(mp.mpf(5), mp.mpf(7)), (mp.mpf(9), mp.mpf(11))]

    def point_exposure(f):
        return sum(f(x) for x in points) / 2

    def bin_exposure(f):
        return sum(mp.quad(f, b) for b in bins) / 4

    # Independent implementation. Integrate the product of the transfer
    # denominators over the scalar-pair orientation; use Gaussian energy
    # nodes and directly sum squares around the weighted centroid.
    nodes, weights = mp.gauss_quadrature(52, "legendre")
    angular_rule = list(zip(nodes, weights))

    def transfer_gap(xi):
        beta = mp.sqrt(1 - 4 / xi)
        h = xi * beta / 2
        return beta / 2 * sum(
            weight * 15 / ((xi / 2 - h*x) * (15 + xi / 2 - h*x))
            for x, weight in angular_rule
        )

    bin_rule = [
        ((left + right + (right - left)*x) / 2, (right - left)*w / 8)
        for left, right in bins for x, w in angular_rule
    ]
    results = {}
    cases = [
        ("two_points", point_exposure, [(x, mp.mpf("0.5")) for x in points]),
        ("two_bins", bin_exposure, bin_rule),
    ]
    for name, exposure, rule in cases:
        close("total exposure", sum(p for x, p in rule), 1)
        moments = scan_moments(exposure)
        centroid, floor = profile(moments)
        weighted = []
        for xi, p in rule:
            j = transfer_gap(xi)
            close("angular transfer integral versus gap logarithm", j, log_gap(xi))
            if j <= 0:
                raise ArithmeticError("nonpositive gap above threshold")
            weighted.append((n(xi), p * rate_kernel(xi, j)))
        s0 = sum(v for _, v in weighted)
        independent_centroid = sum(t*v for t, v in weighted) / s0
        independent_floor = sum(v*(t - independent_centroid)**2 for t, v in weighted)
        close("weighted centroid", independent_centroid, centroid)
        close("centered sum floor", independent_floor, floor)
        if floor <= 0:
            raise ArithmeticError("nonpositive profiled floor")
        # Evaluate the unprofiled functional with genuinely complex phases
        # directly, both at and away from the minimizing complex coefficient.
        g7 = mp.mpc("0.7", "-0.4")
        for g5 in (centroid*g7, mp.mpc("-0.2", "1.4")):
            direct = sum(v * abs(t*g7 - g5)**2 for t, v in weighted)
            projected = abs(g7)**2 * floor + s0 * abs(g5 - centroid*g7)**2
            close("complex square completion", direct, projected)
        if name == "two_points":
            (t1, v1), (t2, v2) = weighted
            close("two-point pairwise variance", v1*v2/(v1 + v2)*(t1 - t2)**2, floor)
            tuned_g5 = n(points[0])*g7
            close("tuned single energy", v1*abs(t1*g7 - tuned_g5)**2, 0)
            close("nonzero cut at second energy after first tuned", abs(t2*g7 - t1*g7), 2*abs(g7))
        results[name] = {
            "moments_S0_S1_S2": [mp.nstr(s, 28) for s in moments],
            "optimal_g5_over_g7": mp.nstr(centroid, 28),
            "a_times_profiled_reciprocal_mean_floor_over_abs_g7_squared": mp.nstr(floor, 28),
            "difference_between_floor_quadratures": mp.nstr(abs(floor - independent_floor), 6),
        }
    xi = mp.mpf(8)
    close("zero-Weinberg limit versus checkpoint at xi=8", rate_kernel(xi, log_gap(xi))*n(xi)**2,
          mp.mpf("0.0000379084757248466264"), relative="3e-18")
    return {
        "status": "PASS",
        "point_exposure": {"xi": [6, 10], "probabilities": [0.5, 0.5]},
        "bin_exposure": {"xi_intervals": [[5, 7], [9, 11]], "density_per_xi_on_union": 0.25},
        "inputs": {"b_over_a": [1, 16], "w_over_w0": [1, -1]},
        "numerics": "58 decimal digits, log gap with tanh-sinh energy integration versus 52-point Gauss–Legendre transfer and energy rules",
        "results": results,
    }


if __name__ == "__main__":
    data = json.dumps(verify(), indent=2, ensure_ascii=False) + "\n"
    Path(__file__).with_name("RATE_SCAN_RESULTS.json").write_text(data, encoding="utf-8")
    print(data, end="")
