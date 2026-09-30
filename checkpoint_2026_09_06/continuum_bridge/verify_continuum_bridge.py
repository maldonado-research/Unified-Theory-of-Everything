"""Independent standard-library checks of a conditional spectral toy model.

No observational data, field equations, or claims of microscopic equivalence.
The numerical integral and elementary closed form are evaluated independently.
"""

from __future__ import annotations

import csv
from fractions import Fraction
import json
import math
from pathlib import Path
import random


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def integrate(f, a: float, b: float, tol: float = 2e-13) -> float:
    """Adaptive Simpson quadrature, with failure on exhausted recursion."""
    fa, fb, fm = f(a), f(b), f((a + b) / 2)
    initial = (b - a) * (fa + 4 * fm + fb) / 6

    def recurse(left, right, fl, fr, fc, old, budget, depth):
        center = (left + right) / 2
        flc = f((left + center) / 2)
        frc = f((center + right) / 2)
        lower = (center - left) * (fl + 4 * flc + fc) / 6
        upper = (right - center) * (fc + 4 * frc + fr) / 6
        delta = lower + upper - old
        if abs(delta) <= 15 * budget:
            return lower + upper + delta / 15
        require(depth > 0, "Quadrature recursion exhausted")
        return recurse(left, center, fl, fc, flc, lower, budget / 2, depth - 1) + recurse(
            center, right, fc, fr, frc, upper, budget / 2, depth - 1
        )

    return recurse(a, b, fa, fb, fm, initial, tol, 28)


def toy_closed(s: float) -> float:
    return 2 / 3 - 2 * s + 2 * s ** 1.5 * math.atan(1 / math.sqrt(s))


def toy_integral(s: float) -> float:
    # lambda=x^2 removes the endpoint fractional power.
    return integrate(lambda x: 2 * x**4 / (x * x + s), 0, 1)


def pole(s: float) -> float:
    return (2 / 9) / (1 / 3 + s)


def regulated_continuum_imag(x: float, eta: float) -> float:
    # lambda=t^2 removes the endpoint fractional power. Explicit landmarks
    # around x ensure adaptive quadrature samples the Lorentzian peak.
    points = sorted({0.0, 1.0} | {
        math.sqrt(max(0.0, min(1.0, x + width * eta)))
        for width in [-1000, -100, -10, -1, 0, 1, 10, 100, 1000]
    })
    f = lambda t: 2 * t**4 * eta / ((t * t - x) ** 2 + eta * eta)
    return math.fsum(
        integrate(f, left, right, tol=2e-11 / (len(points) - 1))
        for left, right in zip(points[:-1], points[1:])
    )


def main() -> None:
    out = Path(__file__).resolve().parent
    mu0, mu1 = Fraction(2, 3), Fraction(2)
    pole_location = mu0 / mu1
    pole_weight = mu0 * mu0 / mu1
    require(pole_location == Fraction(1, 3), "Pole location mismatch")
    require(pole_weight == Fraction(2, 9), "Pole weight mismatch")
    require(pole_weight / pole_location == mu0, "Static moment mismatch")
    require(pole_weight / pole_location**2 == mu1, "First derivative mismatch")

    euclidean = []
    for s in [1e-4, 1e-3, 1e-2, 0.1, 0.3, 1.0, 4.0]:
        numerical, closed, fitted = toy_integral(s), toy_closed(s), pole(s)
        error = abs(numerical - closed)
        require(error < 4e-12, f"Closed form versus integral failed at {s}")
        require(numerical > fitted, f"Strict Jensen inequality failed at {s}")
        euclidean.append(
            dict(s=s, integral=numerical, closed_form=closed, pole=fitted, abs_error=error)
        )

    # These two moments are independently integrated in x=sqrt(lambda).
    require(abs(integrate(lambda x: 2 * x * x, 0, 1) - float(mu0)) < 1e-13, "mu0")
    require(abs(integrate(lambda x: 2.0, 0, 1) - float(mu1)) < 1e-13, "mu1")

    asymptotics = []
    for s in [1e-2, 1e-4, 1e-6, 1e-8]:
        # Stable exact remainders avoid subtracting nearly equal doubles.
        continuum_remainder = 2 * s**1.5 * math.atan(1 / math.sqrt(s))
        pole_remainder = 6 * s * s / (1 + 3 * s)
        susceptibility_deficit = 2 * math.sqrt(s) * math.atan(1 / math.sqrt(s))
        asymptotics.append(
            dict(
                s=s,
                continuum_remainder_over_s_1p5=continuum_remainder / s**1.5,
                pole_remainder_over_s_1p5=pole_remainder / s**1.5,
                susceptibility_deficit_over_sqrt_s=susceptibility_deficit / math.sqrt(s),
            )
        )
    require(
        abs(asymptotics[-1]["continuum_remainder_over_s_1p5"] - math.pi) < 2.01e-4,
        "Nonanalytic coefficient did not converge to pi",
    )
    require(asymptotics[-1]["pole_remainder_over_s_1p5"] < 0.000601, "Pole is not analytic")

    integer_log = []
    for s in [1e-3, 1e-5, 1e-7]:
        # Density=lambda^2 on [0,1]: beta=2, the v294 a=k borderline.
        remainder = s * s * math.log((1 + s) / s)
        value = 0.5 - s + remainder
        numerical = integrate(lambda lam: lam * lam / (lam + s), 0, 1)
        require(abs(value - numerical) < 2e-11, f"Integer logarithmic fixture failed at {s}")
        integer_log.append(dict(s=s, ratio=remainder / (s * s * math.log(1 / s))))
    require(abs(integer_log[-1]["ratio"] - 1) < 1e-7, "Integer log coefficient")

    # Independently exercise the general two-moment Jensen bound on 1,000 mixtures.
    rng = random.Random(20260906)
    jensen_min_gap = float("inf")
    for _ in range(1000):
        locations = [10 ** rng.uniform(-3, 3) for _ in range(rng.randrange(2, 10))]
        weights = [10 ** rng.uniform(-3, 3) for _ in locations]
        s = 10 ** rng.uniform(-4, 4)
        m0 = math.fsum(w / loc for w, loc in zip(weights, locations))
        m1 = math.fsum(w / loc**2 for w, loc in zip(weights, locations))
        exact = math.fsum(w / (loc + s) for w, loc in zip(weights, locations))
        lower = m0 * m0 / (m0 + s * m1)
        normalized_gap = (exact - lower) / max(exact, lower, 1e-300)
        require(normalized_gap >= -2e-13, "Positive-measure Jensen bound failed")
        jensen_min_gap = min(jensen_min_gap, normalized_gap)

    spectral = []
    for x in [0.1, 0.2, 0.5, 0.8]:
        for eta in [1e-2, 1e-3, 1e-4]:
            continuum = regulated_continuum_imag(x, eta)
            target = math.pi * x**1.5
            fitted = float(pole_weight) * eta / ((float(pole_location) - x) ** 2 + eta**2)
            spectral.append(
                dict(
                    x=x,
                    eta=eta,
                    continuum_imag=continuum,
                    continuum_limit=target,
                    continuum_abs_error=abs(continuum - target),
                    pole_imag=fitted,
                    pole_limit_away_from_atom=0.0,
                )
            )
        local = [row for row in spectral if row["x"] == x]
        require(local[-1]["continuum_abs_error"] < local[0]["continuum_abs_error"], "Cut limit")
        require(local[-1]["continuum_abs_error"] < 0.001, "Cut density convergence")
        require(local[-1]["pole_imag"] < 0.011 * local[0]["pole_imag"], "Pole off-atom limit")

    for name, rows in [("euclidean_comparison.csv", euclidean), ("spectral_comparison.csv", spectral)]:
        with (out / name).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    result = {
        "claim_scope": "Conditional dimensionless spectral toy; no observational or microscopic theory claim.",
        "status": "PASS",
        "exact_moments": {"mu0": str(mu0), "mu1": str(mu1)},
        "exact_matched_pole": {"lambda0": str(pole_location), "weight": str(pole_weight)},
        "third_inverse_moment_cutoff": "mu2(epsilon)=2*(epsilon**(-1/2)-1), divergent as epsilon->0",
        "independent_integral_comparisons": euclidean,
        "maximum_closed_form_integral_error": max(row["abs_error"] for row in euclidean),
        "small_s_asymptotics": asymptotics,
        "integer_beta_2_log_fixture": integer_log,
        "general_positive_measure_jensen_cases": 1000,
        "minimum_normalized_jensen_gap": jensen_min_gap,
        "finite_regulator_spectral_cases": 12,
        "maximum_final_spectral_abs_error": max(
            row["continuum_abs_error"] for row in spectral if row["eta"] == 1e-4
        ),
        "limits": [
            "No physical source coupling or particle-production rate has been derived.",
            "Positivity is an explicit scalar spectral assumption, not a Majorana-kernel property.",
            "A finite regulator broadens an atom; only the eta->0 distributional limit distinguishes it exactly.",
            "This script does not rederive the v294 operator or verify its claimed spectral normalization.",
        ],
    }
    (out / "verification_results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in [
        "status", "maximum_closed_form_integral_error", "general_positive_measure_jensen_cases",
        "minimum_normalized_jensen_gap", "maximum_final_spectral_abs_error"
    ]}, indent=2))


if __name__ == "__main__":
    main()
