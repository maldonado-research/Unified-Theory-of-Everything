"""Exact rational guard for a nonzero weak mass-pulse transition bound.

The analytic ingredients (Fourier integral and unitary Duhamel estimate)
are proved in dynamic_eft_review.md. This checks the conservative numerical
inequality without floating arithmetic; it is not an ODE or QFT replay.
"""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import json


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


x = Q(9, 2)
order = 30
partial = sum((x**j / factorial(j) for j in range(order + 1)), Q(0))
first_tail = x**(order + 1) / factorial(order + 1)
ratio = x / (order + 2)
exp_upper = partial + first_tail / (1 - ratio)
require(ratio < 1, "Exponential tail requires a geometric majorant")
require(exp_upper < 91, "Need exp(9/2) < 91")
require(Q(22, 7) * Q(10, 7) < x, "Need pi sqrt(2) < 9/2")
require(Q(10, 7)**2 > 2, "Need sqrt(2) < 10/7")

# M=k=tau=1, a=1/100. b1=(pi/50)/sinh(pi sqrt(2)).
# pi>3 and sinh(pi sqrt(2)) < exp(9/2)/2 < 91/2.
born_lower = Q(3, 2275)
duhamel_remainder_upper = Q(1, 5000)
amplitude_lower = born_lower - duhamel_remainder_upper
probability_lower = amplitude_lower**2
require(amplitude_lower > 0, "Exact transition bound must be nonzero")
require(probability_lower > Q(1, 800000), "Need probability > 1.25e-6")

result = {
    "status": "PASS_EXACT_RATIONAL_NONZERO_TOY_TRANSITION_BOUND",
    "parameters": {"M": "1", "k": "1", "tau": "1", "a": "1/100"},
    "analytic_premises": [
        "Unitary two-level evolution with H=m(t)sigma_z+k sigma_x",
        "m(t)=M[1-a sech^2(t/tau)] and in/out mass M",
        "Exact sech-squared Fourier transform as derived in review",
        "Twice-iterated unitary Duhamel remainder <= L1(V)^2/2",
        "Classical bounds 3 < pi < 22/7"
    ],
    "born_amplitude_strict_lower": str(born_lower),
    "amplitude_remainder_upper": str(duhamel_remainder_upper),
    "true_amplitude_strict_lower": str(amplitude_lower),
    "true_transition_probability_strict_lower": str(probability_lower),
    "displayed_conservative_probability_lower": "0.00000125",
    "exp_4p5_rational_upper": str(exp_upper),
    "scope": "One mode of a toy prescribed background; not cosmological calibration, reheating, lepton asymmetry, or a new external novelty claim."
}
out = Path(__file__).with_name('DYNAMIC_EFT_EXACT_NONZERO_BOUND.json')
out.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: v for k, v in result.items() if k != 'exp_4p5_rational_upper'}, indent=2))
