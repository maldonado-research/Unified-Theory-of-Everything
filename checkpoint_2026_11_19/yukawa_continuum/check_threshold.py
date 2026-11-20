"""Rational certificate for the exclusive three-particle threshold bound."""

import argparse
import json
from fractions import Fraction as F
from pathlib import Path


def require(condition, description):
    if not condition:
        raise ArithmeticError(description)


def verify():
    # m=1. All inequalities used for the spin-vector norms are rational
    # after squaring the (positive) bounds, including the CM energy cap.
    delta = F(1, 100)
    energy = F(501, 500)
    momentum = F(1, 20)
    require(energy**2 > 1+delta/4, "beam and final massive energy bound")
    require(momentum**2 == delta/4, "maximum massive recoil momentum")
    pair_energy = energy+momentum
    qmin, qmax = F(89, 100), F(112, 100)
    require(2*(1-momentum)-1 > qmin, "spacelike transfer lower bound")
    require(2*energy*pair_energy-1 < qmax, "spacelike transfer upper bound")
    matrix_norm = (energy-1)+energy+momentum
    scalar_gap = 2*(1-momentum)
    require(matrix_norm < F(106, 100), "Weyl matrix norm")
    require(scalar_gap == F(19, 10), "scalar propagator gap")
    anomalous_lower = 2*15 / ((1+qmax)*(16+qmax))
    normal_upper = 1/(1+qmin)+1/(4*(16+qmin))
    direct_sq = 4*anomalous_lower**2/(2*pair_energy)
    exchange_sq = 2*energy*pair_energy*(4*energy)*matrix_norm**2*normal_upper**2/scalar_gap**2
    require(direct_sq > F(113, 100)**2, "decay amplitude lower bound")
    require(exchange_sq < F(89, 100)**2, "one exchange amplitude upper bound")
    # Two exchanges, and omega <= Delta/4.
    coefficient_sq = (F(2, 4))**2*exchange_sq/direct_sq
    require(coefficient_sq < F(2, 5)**2, "relative cut coefficient bound")
    return {
        "status": "PASS",
        "domain": "0 < Delta = s/m^2 - 4 <= 1/100",
        "final_state": "charged scalar, muon, N1",
        "outward_transfer_bounds": [float(qmin), float(qmax)],
        "scalar_propagator_gap_lower_bound": str(scalar_gap),
        "weyl_matrix_norm_bound": str(matrix_norm),
        "direct_norm_squared_times_omega_lower_bound": str(direct_sq),
        "one_exchange_norm_squared_divided_by_omega_upper_bound": str(exchange_sq),
        "relative_cut_over_Delta_bound_squared": str(coefficient_sq),
        "outward_relative_cut_bound": "(2/5) Delta",
        "proof_arithmetic": "exact rational inequalities between squared positive quantities",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    options = parser.parse_args()
    report = json.dumps(verify(), indent=2) + "\n"
    if options.output:
        options.output.write_text(report, encoding="utf-8")
    print(report, end="")
