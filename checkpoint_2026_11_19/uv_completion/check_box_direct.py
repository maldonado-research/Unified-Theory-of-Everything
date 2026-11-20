"""Evaluate the physical upper-rim boxes directly in Feynman parameters."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre


def unit_rule(order):
    nodes, weights = roots_legendre(order)
    return (1 + nodes) / 2, weights / 2


def integrated_y(x, s, t, b, c, nodes, weights):
    """Integrate over x3=y after integrating x1 out of the simplex."""
    length = 1 - x
    slope = b - 1 + s*x
    intercept = 1 - s*x*length
    pole = -intercept/slope
    y = length[:, None]*nodes
    a0 = x + c*length
    a1 = b - c - t*length
    a_y = a0[:, None] + y*(a1[:, None] + t*y)
    values = np.empty(len(x), dtype=complex)

    # Plain quadrature is smooth when the pole lies well outside the interval.
    far = pole < -length/4
    yf = y[far]
    integrand = (length[far, None]-yf)/(a_y[far]*(intercept[far, None]+slope[far, None]*yf))
    values[far] = length[far]*np.sum(weights*integrand, axis=1)

    # Subtract H(pole) for physical poles and nearby poles just outside.
    near = ~far
    p = pole[near]
    upper = length[near]
    a_p = a0[near] + p*(a1[near] + t*p)
    h_p = (upper-p)/a_p
    divided = (-a_p[:, None] - (upper-p)[:, None]*(a1[near, None]+t*(y[near]+p[:, None])))
    divided /= a_p[:, None]*a_y[near]
    analytic = np.log(np.abs((upper-p)/p)) + 1j*np.pi*((p > 0) & (p < upper))
    values[near] = (upper*np.sum(weights*divided, axis=1) + h_p*analytic)/slope[near]
    return values


def box(s, t, b, c, rule):
    nodes, weights = rule
    beta = math.sqrt(1-4/s)
    ends = (0., (1-beta)/2, (1+beta)/2, 1.)
    total = 0j
    for left, right in zip(ends, ends[1:]):
        # Both logarithmic endpoints occur at a boundary of one of these panels.
        x = left + (right-left)*np.sin(np.pi*nodes/2)**2
        jacobian = (right-left)*np.pi/2*np.sin(np.pi*nodes)
        total += np.sum(weights*jacobian*integrated_y(x, s, t, b, c, nodes, weights))
    return total


def amplitude(s, cosine, rule):
    transfers = (-s*(1-cosine)/2, -s*(1+cosine)/2)
    value = 0j
    for b, muon_weight in ((1, 1), (16, -1)):
        for c, electron_weight in ((100, 10), (400, -40)):
            value += muon_weight*electron_weight*sum(box(s, t, b, c, rule) for t in transfers)
    return s*value/(16*np.pi**2)


def observables(s, parameter_order, angle_order):
    rule = unit_rule(parameter_order)
    # The two exchanged boxes make the spinor-stripped amplitude even in cos(theta).
    cosines, weights = unit_rule(angle_order)
    amplitudes = np.array([amplitude(s, z, rule) for z in cosines])
    central = amplitude(s, 0, rule)
    mean = np.sum(weights*amplitudes)
    return dict(central_re=central.real, central_im=central.imag, mean_re=mean.real,
                mean_im=mean.imag, m2_sigma=np.sum(weights*abs(amplitudes)**2)/(32*np.pi*s),
                m2_sigma_from_imaginary_part=np.sum(weights*amplitudes.imag**2)/(32*np.pi*s))


def projected_cut(s):
    beta = math.sqrt(1-4/s)
    h = s/2-1
    v = s*beta/2

    def transfer_average(b):
        return math.atanh(v/(b+h))/v

    muon = transfer_average(1)-transfer_average(16)
    electron = 10*transfer_average(100)-40*transfer_average(400)
    return s*beta*muon*electron/(8*math.pi)


def verify():
    refinements = ((64, 6), (128, 10), (224, 16))
    rows = []
    for s in (6., 10.):
        series = [observables(s, param, angle) for param, angle in refinements]
        reference = series[-1]
        for value, tolerance in zip(series[:-1], (3e-6, 2e-8)):
            for name, estimate in value.items():
                if not math.isclose(estimate, reference[name], rel_tol=tolerance):
                    raise ArithmeticError(f"Quadrature refinement failed for s={s}, {name}")
        exact_cut = projected_cut(s)
        if not math.isclose(exact_cut, reference["mean_im"], rel_tol=2e-12):
            raise ArithmeticError(f"Independent cut normalization failed for s={s}")
        floor = exact_cut**2/(32*math.pi*s)
        if not reference["m2_sigma"] > reference["m2_sigma_from_imaginary_part"] > floor:
            raise ArithmeticError(f"Integrated rate or cut floor failed for s={s}")
        rows.append(dict(s_over_m2=int(s), refinement=[
            {name: float(f"{number:.9g}") for name, number in row.items()} for row in series
        ], m2_sigma_floor_J0=float(f"{floor:.9g}")))
    return dict(status="PASS", scalar_mass_unit=1, M_over_m=10, y=1, g=1,
                B34_conjugate_times_B12_over_s=1, central_cos_theta=0,
                orders_parameter_and_angle=refinements, rows=rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write the JSON receipt to this file")
    args = parser.parse_args()
    payload = json.dumps(verify(), indent=2, allow_nan=False) + "\n"
    if args.output is None:
        print(payload, end="")
    else:
        args.output.write_text(payload, encoding="utf-8")
