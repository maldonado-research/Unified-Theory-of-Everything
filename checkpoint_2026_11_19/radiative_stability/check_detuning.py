"""Sensitivity of the complete one-loop rates to a complex Yukawa detuning."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre


PREVIOUS = Path(__file__).resolve().parents[1] / "uv_completion"
RADIUS = 0.01


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def direct_box():
    spec = importlib.util.spec_from_file_location("direct_box_reference", PREVIOUS / "check_box_direct.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.box


def unit_rule(order):
    nodes, weights = roots_legendre(order)
    return (1 + nodes) / 2, weights / 2


def dispersion_box(s, t, b, c, rule):
    """Finite massive box on the upper rim, in units of the scalar mass."""
    nodes, weights = rule
    log_span = math.log(c / b)
    mass = b * np.exp(nodes * log_span)
    fraction = (c - mass) / (c - b)
    jacobian = mass * log_span / (c - b)
    beta = math.sqrt(1 - 4 / s)
    pole = beta * (-math.atanh(beta) + 0.5j * math.pi)
    leading = (1 + mass)**2
    ratio = 4 * (mass - t * fraction * (1 - fraction)) / leading - 1
    # Integral of (1+ratio*q^2)^(-1) from 0 to 1, including ratio=0.
    smooth = np.ones_like(ratio)
    for mask, primitive in ((ratio < 0, np.arctanh), (ratio > 0, np.arctan)):
        root = np.sqrt(np.abs(ratio[mask]))
        smooth[mask] = primitive(root) / root
    return (8 / s) * np.sum(weights * jacobian * (pole + smooth)
                            / (leading * (1 + ratio * beta**2)))


def angular_modes(s, box, mass_order, angular_order):
    mass_rule = unit_rule(mass_order)
    angles, weights = unit_rule(angular_order)
    modes = []
    for z in angles:
        transfers = (-s * (1 - z)/2, -s * (1 + z)/2)
        electron_terms = []
        for c, e_weight in ((100, 10), (400, -40)):
            integral = sum(mu_weight * box(s, t, b, c, mass_rule)
                           for t in transfers for b, mu_weight in ((1, 1), (16, -1)))
            electron_terms.append(s * e_weight * integral / (16 * math.pi**2))
        modes.append((sum(electron_terms), electron_terms[1]))
    return np.array(modes), weights / (32 * math.pi * s)


def gram(modes, weights):
    return (modes.conj().T * weights) @ modes


def at(gram_matrix, delta):
    v = np.array([1, delta])
    return float(np.vdot(v, gram_matrix @ v).real)


def ratio_bounds(matrices, radius):
    low, high = matrices
    real_diag = lambda matrix: matrix[0, 0].real + radius**2 * matrix[1, 1].real
    L, H = real_diag(low), real_diag(high)
    b, e = low[0, 1], high[0, 1]
    require(math.sqrt(low[0, 0].real) > radius*math.sqrt(low[1, 1].real),
            "cannot certify strictly positive denominator on disk")
    # On its circle, a particular ratio is possible iff the constant part is
    # no larger than the magnitude of the variable-phase interference.
    coefficients = (L**2 - 4*radius**2*abs(b)**2,
                    -2*L*H + 8*radius**2*(b.conjugate()*e).real,
                    H**2 - 4*radius**2*abs(e)**2)
    roots = np.roots(np.array(coefficients) / max(map(abs, coefficients)))
    require(np.isreal(roots).all(), "non-real ratio endpoints")
    extrema = []
    for r in sorted(roots.real):
        cross = e - r*b
        curvature = high[1, 1].real - r*low[1, 1].real
        require(curvature < 0, "disk argument requires concave rate difference")
        d = -np.sign(H-r*L)*radius*cross.conjugate()/abs(cross)
        require(abs(at(high, d) / at(low, d) - r) < 1e-12, "circle endpoint witness")
        extrema.append((r, d, abs(cross / curvature)))
    # Minimum of a concave function occurs on boundary. At the upper ratio
    # its unconstrained maximum also lies outside the disk.
    require(extrema[1][2] > radius, "upper extremum might be in disk interior")
    return extrema


def rounded(value):
    return float(f"{value:.10g}")


def run():
    independent_box = direct_box()
    external_rows = json.loads((PREVIOUS / "QCDLOOP.json").read_text())["rows"]
    results, matrices = {}, []
    external = {row["s_over_m2"]: row for row in external_rows}
    for s in (6, 8, 10):
        modes, weights = angular_modes(s, dispersion_box, 36, 16)
        matrix = gram(modes, weights)
        coarse = gram(*angular_modes(s, dispersion_box, 24, 10))
        direct = gram(*angular_modes(s, independent_box, 160, 12))
        require(np.allclose(matrix, coarse, rtol=1e-12, atol=0), "dispersion grid refinement")
        require(np.allclose(matrix, direct, rtol=2e-8, atol=0), "direct simplex cross-check")
        require(np.linalg.eigvalsh(matrix)[0] > 0, "Gram form is not positive definite")
        if s in external:
            require(abs(at(matrix, 0) / external[s]["angle_order_20"]["full"] - 1) < 2e-12,
                    "tuned QCDLoop rate")
        examples = {}
        for name, d in (("minus_one_percent", -RADIUS), ("tuned", 0), ("plus_one_percent", RADIUS),
                        ("plus_one_percent_i", 1j*RADIUS), ("minus_one_percent_i", -1j*RADIUS)):
            independent_integral = np.sum(weights * abs(modes[:, 0] + d*modes[:, 1])**2)
            require(abs(at(matrix, d)/independent_integral-1) < 1e-13, "matrix phase convention")
            examples[name] = rounded(at(matrix, d))
        asymmetry = (at(matrix, 1j*RADIUS)-at(matrix, -1j*RADIUS)) / (
                    at(matrix, 1j*RADIUS)+at(matrix, -1j*RADIUS))
        require(abs(asymmetry + 2*RADIUS*matrix[0, 1].imag /
                    (matrix[0, 0].real + RADIUS**2*matrix[1, 1].real)) < 1e-12,
                "reciprocal sign or phase")
        results[str(s)] = {"a": rounded(matrix[0, 0].real), "b_real": rounded(matrix[0, 1].real),
                           "b_imag": rounded(matrix[0, 1].imag), "c": rounded(matrix[1, 1].real),
                           "scaled_rates": examples,
                           "dlog_rate_dRe_delta_at_zero": float(f"{2*matrix[0, 1].real/matrix[0, 0].real:.15g}"),
                           "exclusive_reciprocal_asymmetry_at_plus_one_percent_i": rounded(asymmetry)}
        matrices.append(matrix)
    pair = (matrices[0], matrices[2])
    ratios = {name: rounded(at(pair[1], value)/at(pair[0], value))
              for name, value in (("minus_one_percent", -RADIUS), ("tuned", 0), ("plus_one_percent", RADIUS))}
    bounds = ratio_bounds(pair, RADIUS)
    return {"status": "PASS", "rate_units": "m^2 sigma / |y_e y_mu|^4", "M_over_m": 10,
            "detuning": "Y_e2^2 = -2 y_e^2 (1+delta), masses and base couplings fixed",
            "full_direct_rate": "a + 2 Re(delta*b) + |delta|^2*c",
            "s_over_m2": results, "rate_ratio_10_over_6": ratios,
            "one_percent_complex_disk_ratio_extrema": [
                {"ratio": rounded(r), "delta_real": rounded(z.real), "delta_imag": rounded(z.imag),
                 "unconstrained_vertex_radius": rounded(vertex)} for r, z, vertex in bounds]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    body = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(body)
    else:
        print(body, end="")
