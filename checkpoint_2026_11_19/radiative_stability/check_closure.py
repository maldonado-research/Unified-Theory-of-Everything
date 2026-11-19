"""Close the finite complex-detuning model using three energies and one inverse rate."""
import argparse
import importlib.util
import json
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp


ROOT = Path(__file__).resolve().parent
ENERGIES = (6, 8, 10)


def require(condition, reason):
    if not condition:
        raise ArithmeticError(reason)


def read(name):
    return json.loads((ROOT/name).read_text(), parse_float=Decimal)


def interval(x):
    """Independent one-part-per-million relative perturbations to a reported input."""
    dx = abs(x)/10**6
    return x-dx, x+dx


def subtract(a, b):
    return a[0]-b[1], a[1]-b[0]


def multiply(a, b):
    endpoints = [x*y for x in a for y in b]
    return min(endpoints), max(endpoints)


def determinant_interval(d, h):
    dx, dy = [subtract(interval(k), interval(d[0])) for k in d[1:]]
    hx, hy = [subtract(interval(k), interval(h[0])) for k in h[1:]]
    return subtract(multiply(dx, hy), multiply(dy, hx))


def finite_box_test(matrix, a, f):
    # Build the data with the earlier direct simplex integral. The closure
    # coefficients came from the separate dispersion calculation.
    spec = importlib.util.spec_from_file_location("detuning", ROOT/"check_detuning.py")
    reference = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reference)
    direct = reference.direct_box()
    samples = [reference.angular_modes(s, direct, 160, 12) for s in ENERGIES]
    rows = []
    for scale, delta in ((0.37, 0.01 + 0j), (2.4, -0.007+0.006j), (7.1, 0.018j)):
        plus, minus = [], []
        for aa, (modes, w) in zip(a, samples, strict=True):
            fold = lambda z: scale*np.sum(w*abs(modes[:, 0]+z*modes[:, 1])**2)/float(aa)
            plus.append(fold(delta))
            minus.append(fold(delta.conjugate()))
        plus, minus = np.array(plus), np.array(minus)
        T = (plus[0]-minus[0])/(2*f[0])
        reconstructed_means = plus - f*T
        require(np.allclose(reconstructed_means, (plus+minus)/2, rtol=1e-8, atol=0),
                "unmeasured inverse rate fails cross-check")
        U, V, W = np.linalg.solve(matrix, reconstructed_means)
        residual = W/U - (V/U)**2 - (T/U)**2
        require(U>0 and W>=0 and abs(U/scale-1)<1e-7, "normalization or positivity")
        require(abs(V/U-delta.real)<2e-8 and abs(T/U-delta.imag)<2e-8, "detuning recovery")
        require(abs(residual)<1e-9, "complex closure against independent boxes")
        rows.append({"normalization": scale, "delta_real": delta.real, "delta_imag": delta.imag,
                     "reconstructed_real": round(V/U, 8)+0.0, "reconstructed_imag": round(T/U, 8)+0.0,
                     "closure_absolute_below_1e_minus_9": True})
    return rows


def main():
    response = read("RESCATTERING.json")["rows"]
    detuning = read("DETUNING.json")["s_over_m2"]
    rows = [detuning[str(s)] for s in ENERGIES]
    require(tuple(item["s_over_m2"] for item in response)==ENERGIES, "unexpected energy list")
    a = [Fraction(item["a"]) for item in rows]
    e = [Fraction(item["c"])/aa for item, aa in zip(rows,a,strict=True)]
    f = [-2*Fraction(item["b_imag"])/aa for item, aa in zip(rows,a,strict=True)]
    d = [Fraction(item["dlog_rate_dRe_delta_at_zero"]) for item in rows]
    ell = [Fraction(item["fractional_rate_derivative_at_lambda_zero"]) for item in response]
    matrix = sp.Matrix([[1, di, ei] for di, ei in zip(d,e,strict=True)])
    require(matrix.det()!=0 and f[0]!=0, "insufficient independent rate information")
    gain = matrix.inv()[2,:]
    require(sum(gain)==0, "overall normalization leaks into closure")
    sensitivity = matrix.inv()*sp.Matrix(ell)

    # Exact arithmetic on rational enclosures for the three inferred columns.
    numerator = determinant_interval(d, ell)
    denominator = determinant_interval(d, e)
    require(numerator[1]<0 and denominator[1]<0, "sign fails under input perturbation")
    bounds = multiply(numerator, (1/denominator[1], 1/denominator[0]))
    require(bounds[0]>0 and bounds[0]<sensitivity[2]<bounds[1], "invalid derivative enclosure")

    # Any complex detuning gives UW=V^2+T^2 with nonnegative W and positive U.
    normal, real, imag = sp.symbols("N x y", real=True)
    vector = sp.Matrix((normal, normal*real, normal*(real**2+imag**2)))
    measured = matrix*vector + sp.Matrix(f)*(normal*imag)
    reversed_rates = matrix*vector - sp.Matrix(f)*(normal*imag)
    T = sp.cancel((measured[0]-reversed_rates[0])/(2*f[0]))
    U, V, W = (matrix.inv()*(measured-sp.Matrix(f)*T)).applyfunc(sp.expand)
    require((U,V,W,T)==(normal,normal*real,normal*real**2+normal*imag**2,normal*imag),
            "symbolic reconstruction")
    require(sp.expand(U*W-V**2-T**2)==0, "symbolic complex closure")
    samples = finite_box_test(np.array(matrix,dtype=float), a, np.array(f,dtype=float))
    gain_norm = sum(abs(v) for v in gain)
    return {"status":"PASS", "energies_s_over_m2":ENERGIES,
            "matrix_determinant":float(f"{float(matrix.det()):.11g}"),
            "quartic_solution_derivatives_U_V_W_over_normalization":[float(f"{float(v):.11g}") for v in sensitivity],
            "quartic_closure_slope_one_ppm_input_enclosure":[float(bounds[0]),float(bounds[1])],
            "reciprocal_mean_first_order_absolute_error_gain":float(f"{float(gain_norm):.11g}"),
            "signed_half_difference_coefficients":[float(f"{float(v):.10g}") for v in f],
            "exact_complex_identity":True, "independent_box_examples":samples}


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args = parser.parse_args()
    data=json.dumps(main(),indent=2,allow_nan=False)+"\n"
    if args.output:
        args.output.write_text(data)
    else:
        print(data,end="")
