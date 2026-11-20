"""Calibrate closed-universe Majorana occupations against soluble conformal-time histories."""
import argparse
import importlib.util
import json
from math import cosh, exp, hypot, pi, sin
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.linalg import expm
from scipy.optimize import brentq


def check(ok, label):
    if not ok:
        raise ArithmeticError(label)


def show(x):
    return float(f"{x:.11g}")


def occupation(k, mi, mf, width):
    """Exact normalized tanh-quench probability for nonnegative endpoint masses (k > 0)."""
    k = np.asarray(k, dtype=float)
    if mi == mf:
        return np.zeros_like(k)
    wi, wf = np.hypot(k, mi), np.hypot(k, mf)
    dm = abs(mf - mi)
    dw = dm * abs(mf + mi) / (wf + wi)
    # Rationalize dm-dw at small momentum for positive endpoints.
    small = dm * (k*k/(wi+abs(mi)) + k*k/(wf+abs(mf))) / (wi+wf)
    x = 0.5*pi*width*small
    y = 0.5*pi*width*(dm+dw)
    # sinh(x) sinh(y) / [sinh(pi*T*wi) sinh(pi*T*wf)] without overflow.
    return (np.exp(pi*width*(dm-wi-wf)) * np.expm1(-2*x)*np.expm1(-2*y)
            / (np.expm1(-2*pi*width*wi)*np.expm1(-2*pi*width*wf)))


def numerical(k, mi, mf, width):
    """Evolve the unrotated spinor from its normalized incoming positive eigenvector."""
    angle_i = np.arctan2(k, mi)/2
    initial = np.r_[np.cos(angle_i), np.sin(angle_i)].astype(complex)

    def rhs(x, spinor):
        m = 0.5*(mf+mi + (mf-mi)*np.tanh(x))
        upper, lower = np.split(spinor, 2)
        return -1j * width * np.r_[m*upper + k*lower, k*upper - m*lower]

    result = solve_ivp(rhs, (-18, 18), initial, t_eval=[18], method="DOP853",
                       rtol=1e-12, atol=2e-14)
    check(result.success, result.message)
    upper, lower = np.split(result.y[:, -1], 2)
    check(np.max(abs(abs(upper)**2 + abs(lower)**2 - 1)) < 1e-10, "unitary norm")
    angle_f = np.arctan2(k, mf)/2
    return abs(-np.sin(angle_f)*upper + np.cos(angle_f)*lower)**2


def weighted_sums(mi, mf, width, nmax):
    n = np.arange(nmax+1)
    k = n + 1.5
    counts = 2*(n+1)*(n+2)*occupation(k, mi, mf, width)
    # Upper envelope P(k) <= C exp(-2*pi*T*k), for the omitted k >= k0.
    j = nmax+1
    k0 = j+1.5
    c = 2*(cosh(pi*width*(mf-mi))-1) / (-np.expm1(-2*pi*width*k0))**2
    first = 2*(j+1)*(j+2)*c*exp(-2*pi*width*k0)
    q = exp(-2*pi*width)*(j+3)/(j+1)
    qe = q*(k0+1)/k0
    check(qe < 1, "energy tail ratios must be below one")
    return {"number": show(sum(counts)), "conformal_energy": show(counts @ np.hypot(k, mf)),
            "number_tail_upper": show(1.001*first/(1-q)),
            "conformal_energy_tail_upper": show(1.001*first*hypot(k0, mf)/(1-qe))}


def calibrate_steps():
    k = np.arange(9)+1.5
    results = {}
    for name, mi, mf, width in [
            ("rapid_expansion", 0.5, 2.0, 0.12), ("slow_expansion", 0.5, 2.0, 0.6),
            ("rapid_contraction", 2.0, 0.5, 0.12), ("massive_static", 1.5, 1.5, 0.4),
            ("massless", 0.0, 0.0, 0.4)]:
        exact = occupation(k, mi, mf, width)
        observed = numerical(k, mi, mf, width)
        check(np.allclose(observed, exact, rtol=3e-6, atol=2e-15), name+" analytic step")
        # Changing the sign of the intrinsic eigenvalue is a constant spinor rotation.
        check(np.allclose(numerical(-k, mi, mf, width), exact, rtol=3e-6, atol=2e-15),
              name+" opposite intrinsic spectrum")
        angle_change = np.arctan2(k, mi)-np.arctan2(k, mf)
        cap = np.sin(angle_change/2)**2
        check(np.all(exact <= cap+2e-15), name+" monotone rotation cap")
        limit = occupation(k, mi, mf, 1e-8)
        check(np.allclose(limit, cap, rtol=1e-12, atol=1e-16), name+" sudden saturation")
        sums = weighted_sums(mi, mf, width, 128)
        check(sums["conformal_energy_tail_upper"] < 1.e-32, name+" integrated tail")
        results[name] = {"M_in": mi, "M_out": mf, "T": width,
                         "P_n_0_through_8": [show(p) for p in exact],
                         "sums_through_n_128": sums}
    return results


def production_calibration():
    """Use the archive solver with a proper-time clock on a soluble geometry."""
    path = Path(__file__).with_name("check_modes.py")
    spec = importlib.util.spec_from_file_location("wall_modes", path)
    modes = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(modes)
    width, ai, af = 0.12, 0.5, 2.
    mid, half = (ai+af)/2, (af-ai)/2

    def clock(x):
        return width*(mid*(x+18) + half*(np.log(np.cosh(x))-np.log(np.cosh(18))))

    end = clock(18.)
    def radius(t):
        x = brentq(lambda x: clock(x)-t, -18., 18.)
        return mid+half*np.tanh(x)

    masses = np.array([1., 4.])
    solve = modes.integrate(radius, end, masses)
    direct = modes.probabilities(solve(end), radius(end), masses)
    exact = np.array([occupation(modes.K, m*ai, m*af, width) for m in masses])
    check(np.allclose(direct, exact, rtol=3e-7, atol=2e-14),
          "published proper-time propagator vs closed-form quench, modes 0..32")
    return {"a_in": ai, "a_out": af, "conformal_T": width, "proper_masses": masses.tolist(),
            "shells": [0, 32], "probability_rtol": 3e-7, "probability_atol": 2e-14,
            "exact_band_numbers": [show(row @ modes.WEIGHTS) for row in exact]}


def cyclic_pulse():
    # Initial/final mass m0, constant intervening mass m1, duration L.
    k, m0, m1 = 1.5, 0.5, 2.0
    h0 = np.array([[m0, k], [k, -m0]])
    h1 = np.array([[m1, k], [k, -m1]])
    _, states = np.linalg.eigh(h0)
    out = []
    for phase in (pi/2, pi):
        duration = phase/hypot(k, m1)
        beta = states[:, 0] @ expm(-1j*duration*h1) @ states[:, 1]
        exact = (k*(m1-m0)/(hypot(k, m0)*hypot(k, m1))*sin(phase))**2
        check(abs(abs(beta)**2-exact) < 2.e-14, "two sudden jumps")
        out.append({"plateau_phase_over_pi": show(phase/pi),
                    "P": 0.0 if phase == pi else show(exact)})
    return out


def flat_limit():
    # Hold the masses and time scale in final physical coordinates fixed as R grows.
    mi, mf, width = 0.4, 1.0, 0.35

    def density_integrand(p):
        return occupation(p, mi, mf, width)*p*p/pi**2

    flat, error = quad(density_integrand, 0, np.inf, epsabs=1e-15, epsrel=1e-12)
    check(error < 1e-13, "flat continuum reference quadrature")
    correction = quad(lambda p: occupation(p, mi, mf, width)/(4*pi**2), 0, np.inf,
                      epsabs=1e-15, epsrel=1e-12)[0]
    rows = []
    for radius in (4, 8, 16, 32, 64):
        n = np.arange(40*radius)
        p = (n+1.5)/radius
        density = np.dot((n+1)*(n+2), occupation(p, mi, mf, width))/(pi**2*radius**3)
        rows.append({"R": radius, "number_per_final_physical_volume": show(density),
                     "relative_difference": show(density/flat-1)})
    check(abs(rows[-1]["relative_difference"]) < 2e-4, "flat two-spin measure")
    gap = flat-rows[-1]["number_per_final_physical_volume"]
    check(abs(gap*rows[-1]["R"]**2/correction-1) < 1.e-5, "finite-radius correction")
    return {"physical_m_in_effective": mi, "physical_m_out": mf, "physical_width": width,
            "flat_number_density": show(flat), "leading_curvature_coefficient": show(correction),
            "sphere_densities": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    record = {"tanh_steps": calibrate_steps(), "cyclic_pulse": cyclic_pulse(),
              "flat_limit": flat_limit(), "production_solver": production_calibration(),
              "checks": {"direct_spinor_ODE": "both signed spatial eigenvalues, n=0..8",
                         "ODE_interval_over_T": [-18, 18], "ODE_method": "DOP853",
                         "probability_rtol": 3e-6, "probability_atol": 2e-15,
                         "sudden_width": 1e-8, "closed_sums_n_max": 128}}
    output = json.dumps(record, indent=2, allow_nan=False)+"\n"
    if args.output:
        args.output.write_text(output, encoding="utf-8")
    else:
        print(output, end="")


if __name__ == "__main__":
    main()
