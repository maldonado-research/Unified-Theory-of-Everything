"""Independent controls in a rotated basis with an implicit real ODE solver.

Does not import or alter benchmark.py. It reads frozen benchmark results.
Not a certified floating-point enclosure.
"""
from pathlib import Path
import json
import hashlib
import math
import sys

if sys.flags.optimize:
    raise SystemExit("Run this review without -O; its scientific assertions must remain enabled.")
import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent


def radau_probability(k, a, tau, cutoff=14, tol=5e-12):
    # Hadamard rotation of the benchmark Hamiltonian:
    # H'=k sigma_z + m sigma_x. Get states from a general eigensolver.
    asymptotic = np.array([[k, 1.0], [1.0, -k]])
    _, eigvecs = np.linalg.eigh(asymptotic)
    negative, positive = eigvecs[:, 0], eigvecs[:, 1]
    y0 = np.r_[negative, [0.0, 0.0]]

    def ham(t):
        m = 1.0 - a / math.cosh(t / tau) ** 2
        return np.array([[k, m], [m, -k]])

    def rhs(t, y):
        h = ham(t)
        return np.r_[h @ y[2:], -h @ y[:2]]

    def jac(t, y):
        h = ham(t)
        z = np.zeros((2, 2))
        return np.block([[z, h], [-h, z]])

    sol = solve_ivp(
        rhs, (-cutoff*tau, cutoff*tau), y0, method="Radau",
        jac=jac, rtol=tol, atol=tol/10)
    if not sol.success:
        raise RuntimeError(sol.message)
    state = sol.y[:2, -1] + 1j*sol.y[2:, -1]
    return {
        "probability": float(abs(positive @ state)**2),
        "norm_error": float(abs(np.vdot(state, state)-1)),
        "function_evaluations": sol.nfev,
    }


def main():
    frozen = HERE / "BENCHMARK_RESULTS.json"
    source_bytes = frozen.read_bytes()
    benchmark = json.loads(source_bytes)
    controls = []
    selected = [
        ("constant", 0.6), ("shallow_fast", 0.6),
        ("crossing_fast", 0.55), ("crossing_slow", 0.25),
        ("shallow_slow", 0.275),
    ]
    for label, k in selected:
        item = next(c for c in benchmark["cases"] if c["name"] == label)
        index = int(np.argmin(abs(np.array(item["k_over_M"])-k)))
        result = radau_probability(k, item["pulse_depth_a"], item["width_M_tau"])
        expected = item["occupation_per_helicity"][index]
        gap = abs(result["probability"]-expected)
        assert gap < 2e-7, (label, result, expected)
        assert result["norm_error"] < 2e-9
        controls.append({"case": label, "k": k, **result,
                         "benchmark_probability": expected, "gap": gap})

    plus = radau_probability(0.55, 2, 1)
    minus = radau_probability(-0.55, 2, 1)
    zero = radau_probability(0, 2, 1)
    assert abs(plus["probability"]-minus["probability"]) < 1e-10
    assert zero["probability"] < 1e-20

    # Analytic first Born term for m=1-a sech^2(t/tau), E=sqrt(1+k^2):
    # |A1| = 2*pi*a*k*tau^2 / sinh(pi*E*tau).
    # The exact unitary Duhamel remainder is <= K^2/2, K=2*a*tau.
    a, tau, k = 0.01, 1.0, 0.5
    energy = math.sqrt(1+k*k)
    born_amplitude = 2*math.pi*a*k*tau*tau/math.sinh(math.pi*energy*tau)
    remainder_bound = 2*a*a*tau*tau
    exact_tail_state_bound = 4*a*tau*math.exp(-28)
    predicted_lower = max(0, born_amplitude-remainder_bound)**2
    predicted_upper = (born_amplitude+remainder_bound)**2
    weak = radau_probability(k, a, tau)
    assert predicted_lower > 0
    assert predicted_lower < weak["probability"] < predicted_upper
    weak_control = {
        "a": a, "tau": tau, "k": k,
        "first_Born_amplitude_magnitude": born_amplitude,
        "unitary_Duhamel_amplitude_remainder_bound": remainder_bound,
        "analytic_infinite_time_probability_lower_bound": predicted_lower,
        "analytic_infinite_time_probability_upper_bound": predicted_upper,
        "finite_window_tail_state_bound": exact_tail_state_bound,
        "rotated_Radau": weak,
    }

    result = {
        "status": "INDEPENDENT_CONTROLS_PASS",
        "benchmark_results_sha256": hashlib.sha256(source_bytes).hexdigest(),
        "solver": "scipy Radau, real four-component state, Hadamard-rotated Hamiltonian",
        "not_certified_float_enclosure": True,
        "selected_modes": controls,
        "helicity_reflection_control": {
            "plus": plus, "minus": minus,
            "absolute_difference": abs(plus["probability"]-minus["probability"]),
        },
        "zero_momentum_crossing_control": zero,
        "weak_pulse_Born_control": weak_control,
    }
    (HERE / "INDEPENDENT_REVIEW_RESULTS.json").write_text(
        json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
