"""Finite-band Majorana occupations on the archived closed-wall geometry."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator
from scipy.optimize import brentq, minimize_scalar

HERE = Path(__file__).resolve().parent
TRACE = HERE / "trace_geometry.csv.gz"
GEOMETRY_SHA256 = "7b89924624d079399e8c260b7ee41a3a9c84cd5d28ce4ec7f0389e2b38412f83"
MASSES = np.array([1., 4.])
SHELL = np.arange(33)
K = SHELL + 1.5
WEIGHTS = 2 * (SHELL + 1) * (SHELL + 2)


def read_geometry(trace_path):
    table = np.loadtxt(trace_path, delimiter=",", skiprows=1, usecols=(0, 1, 2))
    if hashlib.sha256(table.astype("<f8").tobytes()).hexdigest() != GEOMETRY_SHA256:
        raise ValueError("Unexpected numerical t,a,H trace")
    return table.T


def make_geometry_nodes(trace_path=TRACE):
    t, a, hubble = read_geometry(trace_path)
    smooth_h = PchipInterpolator(t, hubble)
    base = CubicHermiteSpline(t, a, a * hubble)
    crossings = np.flatnonzero(hubble[:-1] * hubble[1:] < 0)
    if len(crossings) != 2 or hubble[crossings[0]] >= 0 or hubble[crossings[1]] <= 0:
        raise ValueError("Expected a bounce followed by a turnaround")
    roots = np.array([brentq(smooth_h, t[j], t[j+1]) for j in crossings])
    before_turn = t < roots[-1]
    knots = np.r_[t[before_turn], roots]
    radii = np.r_[a[before_turn], base(roots)]
    slopes = np.r_[(a * hubble)[before_turn], 0., 0.]
    order = np.argsort(knots)
    events = dict(bounce=roots[0], peak_expansion=t[np.argmax(hubble)],
                  expanding_slice=11.03911, turnaround=roots[1])
    return knots[order], radii[order], slopes[order], events


def make_geometry(trace_path=TRACE):
    t, a, slopes, events = make_geometry_nodes(trace_path)
    return CubicHermiteSpline(t, a, slopes, extrapolate=False), events


def angles(radius, masses):
    return 0.5 * np.arctan2(K / float(radius), masses[:, None])


def negative_eigenstate(radius, masses):
    angle = angles(radius, masses)
    return np.stack((-np.sin(angle), np.cos(angle))).astype(complex)


def probabilities(spinors, radius, masses):
    angle = angles(radius, masses)
    return np.abs(np.cos(angle) * spinors[0] + np.sin(angle) * spinors[1])**2


def integrate(radius, end, masses=MASSES, start=0.):
    initial = negative_eigenstate(radius(start), masses)
    mass = masses[:, None]

    def rhs(t, state):
        v, w = state.reshape(initial.shape)
        p = K / radius(t)
        return (-1j * np.stack((mass * v + p * w, p * v - mass * w))).ravel()

    sol = solve_ivp(rhs, (start, end), initial.ravel(), method="DOP853",
                    rtol=1e-11, atol=2e-13, max_step=.018, dense_output=True)
    if not sol.success:
        raise ArithmeticError(sol.message)
    return lambda t: sol.sol(t).reshape(initial.shape + np.shape(t))


def midpoint_probabilities(radius, events, halvings):
    # Unitary exponential of the frozen two-level Hamiltonian at each midpoint.
    knots = np.unique(np.r_[radius.x, list(events.values())])
    for _ in range(halvings):
        knots = np.sort(np.r_[knots, (knots[:-1] + knots[1:]) / 2])
    dt = np.diff(knots)
    mid_a = radius((knots[:-1] + knots[1:]) / 2)
    v, w = negative_eigenstate(radius(0), MASSES)
    mass = MASSES[:, None]
    result = {}
    start = 0
    for name, time in events.items():
        stop = np.searchsorted(knots, time)
        # Batches limit memory as the shell count or trace resolution increases.
        for pos in range(start, stop, 512):
            span = slice(pos, min(pos + 512, stop))
            p = K / mid_a[span, None, None]
            omega = np.hypot(p, mass)
            phase = omega * dt[span, None, None]
            cosine = np.cos(phase)
            sine = -1j * np.sin(phase) / omega
            for c, s, momentum in zip(cosine, sine, p):
                v, w = (c * v + s * (mass * v + momentum * w),
                        c * w + s * (momentum * v - mass * w))
        if np.max(np.abs(np.abs(v)**2 + np.abs(w)**2 - 1)) > 2e-10:
            raise ArithmeticError("Midpoint spinor lost normalization")
        result[name] = probabilities((v, w), radius(time), MASSES)
        start = stop
    return result


def rounded(x, digits=8):
    return float(f"{x:.{digits}g}")


def summarize(p, m, a):
    particles = WEIGHTS * p
    volume = 2 * np.pi**2 * a**3
    number = np.sum(particles)
    energy = np.dot(particles, np.hypot(m, K / a))
    return dict(probabilities=[rounded(x) for x in p],
                number=rounded(number), energy=rounded(energy),
                density=rounded(number / volume), energy_density=rounded(energy / volume),
                shell_of_largest_probability=int(np.argmax(p)),
                shell_of_smallest_probability=int(np.argmin(p)),
                shell_of_largest_number_contribution=int(np.argmax(particles)))


def peak_numbers(solution, radius, stop):
    # Locate all sampled local peaks and refine them in wall proper time.
    grid = np.linspace(0., stop, 4097)
    v, w = solution(grid)
    theta = .5 * np.arctan2(K[:, None] / radius(grid), MASSES[:, None, None])
    p = np.abs(np.cos(theta) * v + np.sin(theta) * w)**2
    numbers = np.einsum("n,mnt->mt", WEIGHTS, p)
    peaks = {}
    for i, sampled in enumerate(numbers):
        def number(t):
            state = solution(t)[:, i:i+1]
            return np.dot(WEIGHTS, probabilities(state, radius(t), MASSES[i:i+1])[0])
        candidates = [(0., sampled[0]), (stop, sampled[-1])]
        indices = 1 + np.flatnonzero((sampled[1:-1] >= sampled[:-2]) &
                                    (sampled[1:-1] > sampled[2:]))
        for j in indices:
            opt = minimize_scalar(lambda t: -number(t), bounds=(grid[j-1], grid[j+1]),
                                  method="bounded", options={"xatol": 1e-10})
            candidates.append((opt.x, -opt.fun))
        t, total = max(candidates, key=lambda point: point[1])
        peaks[str(int(MASSES[i]))] = dict(t=rounded(t, 6), number=rounded(total))
    return peaks


def verify(trace_path=TRACE):
    t, a, hubble = read_geometry(trace_path)
    radius, events = make_geometry(trace_path)
    stop = events["turnaround"]
    solution = integrate(radius, stop)
    coarse = midpoint_probabilities(radius, events, 0)
    medium = midpoint_probabilities(radius, events, 1)
    fine = midpoint_probabilities(radius, events, 2)
    result = {}

    for name, time in events.items():
        spinors = solution(time)
        adaptive = probabilities(spinors, radius(time), MASSES)
        extrapolated = medium[name] + (medium[name] - coarse[name]) / 3
        p = fine[name] + (fine[name] - medium[name]) / 3
        if np.max(np.abs(np.sum(np.abs(spinors)**2, axis=0) - 1)) > 1e-10:
            raise ArithmeticError("Adaptive spinor lost normalization")
        # Compare the probabilities and the degeneracy-weighted totals.
        if (np.min(p) < 0 or np.max(p) > 1 or
                np.max(abs(p - adaptive)) > 8e-12 or
                np.max(abs((p - adaptive) @ WEIGHTS)) > 2e-8 or
                np.max(abs(p - extrapolated)) > 2e-13 or
                np.max(abs((p - extrapolated) @ WEIGHTS)) > 5e-10):
            raise ArithmeticError(f"Independent check of occupations failed at {name}")
        result[name] = dict(t=rounded(time, 13), a=rounded(radius(time), 13),
                            H=rounded(PchipInterpolator(t, hubble)(time), 9)
                            if name not in ("bounce", "turnaround") else 0.,
                            species={str(int(m)): summarize(row, m, float(radius(time)))
                                     for row, m in zip(p, MASSES)})

    # Solve the zero-mass and frozen-geometry problems with the same ODE and
    # projections. Both have a time-independent eigenbasis for every shell.
    for control_radius, masses in ((radius, np.array([0.])),
                                 (lambda t: float(a[0]), MASSES)):
        evolution = integrate(control_radius, stop, masses)
        for when in events.values():
            if np.max(probabilities(evolution(when), control_radius(when), masses)) > 1e-22:
                raise ArithmeticError("A zero-mixing control produced particles")

    # Whittaker-function reference for the exponential-radius comparison in SOURCES.md.
    # The fixed momentum endpoints x=k/a eliminate the shell index from the exact result.
    test_masses = np.array([.2, .55, 1.])
    exact = np.array([.04085885340664992, .04783619338433926, .009110979486389349])
    for n in (0, 32):
        start, end = np.log(K[n] / np.array([8., .25]))
        evolved = integrate(np.exp, end, test_masses, start=start)(end)
        measured = probabilities(evolved, np.exp(end), test_masses)[:, n]
        if np.max(np.abs(measured - exact)) > 3e-11:
            raise ArithmeticError("Exponential-radius analytical reference failed")

    sensitivity = {}
    for label, alternate in (("pchip_radius", PchipInterpolator(t, a)),
                             ("excursion_times_0.99", lambda x: a[0] + .99 * (radius(x) - a[0])),
                             ("excursion_times_1.01", lambda x: a[0] + 1.01 * (radius(x) - a[0]))):
        alternate_solution = integrate(alternate, stop)
        sensitivity[label] = {}
        for name in ("expanding_slice", "turnaround"):
            when = events[name]
            p = probabilities(alternate_solution(when), alternate(when), MASSES)
            totals = p @ WEIGHTS
            sensitivity[label][name] = {str(int(m)): rounded(n) for m, n in zip(MASSES, totals)}
            if label == "pchip_radius" and np.max(abs(totals - fine[name] @ WEIGHTS)) > 1e-7:
                raise ArithmeticError("Radius interpolation sensitivity exceeds tolerance")

    return dict(shells=list(range(len(SHELL))), multiplicity="2*(n+1)*(n+2)",
                t0=0, a0=2, geometry_float64_le_sha256=GEOMETRY_SHA256,
                events=result, sampled_and_refined_peak_number=peak_numbers(solution, radius, stop),
                geometric_sensitivity_number=sensitivity)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--trace", type=Path, default=TRACE)
    args = parser.parse_args()
    data = json.dumps(verify(args.trace), indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(data, encoding="utf8")
    else:
        print(data, end="")
