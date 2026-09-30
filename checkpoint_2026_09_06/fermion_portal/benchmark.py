"""Illustrative, externally prescribed Majorana mass pulse; not an HDBLAST fit.

Units: asymptotic mass M=1.  Each helicity reduces (up to a constant basis
rotation) to H=m(t)*sigma_z+k*sigma_x.  In/out states use the same constant
asymptotic Hamiltonian.  Independent solvers: unitary midpoint products and
SciPy DOP853.  The analytic discretization bound excludes floating rounding.
"""
from pathlib import Path
import json
import math
import platform
import numpy as np
import scipy
from scipy.integrate import solve_ivp, simpson

HERE = Path(__file__).resolve().parent
CASES = [("constant", 0.0, 1.0), ("weak_certified", 0.01, 1.0), ("shallow_fast", 0.5, 1.0),
         ("crossing_fast", 2.0, 1.0), ("crossing_slow", 2.0, 4.0),
         ("shallow_slow", 0.5, 4.0)]


def mass(t, depth, width):
    return 1.0 - depth / np.cosh(t / width) ** 2


def eigenspinors(k):
    k = np.asarray(k, dtype=float)
    angle = np.arctan2(k, 1.0)
    # Columns: negative energy (-sin, cos); positive energy (cos, sin).
    negative = np.stack([-np.sin(angle/2), np.cos(angle/2)], axis=-1)
    positive = np.stack([np.cos(angle/2), np.sin(angle/2)], axis=-1)
    return negative.astype(complex), positive.astype(complex)


def unitary_midpoint(k, depth, width, steps=65536, cutoff=12.0):
    k = np.asarray(k, dtype=float)
    psi, positive = eigenspinors(k)
    dt = 2*cutoff*width/steps
    for j in range(steps):
        t = -cutoff*width+(j+0.5)*dt
        m = mass(t, depth, width)
        energy = np.sqrt(m*m+k*k)
        co = np.cos(energy*dt)
        # sinc includes the continuous energy=0 limit.
        si = dt*np.sinc(energy*dt/np.pi)
        old0, old1 = psi[:, 0].copy(), psi[:, 1].copy()
        psi[:, 0] = co*old0-1j*si*(m*old0+k*old1)
        psi[:, 1] = co*old1-1j*si*(k*old0-m*old1)
    probability = abs(np.sum(positive.conj()*psi, axis=-1))**2
    norm_error = max(abs(np.sum(abs(psi)**2, axis=-1)-1.0))
    # Duhamel + total variation(m)<=2*depth, with uniform bin width dt.
    state_discretization_bound = depth*dt
    state_tail_bound = 4*depth*width*math.exp(-2*cutoff)
    state_bound = state_discretization_bound+state_tail_bound
    return probability, {
        "steps": steps, "time_cutoff_in_widths": cutoff,
        "max_float_norm_error": float(norm_error),
        "state_discretization_bound_exact_arithmetic": state_discretization_bound,
        "state_in_out_tail_bound": state_tail_bound,
        "probability_error_bound_excluding_float_rounding": 2*state_bound+state_bound**2,
    }


def adaptive_probability(k, depth, width, cutoff=12.0, tolerance=2e-12):
    initial, positive = eigenspinors([k])
    def rhs(t, psi):
        m = mass(t, depth, width)
        return -1j*np.array([m*psi[0]+k*psi[1], k*psi[0]-m*psi[1]])
    sol = solve_ivp(rhs, (-cutoff*width, cutoff*width), initial[0],
                    method="DOP853", rtol=tolerance, atol=tolerance/10)
    if not sol.success:
        raise RuntimeError(sol.message)
    state = sol.y[:, -1]
    return float(abs(np.vdot(positive[0], state))**2), float(abs(np.vdot(state, state)-1))


def main():
    ks = np.linspace(0, 8, 321)
    checks = []
    records = []
    selected_indices = [0, 4, 10, 20, 40, 80, 160, 320]
    for name, depth, width in CASES:
        low, low_meta = unitary_midpoint(ks, depth, width, 32768)
        high, high_meta = unitary_midpoint(ks, depth, width, 65536)
        comparisons = []
        for index in selected_indices:
            k = float(ks[index])
            adaptive, norm_error = adaptive_probability(k, depth, width)
            gap = abs(adaptive-high[index])
            comparisons.append({"k_over_M": k, "unitary_midpoint": float(high[index]),
                                "DOP853": adaptive, "absolute_difference": float(gap),
                                "DOP853_norm_error": norm_error})
            if gap > 2e-6 or norm_error > 2e-8:
                raise AssertionError((name, k, gap, norm_error))
        convergence = float(max(abs(low-high)))
        if convergence > 2e-6:
            raise AssertionError((name, "midpoint refinement", convergence))
        if min(high) < -1e-12 or max(high) > 1+1e-9:
            raise AssertionError("Pauli/unitarity failure")
        if high_meta['max_float_norm_error'] > 1e-8 or abs(high[0]) > 1e-12:
            raise AssertionError("norm or exact zero-momentum control failed")
        if depth == 0 and max(high) > 1e-12:
            raise AssertionError("constant-mass control failed")
        energy = np.sqrt(1+ks**2)
        # Two equal helicities of one Majorana species. Finite band only.
        number_band = float(simpson(ks**2*high, x=ks)/np.pi**2)
        energy_band = float(simpson(ks**2*energy*high, x=ks)/np.pi**2)
        coarse_number = float(simpson(ks[::2]**2*high[::2], x=ks[::2])/np.pi**2)
        records.append({"name": name, "pulse_depth_a": depth, "width_M_tau": width,
                        "minimum_signed_mass_over_M": 1-depth,
                        "zero_mass_crossings": 2 if depth>1 else (1 if depth==1 else 0),
                        "k_over_M": ks.tolist(), "occupation_per_helicity": high.tolist(),
                        "solver_metadata": high_meta,
                        "max_probability_change_when_doubling_steps": convergence,
                        "independent_solver_comparisons": comparisons,
                        "max_sampled_occupation": float(max(high)),
                        "sampled_peak_k_over_M": float(ks[np.argmax(high)]),
                        "finite_band_number_density_over_M_cubed": number_band,
                        "finite_band_energy_density_over_M_fourth": energy_band,
                        "number_integral_change_when_halving_momentum_grid": abs(number_band-coarse_number),
                        "integration_domain_k_over_M": [0, 8],
                        "integral_scope": "Floating Simpson diagnostic; not a certified infinite-momentum abundance."})
        checks.append({"case": name, "status": "PASS"})
        print(name, "max n=", max(high), "density band=", number_band, flush=True)
    # Tail-cutoff variation and finer adaptive tolerance on populated modes.
    tail_controls=[]
    for name, depth, width in CASES[1:]:
        for k in [0.25, 1.0]:
            short, _ = adaptive_probability(k, depth, width, 12, 2e-12)
            long, _ = adaptive_probability(k, depth, width, 16, 3e-13)
            change=abs(short-long)
            if change > 1e-8:
                raise AssertionError("tail/tolerance check failed")
            tail_controls.append({"case":name, "k_over_M":k, "change":change})
    output={"status":"NUMERICAL_BENCHMARK_CHECKS_PASS", "physical_model":"prescribed real mass in flat spacetime",
            "not_established":["HDBLAST scalar/trajectory identification", "physical energy normalization", "backreaction", "thermalization", "CP asymmetry", "cosmological abundance", "new fundamental mathematics"],
            "conventions":{"mass":"m/M=1-a*sech^2(t/tau)","hamiltonian":"H/M=(m/M)*sigma_z+(k/M)*sigma_x", "vacuum":"negative-energy in mode; final positive-energy projection", "helicity_degeneracy":2},
            "software":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__},
            "checks":checks,"tail_and_tolerance_checks":tail_controls,"cases":records}
    (HERE/'BENCHMARK_RESULTS.json').write_text(json.dumps(output,indent=2)+'\n')
    print('All benchmark checks pass.',flush=True)


if __name__ == '__main__':
    main()
