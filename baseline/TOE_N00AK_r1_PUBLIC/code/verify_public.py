#!/usr/bin/env python3
"""Replay numerical identities with explicit checks active under python -O.

Results certify tree identities and conditional matrix consistency only.
They do not validate a derivative-vertex UV pole or any physical beta function.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import sys

import numpy as np

import numerical_core as core

ROOT = Path(__file__).resolve().parent
SEED = 20260904
TOLERANCE = 2e-11


def fro(matrix: np.ndarray) -> float:
    return float(np.linalg.norm(matrix, "fro"))


def relative(matrix: np.ndarray, scale: float) -> float:
    return fro(matrix) / max(scale, np.finfo(float).tiny)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def random_y(rng: np.random.Generator, rows: int, columns: int) -> np.ndarray:
    return (rng.normal(size=(rows, columns)) + 1j * rng.normal(size=(rows, columns))) / np.sqrt(columns)


def tree_checks(rng: np.random.Generator) -> dict:
    worst = {name: 0.0 for name in ("direct_spectral_sum", "remainder_identity", "known_node_reconstruction", "moment_recurrence", "complex_symmetry")}
    largest_bound_ratio = 0.0
    repeated_mass_cases = 0
    for trial in range(700):
        rows, columns = int(rng.integers(2, 7)), int(rng.integers(1, 6))
        y = random_y(rng, rows, columns)
        masses = np.sort(rng.uniform(1.5, 6.0, columns))
        if trial % 13 == 0 and columns > 1:
            masses[1] = masses[0]
            repeated_mass_cases += 1
        rho = rng.uniform(0.05, 0.7)
        z = min(masses)**2 * rho * np.exp(1j * rng.uniform(-np.pi, np.pi))
        order = trial % 6
        exact = core.tree_resolvent(y, masses, z)
        direct = np.zeros_like(exact)
        for j, mass in enumerate(masses):
            direct += np.outer(y[:, j], y[:, j]) / mass / (1 - z / mass**2)
        remainder = core.tree_remainder(y, masses, z, order)
        truncation = core.tree_truncation(y, masses, z, order)
        nodes = core.effective_nodes(y, masses)
        moments = [core.tree_moment(y, masses, n) for n in range(len(nodes) + 2)]
        reconstructed = core.tree_reconstruction(moments[:len(nodes)], nodes, z)
        polynomial = core.pole_polynomial(nodes)
        n = len(nodes) + 1
        recurrence = sum((polynomial[k] * moments[n-k] for k in range(len(nodes) + 1)), np.zeros_like(exact))
        recurrence_scale = sum(abs(polynomial[k]) * fro(moments[n-k]) for k in range(len(nodes) + 1))
        errors = {
            "direct_spectral_sum": relative(exact-direct, fro(exact)),
            "remainder_identity": relative(exact-truncation-remainder, fro(exact)+fro(truncation)+fro(remainder)),
            "known_node_reconstruction": relative(exact-reconstructed, fro(exact)),
            "moment_recurrence": relative(recurrence, recurrence_scale),
            "complex_symmetry": relative(exact-exact.T, 2*fro(exact)),
        }
        for name, value in errors.items():
            require(value < TOLERANCE, f"tree trial {trial}: {name}")
            worst[name] = max(worst[name], value)
        bound = core.remainder_bound(y, masses, z, order)
        ratio = float(np.linalg.norm(remainder, 2)) / bound
        require(ratio <= 1 + TOLERANCE, f"tree trial {trial}: spectral bound")
        largest_bound_ratio = max(largest_bound_ratio, ratio)
    # Exact degenerate cancellation and domain controls are deterministic.
    paired = np.array([[1, 1j], [2, 2j]], dtype=complex)
    require(len(core.effective_nodes(paired, np.array([2.0, 2.0]))) == 0, "zero residue deletion")
    require(fro(core.tree_resolvent(paired, np.array([2.0, 2.0]), 0.3)) == 0, "zero combined residue")
    rejected_pole = rejected_bound = False
    try:
        core.tree_resolvent(np.ones((2, 1)), np.array([2.0]), 4.0)
    except ZeroDivisionError:
        rejected_pole = True
    try:
        core.remainder_bound(np.ones((2, 1)), np.array([2.0]), 4.0, 1)
    except ValueError:
        rejected_bound = True
    require(rejected_pole and rejected_bound, "domain rejection")
    return {"pass": True, "trials": 700, "worst_relative_residuals": worst,
            "largest_spectral_remainder_to_bound_ratio": largest_bound_ratio,
            "repeated_mass_cases": repeated_mass_cases,
            "exact_zero_residue_deletion": True, "pole_and_bound_domain_rejection": True}


def conditional_checks(rng: np.random.Generator) -> dict:
    coefficients = [0, 1, 3, 7, 1+2j]
    worst = {name: 0.0 for name in ("local_identity", "rational_identity", "local_to_rational_contact", "local_to_rational_delta_y")}
    mismatched_control_minimum = float("inf")
    mismatch_trials = 0
    for trial in range(500):
        rows, active, removed = int(rng.integers(2, 7)), int(rng.integers(1, 5)), int(rng.integers(1, 5))
        ya, yd = random_y(rng, rows, active), random_y(rng, rows, removed)
        ma, md = np.sort(rng.uniform(0.12, 0.96, active)), np.sort(rng.uniform(2.0, 6.5, removed))
        order, coefficient = trial % 6, coefficients[trial % len(coefficients)]
        _, contact, ds = core.conditional_local(ya, ma, core.tree_moment(yd, md, order), order, coefficient)
        delta_y_exact, contact_exact, ds_exact = core.conditional_rational(ya, ma, yd, md, coefficient)
        contact_sum, delta_y_sum = np.zeros_like(contact), np.zeros_like(ya)
        for n in range(25):
            dy_n, dc_n, _ = core.conditional_local(ya, ma, core.tree_moment(yd, md, n), n, coefficient)
            contact_sum += dc_n
            delta_y_sum += dy_n
        errors = {
            "local_identity": relative(contact+ds, fro(contact)+fro(ds)),
            "rational_identity": relative(contact_exact+ds_exact, fro(contact_exact)+fro(ds_exact)),
            "local_to_rational_contact": relative(contact_sum-contact_exact, fro(contact_exact)),
            "local_to_rational_delta_y": relative(delta_y_sum-delta_y_exact, fro(delta_y_exact)),
        }
        for name, value in errors.items():
            require(value < TOLERANCE, f"conditional trial {trial}: {name}")
            worst[name] = max(worst[name], value)
        if coefficient != 0:
            wrong_contact = 2*contact
            mismatch = relative(ds+wrong_contact, fro(ds)+fro(wrong_contact))
            require(mismatch > 0.25, f"conditional trial {trial}: mismatched contact control")
            mismatched_control_minimum = min(mismatched_control_minimum, mismatch)
            mismatch_trials += 1
    return {"pass": True, "trials": 500, "orders_tested": list(range(6)),
            "coefficient_inputs": [{"real": float(complex(c).real), "imag": float(complex(c).imag)} for c in coefficients],
            "cases_per_coefficient": 100, "worst_relative_residuals": worst,
            "mismatched_contact_control_trials": mismatch_trials,
            "mismatched_contact_control_minimum_relative_residual": mismatched_control_minimum,
            "interpretation": "The assumed paired matrices cancel for every tested common coefficient, including complex and zero values. This identity does not identify a physical coefficient. Changing only the contact coefficient breaks the assumed relation."}


def frozen_checks(y: np.ndarray) -> dict:
    rows = []
    for removed, probe, expected in [((2,), 1, (0.25, 0.0625)), ((2,), 0, (0.0625, 0.00390625)), ((1,), 0, (0.25, 0.0625)), ((1,2), 0, (0.221643102559, 0.0562145822425))]:
        yd, md, z = y[:, removed], core.MASSES_GEV[list(removed)], core.MASSES_GEV[probe]**2
        exact = core.tree_resolvent(yd, md, z)
        errors = [relative(core.tree_truncation(yd, md, z, n)-exact, fro(exact)) for n in (0,1)]
        require(all(abs(a-b) < 2e-11 for a,b in zip(errors,expected)), "frozen benchmark")
        nodes = core.effective_nodes(yd, md)
        reconstructed = core.tree_reconstruction([core.tree_moment(yd, md,n) for n in range(len(nodes))], nodes,z)
        residual = relative(reconstructed-exact, fro(exact))
        require(residual < TOLERANCE, "frozen reconstruction")
        rows.append({"removed_one_based": [j+1 for j in removed], "probe_one_based": probe+1,
                     "tree_d5_relative_error": errors[0], "tree_d5_plus_d7_relative_error": errors[1],
                     "known_node_reconstruction_relative_residual": residual})
    ya, ma, yd, md = y[:,:2], core.MASSES_GEV[:2], y[:,2:], core.MASSES_GEV[2:]
    _, dc, ds = core.conditional_local(ya,ma,core.tree_moment(yd,md,1),1,3)
    require(relative(dc+ds,fro(dc)+fro(ds)) < TOLERANCE, "frozen conditional identity")
    return {"pass": True, "tree_cases": rows,
            "conditional_c_equals_3_d7_contact_norm_GeV_inverse": fro(dc),
            "conditional_c_equals_3_d7_schur_variation_norm_GeV_inverse": fro(ds),
            "conditional_norm_interpretation": "Values of supplied matrices, not established physical beta components or measured quantities."}


def coefficient_nonidentification(y: np.ndarray) -> dict:
    # Hold the same W1,Y,M fixed while varying the coefficient itself.
    ya, ma, yd, md = y[:,:2], core.MASSES_GEV[:2], y[:,2:], core.MASSES_GEV[2:]
    w1 = core.tree_moment(yd,md,1)
    cases = []
    for coefficient in [0,1,3,7,1+2j]:
        _, contact, ds = core.conditional_local(ya,ma,w1,1,coefficient)
        residual = relative(contact+ds,fro(contact)+fro(ds))
        require(residual < TOLERANCE,"common coefficient diagnostic")
        cases.append({"coefficient": {"real": float(complex(coefficient).real), "imag": float(complex(coefficient).imag)},
                      "contact_norm_GeV_inverse": fro(contact), "relative_cancellation_residual": residual})
    return {"pass": True,"same_frozen_matrices": True,"cases": cases,
            "conclusion": "Cancellation alone cannot determine the value 3 (or any common coefficient)."}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--compare",type=Path,help="Compare deterministic results with an earlier normal or optimized run")
    args = parser.parse_args()
    require(sys.version_info >= (3,10),"Python 3.10 or newer is required")
    y = core.load_frozen_yukawa(ROOT/"source"/"TOE_N00AC_YNU_STAR_FULL_PRECISION.csv")
    rng = np.random.default_rng(SEED)
    semantic = {"release":"TOE-N00AK-r1", "all_pass":True,"seed":SEED,"relative_tolerance":TOLERANCE,
                "claim_scope":{"tree_identities":"numerically checked under stated assumptions",
                               "d7_and_all_odd_matrices":"conditional algebra only",
                               "physical_loop_coefficients_established":False,
                               "independent_component_loop_calculation_performed":False},
                "tree":tree_checks(rng),"conditional":conditional_checks(rng),
                "frozen":frozen_checks(y),"common_coefficient_diagnostic":coefficient_nonidentification(y)}
    report = {"runtime":{"python":platform.python_version(),"numpy":np.__version__,
                         "implementation":platform.python_implementation(),"platform":platform.system()+" "+platform.machine(),
                         "optimization_level":sys.flags.optimize,"utc_completed":datetime.now(timezone.utc).isoformat()},
              "input_hashes":{path.as_posix():core.sha256_file(ROOT/path) for path in [Path("numerical_core.py"),Path("verify_public.py"),Path("source/TOE_N00AC_YNU_STAR_FULL_PRECISION.csv")]},
              "semantic":semantic}
    if args.compare:
        previous = json.loads(args.compare.read_text(encoding="utf-8"))
        equal = previous["semantic"] == semantic and previous["input_hashes"] == report["input_hashes"]
        require(equal,"normal/optimized deterministic results differ")
        report["comparison"]={"input":args.compare.name,"semantic_and_input_hashes_equal":True,
                              "runtime_fields_compared":False}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({"pass":True,"tree_trials":700,"conditional_trials":500,"optimization_level":sys.flags.optimize,
                      "output":args.output.name,"physical_loop_coefficients_established":False}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
