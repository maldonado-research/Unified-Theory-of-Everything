"""Resolve shells 33 to 2048 and compare independent bases and mesh refinements."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from check_modes import make_geometry, rounded

SAMPLES = [32, 33, 128, 512, 1024, 2048]
EVENTS = ("bounce", "expanding_slice", "turnaround")
BANDS = ((33, 128), (129, 512), (513, 2048))


def propagate(radius, times, work, executable, refinements, samples_only=False):
    knots = np.unique(np.r_[radius.x, [times[event] for event in EVENTS]])
    for _ in range(refinements):
        knots = np.sort(np.r_[knots, (knots[:-1] + knots[1:]) / 2])
    dt = np.diff(knots)
    mid = (knots[:-1] + knots[1:]) / 2
    offset = np.sqrt(3) / 6 * dt
    lo, hi = mid - offset, mid + offset
    grid = np.column_stack((dt, radius(lo), radius(lo, 1), radius(hi), radius(hi, 1)))
    path = work / "gauss_grid.bin"
    grid.astype(np.double).tofile(path)
    indices = np.searchsorted(knots, [times[event] for event in EVENTS[:-1]])
    args = [executable, path, *(str(i) for i in indices), str(float(radius(0.))),
            *(str(float(radius(times[event]))) for event in EVENTS),
            "samples" if samples_only else "all"]
    result = np.loadtxt(subprocess.check_output(args).splitlines())
    if not np.all(np.isfinite(result)) or np.max(result[:, -1]) > 1e-10:
        raise ArithmeticError("Lost high-band spinor normalization")
    return result


def compare(grids, low_band):
    coarse, fine, sampled = grids
    result = {e: {"species": {}} for e in EVENTS}
    for mass in (1, 4):
        rows = [g[(g[:, 1] == mass) & (g[:, 2] == 1)] for g in grids]
        n = rows[1][:, 0].astype(int)
        if not np.array_equal(n, np.arange(32, 2049)):
            raise ArithmeticError("Incomplete numerical band")
        g_n = 2 * (n + 1) * (n + 2)
        for j, event in enumerate(EVENTS, start=3):
            a = low_band["events"][event]["a"]
            omega = np.hypot(mass, (n + 1.5) / a)
            p = rows[1][:, j]
            diff = np.abs(p - rows[0][:, j])
            # Check the sum of absolute errors, avoiding cancellations between shells.
            if np.dot(g_n, diff) > 3e-9 or np.dot(g_n * omega, diff) > 2e-7:
                raise ArithmeticError("High-band mesh refinement failed")
            if abs(p[0] - low_band["events"][event]["species"][str(mass)]["probabilities"][32]) > 8e-14:
                raise ArithmeticError("Overlap with the 33-shell solution failed")
            spot = np.isin(n, SAMPLES)
            p_refined = rows[2][:, j]
            p_fixed = sampled[(sampled[:, 1] == mass) & (sampled[:, 2] == 0), j]
            if (not np.allclose(p[spot], p_refined, rtol=3e-5, atol=2e-23) or
                    not np.allclose(p_refined, p_fixed, rtol=3e-5, atol=2e-23)):
                raise ArithmeticError("Representative shell mesh/basis check failed")
            group = {}
            for start, end in BANDS:
                sel = (n >= start) & (n <= end)
                occup = g_n[sel] * p[sel]
                group[f"{start}_{end}"] = dict(number=rounded(np.sum(occup)),
                                               energy=rounded(np.dot(occup, omega[sel])))
            numbers = g_n[1:] * p[1:]
            result[event]["species"][str(mass)] = dict(
                number_33_2048=rounded(np.sum(numbers)),
                energy_33_2048=rounded(np.dot(numbers, omega[1:])), bands=group,
                representative_probabilities={str(n_i): rounded(p_i, 6)
                                              for n_i, p_i in zip(n[spot], p[spot])})
    return result


def verify():
    radius, times = make_geometry()
    low_band = json.loads((HERE / "MODES.json").read_text())
    with tempfile.TemporaryDirectory() as folder:
        work = Path(folder)
        executable = work / "band_evolution"
        subprocess.run([os.getenv("CXX", "g++"), "-std=c++17", "-O3", "-pthread",
                        "-ffp-contract=off", str(HERE / "band_evolution.cpp"), "-o", executable],
                       check=True)
        grids = [propagate(radius, times, work, executable, levels, samples_only=levels == 2)
                 for levels in (0, 1, 2)]
    return dict(multiplicity="2*(n+1)*(n+2)", events=compare(grids, low_band))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = json.dumps(verify(), indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(data, encoding="utf8")
    else:
        print(data, end="")
