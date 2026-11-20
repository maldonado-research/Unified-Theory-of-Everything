"""Check local on-shell N2 decay clocks for two independently prepared full angular towers."""

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid


BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE))
from check_modes import make_geometry

PARENT = 4.0
SCALAR = 1.0
YUKAWA = 0.5  # |Y_mu2| at |g| = 1
CHARGED_BRANCH = 0.5
SEED_PER_ACTION = 2 * math.pi**2 * 2**3 * 1e-4 / SCALAR


def require(condition, message):
    if not condition:
        raise ValueError(message)


def clock_grid(times, radius, momenta, width):
    """Integrate the proper clock and energy at birth in units hbar=c=L=1."""
    energy = np.hypot(PARENT, momenta / radius(times)[:, None])
    clock = cumulative_trapezoid(PARENT / energy, times, axis=0, initial=0)
    survival = np.exp(-width * clock)
    probability = -np.expm1(-width * clock[-1])
    # E times the wall-time hazard equals Gamma times the parent mass.
    birth_energy = width * PARENT * np.trapezoid(survival, times, axis=0)
    return clock[-1], probability, birth_energy


def trial(start, stop_time, stop_radius, radius, sample_times, width):
    shell_probability = np.asarray(start["species"]["4"]["probabilities"])
    labels = np.arange(len(shell_probability))
    require(len(labels) == 33, "registered band must be 0..32")
    require(np.all((0 <= shell_probability) & (shell_probability <= 1)), "Pauli bound")
    degeneracy = 2 * (labels + 1) * (labels + 2)
    population = degeneracy * shell_probability
    total = float(population.sum())
    require(math.isclose(total, start["species"]["4"]["number"], rel_tol=1e-7),
            "one-Majorana-species normalization")

    start_time = start["t"]
    delta = stop_time - start_time
    start_radius = float(radius(start_time))
    require(0 < delta and start_radius < stop_radius, "expected an expanding interval")
    t = np.r_[start_time, sample_times[(start_time < sample_times) & (sample_times < stop_time)], stop_time]
    require(np.min(np.diff(radius(t))) > -1e-9, "radius must grow between readouts")
    k = labels + 1.5
    clock, probabilities, energy = clock_grid(t, radius, k, width)
    refined = np.sort(np.r_[t, (t[1:] + t[:-1]) / 2])
    fine_clock, fine_probabilities, fine_energy = clock_grid(refined, radius, k, width)
    require(np.max(np.abs(fine_clock - clock)) < 1e-8 and
            np.max(np.abs(fine_probabilities - probabilities)) < 1e-9 and
            np.max(np.abs(fine_energy - energy)) < 1e-9, "decay clock quadrature")
    lo = delta * PARENT / np.hypot(PARENT, k / start_radius)
    hi = delta * PARENT / np.hypot(PARENT, k / stop_radius)
    require(np.all(lo <= fine_clock + 1e-12) and np.all(fine_clock <= hi + 1e-12),
            "monotone endpoint clock bounds")
    rest_probability = -math.expm1(-width * delta)
    require(np.all(fine_probabilities <= rest_probability), "rest-clock ceiling")

    charged = CHARGED_BRANCH * (population @ fine_probabilities)
    charged_energy = CHARGED_BRANCH * (population @ fine_energy)
    rest_ceiling = CHARGED_BRANCH * rest_probability
    energy_fractions = [(1 + (SCALAR / PARENT)**2) / 2, (1 - (SCALAR / PARENT)**2) / 2]
    energies = [charged_energy * fraction for fraction in energy_fractions]
    require(np.isclose(sum(energies), charged_energy), "two-body energy partition")
    require(np.all((np.hypot(PARENT, k / stop_radius) * fine_probabilities <= fine_energy + 4e-9)
                   & (fine_energy <= np.hypot(PARENT, k / start_radius) * fine_probabilities + 4e-9)),
            "birth energy is bracketed by endpoint parent energies")
    recoil = (PARENT**2 - SCALAR**2) / (2 * PARENT)
    return dict(start_time=start_time, stop_time=stop_time, start_radius=start_radius, stop_radius=stop_radius,
                initial_individual_N2=total, all_channel_decay_probability_at_rest=rest_probability,
                charged_branch_events_per_prepared_N2_upper=rest_ceiling,
                shell_0_and_32_proper_times=fine_clock[[0, -1]].tolist(),
                shell_0_and_32_all_channel_probabilities=fine_probabilities[[0, -1]].tolist(),
                mean_prompt_charged_scalar_births=charged, mean_prompt_charged_lepton_births=charged,
                mean_charged_scalar_births_per_N2=charged / total,
                mean_charged_scalar_birth_energy_sum=energies[0],
                mean_charged_lepton_birth_energy_sum=energies[1],
                mean_scalar_readout_energy_over_rest_lower=np.hypot(1., start_radius / stop_radius * recoil / SCALAR),
                max_secondary_two_body_probability=-math.expm1(-delta / (1024 * math.pi**5)),
                S0_upper_if_all_seed_carriers_are_these_scalar_daughters=charged / SEED_PER_ACTION,
                S0_upper_even_for_filled_N2_band=degeneracy.sum() * rest_ceiling / SEED_PER_ACTION)


def all_shells(low, bands, tail_number, radius, width, minimum_rest_count):
    """Bracket every higher-shell contribution using opposite ends of each number band."""
    first, last = low["start_time"], low["stop_time"]
    t = np.r_[first, radius.x[(first < radius.x) & (radius.x < last)], last]
    refined = np.sort(np.r_[t, (t[1:] + t[:-1]) / 2])
    labels = np.array([33, 128, 129, 512, 513, 2048, 2049])
    _, p, energy = clock_grid(refined, radius, labels + 1.5, width)
    _, p_coarse, energy_coarse = clock_grid(t, radius, labels + 1.5, width)
    require(np.max(np.abs(p - p_coarse)) < 1e-9 and np.max(np.abs(energy - energy_coarse)) < 1e-9,
            "higher-shell clock quadrature")
    p *= CHARGED_BRANCH
    energy *= CHARGED_BRANCH
    require(np.all(np.diff(p) < 0) and np.all(np.diff(energy) > 0), "shell-ordered number and energy")

    charged_low = low["mean_prompt_charged_scalar_births"]
    e_low = low["mean_charged_scalar_birth_energy_sum"] + low["mean_charged_lepton_birth_energy_sum"]
    min_count = max_count = charged_low
    min_energy = max_energy = e_low
    parts = {}
    for i in (0, 2, 4):
        key = f"{labels[i]}_{labels[i+1]}"
        number = bands[key]["number"]
        parts[key] = [number * p[i+1], number * p[i]]
        min_count += number * p[i+1]
        max_count += number * p[i]
        min_energy += number * energy[i]
        max_energy += number * energy[i+1]
    tail_births = tail_number * p[-1]
    max_count += tail_births
    max_energy += tail_number * CHARGED_BRANCH * width * PARENT * (last - first)
    scalar_fraction = (1 + (SCALAR / PARENT)**2) / 2
    return dict(charged_scalar_births_by_band_brackets=parts,
                charged_scalar_births_tail_from_2049_upper=tail_births,
                full_tower_charged_scalar_births_bracket=[min_count, max_count],
                full_tower_charged_lepton_births_bracket=[min_count, max_count],
                full_tower_charged_scalar_birth_energy_bracket=[min_energy * scalar_fraction, max_energy * scalar_fraction],
                full_tower_charged_lepton_birth_energy_bracket=[min_energy * (1-scalar_fraction), max_energy * (1-scalar_fraction)],
                S0_upper_for_massive_scalar_seed_count=max_count / SEED_PER_ACTION,
                minimum_required_over_maximum_supplied_at_late_perturbativity_threshold=minimum_rest_count / max_count)


def calculate():
    modes = json.loads((BASE / "MODES.json").read_text())
    bands = json.loads((BASE / "BAND.json").read_text())["events"]
    tails = json.loads((BASE / "TAIL.json").read_text())
    radius, clock_events = make_geometry()
    # Spin-averaged trace for one sign and one weak component: |h|^2 p_N . p_l.
    momentum = (PARENT**2 - SCALAR**2) / (2 * PARENT)
    matrix_squared = YUKAWA**2 * PARENT * momentum
    partial = momentum * matrix_squared / (8 * math.pi * PARENT**2)
    width = 4 * partial
    closed_form = YUKAWA**2 * PARENT * (1 - (SCALAR / PARENT)**2)**2 / (8 * math.pi)
    require(math.isclose(width, closed_form), "four equally weighted two-body channels")
    events = modes["events"]
    q = math.sqrt(3) * SCALAR * 1e8
    minimum_rest_count = SEED_PER_ACTION * q**2 / (4 * math.pi * 13200.04)
    starts = [("bounce", 10.5017, float(tails["bounce"]["4"]["asymptotic_number_upper"])),
              ("expanding_slice", clock_events["turnaround"],
               float(tails["expanding_slice"]["4"]["number_upper_from_shell"]["2049"]))]
    trials = []
    for event, stop, tail in starts:
        start = dict(events[event], t=clock_events[event])
        low = trial(start, stop, float(radius(stop)), radius, radius.x, width)
        low["all_shells"] = all_shells(low, bands[event]["species"]["4"]["bands"], tail,
                                        radius, width, minimum_rest_count)
        trials.append(low)
    stage_one, stage_two = trials
    return dict(inputs=dict(mL=SCALAR, M2L=PARENT, abs_g=1., Q=1e8, K_initial=.001, K_late=13200.04,
                            initial_radius=2., initial_total_cold_rest_density=1e-4, individual_shells_through=32,
                            grouped_shells_through=2048, tail_from_shell=2049),
                partial_width_L_per_sign_and_weak_component=partial, total_N2_width_L=width,
                secondary_phi_two_body_width_L_leading=1 / (1024 * math.pi**5),
                parent_rest_frame=dict(momentum=momentum, charged_scalar_energy=np.hypot(momentum, SCALAR),
                                       charged_lepton_energy=momentum, scalar_speed=15 / 17,
                                       scalar_isotropic_pressure_over_energy=75 / 289),
                two_independent_trials=dict(bounce_to_reference_onset_time=stage_one,
                                            expanding_slice_to_turnaround=stage_two),
                massive_scalar_seed_match=dict(required_individual_charged_scalars_per_S0=SEED_PER_ACTION,
                                               dimensionless_single_field_charge=q,
                                               S0_lower_for_local_alpha_below_one_initial=q**2 / (4 * math.pi * .001),
                                               S0_lower_for_local_alpha_below_one_late=q**2 / (4 * math.pi * 13200.04)))


def rounded_values(value):
    if isinstance(value, dict):
        return {key: rounded_values(item) for key, item in value.items()}
    if isinstance(value, list):
        return [rounded_values(item) for item in value]
    return float(f"{value:.8g}") if isinstance(value, float) else value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = json.dumps(rounded_values(calculate()), indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(report)
    else:
        print(report, end="")


if __name__ == "__main__":
    main()
