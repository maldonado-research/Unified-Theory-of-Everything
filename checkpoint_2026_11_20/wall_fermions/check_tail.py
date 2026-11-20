"""Enclose the omitted Majorana occupations on the stationary Hermite wall."""
import argparse
import importlib.util
import json
from bisect import bisect_left
from dataclasses import dataclass
from fractions import Fraction as Q
from math import ceil, isqrt
from pathlib import Path

spec = importlib.util.spec_from_file_location("wall_modes", Path(__file__).with_name("check_modes.py"))
wall_modes = importlib.util.module_from_spec(spec)
spec.loader.exec_module(wall_modes)

# Positive contributions are rounded UP to this rational grid before adding.
# This prevents rational denominators at different knots from growing together.
SCALE = 10**12


def upward(x):
    return Q(ceil(x * SCALE), SCALE)


def sqrt_upper(x):
    return Q(isqrt((x.numerator * SCALE**2) // x.denominator) + 1, SCALE)


def report(x):
    integer = ceil(x*10**8)
    return f"{integer//10**8}.{integer%10**8:08d}"


@dataclass
class Block:
    sweep: Q = Q(0)
    acceleration: Q = Q(0)
    a_vel2: Q = Q(0)
    vel2: Q = Q(0)
    jerk_and_jumps: Q = Q(0)
    jerk_correction: Q = Q(0)
    jerk_quintic: Q = Q(0)
    vel_acc: Q = Q(0)
    a_vel3: Q = Q(0)
    boundary_acc: Q = Q(0)
    boundary_cubic: Q = Q(0)


def geometry_bounds(nodes, radii, slopes):
    """Bernstein hulls enclose every point of each exact cubic Hermite cell."""
    t, a, v = ([Q(float(value)) for value in array] for array in (nodes, radii, slopes))
    if v[0] != 0 or v[-1] != 0:
        raise ValueError("Expected two stationary endpoint radii")
    # Split off the short early acceleration pulse. The two blocks can
    # independently use the first or second integration-by-parts inequality.
    split = next(i for i, time in enumerate(t) if time >= Q(float(.01)))
    blocks = [Block(), Block()]
    last_acc = None
    split_jump = Q(0)
    for i, (tl, tr) in enumerate(zip(t, t[1:])):
        block = blocks[int(i >= split)]
        al, ar, vl, vr = a[i], a[i+1], v[i], v[i+1]
        dt = tr - tl
        if dt <= 0:
            raise ValueError("Non-increasing Hermite knots")
        controls = (al, al + dt*vl/3, ar - dt*vr/3, ar)
        amin, amax = min(controls), max(controls)
        if amin <= 0:
            raise ValueError("Positive-radius enclosure failed")
        vmid = 3*(ar-al)/dt - vl - vr
        vmax = max(abs(vl), abs(vmid), abs(vr))
        wl, wr = 2*(vmid-vl)/dt, 2*(vr-vmid)/dt
        wmax = max(abs(wl), abs(wr))
        umax = abs(wr-wl)/dt
        vel = amax*vmax
        acc = amax*vmax**2 + amax**2*wmax
        jerk = amax*(vmax**3 + 4*amax*vmax*wmax + amax**2*umax)
        deta = dt/amin
        block.sweep += upward(deta*vel)
        block.acceleration += upward(deta*acc)
        block.a_vel2 += upward(deta*amax*vel**2)
        block.vel2 += upward(deta*vel**2)
        block.jerk_and_jumps += upward(deta*jerk)
        block.jerk_correction += upward(deta*(10*amax*vel*acc + 3*vel**3))
        block.jerk_quintic += upward(deta*18*amax**2*vel**3)
        block.vel_acc += upward(deta*vel*acc)
        block.a_vel3 += upward(deta*amax*vel**3)

        left_acc = al*vl**2 + al**2*wl
        right_acc = ar*vr**2 + ar**2*wr
        if last_acc is not None:
            jump = upward(abs(left_acc-last_acc))
            if i == split:
                split_jump = jump
            else:
                block.jerk_and_jumps += jump
        if i in (0, split):
            block.boundary_acc += upward(abs(left_acc))
            block.boundary_cubic += upward(3*al*(al*vl)**2)
        if i+1 in (split, len(t)-1):
            block.boundary_acc += upward(abs(right_acc))
            block.boundary_cubic += upward(3*ar*(ar*vr)**2)
        last_acc = right_acc
        if i == 0:
            initial_acc = left_acc

    full = Block(**{key: sum(getattr(b, key) for b in blocks) for key in vars(blocks[0])})
    full.jerk_and_jumps += split_jump
    full.boundary_acc = upward(abs(initial_acc)) + upward(abs(last_acc))
    full.boundary_cubic = Q(0)
    return full, blocks, a[-1], t[split]


def scaled_bounds(block, m, k, beta):
    """Upper bounds after multiplication by k^3, at a given k."""
    common = beta*m**2*block.vel2/8
    first = common + m*k/4*(block.acceleration + 3*m**2*block.a_vel2/k**2)
    second = (common + m/8*(block.jerk_and_jumps + block.boundary_acc
                            + m**2*(block.jerk_correction + block.boundary_cubic)/k**2
                            + m**4*block.jerk_quintic/k**4)
              + beta*m**2/(16*k)*(block.vel_acc + 3*m**2*block.a_vel3/k**2))
    return first, second


def tails(full, blocks, radius, mass):
    m = Q(mass)
    def estimates(k):
        beta = min(Q(1), m*full.sweep/(2*k))
        whole = scaled_bounds(full, m, k, beta)
        pieces = sum(min(scaled_bounds(part, m, k, beta)) for part in blocks)
        return min(*whole, pieces)/k**3, whole[1]

    # Below this cut use a rational bound individually for each omitted shell.
    # Above it the whole-interval second bound on k^3 |beta| is decreasing.
    first_asymptotic_shell = 2049
    num = energy = Q(0)
    for n in range(33, first_asymptotic_shell):
        k = Q(2*n+3, 2)
        bound, _ = estimates(k)
        occupation = 2*(n+1)*(n+2) * bound**2
        num += upward(occupation)
        energy += upward(occupation * sqrt_upper(m**2 + (k/radius)**2))
    k = Q(2*first_asymptotic_shell+3, 2)
    _, coeff = estimates(k)
    far_num = 2*coeff**2/(3*(k-1)**3)
    far_energy = coeff**2*sqrt_upper(1/radius**2 + (m/k)**2)/(k-1)**2
    return dict(number_upper=report(num+far_num), energy_L_upper=report(energy+far_energy),
                first_asymptotic_shell=first_asymptotic_shell,
                asymptotic_amplitude_coefficient_upper=report(coeff),
                asymptotic_number_upper=report(far_num),
                asymptotic_energy_L_upper=report(far_energy))


def expanding_readout(nodes, radii, slopes, time, full):
    """Exact partial-cell value and derivatives from the Hermite node data."""
    i = bisect_left(nodes, time)-1
    t0, t1, al, ar, vl, vr, t = map(
        lambda v: Q(float(v)), (nodes[i], nodes[i+1], radii[i], radii[i+1],
                                slopes[i], slopes[i+1], time))
    h, x = t1-t0, t-t0
    secant = (ar-al)/h
    quadratic = (3*secant-2*vl-vr)/h
    cubic = (vl+vr-2*secant)/h**2
    a = al+x*(vl+x*(quadratic+x*cubic))
    v = vl+x*(2*quadratic+3*x*cubic)
    w = 2*quadratic+6*x*cubic
    da, dda = a*v, a*v**2+a**2*w
    if da <= 0:
        raise ValueError("Expected an expanding slice")
    result = {}
    for mass in (1, 4):
        m = Q(mass)
        d = m*abs(da)/4
        bounds = {}
        for first_shell in (33, 2049):
            k = Q(2*first_shell+3, 2)
            b = min(Q(1), m*full.sweep/(2*k))
            _, coeff = scaled_bounds(full, m, k, b)
            # Using the full history bounds all interior total variations
            # on a prefix. The new right endpoint contributes its own |G|.
            coeff += m/8*(abs(dda) + 3*m**2*a*da**2/k**2)
            z = k-1
            upper = 2*(d**2/z + d*coeff/z**2 + coeff**2/(3*z**3))
            bounds[str(first_shell)] = report(upper)
        result[str(mass)] = dict(number_upper_from_shell=bounds,
                                 energy_log_coefficient_upper=report(m**2*da**2/(8*a)))
    return result


def verify(trace=wall_modes.TRACE):
    times, radii, slopes, events = wall_modes.make_geometry_nodes(trace)
    full, blocks, radius, split = geometry_bounds(times, radii, slopes)
    result = {str(m): tails(full, blocks, radius, m) for m in (1, 4)}
    i = bisect_left(times, events["bounce"])
    if times[i] != events["bounce"]:
        raise ValueError("The bounce must be an exact knot")
    prefix, early, bounce_radius, _ = geometry_bounds(times[:i+1], radii[:i+1], slopes[:i+1])
    bounce = {str(m): tails(prefix, early, bounce_radius, m) for m in (1, 4)}
    expanding = expanding_readout(times, radii, slopes, events["expanding_slice"], full)
    return dict(geometry_float64_le_sha256=wall_modes.GEOMETRY_SHA256, cells=len(times)-1,
                endpoint="turnaround", first_omitted_shell=33,
                multiplicity="2*(n+1)*(n+2)", split_time=float(split), species=result,
                bounce=bounce, expanding_slice=expanding)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trace", type=Path, default=wall_modes.TRACE)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = json.dumps(verify(args.trace), indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(data, encoding="utf8")
    else:
        print(data, end="")
