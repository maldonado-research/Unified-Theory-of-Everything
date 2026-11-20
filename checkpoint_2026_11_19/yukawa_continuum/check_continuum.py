"""Integrate a pole-subtracted scalar self-energy directly in the box simplex."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

HERE = Path(__file__).resolve().parent
LOOP = 16*np.pi**2
spec = importlib.util.spec_from_file_location("leading_box", HERE.parent / "radiative_stability/check_detuning.py")
leading = importlib.util.module_from_spec(spec)
spec.loader.exec_module(leading)


def rule(order):
    nodes, weights = roots_legendre(order)
    return (1+nodes)/2, weights/2


def smoothed_panel(low, high, nodes, weights):
    length = high-low
    angle = np.pi*nodes
    return low+length*np.sin(angle/2)**2, weights*length*np.pi*np.sin(angle)/2


def continuum_transform(p, threshold):
    """Stieltjes transform of (U-B)^2/[U(U-1)^2], with B>1 and p<B."""
    k = 1/threshold
    z = p/threshold
    output = np.empty_like(p)
    local = z >= -1
    # All physical positive arguments here obey p/B < 0.3.
    # The compact v integral avoids coalescing poles when p is near 1.
    nodes, weights = rule(18)
    output[local] = sum(w*((1-v)/(1-k*v))**2/(1-z[local]*v)
                        for v,w in zip(nodes, weights))
    z = z[~local]
    # Partial fractions at large spacelike arguments preserve the log tail.
    major = ((z-1)/(z-k))**2
    rest = (1-k)*(2*z-1-k)/(z-k)**2
    double = (1-k)**2/(k*(k-z))
    log_k = -np.log1p(-k)/k
    output[~local] = (major*(-np.log1p(-z)/z) + rest*log_k
                      + double*(1/(1-k)-log_k))
    return output/threshold


def spectral_grids(s, b, order):
    """Integrate out the variable scalar Feynman parameter and continuum mass."""
    nodes, weights = rule(order)
    beta = np.sqrt(1-4/s)
    edges = [0., (1-beta)/2, (1+beta)/2, 1.]
    panels = []
    for low, high in zip(edges, edges[1:]):
        xs, wxs = smoothed_panel(low, high, nodes, weights)
        for x, wx in zip(xs, wxs):
            length = 1-x
            # p(y)=1 at root. The factorized numerator avoids subtracting
            # nearly equal numbers where the three-body cut opens or closes.
            root = -s*(x-edges[1])*(x-edges[2])/(s*x+b-1)
            split = [0., length]
            if 0 < root < length:
                split.insert(1, root)
            for lo, hi in zip(split, split[1:]):
                y, wy = smoothed_panel(lo, hi, nodes, weights)
                remaining = length-y
                p = s*x - (x+b*y)/remaining
                near_zero = abs(p)<0.25
                one = np.zeros_like(p, dtype=complex)
                one[near_zero] = 1.
                mask = near_zero & (p != 0)
                one[mask] = -np.log1p(-p[mask])/p[mask]
                rest = ~near_zero
                one[rest] = (-np.log((s*x+b-1)*abs(y[rest]-root)/remaining[rest])
                             + 1j*np.pi*(y[rest]<root))/p[rest]
                thresholds = [continuum_transform(p, B) for B in (16,100,400)]
                sources = (one+thresholds[0]/4, thresholds[1]+2*thresholds[2])
                panels.append((x+b*y,y*remaining,remaining,wx*wy*sources[0],wx*wy*sources[1]))
    return tuple(np.concatenate(column) for column in zip(*panels))


def amplitudes(s, cosines, order):
    answers = np.zeros((len(cosines),2), dtype=complex)
    for b, mu_weight in ((1,1),(16,-1)):
        constant, cross, remaining, g, y = spectral_grids(s,b,order)
        for i,z in enumerate(cosines):
            for t in (-s*(1-z)/2, -s*(1+z)/2):
                common = constant-t*cross
                # Combine the two electron residues before quadrature.
                kernel = -30*common / ((common+100*remaining)*(common+400*remaining))
                answers[i] += mu_weight*np.array([np.sum(kernel*g), np.sum(kernel*y)])
    return 2*s*answers/LOOP**2


def energy_row(s, order, angle_order):
    z,w = rule(angle_order)
    corrections = amplitudes(s, np.append(z,0.), order)
    mean = w @ corrections[:-1]
    base, flux = leading.angular_modes(s, leading.dispersion_box, 36, angle_order)
    rate = np.sum(flux*abs(base[:,0])**2)
    slopes = 2*np.real((flux*base[:,0].conj()) @ corrections[:-1])
    return {"delta_g_central":corrections[-1,0], "delta_y_central":corrections[-1,1],
            "delta_g_mean":mean[0], "delta_y_mean":mean[1], "m2_sigma0":rate,
            "m2_sigma_slope_g":slopes[0], "m2_sigma_slope_y":slopes[1]}


def format_value(value):
    if np.iscomplexobj(value):
        return [float(f"{v:.11g}") for v in (value.real,value.imag)]
    return float(f"{value:.11g}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    reference = json.loads((HERE / "QCDLOOP_CONTINUUM.json").read_text())
    fine_reference = {row["s_over_m2"]:row for row in reference["rows"] if row["mass_order"]==28}
    rows = {}
    max_error = 0.
    for s in (6,8,10):
        coarse = energy_row(s, 80, 10)
        fine = energy_row(s, 150, 14)
        if any(abs(coarse[k]-v) > 5e-7*abs(v) for k,v in fine.items()):
            raise ArithmeticError(f"simplex refinement at {s}")
        if fine["delta_y_central"].imag != 0 or fine["delta_y_mean"].imag != 0:
            raise ArithmeticError("closed electron continuum")
        for key,value in fine.items():
            expected = fine_reference[s][key]
            if isinstance(expected,list):
                expected = complex(*expected)
            max_error = max(max_error,abs(value-expected)/abs(expected))
        rows[str(s)] = {k:format_value(v) for k,v in fine.items()}
    if max_error > 3e-8:
        raise ArithmeticError("external varying-mass scalar boxes")
    report = {"s_over_m2": rows, "QCDLoop_relative_difference_below":3e-8, "normalization":"m=1; amplitudes divided by gamma and rates by |gamma|^2",
              "relative_endpoint_ratio_slopes": {
                  block:format_value(rows["10"]["m2_sigma_slope_"+block]/rows["10"]["m2_sigma0"]
                                     -rows["6"]["m2_sigma_slope_"+block]/rows["6"]["m2_sigma0"])
                  for block in ("g","y")}}
    content = json.dumps(report,indent=2,allow_nan=False)+"\n"
    if args.output:
        args.output.write_text(content)
    else:
        print(content,end="")


if __name__ == "__main__":
    main()
