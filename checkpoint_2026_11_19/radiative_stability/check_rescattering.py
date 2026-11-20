"""Linear scalar-quartic response of the four-singlet scattering benchmark."""
import argparse
import importlib.util
import json
from itertools import permutations
from pathlib import Path

from mpmath import mp

HERE = Path(__file__).resolve().parent
PREVIOUS = HERE.parent / "uv_completion"


def load_calculation(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


vertex = load_calculation("upper_rim_triangle", HERE.parent / "energy_scan" / "verify_lower_operator.py")
leading = load_calculation("upper_rim_box", PREVIOUS / "check_bounds.py")


def require_close(actual, reference, tolerance, label):
    if abs(actual-reference) > mp.mpf(tolerance) * max(mp.mpf("1e-28"), abs(reference)):
        raise ArithmeticError(label)


def encoded(value):
    return mp.nstr(value, 16)


def complex_pair(value):
    return {"real": encoded(value.real), "imag": encoded(value.imag)}


def check_triangles():
    external = json.loads((HERE / "QCDLOOP_TRIANGLES.json").read_text())
    if {(row["s"], row["b"]) for row in external["rows"]} != {
        (s,b) for s in (-1,6,8,10) for b in (1,16,100,400)
    }:
        raise ArithmeticError("incomplete external triangle set")
    triangles = {}
    largest = mp.mpf(0)
    for row in external["rows"]:
        s, b = mp.mpf(row["s"]), mp.mpf(row["b"])
        value = vertex.triangle(s, mp.mpf(1), b)
        other = vertex.dispersion(s, mp.mpf(1), b)
        require_close(value, other, "1e-30", "triangle simplex vs spectral integral")
        reference = mp.mpc(str(row["real"]), str(row["imag"]))
        require_close(value, reference, "3e-14", "QCDLoop normalization and ordering")
        largest = max(largest, abs(value-reference)/abs(value))
        triangles[int(s),int(b)] = value
    return triangles, largest


def scalar_contractions():
    # Left endpoints muon, muon, quartic, quartic; right endpoints
    # electron, electron, quartic, quartic, all in the charged component.
    counts = [0, 0, 0]
    for assignment in permutations(range(4)):
        self_pairs = sum(endpoint >= 2 for endpoint in assignment[2:])
        counts[self_pairs] += 1
    if counts != [4, 16, 4]:
        raise ArithmeticError("quartic scalar Wick contractions")
    return counts


def rate_row(energy, triangles, direct, previous):
    s = mp.mpf(energy)
    mu = triangles[energy,1] - triangles[energy,16]
    electron = 10*triangles[energy,100] - 40*triangles[energy,400]
    loop = 16*mp.pi**2
    amplitude = -4*s*mu*electron / loop**2
    beta = mp.sqrt(1-4/s)
    n, h = s/2-1, s*beta/2

    def tree(q, masses, weights):
        return sum(w/(b+q) for b,w in zip(masses, weights))

    def two_tree_cut(x):
        transfers = (n-h*x, n+h*x)
        mu_tree = sum(tree(q, (1,16), (1,-1)) for q in transfers)
        el_tree = sum(tree(q, (100,400), (10,-40)) for q in transfers)
        # One corrected scalar-pair subamplitude on either side. Direction
        # independence of each triangle permits independent angle averages.
        return (-4/loop)*(electron*mu_tree + mu.conjugate()*el_tree)

    sewn = s*beta/(32*mp.pi)*mp.quad(two_tree_cut, [0,1])
    require_close(sewn, amplitude.imag, "1e-35", "unitarity with upper/lower rims")
    samples = [(w, leading.reduced_amplitude(s, x, 32)) for x,w in leading.unit_rule(14)]
    mean = mp.fsum(w*value for w,value in samples)
    rate = mp.fsum(w*abs(value)**2 for w,value in samples)/(32*mp.pi*s)
    if direct is not None:
        require_close(mean, mp.mpc(str(direct["mean_re"]), str(direct["mean_im"])),
                      "6e-9", "direct-box mean")
        require_close(rate, mp.mpf(previous["scaled_full_rate"]), "1e-14", "previous full rate")
    slope = 2*mp.re(mean.conjugate()*amplitude)/(32*mp.pi*s)
    return {"s_over_m2":energy, "muon_triangle_sum":complex_pair(mu),
            "electron_triangle_sum":complex_pair(electron),
            "leading_angle_mean":complex_pair(mean),
            "two_loop_amplitude_per_lambda":complex_pair(amplitude),
            "scaled_rate_at_leading_order":encoded(rate),
            "scaled_rate_derivative_at_lambda_zero":encoded(slope),
            "fractional_rate_derivative_at_lambda_zero":encoded(slope/rate)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    mp.dps = 47
    triangles, delta = check_triangles()
    direct = json.loads((PREVIOUS / "BOX_DIRECT.json").read_text())["rows"]
    previous = json.loads((PREVIOUS / "BOUNDS.json").read_text())["upper_rim"]
    direct_at = {row["s_over_m2"]:row["refinement"][-1] for row in direct}
    rows = [rate_row(s, triangles, direct_at.get(s), previous.get(str(s))) for s in (6,8,10)]
    r6,r8,r10 = rows
    ratio = mp.mpf(r10["scaled_rate_at_leading_order"])/mp.mpf(r6["scaled_rate_at_leading_order"])
    log_slope = (mp.mpf(r10["fractional_rate_derivative_at_lambda_zero"])
                 -mp.mpf(r6["fractional_rate_derivative_at_lambda_zero"]))
    detuned = json.loads((HERE / "DETUNING.json").read_text())["s_over_m2"]
    d = {}
    ell = {}
    for row in rows:
        s = row["s_over_m2"]
        coefficients = detuned[str(s)]
        require_close(mp.mpf(str(coefficients["a"])),
                      mp.mpf(row["scaled_rate_at_leading_order"]), "2e-9", "detuning baseline")
        d[s] = 2*mp.mpf(str(coefficients["b_real"]))/mp.mpf(str(coefficients["a"]))
        ell[s] = mp.mpf(row["fractional_rate_derivative_at_lambda_zero"])
    weight = (d[8]-d[6])/(d[10]-d[6])
    contrast = ell[8]-ell[6]-weight*(ell[10]-ell[6])
    if abs(contrast) < mp.mpf("1e-5"):
        raise ArithmeticError("rank of the three-energy benchmark response")
    def average(s, shift):
        data = detuned[str(s)]
        return (mp.mpf(str(data["a"])) + 2*mp.mpf(str(data["b_real"]))*shift.real
                + mp.mpf(str(data["c"]))*abs(shift)**2)
    def change(shift):
        values = {s:mp.log(average(s, shift)/average(s, 0)) for s in (6,8,10)}
        return values[8]-values[6]-weight*(values[10]-values[6])
    for direction in (1, mp.j, 1+mp.j):
        # Both real and imaginary nuisance directions vanish to first order.
        if abs(mp.diff(lambda h: change(h*direction), 0)) > mp.mpf("1e-35"):
            raise ArithmeticError("linear complex detuning response")
    report = {"charged_scalar_pairings_by_number_of_quartic_self_contractions":scalar_contractions(),
              "three_energy_reciprocal_log_contrast":{
                  "weight_on_10_over_6":encoded(weight),
                  "derivative_in_lambda":encoded(contrast),
                  "real_detuning_per_lambda_mimicking_endpoint_ratio":encoded(log_slope/(d[10]-d[6])),
                  "actual_contrast_at_real_plus_0p01":encoded(change(mp.mpf("0.01"))),
                  "actual_contrast_at_real_minus_0p01":encoded(change(mp.mpf("-0.01"))),
                  "real_detuning_log_slopes":{str(s):encoded(d[s]) for s in d}},
              "mass_units":"m=1, physical scalar pole mass fixed",
              "coupling_convention":"V contains +lambda (H dagger H)^2",
              "normalization":"amplitudes divided by gamma, rates m^2 sigma / |gamma|^2",
              "max_relative_complex_triangle_difference_from_QCDLoop":encoded(delta),
              "rows":rows, "leading_rate_ratio_10_over_6":encoded(ratio),
              "rate_ratio_derivative_at_lambda_zero":encoded(ratio*log_slope),
              "log_rate_ratio_derivative_at_lambda_zero":encoded(log_slope)}
    data = json.dumps(report,indent=2)+"\n"
    if args.output:
        args.output.write_text(data)
    else:
        print(data,end="")


if __name__ == "__main__":
    main()
