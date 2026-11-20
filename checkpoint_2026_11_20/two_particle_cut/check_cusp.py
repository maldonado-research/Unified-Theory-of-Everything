"""Threshold expansion of the complete J=0 Hermitian cut and its local continuation."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from mpmath import mp

HERE = Path(__file__).resolve().parent
BLOCKS = {"E": ((100, 10, 1), (400, -40, 2)),
          "G": ((1, 1, 1), (16, -1, .25))}


def require(actual, expected, relative, label):
    if not abs(actual-expected) <= relative * abs(expected):
        raise ArithmeticError(f"{label}: {actual} vs {expected}")


def vertex_jets(b):
    """Threshold Taylor coefficients of I(q,b) from the sterile-invariant cut."""
    def integrand(t, n):
        if t in (0, 1):
            return mp.mpf(0)
        u = 1/t**2
        start = b - 1 + (u-1)
        span = (u-1)**2/u
        fraction = span/start
        if fraction < mp.mpf(".001"):
            angular = mp.hyp2f1(1, 2, 3, fraction)/(2*start)
        else:
            # start - span = b - 1/u remains well resolved at large u.
            angular = (mp.log(start/(b-1/u))/fraction - 1)/span
        density = ((u-1)/u)**2 * angular
        return 2*density*(-1)**n/(t*(1+t*t)*(u+1)**n)
    return [mp.quad(lambda t: integrand(t, n), [0, mp.mpf(".5"), 1])
            for n in (0, 1, 2)]


def residue(b):
    lepton = mp.quad(lambda x: (1-x)*mp.log1p((b-1)*x), [0, 1])
    scalar = (mp.log(b) - 1 - b - (b**2-1)*mp.log1p(-mp.mpf(1)/b)
              if b != 1 else -2)
    return lepton + 2*scalar


def means_at_threshold(flavor):
    block = BLOCKS[flavor]
    form = {b: vertex_jets(b) for b, _, _ in block}
    if flavor == "G":
        require(form[1][0], mp.pi**2/8-mp.log(2), mp.mpf("1e-32"), "light triangle")
    def tree(q):
        return sum(w/(b+q) for b,w,_ in block)
    def proper(q):
        h = sum(a/(b+q) for b,_,a in block)
        triangle = sum(w*mp.polyval(form[b][::-1],q-1) for b,w,_ in block)
        bubble = 2+1/q-((1+q)/q)**2*mp.log1p(q)
        return -2*q*h*(tree(q)*bubble+triangle)/(16*mp.pi**2)
    def expand(f):
        f0, f1, f2 = mp.taylor(f, 1, 2)
        return (f0, 2*f1+mp.mpf(4)*f2/3)
    legs = sum(a*residue(b) for b,_,a in block)/(16*mp.pi**2)
    return expand(tree), expand(proper), legs


def calculate():
    el, eproper, efields = means_at_threshold("E")
    mu, mproper, mfields = means_at_threshold("G")
    def sew(a, b):
        first = a[0]*b[0]/(2*mp.pi)
        third = first+(a[0]*b[1]+a[1]*b[0])/(2*mp.pi)
        return mp.matrix([first, third])
    lo = sew(el, mu)
    two_e = sew(eproper, mu) + efields*lo
    two_g = sew(el, mproper) + mfields*lo
    three_g = mp.matrix([0, el[0]*mu[0]/(12*mp.pi**3)])
    return {"LO": lo, "E": two_e, "G": two_g+three_g}, two_g, three_g


def load(path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check_open_cut(coeffs):
    sewing = load(HERE / "check_absorptive.py")
    betas = np.array([.08, .06, .04, .02])
    samples = []
    leading_j2 = []
    for beta in betas:
        projected = sewing.two_and_three(4/(1-beta**2), 48, 2)
        samples.append(np.r_[projected["leading"][0], projected["total"][:,0]]/beta)
        leading_j2.append(projected["leading"][1]/beta**5)
    fit = np.polynomial.polynomial.polyfit(betas**2, samples, 3)
    for j, name in enumerate(("LO", "E", "G")):
        require(fit[0,j], float(coeffs[name][0]), 1.e-9, "complete sewn cut, linear")
        require(fit[1,j], float(coeffs[name][1]), 1.e-7, "complete sewn cut, cubic")
    def second(flavor):
        return sum(2*w/(1+b)**3 for b,w,_ in BLOCKS[flavor])
    expected_j2 = (4/15)**2*second("E")*second("G")/(2*np.pi)
    extrapolated = np.polynomial.polynomial.polyfit(betas**2, leading_j2, 3)[0]
    require(extrapolated, expected_j2, 2.e-7, "leading J2 threshold power and coefficient")
    return float(expected_j2)


def check_local_dispersion(coeffs):
    # The measure ds'/(s'-s) differs from dw'/(w'-w), w=beta^2,
    # by a kernel independent of s, which belongs to the analytic part.
    for cutoff in (mp.mpf(".2"), mp.mpf(".4")):
        for eta in (mp.mpf(".01"), mp.mpf(".04")):
            for name, (linear, cubic) in coeffs.items():
                def integrand(t):
                    return 2*t**2*(linear+cubic*t**2)/(t**2+eta**2)
                integral = mp.quad(integrand, [0, eta, cutoff])/mp.pi
                even = 2/mp.pi*(
                    linear*(cutoff+eta*mp.atan(eta/cutoff)) +
                    cubic*(cutoff**3/3-cutoff*eta**2-eta**3*mp.atan(eta/cutoff)))
                require(integral-even, -linear*eta+cubic*eta**3,
                        mp.mpf("1e-36"), "local dispersion, "+name)


def check_real_leading_box(coefficients):
    # Real physical box below threshold; its analytic background is fitted freely.
    path = HERE.parents[1] / "checkpoint_2026_11_19/uv_completion/check_bounds.py"
    with mp.workdps(44):
        box = load(path)
        step = mp.mpf(".0005")
        samples = []
        for i in range(8):
            eta = step*i
            s = 4/(1+eta**2)
            samples.append(sum(w*box.reduced_amplitude(s, c, 24)
                               for c,w in box.unit_rule(6)))
        design = mp.matrix([[mp.mpf(i)**p for p in range(8)] for i in range(8)])
        fit = mp.lu_solve(design, mp.matrix(samples))
        require(fit[1]/step, -coefficients[0], mp.mpf("3e-13"), "real LO box, linear cusp")
        require(fit[3]/step**3, coefficients[1], mp.mpf("2e-10"), "real LO box, cubic cusp")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    mp.dps = 44
    coeffs, two_g, three_g = calculate()
    j2 = check_open_cut(coeffs)
    check_local_dispersion(coeffs)
    check_real_leading_box(coeffs["LO"])
    out = {k: [mp.nstr(x, 16) for x in v] for k,v in coeffs.items()}
    record = {"coordinate": "beta=sqrt(1-4/s), s in fixed scalar pole mass units",
              "coefficient_order": ["kappa_1","kappa_3"],
              "complete_cut_over_gamma": out,
              "muon_two_scalar": [mp.nstr(x,16) for x in two_g],
              "muon_three_particle": [mp.nstr(x,16) for x in three_g],
              "below_threshold_nonanalytic_terms": "-kappa_1*eta + kappa_3*eta^3, eta=sqrt(4/s-1)",
              "leading_J2_over_beta5": float(f"{j2:.13g}"),
              "checks": {"independent_unexpanded_cut_beta": [.08,.06,.04,.02],
                         "relative_tolerances_kappa1_kappa3": [1e-9,1e-7],
                         "below_threshold_real_LO_box_eta": "0, 0.0005, ..., 0.0035",
                         "real_LO_box_relative_tolerances_linear_cubic": [3e-13,2e-10],
                         "local_dispersion_eta": [.01,.04], "local_cutoff_beta": [.2,.4],
                         "light_triangle_at_threshold": "pi^2/8 - ln(2)"}}
    data = json.dumps(record, indent=2, allow_nan=False)+"\n"
    if args.output:
        args.output.write_text(data, encoding="utf-8")
    else:
        print(data, end="")


if __name__ == "__main__":
    main()
