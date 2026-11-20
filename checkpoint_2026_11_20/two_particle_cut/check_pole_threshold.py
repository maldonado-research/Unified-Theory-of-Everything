"""Lightest sterile pole and charged-scalar partial width at the MSbar benchmark."""
import argparse
import json
from pathlib import Path

from mpmath import mp


def require(actual, expected, tolerance, label):
    if not abs(actual-expected) <= tolerance*max(abs(expected), mp.mpf("1e-40")):
        raise ArithmeticError(f"{label}: {actual} vs {expected}")


def bubble(z, scale=1):
    """Real side of the vector bubble in scalar pole mass units; 0 < z <= 1."""
    if z == 1:
        return 1 + 2*mp.log(scale)
    return 2*mp.log(scale) + 2-1/z-(1-z)**2*mp.log1p(-z)/z**2


def mixed_determinant(energy, alpha, scale=1, run_parameters=False):
    """Full 4x4 rest-frame chiral inverse in the two mixed muon columns."""
    row = mp.matrix([[1, -mp.j/2]])
    H = row.H * row
    M = mp.diag([1, 4])
    if run_parameters:
        # One-loop flow in the fixed original basis. Running H affects the next order.
        M += mp.log(scale)*alpha*(H.T*M+M*H)
    factor = alpha*bubble(energy**2, scale)
    K = mp.zeros(4)
    K[:2, :2] = -M.H
    K[2:, 2:] = -M
    K[:2, 2:] = energy*(mp.eye(2)+factor*H)
    K[2:, :2] = energy*(mp.eye(2)+factor*H.T)
    determinant = mp.det(K)
    if abs(determinant.imag) > mp.mpf("1e-42"):
        raise ArithmeticError("complex determinant below threshold")
    return determinant.real


def root_below_scalar(alpha, scale=1, run_parameters=False):
    lower, upper = 1-2*alpha, mp.mpf(1)
    def det(p):
        return mixed_determinant(p, alpha, scale, run_parameters)
    if not det(lower) > 0 > det(upper):
        raise ArithmeticError("light root not bracketed below the scalar")
    for _ in range(105):
        center = (lower+upper)/2
        if det(center)>0:
            lower = center
        else:
            upper = center
    return (lower+upper)/2


def check_retuning(loop):
    shifts = []
    for energy in (6, 8, 10):
        s = mp.mpf(energy)
        beta = mp.sqrt(1-4/s)
        h, v = s/2-1, s*beta/2
        def mean_propagator(b):
            return mp.atanh(v/(h+b))/v
        mean_e = 10*mean_propagator(100)-40*mean_propagator(400)
        analytic = mean_e*(mp.atanh(beta)-beta)/(64*mp.pi**3)
        def leading(mass):
            mu = mp.quad(lambda x: mass/(mass**2+h-v*x) - 1/(16+h-v*x), [-1, 1])/2
            return s*beta*mean_e*mu/(8*mp.pi)
        step = mp.mpf("1e-12")
        shifted = (leading(1+step)-leading(1-step))/(2*step*loop)
        require(shifted, analytic, mp.mpf("1e-21"), "retuned MS mass in LO cut")
        shifts.append({"s_over_m2": energy, "shift_per_G": mp.nstr(analytic,16)})
    tree_e = -mp.mpf(30)/(101*401)
    cubic = tree_e/(192*mp.pi**3)
    beta = mp.mpf("1e-8")
    s = 4/(1-beta**2)
    h, v = s/2-1, s*beta/2
    mean_e = (10*mp.atanh(v/(h+100))-40*mp.atanh(v/(h+400)))/v
    threshold = mean_e*(mp.atanh(beta)-beta)/(64*mp.pi**3)
    require(threshold/beta**3, cubic, mp.mpf("1e-15"), "cubic retuning shift")
    return {"changed_input": "raise only MSbar M1(m) by m G/(16 pi^2); fixed Yukawas",
            "two_loop_reduced_J0_cut_shifts": shifts, "kappa1_muon_shift": "0",
            "kappa3_muon_shift": mp.nstr(cubic,16)}


def verify():
    mp.dps = 48
    loop = 16*mp.pi**2
    for z in map(mp.mpf, (".1", ".76", ".99999", "1")):
        for mu in map(mp.mpf, (".85", "1", "1.3")):
            integral = -2*mp.quad(lambda x: (1-x)*mp.log((1-x)*(1-x*z)/mu**2),
                                  [0, 1])
            require(bubble(z, mu), integral, mp.mpf("1e-42"), "Feynman parameter")
    # The model has a common odd parity and separately conserves each lepton-block parity.
    # Check these signs on every allowed Yukawa monomial and each Majorana mass.
    for phi, le, lmu, ne, nmu in ((-1,1,1,-1,-1),(1,-1,1,-1,1),(1,1,-1,1,-1)):
        if (phi*le*ne, phi*lmu*nmu, ne**2, nmu**2) != (1,1,1,1):
            raise ArithmeticError("discrete symmetry")
    samples = []
    for alpha in map(mp.mpf, ("1e-4", "2.5e-5", "6.25e-6")):
        G = loop*alpha
        mass = root_below_scalar(alpha)
        # Spin trace G*(1-M²) and |p|/(8*pi*m²), m=1.
        momentum = (1-mass**2)/2
        trace_width = momentum*G*(1-mass**2)/(8*mp.pi)
        leading_width = G**3/(1024*mp.pi**5)
        require((1-mass)/alpha, 1, 5*alpha, "mass root")
        require(trace_width/leading_width, 1, 12*alpha, "on-pole phase space")
        # Vary the scale and run the full mass matrix, including off-diagonal entries.
        for scale in map(mp.mpf, (".85", "1.3")):
            same_pole = root_below_scalar(alpha, scale, run_parameters=True)
            if abs(same_pole-mass) > 10*alpha**2:
                raise ArithmeticError("full-matrix scale cancellation")
        samples.append({"G_over_16pi2": mp.nstr(alpha,5),
                        "mass_gap_over_m_alpha": mp.nstr((1-mass)/alpha,12),
                        "partial_width_over_leading": mp.nstr(trace_width/leading_width,12)})
    fmt = lambda x: mp.nstr(x,16)
    return {"model": "gaugeless symmetric phase, fixed Re scalar pole m=1, MSbar masses, Y at mu=m",
            "G": "|g|^2", "M1pole_over_m_minus_one_per_G": fmt(-1/loop),
            "nominal_m_plus_M1pole_over_m_minus_two_per_G": fmt(-1/loop),
            "stable_two_N1_two_mu_threshold_over_m_minus_two_per_G": fmt(-2/loop),
            "charged_scalar_partial_width_over_m_per_G_cubed": fmt(1/(1024*mp.pi**5)),
            "partial_width_over_gap_per_G_squared": fmt(1/(64*mp.pi**3)),
            "scope": "one-loop pole shift and leading threshold-induced two-body partial width",
            "determinant_checks": samples, "pole_degeneracy_retuning": check_retuning(loop)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = json.dumps(verify(), indent=2, allow_nan=False)+"\n"
    if args.output:
        args.output.write_text(data, encoding="utf-8")
    else:
        print(data, end="")


if __name__ == "__main__":
    main()
