"""Pure-Yukawa coefficient of the physical upper-rim box's explicit RG logarithm."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from mpmath import mp

HERE = Path(__file__).resolve().parent
UV = HERE.parent / "uv_completion"


def load(name):
    spec = importlib.util.spec_from_file_location(name, UV / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


dispersion = load("check_bounds")
simplex = load("check_box_direct")


def close(actual, reference, tolerance, label):
    if abs(actual - reference) > tolerance * abs(reference):
        raise ArithmeticError(f"{label}: {actual} vs {reference}")


def angular_derivative(x, value):
    """Derivative of the [0,1] integral of 1/(1+x*q*q)."""
    if abs(x) < mp.mpf("1e-9"):
        return -mp.mpf(1)/3 + 2*x/5 - 3*x*x/7 + 4*x**3/9
    return (1/(1+x) - value)/(2*x)


def amplitude_slopes(s, cosine, order):
    """LO and both explicit logs per reference gamma. Scalar pole mass m=1."""
    beta2 = 1 - 4/s
    beta = mp.sqrt(beta2)
    pole = beta*(-mp.atanh(beta) + mp.j*mp.pi/2)
    total = electron = muon = 0
    for t in (-s*(1-cosine)/2, -s*(1+cosine)/2):
        for b, w, b_log in ((1, 1, 4), (16, -1, 1)):
            for j, (d, v, d_log) in enumerate(((100, 10, 4), (400, -40, 8))):
                span = mp.log(mp.mpf(d)/b)
                for u, weight in dispersion.unit_rule(order):
                    q = b*mp.exp(u*span)
                    lam = (d-q)/(d-b)
                    jac = q*span/(d-b)
                    ell = (1+q)**2
                    h = 4*(q-t*lam*(1-lam))-ell
                    x = h/ell
                    e = dispersion.angular_primitive(x)
                    den = ell+h*beta2
                    kernel = (pole+e)/den
                    # Differentiate q at fixed physical Feynman parameter
                    # lambda, then integrate using its reference log coordinate.
                    dell = 2*(1+q)
                    dh = 4-dell
                    dx = (dh-x*dell)/ell
                    dk = (angular_derivative(x, e)*dx-kernel*(dell+dh*beta2))/den
                    factor = weight*jac*w*v
                    total += factor*kernel
                    electron += factor*((mp.mpf(65)/3+4*j)*kernel + d_log*d*(1-lam)*dk)
                    muon += factor*(mp.mpf(209)/20*kernel + b_log*b*lam*dk)
    norm = 2*mp.pi**2
    loop = 16*mp.pi**2
    return total/norm, -electron/(norm*loop), -muon/(norm*loop)


def spectral_observables(energy, n_mass, n_angle):
    s = mp.mpf(energy)
    center = amplitude_slopes(s, 0, n_mass)
    rate = ds_e = ds_mu = 0
    for cosine, weight in dispersion.unit_rule(n_angle):
        leading, el, mu = amplitude_slopes(s, cosine, n_mass)
        rate += weight*abs(leading)**2
        ds_e += weight*2*(leading.conjugate()*el).real
        ds_mu += weight*2*(leading.conjugate()*mu).real
    ps = 1/(32*mp.pi*s)
    return center, (rate*ps, ds_e*ps, ds_mu*ps)


def fixed_basis_step(e_strength, mu_strength, step):
    """Advance the fundamental full matrices in D-time, then solve their Takagi problem."""
    y, g = np.sqrt(e_strength), np.sqrt(mu_strength)
    yukawa = np.array([[y, 1j*np.sqrt(2)*y, 0, 0], [0, 0, g, -.5j*g]])
    mass = np.diag([10., 20., 1., 4.])
    gram = yukawa.conj().T @ yukawa
    dy = 1.5*yukawa @ gram + np.trace(gram)*yukawa
    dm = gram.T @ mass + mass @ gram
    varied_mass = mass + step*dm
    _, basis = np.linalg.eigh(varied_mass.conj().T @ varied_mass)
    diagonal = np.diag(basis.T @ varied_mass @ basis)
    basis = basis*np.exp(-.5j*np.angle(diagonal))
    masses = np.abs(diagonal)
    if not np.allclose(basis.T @ varied_mass @ basis, np.diag(masses), atol=1e-13, rtol=1e-13):
        raise ArithmeticError("Takagi rediagonalization")
    varied_y = (yukawa + step*dy) @ basis
    mu_residues = masses*varied_y[1].conj()**2
    e_residues = masses*varied_y[0]**2
    # Sorted singular values have muon support at 0,1 and electron at 2,3.
    return list(zip(masses[:2]**2, mu_residues[:2])), list(zip(masses[2:]**2, e_residues[2:]))


def direct_observables(s, params, gamma, mass_rule, angle_rule):
    muon, electron = params
    def amp(c):
        result = 0j
        for b, w in muon:
            for d, v in electron:
                result += w*v*sum(simplex.box(s, t, b, d, mass_rule)
                                  for t in (-s*(1-c)/2, -s*(1+c)/2))
        return result*s/(16*np.pi**2*gamma)
    cs, weights = angle_rule
    rate = sum(weight*abs(amp(c))**2 for weight,c in zip(weights, cs))/(32*np.pi*s)
    return np.array([amp(0), rate])


def check_direct(energy, center, rates):
    # A one-sided fourth-order stencil keeps the varied light sterile mass
    # above the scalar mass, inside the previous direct integrator's domain.
    step = 0.0001
    mass_rule = simplex.unit_rule(96)
    angle_rule = simplex.unit_rule(8)
    for e_strength, mu_strength in ((1, 1), (1, 2)):
        gamma = e_strength*mu_strength
        derivative = np.zeros(2, dtype=complex)
        for n, weight in enumerate((-25, 48, -36, 16, -3)):
            varied = fixed_basis_step(e_strength, mu_strength, n*step)
            derivative += weight*direct_observables(energy, varied, gamma, mass_rule, angle_rule)
        explicit = -derivative/(12*step*16*np.pi**2)
        for got, coeffs in zip(explicit, (center, rates)):
            target = complex(e_strength*coeffs[1] + mu_strength*coeffs[2])
            close(got, target, 3e-7, "matrix finite difference, independent Feynman simplex")


def fmt(value):
    return mp.nstr(value, 14)


def pair(value):
    return {"real": fmt(value.real), "imaginary": fmt(value.imag)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    mp.dps = 42
    rates = {}
    rows = []
    prior = json.loads((HERE.parent / "radiative_stability" / "RESCATTERING.json").read_text())
    for s, old in zip((6,8,10), prior["rows"], strict=True):
        if old["s_over_m2"] != s:
            raise ArithmeticError("previous energy mismatch")
        coarse = spectral_observables(s, 20, 8)
        center, rate = spectral_observables(s, 32, 14)
        for cvals, fvals in zip(coarse, (center, rate)):
            for a, b in zip(cvals, fvals):
                close(a, b, mp.mpf("3e-12"), "quadrature refinement")
        close(rate[0], mp.mpf(old["scaled_rate_at_leading_order"]), mp.mpf("2e-15"), "previous physical rate")
        check_direct(s, center, rate)
        rates[s] = rate
        rows.append({"s_over_m2": s, "central_amplitude_over_gamma": pair(center[0]),
                     "central_explicit_log_per_electron": pair(center[1]),
                     "central_explicit_log_per_muon": pair(center[2]),
                     "scaled_rate": fmt(rate[0]), "scaled_explicit_log_per_electron": fmt(rate[1]),
                     "scaled_explicit_log_per_muon": fmt(rate[2]),
                     "fractional_explicit_log_per_electron": fmt(rate[1]/rate[0]),
                     "fractional_explicit_log_per_muon": fmt(rate[2]/rate[0])})
    low, high = rates[6], rates[10]
    ratio = high[0]/low[0]
    slopes = [high[k]/high[0] - low[k]/low[0] for k in (1,2)]
    report = {"explicit_log": "ln(mu/mu_ref), opposite to implicit variation of physical leading term",
              "mass_unit": "fixed scalar pole mass m=1; unbroken four-singlet full theory in MSbar",
              "coupling_factors": "amplitudes divided by reference gamma; rates multiplied by m^2/|gamma_ref|^2; electron and muon coefficients multiply |y_e|^2 and |g|^2",
              "checks": {"spectral_refinement_relative_tolerance": "3e-12",
                         "matrix_step_direct_simplex_relative_tolerance": "3e-7",
                         "matrix_step_coupling_strengths": [[1,1], [1,2]]},
              "rows": rows, "endpoint_ratio": {"sigma10_over_sigma6": fmt(ratio),
                         "explicit_log_per_electron": fmt(ratio*slopes[0]),
                         "explicit_log_per_muon": fmt(ratio*slopes[1]),
                         "fractional_explicit_log_per_electron": fmt(slopes[0]),
                         "fractional_explicit_log_per_muon": fmt(slopes[1])}}
    data = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(data, encoding="utf-8")
    else:
        print(data, end="")


if __name__ == "__main__":
    main()
