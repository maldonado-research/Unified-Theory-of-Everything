"""Reconstruct the chiral vertex from crossed-channel unitarity and sew the J=0 cut."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad, quad_vec


LOOP = 16 * math.pi**2
# Each column contains squared mass, anomalous signed residue, and Yukawa modulus squared.
BLOCKS = {"electron": ((100., 10., 1.), (400., -40., 2.)),
          "muon": ((1., 1., 1.), (16., -1., .25))}


def check(actual, reference, rtol, name):
    if not np.allclose(actual, reference, rtol=rtol, atol=1.e-18):
        raise ArithmeticError(f"{name}: {actual} vs {reference}")


def integrate(function):
    return quad(function, 0., 1., epsabs=1.e-14, epsrel=2.e-13, limit=160)[0]


def spectral_density(u, mass2):
    """Im V(u)/pi at m_phi=1 after the crossed angular integral."""
    span = (u - 1)**2 / u
    start = (mass2 - 1) + (u - 1)
    fraction = span / start
    if fraction < .001:
        moment = sum(fraction**k / (k + 2) for k in range(6)) / start
    else:
        log_ratio = (-math.log1p(-fraction) if fraction < .8 else
                     math.log(start / ((mass2 - 1) + (u - 1) / u)))
        moment = (log_ratio / fraction - 1) / span
    return ((u - 1) / u)**2 * moment


def dispersion(z, mass2):
    """Map the whole positive u-axis cut, from its threshold to infinity."""
    def transformed(t):
        u = 1/t**2
        return 2*spectral_density(u, mass2)/(t*(1-z/u))
    return quad_vec(transformed, 0., 1., epsabs=1.e-14, epsrel=3.e-13)[0]


def vertex_simplex(z, mass2):
    def integrand(x):
        low, high = 1 - z*x, mass2 - x
        if high == low:
            return (1 - x) / low
        if abs(high - low) < .01 * low:
            return (1 - x) * math.log1p((high - low) / low) / (high - low)
        return (1 - x) * math.log(high / low) / (high - low)
    return integrate(integrand)


def helicity_check():
    # The outgoing lepton points in the +z direction. Test complex incident spinors.
    epsilon = np.array([[0., 1.], [-1., 0.]])
    out = np.array([0., 1.])
    for u in (1.1, 3., 37.):
        energy = (u - 1) / (2*math.sqrt(u))
        for cosine in (-.79, .12, .96):
            for inc in (np.array([.1+.4j, .8]), np.array([.7, -.1j])):
                values = []
                for azimuth in np.arange(8) * (math.pi/4):
                    phi = np.exp(-1j*azimuth)
                    active = math.sqrt(2*energy) * np.array(
                        [-phi*math.sqrt((1-cosine)/2), math.sqrt((1+cosine)/2)])
                    tree_product = (out @ epsilon @ active) * np.vdot(active, inc)
                    values.append(math.sqrt(u)*tree_product/(out @ epsilon @ inc))
                check(np.mean(values), (u-1)*(1-cosine)/2, 2.e-14, "cut helicity sum")


def external_fields(mass2):
    """Finite pole-residue contributions of a unit Yukawa at mu=m=1."""
    lepton = integrate(lambda x: -(1-x)*math.log1p(x*(mass2-1)))
    scalar = integrate(lambda x: -math.log(x*(mass2-1+x)))
    if mass2 != 1:
        scalar -= (mass2-1)*integrate(lambda x: (1-x)/(mass2-1+x))
    return lepton, scalar


def average(s, order):
    xs, weights = np.polynomial.legendre.leggauss(order)
    weights /= 2
    beta = math.sqrt(1-4/s)
    z = 1 - s/2 + s*beta*xs/2
    # The finite vector bubble integral comes directly from shifting the active loop momentum.
    bubble = quad_vec(lambda x: -(1-x)*np.log((1-x)*(1-x*z)), 0., 1.,
                      epsabs=1.e-14, epsrel=3.e-13)[0]
    means, legs, scalar = {}, {}, {}
    for flavor, block in BLOCKS.items():
        tree = sum(w/(mass2-z) for mass2, w, _ in block)
        kinetic = sum(strength/(mass2-z) for mass2, _, strength in block)
        vertex = sum(w*dispersion(z, mass2) for mass2, w, _ in block)
        means[flavor] = np.array([weights @ tree, weights @ (2*z*kinetic*vertex)/LOOP,
                                  weights @ (4*z*kinetic*tree*bubble)/LOOP])
        fields = [np.array(external_fields(mass2))*strength for mass2, _, strength in block]
        legs[flavor], scalar[flavor] = sum(fields)
    norm = s*beta/(8*math.pi)
    leading = norm*means["electron"][0]*means["muon"][0]
    rows = {"leading": leading}
    for flavor, partner in (("electron", "muon"), ("muon", "electron")):
        vertex, sterile = norm*means[partner][0]*means[flavor][1:]
        external = -leading*(legs[flavor]+2*scalar[flavor])/LOOP
        rows[flavor] = {"proper_vertex": vertex, "sterile": sterile,
                        "external_L_and_phi": external, "total": vertex+sterile+external}
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    helicity_check()
    for block in BLOCKS.values():
        for mass2, _, _ in block:
            for z in (-.013, -1., -8.9, 0., .46):
                direct = vertex_simplex(z, mass2)
                dispersive = dispersion(z, mass2)
                check(dispersive, direct, 3.e-12, "crossed cut vs simplex")
            for u in (1.013, 2., 6., 84.):
                from_simplex = quad(lambda x: (1-x)/(mass2-1+x*(u-1)), 1/u, 1,
                                    epsabs=1.e-14, epsrel=3.e-13)[0]
                check(spectral_density(u, mass2), from_simplex, 2.e-12, "timelike discontinuity")
    rows = []
    for s in (6, 8, 10):
        fine, coarse = average(s, 48), average(s, 32)
        check(fine["leading"], coarse["leading"], 2.e-12, "leading angle refinement")
        for flavor in BLOCKS:
            for name in fine[flavor]:
                check(fine[flavor][name], coarse[flavor][name], 2.e-12, "angle refinement")
        rows.append({"s_over_m2": s, **fine})
    result = {"scheme": "Full-theory MSbar Yukawas and sterile masses at mu=m, scalar pole mass m=1",
              "cut": "J=0 scalar pair, divided by gamma, coefficients per extra coupling modulus squared",
              "relative_dispersion_simplex_tolerance": 3.e-12,
              "checks": {"u_over_m2": [1.013, 2, 6, 84], "z_over_m2": [-.013, -1, -8.9, 0, .46],
                         "vertex_angular_helicity_tests": 18, "polar_orders": [32, 48]},
              "rows": rows}
    # Round physics coefficients at a precision below the quadrature validation.
    for row in rows:
        row["leading"] = float(f'{row["leading"]:.12g}')
        for flavor in BLOCKS:
            row[flavor] = {k:float(f"{v:.12g}") for k,v in row[flavor].items()}
    data = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(data, encoding="utf-8")
    else:
        print(data, end="")


if __name__ == "__main__":
    main()
