#!/usr/bin/env python3
"""Pole residues and scale checks for the on-shell scalar-pair cut."""
import argparse
from functools import lru_cache
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import quad

HERE = Path(__file__).resolve().parent
OLD = HERE.parents[1] / "checkpoint_2026_11_19" / "yukawa_continuum"
LOOP = 16 * math.pi**2
# Each record contains M, |Y_fi|²/strength and M (Y_fi*)²/(Y_f1*)².
SPECIES = {"e": ((10, 1., 10.), (20, 2., -40.)),
           "mu": ((1, 1., 1.), (4, .25, -1.))}
CHANNELS = {"e": "E", "mu": "G"}


def lepton_residue(b, mu):
    """Per-sterile coefficient of the full lepton residue, a=m²=1."""
    value, _ = quad(lambda x: (1-x)*math.log((1+(b-1)*x)/mu**2),
                    0., 1., epsabs=3e-14)
    return value


def scalar_residue(b, mu):
    """Per-sterile full scalar residue at a=1; exact threshold limit."""
    if b == 1:
        return math.log(1/mu**2) - 2
    return math.log(b/mu**2) - 1 - b - (b*b-1)*math.log1p(-1/b)


def field_coefficients(mu):
    if not (math.isfinite(mu) and mu > 0):
        raise ValueError("positive finite MSbar scale mu/m required")
    result = {"scalar": {}, "lepton": {}}
    for flavor, species in SPECIES.items():
        key = CHANNELS[flavor]
        result["scalar"][key] = sum(w * scalar_residue(m*m, mu)
                                    for m,w,_ in species) / LOOP
        result["lepton"][key] = sum(w * lepton_residue(m*m, mu)
                                    for m,w,_ in species) / LOOP
    return result


@lru_cache(maxsize=16)
def external_factors(mu=1.0):
    """Two external scalars and two leptons. Includes 1/(16π²); use E,G."""
    field = field_coefficients(mu)
    return {f: {k: field["scalar"][k] +
                (field["lepton"][k] if k == CHANNELS[f] else 0.)
                for k in ("E", "G")} for f in SPECIES}


def tree_and_log(flavor, q):
    """Single transfer. Entries are tree and explicit slopes per E, G."""
    species = SPECIES[flavor]
    tree = sum(w/(m*m+q) for m,_,w in species)
    normal = sum(w/(m*m+q) for m,w,_ in species)
    # Scalar and lepton factors from the actual on-shell inverse kernels.
    logs = {}
    for f, key in CHANNELS.items():
        norm = sum(w for _,w,_ in SPECIES[f])
        logs[key] = -2*norm*tree
        if flavor == f:
            logs[key] -= norm*tree + 4*q*normal*tree
    return np.array([tree, logs["E"]/LOOP, logs["G"]/LOOP])


def cut_observables(s, order):
    """Evaluate J=0 and c=0 by independent on-shell sphere projections."""
    beta = math.sqrt(1-4/s)
    h, v = s/2-1, s*beta/2
    norm = s*beta/(32*math.pi)
    nodes, weights = np.polynomial.legendre.leggauss(order)
    phi = 2*math.pi*(np.arange(2*order)+.5)/(2*order)

    def both(f, x):
        return tree_and_log(f,h+v*x)+tree_and_log(f,h-v*x)

    def product(mu, e):
        return np.array([mu[0]*e[0], mu[0]*e[1]+mu[1]*e[0],
                         mu[0]*e[2]+mu[2]*e[0]])

    mu_avg = sum(w*both("mu",x)/2 for x,w in zip(nodes,weights))
    e_avg = sum(w*both("e",x)/2 for x,w in zip(nodes,weights))
    central = np.zeros(3)
    for x,w in zip(nodes,weights):
        y = math.sqrt(1-x*x)*np.cos(phi)
        central += w*product(both("mu",x),both("e",y).mean(axis=1))/2
    return norm*product(mu_avg,e_avg), norm*central


def require_close(actual, expected, rtol, label, atol=0.):
    if not np.allclose(actual, expected, rtol=rtol, atol=atol):
        raise ArithmeticError(f"{label}: {actual!r} vs {expected!r}")


def beta_functions(y, mass):
    h = y.conj().T @ y
    return 1.5*y@h + np.trace(h)*y, h.T@mass + mass@h


def tree_matrix(y, mass, function):
    # K=M f(M†M) in fixed coordinates, valid away from the diagonal basis.
    b, unitary = np.linalg.eigh(mass.conj().T @ mass)
    kernel = (mass @ unitary * function(b)) @ unitary.conj().T
    return y.conj() @ kernel @ y.conj().T


def implicit_j0(s, e_strength, g_strength, step):
    ye = math.sqrt(e_strength)*np.exp(.47j)
    g = math.sqrt(g_strength)*np.exp(-.73j)
    y = np.array([[ye,1j*math.sqrt(2)*ye,0,0], [0,0,g,-.5j*g]])
    mass = np.diag([10.,20.,1.,4.])
    dy, dm = beta_functions(y, mass)
    h = s/2-1
    beta = math.sqrt(1-4/s)
    v = s*beta/2
    mean = tree_matrix(y + step*dy, mass + step*dm,
                       lambda b: np.log1p(2*v/(b+h-v))/(2*v))
    gamma = ye**2 * g.conjugate()**2
    return s*beta/(8*math.pi)*mean[1,1]*mean[0,0].conjugate()/gamma


def stencil(fn, h=2e-5):
    return (fn(-2*h)-8*fn(-h)+8*fn(h)-fn(2*h))/(12*h)


def check_general_flavor():
    # Deliberately overlapping active rows with independent physical phases.
    y = np.array([[.8+.4j,.2-.3j,.1+.2j], [.7j,.6+.2j,.2-.8j]])
    m = np.array([1.,3.,7.])
    mass = np.diag(m)
    gram = y.conj().T@y
    p = y@y.conj().T
    dy, dm = beta_functions(y,mass)
    errors = []
    for q in (.08, 1., 9.):
        def amplitude(h):
            return tree_matrix(y+h*dy,mass+h*dm,lambda b: 1/(b+q))
        f = amplitude(0)
        inv = 1/(q+m*m)
        kinetic = (m[:,None]*gram+gram.T*m[None,:])*np.outer(inv,inv)
        # Each external lepton contributes minus one half of its own P.
        explicit = (-.5*(p.conj()@f+f@p)-2*np.trace(gram)*f
                    -2*q*y.conj()@kinetic@y.conj().T)
        run = stencil(amplitude)
        require_close(explicit,-run,4e-10,"complex matrix scale derivative",1e-11)
        errors.append(np.max(abs(run+explicit))/np.max(abs(run)))
    return max(errors)


def check_residues():
    """Differentiate the unexpanded pole kernel; integrate lepton moment independently."""
    with mp.workdps(55):
        for m in (1,4,10,20):
            b = mp.mpf(m*m)
            for mu in (1.,1.7):
                mu2 = mp.mpf(mu)**2
                # Finite B0 with its x=0 endpoint log integrated analytically.
                def inverse(z):
                    bubble = 1 - mp.quad(lambda x: mp.log((b-(1-x)*z)/mu2), [0,1])
                    return (z-b)*bubble - b*(1-mp.log(b/mu2))
                # A left derivative exists also for b=a at the branch point.
                # Use a finite difference explicitly one-sided there.
                h = mp.mpf('1e-19')
                rp = -(inverse(1-h)-inverse(1-2*h))/h
                rl = mp.quad(lambda x: (1-x)*mp.log((1+(b-1)*x)/mu2),[0,1])
                require_close(scalar_residue(m*m,mu),float(rp),3e-13,"scalar residue",5e-14)
                require_close(lepton_residue(m*m,mu),float(rl),2e-13,"lepton residue",5e-14)
    for f, key in CHANNELS.items():
        for k in ("E","G"):
            n = 3 if k=="E" else 1.25
            slope = stencil(lambda l: external_factors(math.exp(l))[f][k])
            require_close(slope, -(2+(key==k))*n/LOOP, 5e-10, "external log")


def verify():
    check_residues()
    complex_error = check_general_flavor()
    archived = json.loads((OLD / "RG_LOG.json").read_text())
    rows = []
    j0_by_energy = {}
    for s, prev in zip((6,8,10), archived["rows"], strict=True):
        j0, central = cut_observables(s,36)
        require_close(cut_observables(s,24),(j0,central),4e-12,"angular refinement")
        reference = [prev["central_amplitude_over_gamma"]["imaginary"],
                     prev["central_explicit_log_per_electron"]["imaginary"],
                     prev["central_explicit_log_per_muon"]["imaginary"]]
        require_close(central,np.array(reference,dtype=float),7e-13,"archived discontinuity")
        for e,g in ((1.,1.),(.7,1.4)):
            run = stencil(lambda h: implicit_j0(s,e,g,h))/LOOP
            require_close(run,-(e*j0[1]+g*j0[2]),1e-8,"sewn cut matrix derivative")
        j0_by_energy[s] = j0
        rows.append({"s_over_m2":s,
                     "J0":dict(zip(("leading","explicit_log_E","explicit_log_G"),j0)),
                     "c0":dict(zip(("leading","explicit_log_E","explicit_log_G"),central))})
    low, high = j0_by_energy[6], j0_by_energy[10]
    ratio = high[0]/low[0]
    slopes = ratio*(high[1:]/high[0]-low[1:]/low[0])
    for e,g in ((1.,1.),(.7,1.4)):
        run = stencil(lambda h: implicit_j0(10,e,g,h)/implicit_j0(6,e,g,h))/LOOP
        require_close(run,-(e*slopes[0]+g*slopes[1]),3e-8,"endpoint quotient flow")
    return {"scheme":"MSbar Y,M; fixed scalar pole m=1; unit external residues",
            "factors":"E=|y_e|²; G=|g|²; 1/(16π²) included; derivatives in ln(mu/m)",
            "field_residue_coefficients_mu_eq_m":field_coefficients(1.),
            "two_lepton_two_scalar_factors_mu_eq_m":external_factors(),
            "rows":rows,
            "J0_cut10_over_cut6":{"leading":ratio,"explicit_log_E":slopes[0],
                                  "explicit_log_G":slopes[1]},
            "checks":{"fixed_basis_complex_flavor_max_relative_residual":float(f"{complex_error:.2g}"),
                      "archive":["central upper-rim imaginary part, leading and both explicit logs"],
                      "angular_orders":[24,36],
                      "fixed_basis_step_E_G":[[1.,1.],[.7,1.4]]}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args = parser.parse_args()
    result = json.dumps(verify(),indent=2,allow_nan=False)+"\n"
    if args.output:
        args.output.write_text(result,encoding="utf-8")
    else:
        print(result,end="")
