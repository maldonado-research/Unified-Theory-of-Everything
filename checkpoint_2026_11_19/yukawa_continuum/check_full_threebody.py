#!/usr/bin/env python3
"""Full N1 three-particle unitarity cut via analytic and unreduced spin/phase-space sewing."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss


def mu(q):
    return 15 / ((1 + q) * (16 + q))


def electron(q):
    return -30 * q / ((100 + q) * (400 + q))


def normal(q):
    return 1 / (1 + q) + 0.25 / (16 + q)


def decay_function(u):
    t = u - 1
    if abs(t) < 1e-3:
        return t * (1/2 + t * (-1/6 + t * (1/12 + t * (-1/20 + t/30))))
    return u * math.log1p(t) / t - 1


def j(b, h, v):
    if v == 0:
        return 1 / (b + h)
    return math.log1p(2*v / (b + h - v)) / (2*v)


def projection(u, b, h, v):
    """Mean of q/[(U+q)(b+q)] with stable evaluation near coincident masses."""
    if abs(u - b) < 1e-5 * b:
        # Form the divided difference of b*J_b using log1p for J_U-J_b.
        d = (b + h + v) * (u + h - v)
        z = -2 * v * (u - b) / d
        log_ratio = 1.0 if z == 0 else math.log1p(z) / z
        return j(u, h, v) - b * log_ratio / d
    return (u * j(u, h, v) - b * j(b, h, v)) / (u - b)


def analytic(s, order=40):
    """Invariant-mass integral, decay-angle reduction, J0 logs and c0 azimuthal mean."""
    if s <= 4:
        return np.zeros((2, 2))
    roots, weights = leggauss(order)
    xs, ws = (1 + roots)/2, weights/2
    span = (s - 4) * math.sqrt(s) / (math.sqrt(s) + 2)
    result = np.zeros((2, 2))  # cut (c0,J0) x (decay,exchange)
    for t, weight in zip(xs, ws):
        gap = span * t**2
        u = 1 + span * (1-t*t)
        h = (s - 1 - u)/2
        v = math.sqrt(gap * (4*math.sqrt(s)+gap))/2
        d = decay_function(u)
        plus, minus = h + v*xs, h - v*xs
        sm = mu(plus) + mu(minus)
        sx = d * (normal(plus)*plus/(u+plus) + normal(minus)*minus/(u+minus))
        transverse = v*v * (1-xs*xs)
        se = 2 * (10 / np.sqrt((100+h)**2 - transverse) - 40 / np.sqrt((400+h)**2 - transverse))
        e_bar = 10*j(100,h,v) - 40*j(400,h,v)
        mu_bar = j(1,h,v) - j(16,h,v)
        x_bar = d * (projection(u,1,h,v) + 0.25*projection(u,16,h,v))
        measure = weight * 2*span*t * (2*v) / (256*math.pi**3*u)
        result += measure * np.array([[ws @ (sm*se), ws @ (sx*se)],
                                      [4*e_bar*mu_bar, 4*e_bar*x_bar]])
    return result


def dot(p, q):
    return p[..., 0]*q[..., 0] - np.sum(p[..., 1:]*q[..., 1:], axis=-1)


def spin_kernel(s, k, lepton):
    """Cut spin sum relative to B12*(U-1), before the L,N angular integration."""
    energy = math.sqrt(s)/2
    total = np.array([2*energy, 0., 0., 0.])
    n = total - k - lepton
    decay, exchange = 0.0, 0.0j
    for sign in (-1, 1):
        pi = np.array([energy, 0., 0., sign*energy])
        pj = total-pi
        tr = pj-k
        q = -dot(tr, tr)
        ds = -2*dot(pi,n)
        # Tr[P_{-j} (l0-σ·l)(q0+σ·q)] for q=p_j-k, after using p_jslash*u_j=0.
        imag = lepton[..., 1]*k[2] - lepton[..., 2]*k[1]
        spin = -dot(lepton,k) - sign*(lepton[..., 0]*k[3]-k[0]*lepton[..., 3]-1j*imag)
        decay += mu(q)
        exchange += normal(q)*spin/ds
    return decay, exchange


def rest_directions(roots):
    phi = (np.arange(2*len(roots))+0.5)*math.pi/len(roots)
    x = roots[:, None]
    radius = np.sqrt(1-x*x)
    return np.stack(np.broadcast_arrays(radius*np.cos(phi),radius*np.sin(phi),x), axis=-1)


def direct(s, order=20):
    """Recoil measure, boosted decay spinors' spin trace and explicit production directions."""
    roots, weights = leggauss(order)
    directions = rest_directions(roots)
    axes = np.r_[0., roots]
    phis = (np.arange(2*order)+0.5)*math.pi/order
    pmax = math.sqrt(s-4)/2
    result = np.zeros((2,2))
    for unit, radial_w in zip((1+roots)/2, weights/2):
        p = unit*pmax
        k0 = math.sqrt(1+p*p)
        q0 = math.sqrt(s)-k0
        u = q0*q0 - p*p
        omega = (u-1)/(2*math.sqrt(u))
        boost = q0/math.sqrt(u)
        boost_velocity = p/math.sqrt(u)
        sphere = np.zeros((len(axes),2))
        for x, polar_w in zip(roots, weights/2):
            transverse = math.sqrt(1-x*x)
            khat = np.array([transverse,0.,x])
            k = np.r_[k0, p*khat]
            axis = -khat
            along = directions @ axis
            lepton0 = omega * (boost + boost_velocity*along)
            lepton_spatial = omega*(directions + ((boost-1)*along + boost_velocity)[...,None]*axis)
            lepton = np.concatenate([lepton0[...,None], lepton_spatial], axis=-1)
            decay, raw_exchange = spin_kernel(s,k,lepton)
            exchange = (weights/2) @ raw_exchange.mean(axis=1)
            if abs(exchange.imag) > 1e-12:
                raise ArithmeticError("reflection symmetry of the spin trace")
            y = x*axes[:,None] + transverse*np.sqrt(1-axes[:,None]**2)*np.cos(phis)
            # Direct time and spatial components of both electron-side transfers.
            qe1 = math.sqrt(s)*(k0-p*y) - 1
            qe2 = math.sqrt(s)*(q0+p*y) - u
            se = (electron(qe1)+electron(qe2)).mean(axis=1)
            sphere += polar_w * se[:,None] * np.array([decay,exchange.real])
        measure = pmax * radial_w * s * p*p/(64*math.pi**3 * k0*u)
        result += measure * np.array([sphere[0], (weights/2) @ sphere[1:]])
    return result


def verify():
    def check(got, expected, error, label):
        if not np.allclose(got, expected, rtol=error, atol=0):
            raise ArithmeticError(f"{label}: {got} vs {expected}")

    def rounded(x):
        return float(f"{x:.10g}")

    rows = []
    for s in (6,8,10):
        exact = analytic(s)
        raw = direct(s)
        check(raw, exact, 6e-10, "direct recoil and decay sphere")
        check(analytic(s,64), exact, 1e-11, "analytic quadrature refinement")
        if not np.all(exact < 0):
            raise ArithmeticError("charged-cut signs")
        row = {"s_over_m2":s}
        for name, vals, numeric in zip(("c0","J0"),exact,raw):
            row[name] = {"decay":rounded(vals[0]),"exchange":rounded(vals[1]),
                         "full":rounded(sum(vals)),"direct_full":rounded(sum(numeric)),
                         "full_over_decay":rounded(sum(vals)/vals[0])}
        rows.append(row)
    mg = json.loads(Path(__file__).with_name("MADGRAPH.json").read_text(encoding="utf-8"))
    mg_s = mg["s_over_m2"]
    overlaps = []
    for point in [mg["fixed_point"], *mg["other_points"]]:
        p1, p2, lepton, n, k = np.array(point["p_mu_or_e_in_p_mu_out_p_N_p_phi"])
        pole, extra = spin_kernel(mg_s, k, lepton)
        ratio = 1 + extra / pole
        u = dot(n + lepton, n + lepton)
        se = sum(electron(-dot(pi - k, pi - k)) for pi in (p1, p2))
        prefactor = mg_s * se / (u - 1)
        check(prefactor * pole, complex(*point["pole_overlap"]), 2e-13, "MadGraph decay overlap")
        check(prefactor * (pole + extra), complex(*point["full_overlap"]), 2e-13, "MadGraph spin overlap")
        check(ratio, complex(*point["overlap_ratio"]), 2e-13, "MadGraph overlap ratio")
        overlaps.append([rounded(ratio.real), rounded(ratio.imag)])
    projected = mg["projected_J0"][-1]
    decay_j0, exchange_j0 = analytic(mg_s)[1]
    full_j0 = decay_j0 + exchange_j0
    check(decay_j0, projected["pole"], 3e-10, "MadGraph projected scalar contribution")
    check(full_j0, projected["full"], 3e-10, "MadGraph projected full contribution")
    threshold_rows = []
    for offset in (1e-3, 1e-4, 1e-5):
        near = analytic(4 + offset)
        slopes = near[:, 1] / near[:, 0] / offset
        check(slopes, 7/120, offset, "exchange onset relative to decay")
        threshold_rows.append(dict(delta=offset, at_c0=rounded(slopes[0]), J0=rounded(slopes[1])))
    return dict(scope="complete phi- mu- N1 contribution to Hermitian cut divided by gamma*abs(g)^2",
                phase="B34_conjugate_times_B12 = s", m=1, M_over_m=10,
                J0="external-cosine mean", threshold_exchange_over_decay_per_s_minus_4="7/120",
                threshold_slope_samples=threshold_rows,
                generator_audit=dict(reference="MADGRAPH.json", s_over_m2=mg_s,
                                     forward_spin_overlap_ratios=overlaps,
                                     projected_full=rounded(projected["full"]),
                                     projected_scalar_insertion=rounded(projected["pole"])), rows=rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    payload=json.dumps(verify(),indent=2,allow_nan=False)+"\n"
    if args.output is None:
        print(payload,end="")
    else:
        args.output.write_text(payload,encoding="utf-8")


if __name__ == "__main__":
    main()
