"""Sew the one-loop production kernels and the three-body trees in angular moments."""
import argparse
from functools import lru_cache
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss, legvander, legval
from scipy.integrate import quad
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parents[1] / "checkpoint_2026_11_19"
LOOP = 16*np.pi**2
BLOCKS = {
    "electron": ((100., 400.), (10., -40.), (1., 2.)),
    "muon": ((1., 16.), (1., -1.), (1., 0.25)),
}


@lru_cache(maxsize=None)
def rule(n):
    return leggauss(n)


def tree_and_normal(flavor, q):
    masses, _, weights = BLOCKS[flavor]
    if flavor == "electron":
        tree = -30*q / ((100+q)*(400+q))
    else:
        tree = 15 / ((1+q)*(16+q))
    normal = sum(w/(b+q) for b,w in zip(masses, weights))
    return tree, normal


def triangle(q, b, order):
    # Integrate the linear denominator along a line of fixed massless-line
    # parameter; u=t^4 resolves the logarithmic endpoint at b=a=1.
    nodes, weights = rule(order)
    t = (1+nodes)/2
    u = t**4
    z = 1-u
    upper = 1 + np.multiply.outer(q.ravel(), z)
    lower = b-1+u
    ratio = (upper-lower)/lower
    same = abs(ratio) < 1e-12
    logmean = np.log1p(ratio)/np.where(same, 1., ratio)
    logmean = np.where(same, 1.-ratio/2, logmean)
    return np.sum((logmean/lower)*(weights*2*t**3*u), axis=-1).reshape(q.shape)


def production(flavor, q, order):
    tree, normal = tree_and_normal(flavor, q)
    masses, residues, _ = BLOCKS[flavor]
    bubble = 2+1/q-((1+q)/q)**2*np.log1p(q)
    vertex = sum(w*triangle(q, b, order) for b,w in zip(masses, residues))
    factor = -2*q*normal/LOOP
    return tree, factor*bubble*tree, factor*vertex


@lru_cache(maxsize=None)
def external(flavor):
    # Combined two sides. Four scalar legs and the two active legs of this flavor.
    masses, _, weights = BLOCKS[flavor]
    total = 0.
    for b,w in zip(masses, weights):
        lepton = quad(lambda x: (1-x)*np.log1p(x*(b-1)), 0,1)[0]
        scalar = quad(lambda x: np.log((1-x)*(b-x)) +
                      ((b-1)*x/(b-x) if b != 1 else 0.), 0,1)[0]
        total += w*(lepton+2*scalar)
    return total/LOOP


def project(order, ell_max):
    x, weights = rule(order)
    return x, legvander(x, ell_max)[:, ::2].T * (weights/2)


def recoil(s, t):
    span = (s-4)*np.sqrt(s)/(np.sqrt(s)+2)
    gap = span*t*t
    u = 1+span*(1-t*t)
    h = (s-1-u)/2
    v = np.sqrt(gap*(4*np.sqrt(s)+gap))/2
    measure = (2*span*t)*(2*v)/(256*np.pi**3*u)
    return u, h, v, measure


def decay_function(u):
    gap = u-1
    if gap < 1e-4:
        return gap*(0.5 + gap*(-1/6 + gap/12))
    return u*np.log1p(gap)/gap-1


def two_and_three(s, order, ell_max):
    if not 4 < s < 25:
        raise ValueError("open first-threshold interval 4 < s/m² < 25 required")
    x, mean = project(order, ell_max)
    h = s/2-1
    v = np.sqrt(s*(s-4))/2
    q = h-v*x
    prefactor = v/(4*np.pi)
    ef, es, ev = [mean @ f for f in production("electron", q, order)]
    mf, ms, mv = [mean @ f for f in production("muon", q, order)]
    # Moments are <A P_l>, hence the common (2l+1) belongs only in synthesis.
    leading = prefactor*ef*mf
    sterile = prefactor*np.array([es*mf, ef*ms])
    vertex = prefactor*np.array([ev*mf, ef*mv])
    legs = np.outer([external("electron"), external("muon")], leading)
    three = np.zeros_like(leading)
    for t,w in zip((1+x)/2, rule(order)[1]/2):
        u, h3, v3, measure = recoil(s,t)
        q3 = h3-v3*x
        m3, normal = tree_and_normal("muon",q3)
        e3, _ = tree_and_normal("electron",q3)
        full = m3+decay_function(u)*normal*q3/(u+q3)
        three += w*4*measure*(mean @ full)*(mean @ e3)
    two = sterile+vertex+legs
    total = two.copy()
    total[1] += three
    return dict(leading=leading, sterile=sterile, vertex=vertex, external=legs,
                two=two, three=three, total=total)


def angular(moments, cosine):
    degree = 2*np.arange(len(moments))
    coefficients = np.zeros(degree[-1]+1)
    coefficients[::2] = (2*degree+1)*moments
    return legval(cosine, coefficients)


def sphere_at_angle(s, cosine, n):
    """Direct sphere integration at fixed external angle, including the continuum."""
    x, wx = rule(n)
    phis = 2*np.pi*(np.arange(2*n)+0.5)/(2*n)
    y = cosine*x[:,None] + np.sqrt((1-cosine**2)*(1-x*x))[:,None]*np.cos(phis)
    h = s/2-1
    v = np.sqrt(s*(s-4))/2

    def boson_sum(flavor, axis):
        a = production(flavor,h-v*axis,n)
        b = production(flavor,h+v*axis,n)
        return [f+g for f,g in zip(a,b)]

    ef, es, ev = [f.mean(axis=1) for f in boson_sum("electron",y)]
    mf, ms, mv = boson_sum("muon",x)
    factor = v/(16*np.pi)
    base = factor*(wx/2)@(mf*ef)
    vals = np.array([factor*(wx/2)@(mf*(es+ev)),factor*(wx/2)@((ms+mv)*ef)])
    vals += base*np.array([external("electron"),external("muon")])
    three = 0.
    for t,w in zip((1+x)/2,wx/2):
        u, h3, v3, measure = recoil(s,t)
        ms = np.zeros_like(x)
        es = np.zeros_like(y)
        for sign in (-1,1):
            q = h3+sign*v3*x
            tree, normal = tree_and_normal("muon",q)
            ms += tree+decay_function(u)*normal*q/(u+q)
            es += tree_and_normal("electron",h3+sign*v3*y)[0]
        three += w*measure*((wx/2)@(ms*es.mean(axis=1)))
    return base, vals, three


def close(actual, reference, relative, label):
    if not np.allclose(actual,reference,rtol=relative,atol=1e-19):
        raise ArithmeticError(f"{label}: {actual} versus {reference}")


def rounded(value):
    return float(f"{value:.11g}")


def pair(vals):
    return {key:rounded(x) for key,x in zip(BLOCKS,vals)}


def energy_shape():
    """Locate a half cancellation and a local maximum of the muon coefficient."""
    def pieces(s, order):
        result = two_and_three(s, order, 0)
        return result["two"][1,0], result["three"][0]

    def locate(order, step):
        def half(s):
            two, three = pieces(s, order)
            return two/2 + three

        def derivative(s):
            lower = sum(pieces(s-step, order))
            upper = sum(pieces(s+step, order))
            return (upper-lower)/(2*step)

        if not (half(8)>0>half(10) and derivative(16)>0>derivative(19)):
            raise ArithmeticError("missing bracket for the cancellation or turning point")
        return brentq(half,8,10), brentq(derivative,16,19,xtol=1e-9)

    coarse = locate(36,0.004)
    half, maximum = locate(56,0.001)
    close(coarse,(half,maximum),5e-8,"half cancellation and peak refinement")
    scan = []
    for s in range(5,25):
        two, three = pieces(s,56)
        scan.append([s,rounded(two),rounded(three),rounded(two+three)])
    return {"scan_columns":["s/m^2","two-particle","three-particle","sum"],
            "muon_J0_scan":scan, "half_cancellation_s_over_m2":rounded(half),
            "local_maximum_s_over_m2":rounded(maximum),
            "local_maximum_combined_J0":rounded(sum(pieces(maximum,56)))}


def verify():
    external_ref = json.loads((HERE/"QCDLOOP_VERTEX.json").read_text())
    q_ref = np.array(external_ref["q_over_m2"])
    for b, values in external_ref["I_per_squared_mass"].items():
        close(triangle(q_ref,float(b),64),values,3e-13,"QCDLoop scalar integral reduction")
    previous = json.loads((PRIOR/"uv_completion"/"CUT.json").read_text())
    leading_refs = {row["s_over_m2"]:row for row in previous["rows"]}
    three_refs = json.loads((PRIOR/"yukawa_continuum"/"FULL_THREE_BODY.json").read_text())["rows"]
    output = []
    for s, prior3 in zip((6,8,10),three_refs, strict=True):
        fine = two_and_three(s,56,16)
        coarse = two_and_three(s,36,12)
        for key in fine:
            close(fine[key][...,0],coarse[key][...,0],1e-9,"projection refinement "+key)
        lo, two, three, total = [fine[k] for k in ("leading","two","three","total")]
        close(lo[0],leading_refs[s]["J0_cut_analytic"],6e-10,"leading J0")
        close(lo[1],leading_refs[s]["J2_moment_from_boxes"],9e-10,"leading J2")
        close(three[0],prior3["J0"]["full"],6e-10,"prior complete three-particle cut")
        for c in (0.,0.37,1.):
            raw_lo, raw_two, raw_three = sphere_at_angle(s,c,36)
            close(raw_lo,angular(lo,c),2e-10,"leading angular convolution")
            close(raw_two,[angular(a,c) for a in two],2e-10,"finite two-particle angular convolution")
            close(raw_three,angular(three,c),2e-10,"three-particle angular convolution")
        ster = json.loads((HERE/"STERILE.json").read_text())
        record = next(r for r in ster["rows"] if r["s_over_m2"]==s)
        close(fine["sterile"][:,0],[float(record[f]) for f in BLOCKS],2e-9,"signed sterile insertion")
        vertex = json.loads((HERE/"VERTEX.json").read_text())["energies_s_over_m2"][str(s)]
        close(fine["vertex"][:,0],[vertex["A2_vertex_e"],vertex["A2_vertex_mu"]],2e-9,"signed vertex")
        fields = json.loads((HERE/"EXTERNAL.json").read_text())["two_lepton_two_scalar_factors_mu_eq_m"]
        close(fine["external"][:,0], [lo[0]*(fields["e"][key]+fields["mu"][key]) for key in ("E","G")],
              1e-11, "external pole residues")
        ell = 2*np.arange(len(lo))
        cut_weight = (2*ell+1)*lo/(32*np.pi*s)
        square = cut_weight @ lo
        slopes = 2*(total @ cut_weight)
        close(square,leading_refs[s]["m2_sigma_floor_full_cut"],1e-9,"leading squared cut")
        output.append(dict(s_over_m2=s, leading_J0=rounded(lo[0]),
            sterile_J0=pair(fine["sterile"][:,0]), vertex_J0=pair(fine["vertex"][:,0]),
            external_J0=pair(fine["external"][:,0]), two_particle_J0=pair(two[:,0]),
            three_particle_J0_per_muon=rounded(three[0]), combined_J0=pair(total[:,0]),
            fractional_J0=pair(total[:,0]/lo[0]),
            combined_J2=pair(total[:,1]), combined_at_cosine_zero=pair([angular(a,0.) for a in total]),
            leading_absorptive_rate_part=rounded(square), absorptive_rate_interference=pair(slopes)))
    low, high = output[0], output[-1]
    response = {f:rounded(high["fractional_J0"][f]-low["fractional_J0"][f]) for f in BLOCKS}
    return dict(
        process="mu_L- mu_L- -> e_L- e_L-, gaugeless unbroken doublet, stated four-singlet benchmark",
        units="scalar pole m=1; MSbar flavor and sterile parameters at mu=m; B34* B12=s",
        normalization="anti-Hermitian s-channel cut divided by gamma; corrections multiply |y|^2 or |g|^2; J_l=<P_l A>",
        observable="rate component m^2 <(A/gamma)^2>/(32 pi s), coefficients linear in the extra Yukawa pair",
        angular_checks="even Legendre moments through 16; direct sphere at c=0,0.37,1; orders 36 and 56",
        rows=output, leading_signed_J0_ratio_10_over_6=rounded(high["leading_J0"]/low["leading_J0"]),
        fractional_correction_to_J0_ratio_10_over_6=response, muon_energy_dependence=energy_shape())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args = parser.parse_args()
    data = json.dumps(verify(),indent=2,allow_nan=False)+"\n"
    if args.output:
        args.output.write_text(data,encoding="utf-8")
    else:
        print(data,end="")
