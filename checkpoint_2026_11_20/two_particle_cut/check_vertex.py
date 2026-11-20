"""Proper Yukawa endpoint insertions on the charged two-scalar cut, in units m=1."""
import argparse
import json
import math
from pathlib import Path

from scipy.integrate import quad

# Entries are (squared Takagi mass, relative squared Yukawa).
BLOCKS = {"muon": ((1, 1), (16, -0.25)), "electron": ((100, 1), (400, -2))}
LOOP = 16 * math.pi**2


def integrate(function):
    return quad(function, 0, 1, epsabs=2e-15, epsrel=2e-12, limit=120)[0]


def triangle(q, b, a=1):
    """Feynman triangle after the exact integration over the relative massive-line weight."""
    def kernel(z):
        left, right = a + z*q, b - z*a
        ratio = right / left
        if abs(ratio - 1) < 1e-6:
            diff = ratio - 1
            log_mean = (1 - diff/2 + diff**2/3) / left
        else:
            log_mean = math.log(ratio) / (right - left)
        return (1-z)*log_mean
    return integrate(kernel)


def cut_density(q, b, u):
    """Positive Im(I)/pi at scalar virtuality u, with lepton and sterile on shell."""
    span = (u - b)/u
    first = 1 + b*q/u
    slope = u + q
    eps = slope*span/first
    if eps < 1e-4:
        primitive = eps**2*(0.5 - eps/6 + eps**2/12 - eps**3/20)
    else:
        primitive = (1+eps)*math.log1p(eps) - eps
    return first*primitive/slope**2


def triangle_dispersion(q, b):
    """Unsubtracted integral in scalar virtuality, mapping u=b/(1-v)."""
    def kernel(v):
        u = b/(1-v)
        return b*cut_density(q, b, u) / ((1-v)*(b-1+v))
    return integrate(kernel)


def tree(q, block):
    return sum(r*math.sqrt(b)/(b+q) for b, r in block)


def vertex(q, block, form=triangle):
    normal = sum(abs(r)/(b+q) for b, r in block)
    anomalous = sum(r*math.sqrt(b)*form(q, b) for b, r in block)
    return -2*q*normal*anomalous/LOOP


def projection(s, block, kernel):
    center = s/2 - 1
    radius = math.sqrt(s*(s-4))/2
    return integrate(lambda x: kernel(center + radius*(2*x-1), block))


def figures(s, form):
    tree_mu, tree_e = (projection(s, block, tree) for block in BLOCKS.values())
    v_mu, v_e = (projection(s, block, lambda q, b: vertex(q, b, form))
                 for block in BLOCKS.values())
    scale = s*math.sqrt(1-4/s)/(8*math.pi)
    return {"mean_F_mu": tree_mu, "mean_F_e": tree_e,
            "mean_V_mu": v_mu, "mean_V_e": v_e,
            "A2_vertex_mu": scale*v_mu*tree_e, "A2_vertex_e": scale*tree_mu*v_e}


def close(a, b, tolerance, label):
    if not math.isclose(a, b, rel_tol=tolerance, abs_tol=1e-17):
        raise ArithmeticError(f"{label}: {a} != {b}")


def result():
    values = {}
    for s in (6, 8, 10):
        params = figures(s, triangle)
        dispersive = figures(s, triangle_dispersion)
        for key, value in params.items():
            close(value, dispersive[key], 5e-11, f"s={s}, {key}")
        values[str(s)] = {key: float(f"{v:.12g}") for key, v in params.items()}

    # Recover the previously cut exchange term from this triangle's scalar-leg cut.
    for u in (1.0001, 1.01, 1.4, 2.25, 4.0):
        exchange = u*math.log1p(u-1)/(u-1) - 1
        for q in (0.2, 1, 3.5, 8):
            prediction = (u-1)/u * exchange/(u+q)
            close(cut_density(q, 1, u), prediction, 1e-10, "light spectral exchange")

    # Buchmüller-Plümacher Eqs. (17)-(18), in their normalization b_k=-I/(16 pi^2).
    for b, t in ((1, 1), (1, 8), (16, 9), (400, 4)):
        imag_b = -integrate(lambda z: (1-z)/(b+t*z))/(16*math.pi)
        reference = (1 - (1+b/t)*math.log1p(t/b))/(16*math.pi*t)
        close(imag_b, reference, 5e-12, "massless scalar, positive sterile invariant")
        close(triangle(b, b, a=0), 1/b, 5e-12, "massless scalar, spacelike diagonal")

    mu, el = BLOCKS.values()
    threshold = {"muon": vertex(1, mu)*tree(1, el)/(2*math.pi),
                 "electron": tree(1, mu)*vertex(1, el)/(2*math.pi)}
    return {"mass_unit": "m=1", "spin_phase": "B34* B12=s",
            "couplings_removed": {"A2_vertex_mu": "gamma |g|^2",
                                  "A2_vertex_e": "gamma |y|^2"},
            "energies_s_over_m2": values,
            "A2_vertex_over_beta_at_threshold":
                {key: float(f"{v:.12g}") for key, v in threshold.items()},
            "checks": ["Feynman parameters vs scalar-leg dispersion at all scan energies",
                       "light scalar-leg spectral density vs Candidate.22 three-particle exchange",
                       "zero scalar mass vs Buchmuller-Plumacher Eqs. (17)-(18) and I(b,b)=1/b"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    receipt = json.dumps(result(), indent=2) + "\n"
    if args.output:
        args.output.write_text(receipt, encoding="utf-8")
    else:
        print(receipt, end="")
