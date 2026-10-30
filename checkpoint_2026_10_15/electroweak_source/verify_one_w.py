#!/usr/bin/env python3
"""One-W seed pole: exact weak contractions and crossed-triangle integration.

Standard library only, no repository imports.  The rank-two/three triangle
and bubble UV moments are explicit inputs; crossed residues are integrated
from their numerators, not obtained by charge substitution in a B result.
All checks remain active with python -O.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
from random import Random
import json


@dataclass(frozen=True)
class Z:
    re: F = F(0)
    im: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, "re", F(self.re))
        object.__setattr__(self, "im", F(self.im))

    @staticmethod
    def cast(x):
        return x if isinstance(x, Z) else Z(x)

    def __add__(self, b):
        b = Z.cast(b)
        return Z(self.re+b.re, self.im+b.im)

    __radd__ = __add__

    def __neg__(self):
        return Z(-self.re, -self.im)

    def __sub__(self, b):
        return self + -Z.cast(b)

    def __mul__(self, b):
        b = Z.cast(b)
        return Z(self.re*b.re-self.im*b.im, self.re*b.im+self.im*b.re)

    __rmul__ = __mul__

    def __truediv__(self, b):
        b = Z.cast(b)
        norm = b.re*b.re+b.im*b.im
        return self*Z(b.re/norm, -b.im/norm)


I = Z(0, 1)
eta = (1, -1, -1, -1)


def mat(rows):
    return tuple(tuple(Z.cast(x) for x in row) for row in rows)


def ident(n):
    return mat([[int(a == b) for b in range(n)] for a in range(n)])


def zeros(n):
    return mat([[0]*n for _ in range(n)])


def add(*ms):
    n = len(ms[0])
    return tuple(tuple(sum((m[a][b] for m in ms), Z())
                       for b in range(n)) for a in range(n))


def scale(c, m):
    return tuple(tuple(Z.cast(c)*v for v in row) for row in m)


def mul(*ms):
    n = len(ms[0])
    out = ident(n)
    for m in ms:
        out = tuple(tuple(sum((out[a][c]*m[c][b] for c in range(n)), Z())
                          for b in range(n)) for a in range(n))
    return out


def trans(m):
    return tuple(zip(*m))


def trace(m):
    return sum((m[a][a] for a in range(len(m))), Z())


def sub(a, b):
    return add(a, scale(-1, b))


gam = (
    mat([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]]),
    mat([[0, 0, 0, 1], [0, 0, 1, 0], [0, -1, 0, 0], [-1, 0, 0, 0]]),
    mat([[0, 0, 0, -I], [0, 0, I, 0], [0, I, 0, 0], [-I, 0, 0, 0]]),
    mat([[0, 0, 1, 0], [0, 0, 0, -1], [-1, 0, 0, 0], [0, 1, 0, 0]]),
)
eps = mat([[0, 1], [-1, 0]])
T = tuple(scale(F(1, 2), m) for m in (
    mat([[0, 1], [1, 0]]), mat([[0, -I], [I, 0]]), mat([[1, 0], [0, -1]])))


def dot(a, b):
    return sum(eta[j]*a[j]*b[j] for j in range(4))


def vs(*v):
    return tuple(sum(x[j] for x in v) for j in range(4))


def vm(c, v):
    return tuple(c*x for x in v)


def slash(v):
    return add(*(scale(eta[j]*v[j], gam[j]) for j in range(4)))


def metric(a, b):
    return eta[a] * int(a == b)


def tri2(a, b):
    return F(metric(a, b), 4)


def tri3(a, b, c, shifts):
    return (metric(a, b)*shifts[c]+metric(a, c)*shifts[b]
            +metric(b, c)*shifts[a])/12


def bubble2(a, b, p, ma2, mb2):
    return p[a]*p[b]/3 + (F(1, 4)*(ma2+mb2)-dot(p, p)/12)*metric(a, b)


def integrate(p, h, k, m2, e):
    """Return per-pair Ldirect,Hdirect,Lcross,Hcross UV residues.

    Crossed numerators: (q-h-k)^2 (/q-/k)/e/q on L, and
    (q-h)^2 /q (2p-2q-k).e on H. Cancel q^2 once, then integrate
    remaining component monomials using rank<=3 tensors.
    """
    t, v, u = vs(p, vm(-1, k)), vs(h, k), vs(vm(2, p), vm(-1, k))
    E, K = slash(e), slash(k)
    lc, hc, ld, hd = (zeros(4) for _ in range(4))
    # Fermion q^2 cancellation: q=l+k; denominators l^2, (l-t)^2-m2.
    for a, b in product(range(4), repeat=2):
        lc = add(lc, scale(eta[a]*eta[b]*bubble2(a, b, t, F(0), m2),
                           mul(gam[a], E, gam[b])))
    lc = add(lc, scale(F(1, 2), mul(slash(t), E, K)))
    # Scalar q^2 cancellation: q=l+p; equal masses, shift -k.
    for a, b in product(range(4), repeat=2):
        hc = add(hc, scale(-2*eta[a]*eta[b]*e[b]
                           *bubble2(a, b, vm(-1, k), m2, m2), gam[a]))
    hc = add(hc, scale(dot(k, e)/2, K))
    # Remaining triangle polynomials, including contractions with the shift sum.
    for a, b in product(range(4), repeat=2):
        lc = add(lc, scale(dot(v, v)*eta[a]*eta[b]*tri2(a, b),
                           mul(gam[a], E, gam[b])))
        lc = add(lc, scale(2*eta[a]*v[a]*eta[b]*tri2(a, b),
                           mul(K, E, gam[b])))
        hc = add(hc, scale(-2*dot(u, e)*eta[a]*h[a]*eta[b]*tri2(a, b),
                           gam[b]))
        hc = add(hc, scale(-2*dot(h, h)*eta[a]*eta[b]*e[b]*tri2(a, b),
                           gam[a]))
        ld = add(ld, scale(dot(t, t)*eta[a]*eta[b]*tri2(a, b),
                           mul(gam[a], E, gam[b])))
        hd = add(hd, scale(-2*dot(t, t)*eta[a]*eta[b]*e[b]*tri2(a, b),
                           gam[a]))
    for a, b, c in product(range(4), repeat=3):
        lc = add(lc, scale(-2*eta[a]*v[a]*eta[b]*eta[c]
                           *tri3(a, b, c, vs(p, k)), mul(gam[b], E, gam[c])))
        hc = add(hc, scale(4*eta[a]*h[a]*eta[b]*eta[c]*e[c]
                           *tri3(a, b, c, u), gam[b]))
    return ld, hd, lc, hc


def crossed_formulas(p, h, k, m2, e):
    t, v, u = vs(p, vm(-1, k)), vs(h, k), vs(vm(2, p), vm(-1, k))
    pk, E = vs(p, k), slash(e)
    chain = lambda a, b: mul(slash(a), E, slash(b))
    lc = add(scale(F(1, 3), chain(t, t)), scale(F(1, 2), chain(t, k)),
             scale(-F(1, 6), add(chain(v, pk), chain(pk, v))),
             scale(F(1, 2), chain(k, v)),
             scale(dot(t, t)/6-dot(v, v)/2+dot(v, pk)/3-m2/2, E))
    hc = add(scale(-dot(k, e)/6, slash(k)), scale(-dot(u, e)/6, slash(h)),
             scale(dot(h, e)/3, slash(u)),
             scale(dot(k, k)/6-m2+dot(h, u)/3-dot(h, h)/2, E))
    return lc, hc


def green(p, h, k, m2, e):
    r = vs(p, vm(-1, h), vm(-1, k))
    return add(scale(-dot(vs(vm(2, r), k), e)/3
                     +dot(vs(vm(2, h), k), e)/6, slash(p)),
               scale(dot(p, p)/6-m2/2, slash(e)))


def transverse(p, h, k, e):
    r, K, E = vs(p, vm(-1, h), vm(-1, k)), slash(k), slash(e)
    return add(scale(F(1, 4), mul(sub(mul(K, E), mul(E, K)), slash(r))),
               scale(dot(k, e)/3, K), scale(-dot(k, k)/3, E))


def bubble(p, h, m2):
    return add(scale(dot(p, p)+m2-F(2, 3)*dot(p, h)+dot(h, h)/2, slash(p)),
               scale(dot(p, p)/6-m2/2, slash(h)))


counts = {}


def check(a, b, group):
    if a != b:
        raise RuntimeError("FAILED: " + group)
    counts[group] = counts.get(group, 0)+1


def reject(a, b, label):
    if a == b:
        raise RuntimeError("Mutation survived: " + label)
    counts["mutation_rejections"] = counts.get("mutation_rejections", 0)+1


def group_factors(generator, pairing, line):
    """Literal phi-basis contractions; keep four internal indices explicit."""
    out = [[Z() for _ in range(2)] for _ in range(2)]
    for a, c, b, d, e, f in product(range(2), repeat=6):
        seed = int(a == c and b == d) if pairing == "direct" else int(a == d and b == c)
        insertion = -generator[e][b]*int(d == f) if line == "L" else generator[d][f]*int(b == e)
        out[a][c] += seed*insertion*int(e == f)
    return mat(out)


def epsilon_factors(generator, pairing, line):
    """Independently contract epsilon seed in original H* scalar basis."""
    out = [[Z() for _ in range(2)] for _ in range(2)]
    anti = scale(-1, trans(generator))
    for a, c, b, d, e, f in product(range(2), repeat=6):
        seed = eps[a][c]*eps[b][d] if pairing == "direct" else eps[a][d]*eps[b][c]
        insertion = anti[b][e]*int(d == f) if line == "L" else anti[d][f]*int(b == e)
        out[a][c] += seed*insertion*eps[e][f]
    return mat(out)


def main():
    for a, b in product(range(4), repeat=2):
        check(add(mul(gam[a], gam[b]), mul(gam[b], gam[a])),
              scale(2*metric(a, b), ident(4)), "clifford")
    for i, j in product(range(3), repeat=2):
        check(trace(mul(T[i], T[j])), Z(F(int(i == j), 2)), "weak_normalization")
    for t in T:
        check(trace(t), Z(), "weak_normalization")
        check(mul(t, eps), scale(-1, mul(eps, trans(t))), "epsilon_intertwiner")
        for pair, line in product(("direct", "crossed"), ("L", "H")):
            expected = zeros(2) if pair == "direct" else scale(-1 if line == "L" else 1, t)
            check(group_factors(t, pair, line), expected, "literal_weak_contractions")
            check(epsilon_factors(t, pair, line), mul(expected, eps), "epsilon_contractions")
    # General scalar matrix generator tests also retain the U(1) trace weights.
    y = -F(1, 2)
    for gen in (scale(y, ident(2)), add(scale(F(2, 5), ident(2)), T[0], scale(F(3, 7), T[1]))):
        for pair, line in product(("direct", "crossed"), ("L", "H")):
            expected = scale(trace(gen), ident(2)) if pair == "direct" else gen
            check(group_factors(gen, pair, line), scale(-1 if line == "L" else 1, expected),
                  "trace_and_crossed_weights")
    reject(mul(T[1], eps), scale(-1, mul(eps, T[1])), "forgetting transpose")
    reject(group_factors(T[1], "crossed", "L"), T[1], "wrong conjugate-line sign")
    # Actual source-chain phases, with one gauge insertion: i*i*(-i)*i*(-i)*i=-1.
    check(I*I*(-I)*I, Z(-1), "action_phases")
    check(I*I*(-I)*I*(-I)*I, Z(-1), "action_phases")
    # Complex symmetric seed times Yukawa*, ordered in flavor space.
    cseed = mat([[1, I], [I, 2]])
    yc = mat([[Z(2, -1), 3], [Z(1, 2), Z(0, -1)]])
    f = mul(cseed, yc)
    explicit = mat([[sum((cseed[a][b]*yc[b][i] for b in range(2)), Z())
                     for i in range(2)] for a in range(2)])
    check(f, explicit, "flavor_order")
    reject(f, mul(yc, cseed), "commuted flavor matrices")
    rng = Random(15102026)
    for sample in range(24):
        vector = lambda: tuple(F(rng.randint(-7, 7), rng.randint(1, 4)) for _ in range(4))
        p, h, k, e = (vector() for _ in range(4))
        m2 = F(sample-9, 3)
        ld, hd, lc, hc = integrate(p, h, k, m2, e)
        cf, cs = crossed_formulas(p, h, k, m2, e)
        check(lc, cf, "crossed_triangle_integrations")
        check(hc, cs, "crossed_triangle_integrations")
        t = vs(p, vm(-1, k))
        check(ld, scale(-dot(t, t)/2, slash(e)), "direct_triangle_integrations")
        check(hd, ld, "direct_triangle_integrations")
        result = sub(hc, lc)
        check(result, add(green(p, h, k, m2, e), transverse(p, h, k, e)), "local_decomposition")
        check(sub(add(hc, scale(2, hd)), add(lc, scale(2, ld))), result, "hypercharge_direct_cancellation")
        for gen in T:
            gdL, gdH = group_factors(gen, "direct", "L"), group_factors(gen, "direct", "H")
            gcL, gcH = group_factors(gen, "crossed", "L"), group_factors(gen, "crossed", "H")
            for a, c in product(range(2), repeat=2):
                contracted = add(scale(gdL[a][c], ld), scale(gdH[a][c], hd),
                                 scale(gcL[a][c], lc), scale(gcH[a][c], hc))
                check(contracted, scale(gen[a][c], result), "full_open_weak_pole")
        _, _, wlc, whc = integrate(p, h, k, m2, k)
        check(sub(whc, wlc), sub(bubble(p, vs(h, k), m2), bubble(p, h, m2)), "ward_identity")
        check(transverse(p, h, k, k), zeros(4), "transversality")
        _, _, mlc, mhc = integrate(p, h, k, m2+5, e)
        check(sub(sub(mhc, mlc), green(p, h, k, m2+5, e)),
              transverse(p, h, k, e), "mass_independent_transverse_part")
        # Choose rational weak mixing s=3/5,c=4/5 and g=e/s,g'=e/c.
        sw, cw, electric = F(3, 5), F(4, 5), F(2, 7)
        weak, hyper = electric/sw, electric/cw
        photon = add(scale(cw*hyper*y, result), scale(sw*weak*T[2][0][0], result))
        check(photon, zeros(4), "neutral_photon_projection")
        zproj = add(scale(-sw*hyper*y, result), scale(cw*weak*T[2][0][0], result))
        check(zproj, scale(electric/(2*sw*cw), result), "neutral_z_projection")
        reject(result, green(p, h, k, m2, e), "dropped transverse pole")
        reject(result, add(hc, lc), "same physical sign on both lines")
        reject(photon, scale(cw*hyper*y, result), "hypercharge treated as photon")
    qphi = add(T[2], scale(y, ident(2)))
    check(qphi, mat([[0, 0], [0, -1]]), "electromagnetic_charge_matrix")
    check(trace(qphi), Z(-1), "nonzero_direct_photon_trace")
    print(json.dumps({"status": "pass", "checks": sum(counts.values()),
                      "breakdown": counts,
                      "scope": "selected simple UV residue; no finite-amplitude claim"}, indent=2))


if __name__ == "__main__":
    main()
