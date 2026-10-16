#!/usr/bin/env python3
"""Exact action/Fourier/charge-conjugation audit; Python standard library only.

The previously evaluated R_L and R_H pole tensors are INPUTS, not independently
integrated here. New checks derive their physical prefactors and local-action
vertices, including the charge-conjugate line and its ordered gamma chain.
Explicit failures remain active under python -O. No repository module is read.
"""
from fractions import Fraction as F
from dataclasses import dataclass
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
    def coerce(x):
        return x if isinstance(x, Z) else Z(x)

    def __add__(self, x):
        x = Z.coerce(x)
        return Z(self.re + x.re, self.im + x.im)

    __radd__ = __add__

    def __neg__(self):
        return Z(-self.re, -self.im)

    def __sub__(self, x):
        return self + -Z.coerce(x)

    def __mul__(self, x):
        x = Z.coerce(x)
        return Z(self.re*x.re - self.im*x.im,
                 self.re*x.im + self.im*x.re)

    __rmul__ = __mul__

    def __truediv__(self, x):
        x = Z.coerce(x)
        norm = x.re*x.re + x.im*x.im
        return self * Z(x.re/norm, -x.im/norm)


I = Z(0, 1)
eta = (1, -1, -1, -1)


def matrix(rows):
    return tuple(tuple(Z.coerce(x) for x in row) for row in rows)


zero = matrix([[0]*4 for _ in range(4)])
one = matrix([[int(a == b) for b in range(4)] for a in range(4)])
gam = (
    matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]]),
    matrix([[0, 0, 0, 1], [0, 0, 1, 0], [0, -1, 0, 0], [-1, 0, 0, 0]]),
    matrix([[0, 0, 0, -I], [0, 0, I, 0], [0, I, 0, 0], [-I, 0, 0, 0]]),
    matrix([[0, 0, 1, 0], [0, 0, 0, -1], [-1, 0, 0, 0], [0, 1, 0, 0]]),
)


def add(*ms):
    return tuple(tuple(sum((m[a][b] for m in ms), Z())
                       for b in range(4)) for a in range(4))


def sc(c, m):
    c = Z.coerce(c)
    return tuple(tuple(c*x for x in row) for row in m)


def mm(*ms):
    out = one
    for m in ms:
        out = tuple(tuple(sum((out[a][c]*m[c][b] for c in range(4)), Z())
                          for b in range(4)) for a in range(4))
    return out


def tr(m):
    return tuple(tuple(m[b][a] for b in range(4)) for a in range(4))


def sub(a, b):
    return add(a, sc(-1, b))


def dot(a, b):
    return sum(eta[j]*a[j]*b[j] for j in range(4))


def vs(*vectors):
    return tuple(sum(v[j] for v in vectors) for j in range(4))


def vm(c, v):
    return tuple(c*x for x in v)


def slash(v):
    return add(*(sc(eta[j]*v[j], gam[j]) for j in range(4)))


g5 = sc(I, mm(*gam))
PL, PR = sc(F(1, 2), sub(one, g5)), sc(F(1, 2), add(one, g5))
C = sc(I, mm(gam[2], gam[0]))
Ci = sc(-1, C)
counts = {}


def check(actual, expected, label, group):
    if actual != expected:
        raise RuntimeError("FAILED: " + label)
    counts[group] = counts.get(group, 0) + 1


def reject(actual, wrong, label, group="mutation_rejections"):
    if actual == wrong:
        raise RuntimeError("Mutation survived: " + label)
    counts[group] = counts.get(group, 0) + 1


def epsilon(a, b):
    return ((0, 1), (-1, 0))[a][b]


def seed_from_action(a, b, c, d, fa, fb, hc, hd):
    """Differentiate (1/2) barX(-Box)Xc in two fermion/two H slots.

    The antisymmetric charge-conjugation matrix cancels the Grassmann swap
    of the two fermion slots. Complementary composite momenta square equally.
    The four terms retain which composite carries the literal -Box.
    """
    return F(1, 2) * (
        epsilon(a, c)*epsilon(b, d)*dot(vs(fb, hd), vs(fb, hd))
        + epsilon(a, d)*epsilon(b, c)*dot(vs(fb, hc), vs(fb, hc))
        + epsilon(b, c)*epsilon(a, d)*dot(vs(fa, hd), vs(fa, hd))
        + epsilon(b, d)*epsilon(a, c)*dot(vs(fa, hc), vs(fa, hc)))


def propagator(q, chirality):
    return sc(I / dot(q, q), mm(chirality, slash(q)))


def vertex(charge, e, chirality, coupling=F(2, 3)):
    return sc(-I*coupling*charge, mm(slash(e), chirality))


def bubble(p, h, m2):
    return add(sc(dot(p, p)+m2-F(2, 3)*dot(p, h)+dot(h, h)/2,
                  slash(p)), sc(dot(p, p)/6-m2/2, slash(h)))


def residues(p, h, k, m2, e, reverse=False):
    """Previously audited pole tensors, contracted with generic e_mu."""
    t, v, u = vs(p, vm(-1, k)), vs(h, k), vs(vm(2, p), vm(-1, k))
    pk = vs(p, k)

    def chain(a, b):
        return mm(slash(b), slash(e), slash(a)) if reverse else mm(slash(a), slash(e), slash(b))

    RL = add(sc(F(1, 3), chain(t, t)), sc(F(1, 2), chain(t, k)),
             sc(-F(1, 6), add(chain(v, pk), chain(pk, v))),
             sc(F(1, 2), chain(k, v)),
             sc(-F(5, 6)*dot(t, t)-dot(v, v)/2+dot(v, pk)/3-m2/2,
                slash(e)))
    RH = add(sc(-dot(k, e)/6, slash(k)), sc(-dot(u, e)/6, slash(h)),
             sc(dot(h, e)/3, slash(u)),
             sc(dot(k, k)/6-m2+dot(h, u)/3-dot(h, h)/2-dot(t, t),
                slash(e)))
    return RL, RH


def green_from_action(p, r, h, k, m2, e):
    """Direct plane-wave expansion of the six placements; strip g'y_L."""
    # -D² row: +i(∂.B)+2iB.∂; -D² scalar: -i(∂.B)-2iB.∂.
    row_laplacian = I*I*(dot(k, e)+2*dot(r, e))
    higgs_laplacian = -I*I*(dot(k, e)+2*dot(h, e))
    row_dirac = -I*(-I)  # K barL = -i D_mu barL gamma^mu
    return add(sc(row_laplacian/3 + higgs_laplacian/6, slash(p)),
               sc(row_dirac*(-dot(p, p)/6+m2/2), slash(e)))


def transverse_from_action(r, k, e):
    """Q_sigma/4 + Q_div/3 with Q_sigma=i(D_rho barL)sigma B gamma^rho."""
    sigma_ke = sc(I/F(2), sub(mm(slash(k), slash(e)), mm(slash(e), slash(k))))
    fieldstrength_sigma = sc(2*I, sigma_ke)
    row_derivative = sc(I, slash(r))
    q_sigma = sc(I, mm(fieldstrength_sigma, row_derivative))
    q_div = add(sc(dot(k, e), slash(k)), sc(-dot(k, k), slash(e)))
    return add(sc(F(1, 4), q_sigma), sc(F(1, 3), q_div))


def main():
    for a in range(4):
        for b in range(4):
            check(add(mm(gam[a], gam[b]), mm(gam[b], gam[a])),
                  sc(2*eta[a]*int(a == b), one), "Clifford", "conventions")
        check(mm(C, tr(gam[a]), Ci), sc(-1, gam[a]),
              "charge conjugation of vector", "conventions")
    check(mm(C, Ci), one, "C inverse", "conventions")
    check(tr(C), sc(-1, C), "C antisymmetry", "conventions")
    check(mm(C, tr(PL), Ci), PL, "charge conjugation of scalar PL", "conventions")
    # RHS barLc P_L Nc = -L^T C^-1 P_L C barN^T. Reordering the
    # Grassmann variables supplies the second minus and transposes the matrix.
    check(tr(mm(Ci, PL, C)), PL, "Yukawa Majorana flip has no minus", "conventions")

    g, yL, yH = F(2, 3), -F(1, 2), F(1, 2)
    ychi, yphi = -yL, -yH
    phase0 = I*I*(-I)*I
    phaseL = I*I*(-I*g*ychi)*I*(-I)*I
    phaseH = I*I*(-I)*I*(-I*g*yphi)*I
    check(phase0, Z(-1), "zero-gauge phase", "phase_count")
    check(phaseL, Z(-g*ychi), "fermion-insertion phase", "phase_count")
    check(phaseH, Z(-g*yphi), "scalar-insertion phase", "phase_count")

    rng = Random(15102026)
    for sample in range(32):
        vec = lambda: tuple(F(rng.randint(-7, 7), rng.randint(1, 4)) for _ in range(4))
        p, h, k, e, q = (vec() for _ in range(5))
        r, t, m2 = vs(p, vm(-1, h), vm(-1, k)), vs(p, vm(-1, k)), F(sample-11, 3)
        if dot(q, q) == 0 or dot(vs(q, vm(-1, k)), vs(q, vm(-1, k))) == 0:
            q = (F(29), F(2), F(3), F(5))

        r0 = vs(p, vm(-1, h))
        incoming_A = I*(-I)
        outgoing_K = -I*I
        negative_box = -(-I)*(-I)
        check(negative_box, Z(1), "literal -Box positive momentum square", "zero_gauge_action")
        OJ0 = add(sc(incoming_A*(F(5, 6)*dot(p, p)+dot(r0, r0)/3+dot(h, h)/6+m2/2), slash(p)),
                  sc(outgoing_K*(-dot(p, p)/6+m2/2), slash(r0)))
        check(OJ0, bubble(p, h, m2), "zero-gauge six-operator map", "zero_gauge_action")
        check(add(sc(I*phase0, bubble(p, h, m2)), sc(I, OJ0)), zero,
              "zero-gauge pole/counterterm cancellation", "zero_gauge_action")

        # Charge-conjugate the actual ordered left-chiral line. Transposition
        # reverses the original chain, fixing seed-side q-k on the left.
        left = mm(propagator(vm(-1, q), PL), vertex(yL, e, PL),
                  propagator(vs(vm(-1, q), k), PL))
        right = mm(propagator(vs(q, vm(-1, k)), PR), vertex(ychi, e, PR),
                   propagator(q, PR))
        check(mm(C, tr(left), Ci), right, "conjugated ordered fermion line", "charge_conjugation")
        check(mm(C, tr(propagator(vm(-1, q), PL)), Ci), propagator(q, PR),
              "charge-conjugate propagator momentum", "charge_conjugation")
        s, sp = vs(p, vm(-1, q)), vs(p, vm(-1, q), vm(-1, k))
        scalar_action = I*g*yphi*(I*dot(sp, e)-(-I)*dot(s, e))
        check(I*scalar_action, -I*g*yphi*dot(vs(vm(2, p), vm(-2, q), vm(-1, k)), e),
              "scalar vertex from kinetic action", "kinetic_vertices")

        # Literal derivative seed at both attachments, with all seed momenta
        # outgoing. This also checks direct/crossed contraction weights 2:1.
        for kind, f_internal, h_internal, crossed in (
            ("fermion", vs(vm(-1, q), k), vs(q, vm(-1, p)), vs(q, vm(-1, h), vm(-1, k))),
            ("scalar", vm(-1, q), vs(q, vm(-1, t)), vs(q, vm(-1, h)))):
            check(vs(r, h, f_internal, h_internal), (F(0),)*4,
                  kind + " seed momentum conservation", "seed_action")
            for a in range(2):
                for c in range(2):
                    result = sum(seed_from_action(a, b, c, d, r, f_internal, h, h_internal)*epsilon(b, d)
                                 for b in range(2) for d in range(2))
                    expected = epsilon(a, c)*(2*dot(t, t)+dot(crossed, crossed))
                    check(result, expected, kind + " seed Fourier differentiation", "seed_action")

        RL, RH = residues(p, h, k, m2, e)
        R = sub(RH, RL)
        G = green_from_action(p, r, h, k, m2, e)
        D = transverse_from_action(r, k, e)
        check(R, add(G, D), "full pole equals action lift", "local_match")
        check(transverse_from_action(r, k, k), zero, "transverse action", "ward_controls")
        check(green_from_action(p, r, h, k, m2, k),
              sub(bubble(p, vs(h, k), m2), bubble(p, h, m2)),
              "covariant representative Ward identity", "ward_controls")

        # Multiply by the loop pole i/(16pi²epsilon) and compare action vertex i.
        pole = sc(I, add(sc(phaseL, RL), sc(phaseH, RH)))
        counterterm = sc(I*g*yL, add(G, D))
        check(add(pole, counterterm), zero, "full one-B pole/counterterm cancellation", "phase_match")
        # These mutations were formerly hidden by reducing each side to a
        # separate Ward normalization; here they must fail before that step.
        wrong_charges = sc(I, add(sc(-g*yL, RL), sc(-g*yH, RH)))
        reject(add(wrong_charges, counterterm), zero, "unconjugated internal charges")
        reject(add(pole, sc(-1, counterterm)), zero, "wrong counterterm sign")
        reject(add(pole, sc(I*g*yL, G)), zero, "missing transverse action")
        RLrev, _ = residues(p, h, k, m2, e, reverse=True)
        reject(RLrev, RL, "reversed ordered line")
        RLk, _ = residues(p, h, k, m2, k)
        RLrevk, _ = residues(p, h, k, m2, k, reverse=True)
        check(RLrevk, RLk, "reversed line remains Ward-blind", "ward_controls")

        # A prediction of this covariant lift, NOT a two-B loop calculation.
        k2, e2 = vec(), vec()
        comm1 = sub(mm(slash(k), slash(e)), mm(slash(e), slash(k)))
        comm2 = sub(mm(slash(k2), slash(e2)), mm(slash(e2), slash(k2)))
        twoB_D = sc(-g*g*yL*yL/4, add(mm(comm1, slash(e2)), mm(comm2, slash(e))))
        from_connection = sc(g*yL*I*(-I*g*yL)/4,
                             add(mm(sc(-1, comm1), slash(e2)), mm(sc(-1, comm2), slash(e))))
        check(twoB_D, from_connection, "predicted two-B transverse contact", "two_B_lift_only")
        laplacian_seagull = -(-I*g*yL)*(-I*g*yL)
        twoB_G = sc(2*(F(1, 3)+F(1, 6))*laplacian_seagull*dot(e, e2), slash(p))
        check(twoB_G, sc(g*g*yL*yL*dot(e, e2), slash(p)),
              "predicted two-B representative contact", "two_B_lift_only")

    report = {
        "status": "PASS", "checks": counts, "total_checks": sum(counts.values()),
        "kinematics": 32,
        "loop_residue_inputs": "Previously audited R_L,R_H; not reintegrated by this checker",
        "full_pole": "-i g_prime y_L (C_seed Y_nu*) (R_H-R_L)/(16 pi^2 epsilon)",
        "counterterm": "+(C_seed Y_nu*) [O_J+g_prime y_L(Q_sigma/4+Q_div/3)]/(16 pi^2 epsilon)+h.c.",
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
