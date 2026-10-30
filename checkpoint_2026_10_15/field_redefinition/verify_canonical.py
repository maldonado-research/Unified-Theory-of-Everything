#!/usr/bin/env python3
"""Exact neutral-projection/canonical-coordinate checks. No repository imports.

Checks concern a declared quadratic sector, not complete EFT matching or RG.
Run: python -I -B -O verify_canonical.py
"""
import json
from pathlib import Path
import sympy as s

I = s.I
R = s.Rational
results = []
rejected = []


def clean(a):
    if isinstance(a, s.MatrixBase):
        return a.applyfunc(s.cancel)
    return s.cancel(a)


def zero(a):
    a = clean(a)
    return a == s.zeros(*a.shape) if isinstance(a, s.MatrixBase) else a == 0


def check(name, value):
    if not value:
        raise RuntimeError(name)
    results.append(name)


def reject(name, residual):
    if zero(residual):
        raise RuntimeError("negative control unexpectedly passed: " + name)
    rejected.append(name)


def det1(a, b):
    """Coefficient of t in det(a+t*b), by determinant multilinearity.

    No inverse is used; valid also when a is singular.
    """
    total = 0
    for col in range(a.cols):
        term = a.copy()
        term[:, col] = b[:, col]
        total += term.det()
    return clean(total)


def block(a, b, c, d):
    return a.row_join(b).col_join(c.row_join(d))


# Literal Fourier projection of six covariant placements, with h=0 and B=0.
sigma = [s.Matrix([[0, 1], [1, 0]]),
         s.Matrix([[0, -I], [I, 0]]), s.diag(1, -1)]
z2 = s.zeros(2)
g = [block(z2, s.eye(2), s.eye(2), z2)]
g += [block(z2, t, -t, z2) for t in sigma]
p = s.symbols("p0:4", real=True)
m2 = s.symbols("m_H_squared", real=True)
slash = sum((g[mu] * p[mu] for mu in range(4)), s.zeros(4))
p2 = p[0]**2 - sum(x**2 for x in p[1:])
check("Clifford square determines the projected derivative order",
      zero(slash * slash - p2 * s.eye(4)))
placements = [R(5, 6) * slash * p2, R(1, 3) * p2 * slash,
              s.zeros(4), -R(1, 6) * slash * p2,
              R(1, 2) * m2 * slash, R(1, 2) * m2 * slash]
kernel = sum(placements, s.zeros(4))
check("all six O_J placements give slash p times (p squared plus m_H squared)",
      zero(kernel - slash * (p2 + m2)))
reject("discarding row-derivative placements changes the neutral kernel",
       sum((placements[i] for i in [0, 2, 4]), s.zeros(4)) - kernel)

# Actual sector: two active left Weyl fields and one retained left sterile.
# q=(nu_L, N_L), mass convention -1/2(q^T epsilon M q + h.c.).
u = R(3, 2)
mh2 = -R(2, 3)
Y = s.Matrix([1 + I, 2 - I])
Cseed = s.Matrix([[1 + 2*I, 2 - I], [2 - I, -1 + I]])
F = Cseed * Y.conjugate()
B = u * Y.conjugate()
Mref = block(s.zeros(2), B, B.T, s.Matrix([[5]]))
G = s.Matrix([[1 + I, 2, 1 - 2*I], [2 - I, -1 + I, 3]])
W = s.Matrix([[1, I, 2, 1 - I, 0, 3], [I, 2, 1, -1, 2 + I, I]])
check("selected symmetric-seed flavor orientation", Cseed.T == Cseed)
check("left-Weyl reference mass is symmetric", Mref.T == Mref)
check("two-active one-sterile reference mass is rank deficient", Mref.rank() == 2)

sample = []
for energy in [s.Integer(2), s.Integer(3), R(5, 2)]:
    z_as = clean(u * (energy**2 + mh2) * F)
    Zactual = block(s.zeros(2), z_as, z_as.H, s.zeros(1))
    # A diagonal spectator is solely an algebra control for nonzero trace.
    for label, Z in [("selected", Zactual),
                     ("generic_hermitian_control", Zactual + s.diag(R(1, 3), -R(2, 5), R(3, 7)))]:
        tag = label + " E=" + str(energy)
        n = Mref.rows
        T1 = -Z / 2
        dm = clean(T1.T * Mref + Mref * T1)
        check(tag + " Hermitian kinetic input", Z.H == Z)
        check(tag + " canonical kinetic coefficient vanishes",
              zero(T1.H + Z + T1))
        check(tag + " canonical Majorana mass stays symmetric", zero(dm.T - dm))

        N1 = s.diag(T1, T1.conjugate())
        D0 = block(energy * s.eye(n), -Mref.H, -Mref, energy * s.eye(n))
        D1 = s.diag(energy * Z, energy * Z.T)
        Dc1 = block(s.zeros(n), -dm.H, -dm, s.zeros(n))
        check(tag + " doubled inverse kernel is a congruence",
              zero(Dc1 - (N1.H * D0 + D1 + D0 * N1)))
        check(tag + " doubled pole determinant has the correct Jacobian",
              zero(det1(D0, Dc1) - det1(D0, D1) + 2 * s.trace(Z) * D0.det()))

        # Independently eliminate the conjugate Weyl equation. This produces
        # Q(E^2)=E^2*K-M^dagger*K^{-T}*M, rather than the doubled kernel.
        Q0 = energy**2 * s.eye(n) - Mref.H * Mref
        Q1 = energy**2 * Z + Mref.H * Z.T * Mref
        Qc1 = -(dm.H * Mref + Mref.H * dm)
        check(tag + " generalized squared-pole kernel transforms covariantly",
              zero(Qc1 - (T1.H * Q0 + Q1 + Q0 * T1)))
        check(tag + " squared-pole determinant has the correct Jacobian",
              zero(det1(Q0, Qc1) - det1(Q0, Q1) + s.trace(Z) * Q0.det()))

        # External-response equality tests both leg/source factors and inverse
        # propagator, rather than accepting a mass-spectrum check alone.
        inv0 = D0.inv()
        ro1 = -W * inv0 * D1 * inv0 * W.H
        W1 = W * N1
        rc1 = W1 * inv0 * W.H - W * inv0 * Dc1 * inv0 * W.H + W * inv0 * W1.H
        check(tag + " source-dressed response is invariant", zero(rc1 - ro1))
        reject(tag + " omitting source transport breaks external response",
               -W * inv0 * Dc1 * inv0 * W.H - ro1)
        reject(tag + " omitting induced mass terms changes pole determinant",
               det1(D0, D1) + (-2 * s.trace(Z) * D0.det()))

        # Flavor-coordinate freedom: add an arbitrary anti-Hermitian Omega.
        Om = s.Matrix([[I, 1 + I, 2], [-1 + I, -2*I, I], [-2, I, 0]]) / 11
        check(tag + " residual flavor generator is anti-Hermitian", Om.H == -Om)
        Trot = T1 + Om
        dmrot = clean(Trot.T * Mref + Mref * Trot)
        check(tag + " residual flavor freedom preserves canonical kinetic term",
              zero(Trot.H + Z + Trot))
        Qrot1 = -(dmrot.H * Mref + Mref.H * dmrot)
        check(tag + " residual flavor coordinates preserve pole determinant",
              zero(det1(Q0, Qrot1) - det1(Q0, Qc1)))
        reject(tag + " replacing transpose with adjoint destroys Majorana symmetry",
               (T1.H * Mref + Mref * T1).T - (T1.H * Mref + Mref * T1))

        if label == "selected":
            check(tag + " selected trace vanishes", s.trace(Z) == 0)
            # Deliberately retain only the one-derivative m_H^2 part of Z.
            Zlow = Zactual * mh2 / (energy**2 + mh2)
            dmlow = -(Zlow.T * Mref + Mref * Zlow) / 2
            Dclow1 = block(s.zeros(n), -dmlow.H, -dmlow, s.zeros(n))
            reject(tag + " discarding p squared induced mass terms changes poles",
                   det1(D0, Dclow1) - det1(D0, Dc1))

        sample.append({"sector": label, "energy": str(energy),
                       "kinetic_trace": str(s.trace(Z)), "mass_rank": Mref.rank()})

# No inverse mass is necessary for the covariance or source-response checks.
# Exact invertible congruence preserves a massless mode even for complex Z.
Z = block(s.zeros(2), u * F, u * F.H, s.zeros(1))
Tfinite = s.eye(3) - Z / 26
v = Mref.nullspace()[0]
Mfinite = clean(Tfinite.T * Mref * Tfinite)
check("finite congruence is invertible", Tfinite.det() != 0)
check("finite congruence preserves rank-deficient mass", Mfinite.rank() == Mref.rank())
check("massless null vector transports without a mass inverse",
      zero(Mfinite * Tfinite.inv() * v))

# Zero-momentum holomorphic Schur complement, explicitly assuming only that
# the retained sterile block is invertible. The active Schur block may vanish.
A = Mref[:2, :2]
B = Mref[:2, 2:]
D = Mref[2:, 2:]
Schur = clean(A - B * D.inv() * B.T)
Tri = block(s.eye(2), s.zeros(2, 1), -D.inv() * B.T, s.eye(1))
check("sterile block invertibility is explicit", D.det() != 0)
check("zero-momentum Schur block diagonalization", zero(Tri.T * Mref * Tri - s.diag(Schur, D)))
check("rank formula remains valid for a singular active Schur block",
      Mref.rank() == D.rank() + Schur.rank() and Schur.det() == 0)
Ga, Gs = G[:, :2], G[:, 2:]
Geff = clean(Ga - Gs * D.inv() * B.T)
check("Schur reduction transports the external vertex", zero((G * Tri)[:, :2] - Geff))

# A declared invertible active spectator permits a full zero-momentum response
# benchmark. It is not an additional coefficient predicted by the packet.
Af = s.Matrix([[2 + I, 1 - I], [1 - I, -1 + 2*I]])
Mf = block(Af, B, B.T, D)
Sf = clean(Af - B * D.inv() * B.T)
contact = clean(Gs * D.inv() * Gs.T)
response = clean(G * Mf.inv() * G.T)
check("invertible-spectator zero-momentum source response",
      zero(response - contact - Geff * Sf.inv() * Geff.T))
reject("omitting the induced Schur source contact changes the response",
       response - Geff * Sf.inv() * Geff.T)
reject("omitting the induced Schur source vertex changes the response",
       response - contact - Ga * Sf.inv() * Ga.T)

out = {
    "status": "PASS", "scope": "Exact quadratic neutral projection and canonical-coordinate algebra; not full EFT matching or physical running",
    "positive_checks": len(results), "rejected_negative_controls": len(rejected),
    "checks": results, "negative_controls": rejected, "samples": sample,
    "input_projection": "phi=(u,0) real constant; gauge backgrounds zero; original selected O_J+O_D before EOM reduction",
    "projected_kinetic_mixing": "Z_as(p^2)=u*(p^2+m_H^2)*F; Z_sa=Z_as^dagger; F=C_seed*Y_nu^*",
    "mass_convention": "-1/2(q^T epsilon Mcal q+h.c.), q=(nu_L,N_L)",
    "rank_deficient_reference": True,
    "novelty": "Application/check certificate of standard field-redefinition and linear-algebra identities; no new fundamental mathematics claimed"
}
Path(__file__).with_name("CANONICAL_RESULTS.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps({"status": out["status"], "positive_checks": len(results), "rejected_negative_controls": len(rejected)}))
