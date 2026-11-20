#!/usr/bin/env python3
"""Exact controls for the selected U(1) source-retaining representative.

Run with Python and SymPy. No repository code is imported. The displayed
analytic residual certificate is the derivation; sampled exact jets are an
implementation audit, not a substitute for a full EFT matching calculation.
"""

import json
import random
from collections import defaultdict

import sympy as s

I = s.I
Q = s.Rational
eta = (1, -1, -1, -1)
zero2 = s.zeros(2)
one2 = s.eye(2)
pauli = (
    s.Matrix([[0, 1], [1, 0]]),
    s.Matrix([[0, -I], [I, 0]]),
    s.diag(1, -1),
)
gamma = [s.diag(1, 1, -1, -1)]
gamma.extend(zero2.row_join(x).col_join((-x).row_join(zero2)) for x in pauli)
one4 = s.eye(4)
sigma = [[I * (gamma[a] * gamma[b] - gamma[b] * gamma[a]) / 2
          for b in range(4)] for a in range(4)]
rng = random.Random(20261015)
counts = defaultdict(int)
mutants = set()


def check(label, left, right):
    diff = left - right
    if isinstance(diff, s.MatrixBase):
        good = all(s.expand(x) == 0 for x in diff)
    else:
        good = s.expand(diff) == 0
    if not good:
        raise RuntimeError(f"{label} failed: {diff}")
    counts[label] += 1


def witness(name, difference):
    if isinstance(difference, s.MatrixBase):
        differs = any(s.expand(x) != 0 for x in difference)
    else:
        differs = s.expand(difference) != 0
    if differs:
        mutants.add(name)


def integer():
    return s.Integer(rng.randint(-3, 3))


def cmatrix(n, m):
    return s.Matrix(n, m, lambda a, b: integer() + I * integer())


def vector():
    return [integer() for _ in range(4)]


def dot(a, b):
    return sum(eta[mu] * a[mu] * b[mu] for mu in range(4))


def slash(a):
    return sum((eta[mu] * a[mu] * gamma[mu] for mu in range(4)), s.zeros(4))


def comm(a, b):
    return a * b - b * a


def v_kernel(a, k, e):
    return dot(a, e) * slash(k) - dot(a, k) * slash(e)


def e_kernel(a, k, e):
    return (slash(k) * slash(e) * slash(a) - dot(k, e) * slash(a)
            - dot(a, e) * slash(k) + dot(k, a) * slash(e))


for mu in range(4):
    for nu in range(4):
        check("Clifford anticommutator", gamma[mu] * gamma[nu] + gamma[nu] * gamma[mu],
              2 * eta[mu] * one4 if mu == nu else s.zeros(4))

# Symmetry of the chiral spinor contraction matrix gives antisymmetry in
# flavor after one Grassmann exchange. This checks the ordered-pair factor.
charge_c = I * gamma[2] * gamma[0]
gamma5 = I * gamma[0] * gamma[1] * gamma[2] * gamma[3]
pr = (one4 + gamma5) / 2
for mu in range(4):
    for nu in range(mu + 1, 4):
        tensor_spin = pr * sigma[mu][nu] * charge_c * pr.T
        check("Grassmann tensor flavor antisymmetry", tensor_spin.T, tensor_spin)

# This verifies the required derivative transfer from explicit Fourier
# kernels, independently of the operator coefficients in the residual audit.
for _ in range(16):
    r, h, k, e = vector(), vector(), vector(), vector()
    p = [r[a] + h[a] + k[a] for a in range(4)]
    d = comm(slash(k), slash(e)) * slash(r) / 4 + v_kernel(k, k, e) / 3
    q_an = -comm(slash(k), slash(e)) * slash(p)
    higgs = -(v_kernel(h, k, e) + e_kernel(h, k, e)) / 2
    reduced = higgs - v_kernel(k, k, e) / 6 - q_an / 4
    check("explicit one-B IBP kernel", d, reduced)
    witness("opposite Q_AN sign", d - (higgs - v_kernel(k, k, e) / 6 + q_an / 4))
    witness("opposite divergence sign", d - (higgs + v_kernel(k, k, e) / 6 - q_an / 4))


required_mutants = {
    "opposite Q_AN sign", "opposite divergence sign", "opposite row curvature sign",
    "isolated -1/4 used for combined dipole", "drop combined sterile tensor source",
    "opposite sourced current sign", "commute retained mass through flavor matrix",
    "replace complex alternating mass tower by M cubed",
}

# Two active and three retained sterile flavors. All sources and residuals
# are nonzero independent jets. Constant M is positive diagonal in this audit.
# Matrix columns label sterile flavor; M acts on them as Z M^T.
for trial in range(8):
    na, ns = 2, 3
    F = cmatrix(na, ns)
    M = s.diag(2, 3, 5)
    gy = Q(-2, 3)
    gp = Q(4, 3)
    yl = Q(-1, 2)
    m2, lam, rho, phi, Sphi, lap_phi = (integer() for _ in range(6))
    phi += 5
    L, SLbar, KSLbar = (cmatrix(na, 4) for _ in range(3))
    NL, NR, ANL, ANR, A2NL, A2NR, A3NL = (cmatrix(4, ns) for _ in range(7))
    SR, SL, ASR, ASL, A2SR = (cmatrix(4, ns) for _ in range(5))
    dL = [cmatrix(na, 4) for _ in range(4)]
    B = s.zeros(4)
    for mu in range(4):
        for nu in range(mu + 1, 4):
            B[mu, nu] = integer()
            B[nu, mu] = -B[mu, nu]
    sigB = sum((sigma[mu][nu] * B[mu, nu] for mu in range(4) for nu in range(4)),
               s.zeros(4))
    symmetric = {}
    for mu in range(4):
        for nu in range(mu, 4):
            symmetric[mu, nu] = cmatrix(na, 4)
    d2L = [[symmetric[min(mu, nu), max(mu, nu)] - I * gy * B[mu, nu] * L / 2
            for nu in range(4)] for mu in range(4)]
    for mu in range(4):
        for nu in range(mu + 1, 4):
            check("row derivative commutator", d2L[mu][nu] - d2L[nu][mu],
                  -I * gy * B[mu, nu] * L)
    KL = sum((-I * dL[mu] * gamma[mu] for mu in range(4)), s.zeros(na, 4))
    lapL = -sum((eta[mu] * d2L[mu][mu] for mu in range(4)), s.zeros(na, 4))
    K2L = -sum((d2L[mu][nu] * gamma[nu] * gamma[mu]
                for mu in range(4) for nu in range(4)), s.zeros(na, 4))
    check("row curvature from derivative jets", K2L, lapL - gy * L * sigB / 2)
    witness("opposite row curvature sign", K2L - (lapL + gy * L * sigB / 2))

    def mass(z):
        return z * M.T

    def b(l, z, middle=one4, factor=1):
        contraction = l * middle * z
        return factor * sum(F[a, i] * contraction[a, i] for a in range(na) for i in range(ns))

    R = mass(NR) + SR
    T = mass(mass(NL)) + mass(SL) + ASR
    U = mass(mass(mass(NR))) + mass(mass(SR)) + mass(ASL) + A2SR
    er = ANL - mass(NR) - SR
    el = ANR - mass(NL) - SL
    Aer = A2NL - mass(ANR) - ASR
    Ael = A2NR - mass(ANL) - ASL
    A2er = A3NL - mass(A2NR) - A2SR
    delta2 = mass(el) + Aer
    delta3 = mass(mass(er)) + mass(Ael) + A2er
    check("sterile first residual", ANL - R, er)
    check("sterile second residual", A2NL - T, delta2)
    check("sterile third residual", A3NL - U, delta3)

    ebar = KL - SLbar
    Kebar = K2L - KSLbar
    ephi = lap_phi - (m2 + 2 * lam * rho) * phi + Sphi
    divB, Jy = vector(), vector()
    eb = [divB[mu] - gp * Jy[mu] for mu in range(4)]
    qpair = integer() + I * integer()
    scalar_coefficient = 2 * m2 / 3 + lam * rho / 3

    # Route A starts at all six original O_J placements and the exact
    # independently checked IBP form of O_D. No EOM source is discarded.
    raw = (Q(5, 6) * b(L, A3NL, factor=phi)
           + Q(1, 3) * b(lapL, ANL, factor=phi)
           + Q(1, 6) * b(L, ANL, factor=lap_phi)
           - Q(1, 6) * b(KL, A2NL, factor=phi)
           + m2 / 2 * b(L, ANL, factor=phi)
           + m2 / 2 * b(KL, NL, factor=phi)
           - qpair / 2 - gy / 6 * b(L, NL, slash(divB), phi)
           - gy / 4 * b(L, ANL, sigB, phi))

    # Route B is the full selected source packet written independently.
    red = (Q(5, 6) * b(L, U, factor=phi)
           + scalar_coefficient * b(L, R, factor=phi)
           + Q(1, 3) * b(KSLbar, R, factor=phi)
           - Q(1, 6) * b(L, R, factor=Sphi)
           - Q(1, 6) * b(SLbar, T, factor=phi)
           + m2 / 2 * b(SLbar, NL, factor=phi)
           - qpair / 2 - gp**2 * yl / 6 * b(L, NL, slash(Jy), phi)
           - gy / 12 * b(L, R, sigB, phi))

    residual = (Q(1, 3) * b(Kebar, ANL, factor=phi)
                + Q(1, 6) * b(L, ANL, factor=ephi)
                - Q(1, 6) * b(ebar, A2NL, factor=phi)
                + m2 / 2 * b(ebar, NL, factor=phi)
                - gy / 6 * b(L, NL, slash(eb), phi)
                + Q(5, 6) * b(L, delta3, factor=phi)
                + scalar_coefficient * b(L, er, factor=phi)
                + Q(1, 3) * b(KSLbar, er, factor=phi)
                - Q(1, 6) * b(L, er, factor=Sphi)
                - Q(1, 6) * b(SLbar, delta2, factor=phi)
                - gy / 12 * b(L, er, sigB, phi))
    check("full off-shell residual certificate", raw - red, residual)

    # Individual source descendants must remain. A nonzero exact witness
    # for each deletion is required; no source-equals-zero specialization.
    descendants = {
        "cubic M squared S_R": Q(5, 6) * b(L, mass(mass(SR)), factor=phi),
        "cubic M A S_L": Q(5, 6) * b(L, mass(ASL), factor=phi),
        "cubic A squared S_R": Q(5, 6) * b(L, A2SR, factor=phi),
        "scalar coefficient S_R": scalar_coefficient * b(L, SR, factor=phi),
        "row source M N_R": Q(1, 3) * b(KSLbar, mass(NR), factor=phi),
        "row source S_R": Q(1, 3) * b(KSLbar, SR, factor=phi),
        "scalar source M N_R": -Q(1, 6) * b(L, mass(NR), factor=Sphi),
        "scalar source S_R": -Q(1, 6) * b(L, SR, factor=Sphi),
        "bar source M squared N_L": -Q(1, 6) * b(SLbar, mass(mass(NL)), factor=phi),
        "bar source M S_L": -Q(1, 6) * b(SLbar, mass(SL), factor=phi),
        "bar source A S_R": -Q(1, 6) * b(SLbar, ASR, factor=phi),
        "scalar mass bar source": m2 / 2 * b(SLbar, NL, factor=phi),
    }
    for name, term in descendants.items():
        required_mutants.add("drop " + name)
        witness("drop " + name, term)
    witness("isolated -1/4 used for combined dipole", gy / 6 * b(L, R, sigB, phi))
    witness("drop combined sterile tensor source", gy / 12 * b(L, SR, sigB, phi))
    witness("opposite sourced current sign", gp**2 * yl / 3 * b(L, NL, slash(Jy), phi))

    # Right multiplication is tested with rectangular F and arbitrary fields.
    check("retained mass flavor order", b(L, mass(NR)),
          sum((F * M)[a, i] * (L * NR)[a, i] for a in range(na) for i in range(ns)))

# A real sterile rotation checks the ordered mass products without staying
# in a diagonal basis. The same identities need no inverse mass.
O = s.Matrix([[Q(3, 5), Q(4, 5), 0], [-Q(4, 5), Q(3, 5), 0], [0, 0, 1]])
M = s.diag(2, 3, 5)
F = cmatrix(2, 3)
Mp, Fp = O.T * M * O, F * O
check("orthogonal sterile basis", O.T * O, s.eye(3))
for power in (1, 2, 3):
    check("ordered mass basis covariance", Fp * Mp**power, F * M**power * O)
Fsquare = cmatrix(3, 3)
witness("commute retained mass through flavor matrix", Fsquare * M - M * Fsquare)

# For a complex symmetric chiral mass convention, alternating products are
# required by the paired equation. This is a convention control, not an
# asserted generalization of the full field-source basis.
Mc = s.Matrix([[1 + I, 2 - I], [2 - I, 3 + 2 * I]])
check("complex Majorana mass symmetric", Mc, Mc.T)
witness("replace complex alternating mass tower by M cubed", Mc * Mc.conjugate().T * Mc - Mc**3)

# Source tensor normalization and antisymmetry, including the actual selected F.
Y = s.Matrix([[1, 0, I], [0, 2, 1]])
C = s.Matrix([[1 + I, 2 - I], [2 - I, 3]])
F = C * Y.conjugate()
H = Y * Y.conjugate().T
K = (F * Y.T - Y * F.T) / 2
check("combined active source flavor order", K, (C * H.T - H * C) / 2)
check("active source antisymmetry", K.T, -K)
t = s.Symbol("t")
tensor = s.Matrix([[0, t], [-t, 0]])
ordered = sum(K[a, b] * tensor[a, b] for a in range(2) for b in range(2))
check("ordered versus independent tensor pairs", ordered, 2 * K[0, 1] * t)

if required_mutants - mutants:
    raise RuntimeError("Unrejected mutations: " + repr(sorted(required_mutants - mutants)))

print(json.dumps({
    "status": "passed",
    "arithmetic": "exact Gaussian rational matrices and polynomial expressions",
    "full_residual_trials": 8,
    "active_flavors": 2,
    "retained_sterile_flavors": 3,
    "kernel_trials": 16,
    "checks_by_category": dict(counts),
    "negative_controls_rejected": sorted(mutants),
    "negative_control_count": len(mutants),
    "scope": "Selected U(1) source-retaining representative; not a physical beta function or minimal basis",
}, indent=2, sort_keys=True))
