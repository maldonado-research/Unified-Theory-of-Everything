#!/usr/bin/env python3
"""Independent exact symbolic audit; imports no repository/producer code.

Uses SymPy for rational polynomial identities in a Weyl gamma representation.
This checks the declared local lift, source identities, and phase conventions;
it does not establish a complete operator basis or perform another loop integral.
"""
import json
import random
import sympy as s

I = s.I
eta = (1, -1, -1, -1)
zero2, one2 = s.zeros(2), s.eye(2)
pauli = [s.Matrix([[0, 1], [1, 0]]),
         s.Matrix([[0, -I], [I, 0]]), s.diag(1, -1)]
gamma = [zero2.row_join(one2).col_join(one2.row_join(zero2))]
gamma += [zero2.row_join(x).col_join((-x).row_join(zero2)) for x in pauli]
one, zero = s.eye(4), s.zeros(4)
g5 = I * gamma[0] * gamma[1] * gamma[2] * gamma[3]
pl = (one - g5) / 2
c = I * gamma[2] * gamma[0]
checks = []


def eq(a, b):
    if isinstance(a, s.MatrixBase):
        return (a - b).applyfunc(s.expand) == s.zeros(*a.shape)
    return s.expand(a - b) == 0


def check(value, label):
    if not value:
        raise RuntimeError(label)
    checks.append(label)


def slash(v):
    return sum((eta[a] * v[a] * gamma[a] for a in range(4)), zero)


def dot(v, w):
    return sum(eta[a] * v[a] * w[a] for a in range(4))


def sigma(a, b):
    return I * (gamma[a] * gamma[b] - gamma[b] * gamma[a]) / 2


def structure(v, k, eps):
    return dot(v, eps) * slash(k) - dot(v, k) * slash(eps)


for a in range(4):
    for b in range(4):
        if not eq(gamma[a] * gamma[b] + gamma[b] * gamma[a],
                  2 * (eta[a] if a == b else 0) * one):
            raise RuntimeError("Clifford algebra")
check(True, "Weyl Clifford algebra")

r = s.symbols("r0:4", real=True)
k = s.symbols("k0:4", real=True)
eps = s.symbols("e0:4", real=True)
h = s.symbols("h0:4", real=True)
p = tuple(r[a] + h[a] + k[a] for a in range(4))
f = [[I * eta[a] * eta[b] * (k[a] * eps[b] - k[b] * eps[a])
      for b in range(4)] for a in range(4)]
comm = slash(k) * slash(eps) - slash(eps) * slash(k)
qsigma = sum((I * (I * eta[m] * r[m]) * sigma(a, b) * gamma[m] * f[a][b]
              for m in range(4) for a in range(4) for b in range(4)), zero)
qdiv = sum((gamma[m] * (I * k[n]) * f[n][m]
            for m in range(4) for n in range(4)), zero)
check(eq(qsigma, comm * slash(r)), "literal Fourier kernel of Q_sigma")
check(eq(qdiv, structure(k, k, eps)), "literal Fourier kernel of Q_div")
target = comm * slash(r) / 4 + structure(k, k, eps) / 3
check(eq(qsigma / 4 + qdiv / 3, target), "exact D = Q_sigma/4 + Q_div/3")
check(eq(target.subs(dict(zip(eps, k))), zero), "symbolic Ward transversality")

vrow = sum(((I * r[n]) * gamma[m] * f[m][n]
            for m in range(4) for n in range(4)), zero)
vh = sum((gamma[m] * (I * h[n]) * f[m][n]
          for m in range(4) for n in range(4)), zero)
vn = sum((gamma[m] * (-I * p[n]) * f[m][n]
          for m in range(4) for n in range(4)), zero)
qsource = sum((slash(r) * sigma(a, b) * f[a][b]
               for a in range(4) for b in range(4)), zero)
check(eq(qsigma, -qsource - 4 * vrow), "ordered row-EOM Clifford identity")
check(eq(vrow + vh + vn - qdiv, zero), "IBP with outgoing B and incoming N")
check(eq(target, vh + vn - s.Rational(2, 3) * qdiv - qsource / 4),
      "combined sourced/IBP reduction")
check(not eq(target, vh + vn + s.Rational(2, 3) * qdiv - qsource / 4),
      "wrong divergence sign mutation rejected")
check(not eq(qsigma, -qsource), "naive rightmost slash-r EOM deletion rejected")

# Convention epsilon^{0123}=+1, hence epsilon_{0123}=-1 for (+---).
# The dual here has both lower indices, as in gamma^mu D^nu H dualB_mu_nu.
dual_f = [[sum((-s.LeviCivita(a, b, u, v) * eta[u] * eta[v] * f[u][v] / 2
                for u in range(4) for v in range(4)), s.Integer(0))
           for b in range(4)] for a in range(4)]
qhdual = sum((gamma[m] * I * h[n] * dual_f[m][n]
              for m in range(4) for n in range(4)), zero)
qan = sum((sigma(a, b) * f[a][b] * slash(p)
           for a in range(4) for b in range(4)), zero)
higgs_form = (vh + I * qhdual) / 2 - qdiv / 6 - qan / 4
check(eq(target * pl, higgs_form * pl),
      "Higgs-derivative/dual/sterile-EOM representation")
check(not eq(target * pl, ((vh - I * qhdual) / 2 - qdiv / 6 - qan / 4) * pl),
      "opposite dual sign mutation rejected")

# Independent stacked matrix evaluations certify that the classifier's
# 19 declared coefficient columns are independent and Ward rank is 14.
# Exact nonzero minors are witnessed by rank: this is not a random proof
# of equality, but a rigorous lower bound saturated by the stated sizes.
def grade3(a, nu, b):
    return (slash(a) * gamma[nu] * slash(b) - a[nu] * slash(b)
            - b[nu] * slash(a) + dot(a, b) * gamma[nu])


def declared_columns(pv, rv, kv, m2):
    cols = []
    for scalar in (dot(pv, pv), dot(rv, rv), dot(kv, kv),
                   dot(pv, rv), dot(pv, kv), dot(rv, kv)):
        cols.append([scalar * gamma[nu] * pl for nu in range(4)])
    for a in (pv, rv, kv):
        for b in (pv, rv, kv):
            cols.append([a[nu] * slash(b) * pl for nu in range(4)])
    for a, b in ((pv, rv), (kv, pv), (kv, rv)):
        cols.append([grade3(a, nu, b) * pl for nu in range(4)])
    cols.append([m2 * gamma[nu] * pl for nu in range(4)])
    return cols


rng = random.Random(161204527)
uncontracted_rows, ward_rows = [], []
for sample in range(8):
    pv, rv, kv = [tuple(rng.randint(-3, 3) for _ in range(4)) for _ in range(3)]
    cols = declared_columns(pv, rv, kv, sample + 1)
    # In this Weyl convention all odd chains times P_L have only lower-left block.
    for nu in range(4):
        for a in (2, 3):
            for b in (0, 1):
                uncontracted_rows.append([column[nu][a, b] for column in cols])
    wcols = [sum((eta[nu] * kv[nu] * column[nu] for nu in range(4)), zero)
             for column in cols]
    for a in (2, 3):
        for b in (0, 1):
            ward_rows.append([column[a, b] for column in wcols])
check(s.Matrix(uncontracted_rows).rank() == 19,
      "independent matrix evaluation certifies all 19 declared columns")
check(s.Matrix(ward_rows).rank() == 14,
      "independent matrix contraction certifies Ward rank 14")

# A finite exact witness of rank two; no full-basis inference is made.
rv, kv = (2, 1, 0, 1), (1, 0, 1, 0)
columns = []
for kind in (0, 1):
    values = []
    for nu in range(4):
        ev = tuple(int(a == nu) for a in range(4))
        m = ((slash(kv) * slash(ev) - slash(ev) * slash(kv)) * slash(rv)
             if kind == 0 else structure(kv, kv, ev))
        values.extend(list(m))
    columns.append(s.Matrix(values))
check(s.Matrix.hstack(*columns).rank() == 2,
      "two declared one-B tensors linearly independent off shell")

# Massless left-row projection preserves the nonzero field-strength term.
rs, ks, es = (1, 0, 0, 1), (s.Rational(1, 2), 0, 0, -s.Rational(1, 2)), (0, 1, 0, 0)
ds = (slash(ks) * slash(es) - slash(es) * slash(ks)) * slash(rs) / 4
check(eq(slash(rs) * (ds - structure(rs, ks, es)), zero),
      "massless-row projected commutator identity")
check(not eq(slash(rs) * ds * pl, zero),
      "nonzero projected remainder survives free-row EOM")

# Charge conjugation fixes the relative gauge charge and preserves S(q)
# when it is accompanied by the required momentum reversal.
q = s.symbols("q0:4", real=True)
g, y = s.symbols("g y", real=True)
for a in range(4):
    if not eq(c * gamma[a].T * c.inv(), -gamma[a]):
        raise RuntimeError("C gamma transpose")
    if not eq(c * (-I * g * y * gamma[a]).T * c.inv(), I * g * y * gamma[a]):
        raise RuntimeError("charge-conjugated gauge vertex")
check(True, "charge-conjugated gauge charge y_chi=-y_L")
check(eq(c * (I * slash(tuple(-a for a in q))).T * c.inv(), I * slash(q)),
      "charge-conjugated propagator with reversed momentum")
check(eq(I * I * (-I) * I, -1), "zero-bubble vertex/propagator phase")
check(eq((-I) * I, 1), "gauge insertion vertex/extra-propagator phase")

# Vary a scalar kinetic term, treating field and its conjugate as independent.
H, Hb, dH, dHb, B, yh = s.symbols("H Hb dH dHb B yh")
dh = dH + I * g * yh * B * H
dhb = dHb - I * g * yh * B * Hb
current_h = I * yh * (Hb * dh - dhb * H)
check(eq(s.diff(dhb * dh, B), -g * current_h),
      "scalar kinetic variation gives deltaL/deltaB=-g*j_H")
check(not eq(s.diff(dhb * dh, B), g * current_h),
      "opposite scalar-current sign mutation rejected")
Jq, Ju, Jd, Jl, Je = s.symbols("Jq Ju Jd Jl Je")
jf = Jq / 6 + 2 * Ju / 3 - Jd / 3 - Jl / 2 - Je
lf = -g * B * jf
check(eq(s.diff(lf, B), -g * jf), "all SM hypercharged chiral fermion currents")
check(s.diff(jf, Jq) != 0 and s.diff(jf, Ju) != 0 and s.diff(jf, Jd) != 0,
      "omitting quark currents changes gauge source")

# SU(2) and Grassmann/charge-conjugation identities for source specializations.
ha, hb, hac, hbc = s.symbols("ha hb hac hbc")
Hrow = s.Matrix([[hac, hbc]])
htilde = s.Matrix([hbc, -hac])
check(eq((Hrow * htilde)[0], 0), "charged-lepton source Hdagger*tildeH vanishes")
for a in range(4):
    for b in range(4):
        m = c * sigma(a, b) * pl
        if not eq(m.T, m):
            raise RuntimeError("Majorana tensor flavor antisymmetry")
check(True, "C*sigma*PL symmetric implies Majorana flavor antisymmetry")
Y = s.Matrix([[1 + I, 2], [3, 1 - I], [-I, 4 + I]])
Cseed = s.Matrix([[1 + I, 2 - I, 3], [2 - I, -1, I], [3, I, 2]])
Fseed = Cseed * Y.conjugate()
Asource = Y.conjugate().T * Fseed
check(eq(Asource, Asource.T), "Ydagger*C_sym*Ystar is symmetric")
Fa = s.Matrix([[I, 2], [3, 1], [0, 4 + I]])
Aa = Y.conjugate().T * Fa
check(not eq(Aa, Aa.T), "arbitrary flavor tensor need not kill Majorana source")

# Sterile EOM instead yields an active-composite tensor, generically nonzero.
# barX sigma Xc has a symmetric spinor matrix and is Grassmann-antisymmetric
# under exchanging its two active flavor labels.
pr = one - pl
for a in range(4):
    for b in range(4):
        m = pr * sigma(a, b) * c * pr.T
        if not eq(m.T, m):
            raise RuntimeError("active-composite tensor antisymmetry")
check(True, "active barX*sigma*Xc bilinear is flavor antisymmetric")
y2 = s.diag(1, 2)
c2 = s.Matrix([[0, 1], [1, 0]])
h2 = y2 * y2.conjugate().T
k2 = (c2 * h2.T - h2 * c2) / 2
check(eq(k2, s.Matrix([[0, s.Rational(3, 2)], [-s.Rational(3, 2), 0]])),
      "nonzero two-flavor active-source counterexample")

r12 = s.Matrix([[s.Rational(3, 5), s.Rational(4, 5), 0],
                [-s.Rational(4, 5), s.Rational(3, 5), 0], [0, 0, 1]])
r23 = s.Matrix([[1, 0, 0], [0, s.Rational(5, 13), s.Rational(12, 13)],
                [0, -s.Rational(12, 13), s.Rational(5, 13)]])
u = s.diag(1, I, (3 + 4 * I) / 5) * r12 * r23
check(eq(u.conjugate().T * u, s.eye(3)), "exact complex unitary flavor rotation")
chat = s.Matrix([[1 + I, 2 - I, 3], [2 - I, 4, 1 + I], [3, 1 + I, 5 - I]])
cf = u * chat * u.T
lam = (0, 1, 4)
yf = u[:, 1:3] * s.diag(1, 2)
hf = yf * yf.conjugate().T
ff = cf * yf.conjugate()
ka = (ff * yf.T - yf * ff.T) / 2
kformula = (cf * hf.T - hf * cf) / 2
check(eq(ka, kformula), "active source FYtranspose yields C*Htranspose-H*C")
check(eq(ka.T, -ka), "active source coefficient is skew symmetric")
kh = u.conjugate().T * ka * u.conjugate()
kh_expected = s.Matrix(3, 3, lambda a, b: (lam[b] - lam[a]) * chat[a, b] / 2)
check(eq(kh, kh_expected), "complex unitary congruence gives spectral differences")
wrong_transpose = (cf * hf - hf * cf) / 2
check(not eq(ka, wrong_transpose), "H rather than Htranspose mutation rejected")
wrong_f = cf * yf
wrong_k = (wrong_f * yf.T - yf * wrong_f.T) / 2
check(not eq(ka, wrong_k), "Y rather than Ystar mutation rejected")
norm = lambda a: s.expand(sum(s.conjugate(x) * x for x in a))
spectral_norm = sum((lam[b] - lam[a]) ** 2 * s.conjugate(chat[a, b]) * chat[a, b] / 4
                    for a in range(3) for b in range(3))
check(eq(norm(ka), spectral_norm), "Frobenius norm equals spectral-gap sum")

def symmetric_to_skew_rank(eigenvalues):
    n = len(eigenvalues)
    pairs = [(a, b) for a in range(n) for b in range(a, n)]
    outpairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    jacobian = s.Matrix(len(outpairs), len(pairs),
                        lambda out, inp: ((eigenvalues[outpairs[out][1]] - eigenvalues[outpairs[out][0]]) / s.Integer(2)
                                          if outpairs[out] == pairs[inp] else 0))
    return jacobian.rank(), len(pairs) - jacobian.rank()

for eigenvalues, multiplicities in [((0, 0, 4), (2, 1)),
                                    ((0, 1, 1), (1, 2)),
                                    ((0, 0, 0), (3,)),
                                    ((2, 2, 2), (3,)),
                                    ((0, 1, 4), (1, 1, 1))]:
    rank_value, nullity = symmetric_to_skew_rank(eigenvalues)
    n = len(eigenvalues)
    check(rank_value == (n * n - sum(d * d for d in multiplicities)) // 2
          and nullity == sum(d * (d + 1) // 2 for d in multiplicities),
          "kernel/rank with eigenvalue blocks " + str(eigenvalues))
deg_h = u * s.diag(0, 0, 4) * u.conjugate().T
allowed_chat = s.Matrix([[1, 2 + I, 0], [2 + I, 3, 0], [0, 0, 5]])
allowed_c = u * allowed_chat * u.T
check(eq(allowed_c * deg_h.T - deg_h * allowed_c, s.zeros(3)),
      "nontrivial null-eigenspace symmetric block belongs to kernel")
forbidden_chat = allowed_chat + s.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]])
forbidden_c = u * forbidden_chat * u.T
check(not eq(forbidden_c * deg_h.T - deg_h * forbidden_c, s.zeros(3)),
      "coupling unequal eigenspaces exits source kernel")

print(json.dumps({"status": "PASS", "checks": checks, "check_count": len(checks),
                  "dependencies": {"sympy": s.__version__},
                  "scope": "exact local-operator and declared-phase checks; no loop recomputation, full basis, beta function, or observable"}, indent=2))
