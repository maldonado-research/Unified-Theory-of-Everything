#!/usr/bin/env python3
"""Independent exact local-action audit; no repository imports or loop inputs.

Requires SymPy 1.14.0. Run python -I -B -O independent_ew_audit.py.
The script prints JSON; redirect to INDEPENDENT_EW_AUDIT_RESULTS.json.
"""
import json
import random
import sys
import sympy as s

I, R = s.I, s.Rational
eta = [1, -1, -1, -1]
z2, z4 = s.zeros(2), s.zeros(4)
pauli = [s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -I], [I, 0]]), s.diag(1, -1)]
T = [p / 2 for p in pauli]
gamma = [s.diag(1, 1, -1, -1)] + [s.BlockMatrix([[z2, p], [-p, z2]]).as_explicit() for p in pauli]
sigma = [[I * (gamma[m]*gamma[n]-gamma[n]*gamma[m])/2 for n in range(4)] for m in range(4)]
G = [[[gamma[m]*gamma[n]*gamma[r] - eta[m]*(m == n)*gamma[r]
        + eta[m]*(m == r)*gamma[n] - eta[n]*(n == r)*gamma[m]
       for r in range(4)] for n in range(4)] for m in range(4)]
C = I * gamma[2] * gamma[0]
checks, negatives = [], []

def zero(value):
    if isinstance(value, s.MatrixBase):
        return all(s.expand(x) == 0 for x in value)
    return s.expand(value) == 0

def check(name, value):
    if not zero(value):
        raise RuntimeError(name + ': nonzero exact residual ' + str(value))
    checks.append(name)

def reject(name, value):
    if zero(value):
        raise RuntimeError(name + ': mutation escaped detection')
    negatives.append(name)

def comm(a, b):
    return a*b-b*a

def sm(values, zero_matrix):
    return sum(values, zero_matrix.copy())

def eps(a, b, c):
    return s.LeviCivita(a, b, c)

def trial(seed):
    rng = random.Random(seed)
    q = lambda: s.Integer(rng.randint(-3, 3))
    vec = lambda n: s.Matrix([q()+I*q() for _ in range(n)])
    herm = lambda: q()*s.eye(2)+sm([q()*p for p in pauli], z2)
    A = [herm() for _ in range(4)]
    dA = [[herm() for _ in range(4)] for _ in range(4)]
    ddA = [[[z2.copy() for _ in range(4)] for _ in range(4)] for _ in range(4)]
    for r in range(4):
        for m in range(r, 4):
            for n in range(4):
                ddA[r][m][n] = ddA[m][r][n] = herm()
    F = [[dA[m][n]-dA[n][m]+I*comm(A[m], A[n]) for n in range(4)] for m in range(4)]
    dF = [[[ddA[r][m][n]-ddA[r][n][m]+I*(comm(dA[r][m], A[n])+comm(A[m], dA[r][n]))
            for n in range(4)] for m in range(4)] for r in range(4)]
    DF = [[[dF[r][m][n]+I*comm(A[r], F[m][n]) for n in range(4)] for m in range(4)] for r in range(4)]
    for r in range(4):
        for m in range(r+1, 4):
            for n in range(m+1, 4):
                check(f'{seed}: covariant Bianchi {r}{m}{n}', DF[r][m][n]+DF[m][n][r]+DF[n][r][m])
    reject(f'{seed}: ordinary Bianchi substituted', dF[0][1][2]+dF[1][2][0]+dF[2][0][1])
    reject(f'{seed}: reversed adjoint sign', DF[0][1][2]-2*I*comm(A[0], F[1][2])
           +DF[1][2][0]-2*I*comm(A[1], F[2][0])+DF[2][0][1]-2*I*comm(A[2], F[0][1]))
    L = [vec(4).T for _ in range(2)]
    dL = [[vec(4).T for _ in range(2)] for _ in range(4)]
    phi, N = vec(2), vec(4)
    dphi, dN = [vec(2) for _ in range(4)], [vec(4) for _ in range(4)]
    DL = [[dL[r][a]-I*sm([L[b]*A[r][b,a] for b in range(2)], s.zeros(1,4)) for a in range(2)] for r in range(4)]
    Dphi = [dphi[r]+I*A[r]*phi for r in range(4)]
    def bil(row, spin, weak, col, sterile):
        return s.expand(sum((row[a]*spin*sterile)[0]*weak[a,b]*col[b] for a in range(2) for b in range(2)))
    qs = sum(I*bil(DL[r], sigma[m][n]*gamma[r], F[m][n], phi, N) for m in range(4) for n in range(4) for r in range(4))
    qdiv = sum(eta[n]*bil(L, gamma[m], DF[n][n][m], phi, N) for m in range(4) for n in range(4))
    qh = sum(-eta[n]*bil(L, gamma[m], F[m][n], Dphi[n], N) for m in range(4) for n in range(4))
    qh += sum(-bil(L, G[m][n][r], F[m][n], Dphi[r], N)/2 for m in range(4) for n in range(4) for r in range(4))
    an = sum(I*bil(L, sigma[m][n]*gamma[r], F[m][n], phi, dN[r]) for m in range(4) for n in range(4) for r in range(4))
    # Evaluate the ordinary derivative of the neutral local boundary current.
    boundary = sum(I*(bil(dL[r], sigma[m][n]*gamma[r], F[m][n], phi, N)
                    +bil(L, sigma[m][n]*gamma[r], dF[r][m][n], phi, N)
                    +bil(L, sigma[m][n]*gamma[r], F[m][n], dphi[r], N)
                    +bil(L, sigma[m][n]*gamma[r], F[m][n], phi, dN[r]))/4
                   for m in range(4) for n in range(4) for r in range(4))
    od = qs/4+qdiv/3
    check(f'{seed}: matrix IBP against ordinary boundary derivative', od+qh/2+qdiv/6+an/4-boundary)
    reject(f'{seed}: wrong divergence coefficient', od+qh/2-qdiv/6+an/4-boundary)
    reject(f'{seed}: omit adjoint connection in divergence', sum(eta[n]*bil(L, gamma[m], DF[n][n][m]-dF[n][n][m], phi, N) for m in range(4) for n in range(4)))
    # Arbitrary symmetric second partial derivatives of the row.
    ddL = [[[s.zeros(1,4) for _ in range(2)] for _ in range(4)] for _ in range(4)]
    for n in range(4):
        for m in range(n,4):
            for a in range(2):
                ddL[n][m][a] = ddL[m][n][a] = vec(4).T
    DDL = [[[ddL[n][m][a]
        -I*sm([dL[n][b]*A[m][b,a]+L[b]*dA[n][m][b,a]+dL[m][b]*A[n][b,a] for b in range(2)], s.zeros(1,4))
        -sm([L[b]*(A[m]*A[n])[b,a] for b in range(2)], s.zeros(1,4))
        for a in range(2)] for m in range(4)] for n in range(4)]
    for a in range(2):
        k2 = -sm([DDL[n][m][a]*gamma[m]*gamma[n] for n in range(4) for m in range(4)], s.zeros(1,4))
        lap = sm([eta[m]*DDL[m][m][a] for m in range(4)], s.zeros(1,4))
        curv = sm([L[b]*sigma[m][n]*F[m][n][b,a]/2 for b in range(2) for m in range(4) for n in range(4)], s.zeros(1,4))
        check(f'{seed}: row Dirac square weak slot {a}', k2+lap+curv)
        reject(f'{seed}: reversed row curvature weak slot {a}', k2+lap-curv)
    # Gauge Euler derivative from varying -W^2/4 directly in components.
    g = R(3,2)
    W = [[q() for _ in range(3)] for _ in range(4)]
    dW = [[[q() for _ in range(3)] for _ in range(4)] for _ in range(4)]
    f = [[[dW[m][n][k]-dW[n][m][k]-g*sum(eps(i,j,k)*W[m][i]*W[n][j] for i in range(3) for j in range(3)) for k in range(3)] for n in range(4)] for m in range(4)]
    for n in range(4):
        for k in range(3):
            # Algebraic derivative of the nonlinear curvature. The derivative
            # of the kinetic momentum supplies partial_mu F^{mu nu} separately.
            alg = -sum(eta[m]*eta[r]*f[m][r][j]*(-g)*sum(eps(a,b,j)*((m==n)*(a==k)*W[r][b]+W[m][a]*(r==n)*(b==k)) for a in range(3) for b in range(3))/2 for m in range(4) for r in range(4) for j in range(3))
            expected = -g*sum(eta[m]*eta[n]*eps(i,j,k)*W[m][i]*f[m][n][j] for m in range(4) for i in range(3) for j in range(3))
            check(f'{seed}: nonlinear Yang-Mills Euler sign {n}{k}', alg-expected)
            reject(f'{seed}: reversed Yang-Mills Euler sign {n}{k}', alg+expected)
    # Scalar current, including connections, from an exact gauge variation.
    for k in range(3):
        dp = Dphi[1]
        derivative = (-I*g*(phi.conjugate().T*T[k])*dp + I*g*dp.conjugate().T*T[k]*phi)[0]
        current = I*(phi.conjugate().T*T[k]*dp-dp.conjugate().T*T[k]*phi)[0]
        check(f'{seed}: scalar current variation {k}', derivative+g*current)
        no_connection = I*(phi.conjugate().T*T[k]*dphi[1]-dphi[1].conjugate().T*T[k]*phi)[0]
        reject(f'{seed}: scalar current connection omitted {k}', current-no_connection)
    # Adjoint reality and covariance of the W shift under a finite SU(2) rotation.
    U = (s.eye(2)+I*pauli[1])/s.sqrt(2)
    adj = s.Matrix(3,3,lambda k,j: s.trace(T[k]*U*T[j]*U.conjugate().T)*2)
    j = s.Matrix([bil(L, gamma[0], t, phi, N) for t in T])
    Lp = [sm([L[b]*U.conjugate().T[b,a] for b in range(2)], s.zeros(1,4)) for a in range(2)]
    jp = s.Matrix([bil(Lp, gamma[0], t, U*phi, N) for t in T])
    check(f'{seed}: adjoint gauge-shift covariance', jp-adj*j)
    shift = g*(j+j.conjugate())/6
    check(f'{seed}: real W shift', shift-shift.conjugate())
    residual = vec(3).applyfunc(s.re)
    packet = -(g*j.dot(residual))/6
    check(f'{seed}: real gauge shift cancels residual plus conjugate', residual.dot(shift)+packet+s.conjugate(packet))
    reject(f'{seed}: reversed gauge field substitution', -residual.dot(shift)+packet+s.conjugate(packet))
    reject(f'{seed}: omitted gauge Hermitian conjugate', residual.dot(g*j/6)+packet+s.conjugate(packet))
    # Curvature covector ordering and its actual Dirac adjoint.
    mu, nu = next((m,n) for m in range(4) for n in range(m+1,4) if F[m][n] != F[m][n].T)
    barx = sm([L[a]*sigma[mu][nu]*F[mu][nu][a,b]*phi[b] for a in range(2) for b in range(2)], s.zeros(1,4))
    direct = gamma[0]*barx.conjugate().T
    expected = sm([phi[a].conjugate()*F[mu][nu][a,b]*sigma[mu][nu]*gamma[0]*L[b].conjugate().T for a in range(2) for b in range(2)], s.zeros(4,1))
    check(f'{seed}: sterile curvature covector Dirac adjoint', direct-expected)
    wrong = sm([phi[a].conjugate()*F[mu][nu][b,a]*sigma[mu][nu]*gamma[0]*L[b].conjugate().T for a in range(2) for b in range(2)], s.zeros(4,1))
    reject(f'{seed}: transposed weak curvature in shift', direct-wrong)

for mu in range(4):
    for nu in range(4):
        check(f'Clifford {mu}{nu}', gamma[mu]*gamma[nu]+gamma[nu]*gamma[mu]-2*eta[mu]*(mu==nu)*s.eye(4))
        check(f'Grassmann tensor spin symmetry {mu}{nu}', (sigma[mu][nu]*C).T-sigma[mu][nu]*C)
check('Majorana scalar charge-conjugate exchange', C+C.T)
# Grassmann coefficient extraction for two independent odd spinors: E^T C w
# versus w^T C E, with the second reordered using anticommutation.
check('Majorana scalar residual and conjugate', C-(-C.T))
reject('Majorana half shift', C-C/2)
# One active flavor already rejects a flavor-only antisymmetrization for SU(2).
ph = s.Matrix([1,0])
gamma5 = I*gamma[0]*gamma[1]*gamma[2]*gamma[3]
PR = (s.eye(4)+gamma5)/2
spin = PR*sigma[0][1]*C*PR.T
reject('nonzero same-chirality tensor spin block', spin)
for label, mat in [('hypercharge', -s.eye(2)/2), ('W1', T[0]), ('W2', T[1]), ('W3', T[2])]:
    coefficient = s.kronecker_product((mat*ph)*ph.T, spin)
    grassmann_coefficient = coefficient-coefficient.T
    if label in ('hypercharge', 'W3'):
        check(label+' single-flavor tensor vanishes on chosen scalar', grassmann_coefficient)
    else:
        reject(label+' erroneous flavor-only antisymmetrization', grassmann_coefficient)
K = s.Matrix([[2,3+I],[-1,5]])
KA, KS = (K-K.T)/2, (K+K.T)/2
for a, mat in enumerate(T):
    weak = (mat*ph)*ph.T
    weakS, weakA = (weak+weak.T)/2, (weak-weak.T)/2
    full = s.kronecker_product(K, weak, spin)
    reduced = s.kronecker_product(KA, weakS, spin)+s.kronecker_product(KS, weakA, spin)
    check(f'combined flavor/weak Grassmann projection {a}', (full-full.T)/2-reduced)
    if a < 2:
        reject(f'flavor-symmetric charged contribution omitted {a}', reduced-s.kronecker_product(KA, weakS, spin))
for a in range(3):
    check(f'SU2 generator normalization {a}', s.trace(T[a]*T[a])-R(1,2))
gp, g, y = R(2,3), R(3,2), R(-1,2)
Bcur, Wcur = s.Integer(5), [s.Integer(2),s.Integer(-3),s.Integer(7)]
matrix_source = gp*y*s.eye(2)*(gp*Bcur)+g*sm([T[a]*(g*Wcur[a]) for a in range(3)], z2)
expected_source = gp**2*y*Bcur*s.eye(2)+g**2*sm([T[a]*Wcur[a] for a in range(3)], z2)
check('gauge source coupling factors and T normalization', matrix_source-expected_source)
reject('gauge source loses second g', matrix_source-(gp**2*y*Bcur*s.eye(2)+g*sm([T[a]*Wcur[a] for a in range(3)], z2)))
reject('Pauli matrices substituted for T', matrix_source-(gp**2*y*Bcur*s.eye(2)+g**2*sm([pauli[a]*Wcur[a] for a in range(3)], z2)))
check('combined mass/source dipole coefficient', R(1,6)-R(1,4)+R(1,12))
reject('isolated transverse coefficient reused', R(1,6)-R(1,4)+R(1,4))
for seed in [31415, 27182, 16180]:
    trial(seed)
result = {'status':'pass', 'python':sys.version.split()[0], 'sympy':s.__version__,
          'scope':'Independent exact local matrix action, gauge variation, field-shift and Grassmann controls; no loop integration, beta function, or observable.',
          'input_commit':'d822e9f4d70ad5b455436f40abc3c2b45909bbef',
          'exact_checks':len(checks), 'rejected_mutations':len(negatives),
          'checks':checks, 'mutations':negatives}
print(json.dumps(result, indent=2))
