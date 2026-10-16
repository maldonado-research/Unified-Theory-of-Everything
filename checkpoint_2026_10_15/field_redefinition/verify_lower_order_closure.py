#!/usr/bin/env python3
"""Exact polynomial action substitution for a declared lower-order benchmark."""
import json
from pathlib import Path
import random
import sympy as s

t, rho, drho = s.symbols('t rho drho', real=True)
checks = 0
mutations = []


def require(condition, label):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(label)


def zero(matrix):
    return all(s.expand(x) == 0 for x in matrix)


def coefficient(matrix, degree):
    return matrix.applyfunc(lambda x: s.expand(x).coeff(t, degree))


def reject(matrix, label):
    require(not zero(matrix), 'Mutation survived: ' + label)
    mutations.append(label)


rng = random.Random(71015)


def complex_matrix(rows, cols):
    return s.Matrix(rows, cols, lambda i, j:
                    s.Rational(rng.randrange(-5, 6), rng.randrange(1, 5))
                    + s.I*s.Rational(rng.randrange(-5, 6), rng.randrange(1, 5)))


for trial, n in enumerate([1, 2, 3, 2, 3, 3]):
    a, b, c, d = [complex_matrix(n, n) for _ in range(4)]
    Z = a + a.conjugate().T
    M = b + b.T
    if trial == 5:
        v = complex_matrix(n, 1)
        M = v*v.T  # Singular symmetric mass; no inverse is used.
        require(M.rank() == 1, 'Rank-deficient mass fixture')
    C, D = c+c.T, d+d.T
    Y = complex_matrix(n, n+1)
    K = s.eye(n)+t**2*rho*Z
    S = s.eye(n)-t**2*rho*Z/2
    kinetic = S.conjugate().T*K*S
    mass = S.T*(M+t*rho*C+t**3*rho**2*D)*S
    source = S.T*Y
    connection = S.conjugate().T*K*S.diff(rho)*drho-S.diff(rho).conjugate().T*drho*K*S
    expected = [M, rho*C, -rho*(Z.T*M+M*Z)/2,
                rho**2*(D-(Z.T*C+C*Z)/2)]
    for order in range(4):
        require(zero(coefficient(kinetic, order)-(s.eye(n) if order == 0 else s.zeros(n))),
                f'kinetic normalization trial {trial} degree {order}')
        require(zero(coefficient(mass, order)-expected[order]),
                f'full action substitution trial {trial} degree {order}')
        require(zero(coefficient(connection, order)),
                f'local kinetic connection trial {trial} degree {order}')
    require(zero(coefficient(source, 0)-Y), 'source leading term')
    require(zero(coefficient(source, 2)+rho*Z.T*Y/2), 'source order-two term')
    delta = -(Z.T*C+C*Z)/2
    require(zero(delta-delta.T), 'symmetric dimension-seven shift')
    # A dense exact unitary matrix (discrete Fourier for n=2,3) with complex phases.
    if n == 1:
        U = s.Matrix([[s.I]])
    elif n == 2:
        U = s.Matrix([[1, s.I], [s.I, 1]])/s.sqrt(2)
    else:
        # Householder reflection, followed by distinct phases, stays exact rational.
        u = s.Matrix([1, 2, 2])
        U = (s.eye(3)-2*u*u.T/(u.T*u)[0])*s.diag(1, s.I, -1)
    require(zero(s.simplify(U.conjugate().T*U)-s.eye(n)), 'unitary fixture')
    Zp = s.simplify(U.conjugate().T*Z*U)
    Cp = s.simplify(U.T*C*U)
    delta_p = -(Zp.T*Cp+Cp*Zp)/2
    require(zero(s.simplify(delta_p-U.T*delta*U)), 'complex basis covariance')
    reject(coefficient(mass, 3)-rho**2*D, f'omit dimension-five cross term {trial}')
    if n > 1:
        reject(delta+(Z.conjugate().T*C+C*Z)/2, f'dagger instead of transpose {trial}')
        reject(coefficient(source, 2)+rho*Z*Y/2, f'reverse holomorphic source adjoint {trial}')
        reject(coefficient(mass, 2)+rho*(Z*M+M*Z)/2, f'wrong mass congruence {trial}')

report = {
    'status': 'PASS',
    'positive_and_mutation_controls': checks,
    'exact_matrix_trials': 6,
    'rejected_mutations': len(mutations),
    'mutation_labels': mutations,
    'scope': 'Declared optional dimension-five sterile-scalar and dimension-six kinetic operators; exact formal action substitution through inverse-scale order three. No loop coefficient, complete EFT closure, physical running or observable is calculated.',
    'uses_mass_inverse': False,
    'checks_active_under_optimization': True,
}
Path(__file__).with_name('LOWER_ORDER_CLOSURE_RESULTS.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
