"""Check the open-line projection and local evanescent-coordinate conversion."""

import itertools
import json
from pathlib import Path

import sympy as sp


positive = []
adverse = []


def require_equal(name, left, right):
    if sp.simplify(left - right) != 0:
        raise ValueError(f"{name}: {left} != {right}")
    positive.append(name)


def reject_equal(name, left, right):
    if sp.simplify(left - right) == 0:
        raise ValueError(f"{name}: an incorrect alternative survived")
    adverse.append(name)


def require_matrix(name, actual, expected):
    if (actual - expected).applyfunc(sp.simplify) != sp.zeros(*actual.shape):
        raise ValueError(f"{name}: unequal matrices")
    positive.append(name)


def reorder_sign(slots):
    inversions = sum(a > b for a, b in itertools.combinations(slots, 2))
    return (-1) ** inversions


# The gamma matrices below test the physical projection. They do not stand in
# for a representation of the formal open chains in non-integer dimension.
i = sp.I
eye = sp.eye(2)
zero = sp.zeros(2)
metric = [1, -1, -1, -1]
pauli = [
    sp.Matrix([[0, 1], [1, 0]]),
    sp.Matrix([[0, -i], [i, 0]]),
    sp.diag(1, -1),
]


def blocks(ul, ur, ll, lr):
    return ul.row_join(ur).col_join(ll.row_join(lr))


gamma = [blocks(zero, eye, eye, zero)]
gamma.extend(blocks(zero, s, -s, zero) for s in pauli)
gamma5 = i * gamma[0] * gamma[1] * gamma[2] * gamma[3]
charge = i * gamma[2] * gamma[0]
spin_epsilon = sp.Matrix([[0, 1], [-1, 0]])
require_matrix("chirality of the selected blocks", gamma5, sp.diag(-1, -1, 1, 1))
require_matrix("charge conjugation blocks", charge, sp.diag(spin_epsilon, -spin_epsilon))
require_matrix("transpose of charge conjugation", charge.T, -charge)
for mu, matrix in enumerate(gamma):
    require_matrix(f"transpose intertwiner {mu}", charge * matrix.T * charge.inv(), -matrix)
require_matrix("chirality transpose", charge * gamma5.T * charge.inv(), gamma5)

# Canonical odd-field order is [bar_alpha, bar_beta, L_gamma, L_delta].
# Current multiplication has order [bar_alpha, L_gamma, bar_beta, L_delta].
# Fix one pair of weak indices. Every pair has the identical spin tensor.
components = list(itertools.product(range(2), repeat=4))
ordered, currents = [], []
for r, t, u, v in components:
    ordered.append(charge[r + 2, t + 2] * charge.inv()[u, v])
    currents.append(reorder_sign([0, 2, 1, 3]) * sum(
        sign * matrix[r + 2, u] * matrix[t + 2, v]
        for sign, matrix in zip(metric, gamma)
    ))
ordered = sp.Matrix(ordered)
currents = sp.Matrix(currents)
weyl = sp.Matrix([spin_epsilon[r, t] * spin_epsilon[u, v] for r, t, u, v in components])
require_matrix("ordered Weyl lift for all sixteen components", ordered, weyl)
require_matrix("physical Fierz projection for all sixteen components", currents, -2 * ordered)
witness = components.index((0, 1, 0, 1))
require_equal("nonzero ordered spin component", ordered[witness], 1)
reject_equal("charge matrix in place of its inverse in the unbarred pair",
             -ordered[witness], weyl[witness])
reject_equal("commuting the middle odd fields without a sign", -currents[witness],
             -2 * ordered[witness])

# A single ordered flavor channel. h includes the action half and 16 pi^2.
eps = sp.Symbol("epsilon", nonzero=True)
a, b, h = sp.symbols("a b h")
p, j = sp.symbols("P_d J_d")
half = sp.Rational(1, 2)


def evanescent(parameter):
    return p + (half + parameter * eps) * j


def pole_part(parameter):
    return h * (j / 2 - evanescent(parameter)) / eps


fixed_counterterm = -h * p / eps
for parameter in (a, b):
    require_equal(f"same ordered counterterm in frame {parameter}",
                  pole_part(parameter) + h * parameter * j, fixed_counterterm)
    require_equal(f"four-dimensional projection of E_{parameter}",
                  evanescent(parameter).subs({p: -j / 2, eps: 0}), 0)
require_equal("translation of evanescent basis",
              evanescent(b) - evanescent(a), (b - a) * eps * j)
require_equal("finite difference of pole-only prescriptions",
              pole_part(b) - pole_part(a), -h * (b - a) * j)
reject_equal("wrong sign of compensating finite counterterm",
             pole_part(a) - h * a * j, fixed_counterterm)
reject_equal("drop the evanescent pole before converting coordinates",
             h * j / (2 * eps) + h * a * j, fixed_counterterm)

# Different O(epsilon) physical lifts with the same four-dimensional limit.
x, y = sp.symbols("X Y")
p_new, j_new = p + eps * x, j + eps * y
e_new = p_new + (half + b * eps) * j_new
require_equal("general local lift transformation",
              (e_new - evanescent(a)).expand().coeff(eps, 1),
              x + y / 2 + (b - a) * j)
require_equal("ordered finite shift under changing its lift",
              -h * p_new / eps - fixed_counterterm, -h * x)

# A reduced form factor has the same external operator stripped at the two
# momenta. No assumption about its analytic dependence on t enters here.
t, local, derivative = sp.symbols("t local derivative")
kernel = sp.Function("K")
difference = (kernel(t) + local) - (kernel(0) + local)
require_equal("constant contact cancels in reduced form-factor subtraction",
              difference, kernel(t) - kernel(0))
reject_equal("extend the same conclusion to a derivative contact",
             difference + derivative * t, kernel(t) - kernel(0))

# A channel with distinct barred and unbarred flavor pairs has an independent
# conjugate slot. Real scheme parameters shift the two slots conjugately.
channel = sp.Matrix([2 + 3*i, 2 - 3*i])
finite_shift = sp.Rational(3, 5) * channel
require_equal("conjugate slot has conjugate shifted coefficient",
              finite_shift[1], sp.conjugate(finite_shift[0]))
reject_equal("copy complex channel coefficient into the conjugate slot",
             finite_shift[1], finite_shift[0])

receipt = {
    "scope": "Physical spin projection and finite local conversion for a single ordered pole.",
    "positive_checks": positive,
    "adverse_controls": adverse,
    "spin_components": len(components),
    "nonzero_component": {"slots": list(components[witness]), "P_4": 1, "J_4": -2},
}
Path(__file__).with_name("FIERZ_SCHEME_RESULTS.json").write_text(
    json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
)
print(f"PASS {len(positive)} checks and {len(adverse)} adverse controls")
