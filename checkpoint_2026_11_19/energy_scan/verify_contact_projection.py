"""Check the two-derivative contact projection using exact Weyl matrices."""
import json
from collections import Counter
from itertools import product
from pathlib import Path

import sympy as sp


E, x, y, z = sp.symbols("E x y z", real=True)
mass2 = E**2 - x**2 - y**2 - z**2
s = 4*E**2
t = mass2 - 2*E*(E-z)
u = mass2 - 2*E*(E+z)
p1 = sp.Matrix([E, 0, 0, E])
p2 = sp.Matrix([E, 0, 0, -E])
k1 = sp.Matrix([E, x, y, z])
k2 = p1 + p2 - k1
momenta = (p1, p2, k1, k2)
spin1, spin2 = sp.Matrix([0, 1]), sp.Matrix([1, 0])
epsilon = sp.Matrix([[0, 1], [-1, 0]])
metric = sp.diag(1, -1, -1, -1)


def dot(q, r):
    return (q.T*metric*r)[0]


def sigma(q, sign):
    energy, qx, qy, qz = q
    return sp.Matrix([[energy + sign*qz, sign*(qx-sp.I*qy)],
                      [sign*(qx+sp.I*qy), energy - sign*qz]])


def det_spin(left, right):
    return (left.T*epsilon*right)[0]


def reduced_chain(q, r, left, right):
    chain = sigma(q, -1)*sigma(r, 1)
    return (left.T*epsilon*chain*right)[0] / det_spin(left, right)


def symmetrized(i, j, contraction):
    """Factoring the spin determinant converts Fermi exchange to a plus."""
    images = []
    for leptons, spinors in [((p1, p2), (spin1, spin2)),
                             ((p2, p1), (spin2, spin1))]:
        for scalars in [(k1, k2), (k2, k1)]:
            slots = leptons + scalars
            if contraction == "metric":
                images.append(dot(slots[i], slots[j]))
            else:
                images.append(reduced_chain(slots[i], slots[j], *spinors))
    return sp.expand(sum(images)/4)


def main():
    checks = Counter()

    def equal(label, expr, target):
        if sp.expand(expr-target) != 0:
            raise ArithmeticError(f"{label}: {expr} != {target}")
        checks[label] += 1

    def different(label, expr, target):
        if sp.expand(expr-target) == 0:
            raise ArithmeticError(f"surviving mutation: {label}")
        checks[label] += 1

    table = ((0, 0, 0, 0), (s, 0, mass2-u, mass2-t),
             (mass2-t, 0, mass2, -t), (mass2-u, 0, -u, mass2))
    equal("on-shell Mandelstam sum", s+t+u, 2*mass2)
    equal("antisymmetric lepton determinant", det_spin(spin2, spin1), -det_spin(spin1, spin2))
    images = []
    monomials = [(2, 0, 0, 0), (0, 2, 0, 0), (0, 0, 2, 0), (0, 0, 0, 2)]
    for i, j in product(range(4), repeat=2):
        q, r = momenta[i], momenta[j]
        chain = reduced_chain(q, r, spin1, spin2)
        reverse = reduced_chain(r, q, spin1, spin2)
        equal("ordered Weyl matrix elements", chain, table[i][j])
        equal("Weyl anticommutator", chain+reverse, 2*dot(q, r))
        for contraction in ("metric", "chain"):
            expr = symmetrized(i, j, contraction)
            coeffs = [sp.Poly(expr, E, x, y, z).coeff_monomial(m) for m in monomials]
            images.append(coeffs)
            # Reconstruct exclusively from s and a using the E^2 and x^2 components.
            reconstructed = ((coeffs[0]+coeffs[1])/4)*s - coeffs[1]*mass2
            equal("external exchanges leave only s and scalar mass squared", expr, reconstructed)
    equal("rank of derivative images as both a and s vary", sp.Matrix(images).rank(), 2)
    equal("scalar derivative contraction", dot(k1, k2), (s-2*mass2)/2)
    equal("lepton derivative contraction", dot(p1, p2), s/2)
    equal("scalar self-contraction", dot(k1, k1), mass2)
    equal("mean literal composite derivative", (t+u)/2, mass2-s/2)
    equal("Bose projection of scalar-pair sigma tensor",
          symmetrized(2, 3, "chain")-symmetrized(3, 2, "chain"), 0)
    q, r = p1-p2, k1+k2
    equal("mixed sigma tensor retains an energy slope",
          (reduced_chain(q, r, spin1, spin2)-reduced_chain(r, q, spin1, spin2))/2, -s)
    different("ordered t transfer differs from its Bose average", t, (t+u)/2)
    different("discarding all tensor representatives loses a slope",
              (reduced_chain(q, r, spin1, spin2)-reduced_chain(r, q, spin1, spin2))/2, 0)
    different("derivative slopes can cancel the literal seed",
              (t+u)/2 + dot(k1, k2), (t+u)/2)
    equal("two derivative representatives cancel for all energies", (t+u)/2 + dot(k1, k2), 0)
    # Fix the scalar mass and s, then change the angle at dimension nine.
    higher = (t-u)**2
    different("four derivatives permit even angular dependence",
              higher.subs({x: 0, y: 0, z: 1}), higher.subs({x: 1, y: 0, z: 0}))

    result = {"status": "PASS", "method": "exact symbolic two-component spinors and both external exchange projectors",
              "ordered_chain_rows": ["p1", "p2", "k1", "k2"], "polynomial_identity_checks": dict(checks),
              "checks": sum(checks.values()), "derivative_rank_varying_s_and_a": 2,
              "fixed_mass_shapes_through_dimension_seven": ["constant", "s"],
              "repository_imports": False, "assertions_required": False}
    data = json.dumps(result, indent=2)+"\n"
    Path(__file__).with_name("CONTACT_PROJECTION_RESULTS.json").write_text(data, encoding="utf-8")
    print(data, end="")


if __name__ == "__main__":
    main()
