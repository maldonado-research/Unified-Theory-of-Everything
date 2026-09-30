"""N00AK-r1 tree mathematics and CONDITIONAL matrix constructions.

No function in this module evaluates a loop integral or establishes a
dimension-seven beta function. Adaptation provenance is in PROVENANCE.json.
"""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import numpy as np

SOURCE_SHA256 = "5806d06fd82c3f1ca28f7a6e60db0383a1160304a090577e5457566eccf6e082"
MASSES_GEV = np.array([2e12, 4e12, 8e12], dtype=float)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_frozen_yukawa(path: Path) -> np.ndarray:
    if sha256_file(path) != SOURCE_SHA256:
        raise ValueError("frozen Yukawa source hash mismatch")
    value = np.zeros((3, 3), dtype=complex)
    seen = set()
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("matrix") not in (None, "", "Ynu_star"):
                continue
            key = (int(row["row"]) - 1, int(row["col"]) - 1)
            if key in seen or not (0 <= key[0] < 3 and 0 <= key[1] < 3):
                raise ValueError("invalid or duplicate Yukawa coordinate")
            value[key] = complex(float(row["real"]), float(row["imag"]))
            seen.add(key)
    if len(seen) != 9:
        raise ValueError("incomplete Yukawa source")
    return value


def tree_moment(y: np.ndarray, masses: np.ndarray, order: int) -> np.ndarray:
    if order < 0:
        raise ValueError("order must be nonnegative")
    return (y * masses ** (-2 * order - 1)) @ y.T


def tree_resolvent(y: np.ndarray, masses: np.ndarray, z: complex) -> np.ndarray:
    denominators = 1 - z / masses**2
    if np.any(denominators == 0):
        raise ZeroDivisionError("probe lies on a supplied pole")
    return (y * (1 / (masses * denominators))) @ y.T


def tree_truncation(y: np.ndarray, masses: np.ndarray, z: complex, order: int) -> np.ndarray:
    result = np.zeros((y.shape[0], y.shape[0]), dtype=complex)
    for n in range(order + 1):
        result += z**n * tree_moment(y, masses, n)
    return result


def tree_remainder(y: np.ndarray, masses: np.ndarray, z: complex, order: int) -> np.ndarray:
    ratio = z / masses**2
    if np.any(ratio == 1):
        raise ZeroDivisionError("probe lies on a supplied pole")
    weights = ratio ** (order + 1) / (masses * (1 - ratio))
    return (y * weights) @ y.T


def remainder_bound(y: np.ndarray, masses: np.ndarray, z: complex, order: int) -> float:
    minimum = float(min(masses))
    rho = abs(z) / minimum**2
    if not rho < 1:
        raise ValueError("bound requires |z| below every supplied pole")
    return float(np.linalg.norm(y, 2)**2 * abs(z)**(order + 1)
                 / (minimum**(2 * order + 3) * (1 - rho)))


def effective_nodes(y: np.ndarray, masses: np.ndarray) -> np.ndarray:
    """Known distinct nodes after exact equal-mass combination/zero deletion.

    No approximate equality or numerical inference of unknown nodes is used.
    """
    nodes = []
    for mass in np.unique(masses):
        selected = y[:, masses == mass]
        residue = selected @ selected.T / mass
        if np.any(residue != 0):
            nodes.append(1 / mass**2)
    return np.array(nodes, dtype=float)


def pole_polynomial(nodes: np.ndarray) -> np.ndarray:
    coefficients = np.array([1.0])
    for node in nodes:
        coefficients = np.convolve(coefficients, np.array([1.0, -node]))
    return coefficients


def tree_reconstruction(moments: list[np.ndarray], nodes: np.ndarray, z: complex) -> np.ndarray:
    """Known-node reconstruction; nodes must be distinct and effective."""
    q = len(nodes)
    if q == 0 or len(moments) != q:
        raise ValueError("one initial moment per nonzero effective node is required")
    polynomial = pole_polynomial(nodes)
    numerator = np.zeros_like(moments[0], dtype=complex)
    for n in range(q):
        for k in range(n + 1):
            numerator += z**n * polynomial[k] * moments[n - k]
    denominator = np.prod(1 - nodes * z)
    if denominator == 0:
        raise ZeroDivisionError("probe lies on an effective pole")
    return numerator / denominator


def even_spurion(y: np.ndarray, masses: np.ndarray, order: int) -> np.ndarray:
    return (y * masses ** (2 * order)) @ y.conj().T


def schur_variation(y: np.ndarray, masses: np.ndarray, delta_y: np.ndarray) -> np.ndarray:
    """First variation of Y M^-1 Y^T at fixed M; uses transpose, not adjoint."""
    return (delta_y / masses) @ y.T + (y / masses) @ delta_y.T


def conditional_local(y_a: np.ndarray, m_a: np.ndarray, moment: np.ndarray,
                      order: int, coefficient: complex = 3) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Assumed deltaY/contact and their Schur variation, for any common c.

    This is an algebraic ansatz. The factor c is an input, not a derived
    loop coefficient, including when c=3 and order=1.
    """
    delta_y = -coefficient * (moment @ y_a.conj()) * m_a**(2 * order + 1)
    spurion = even_spurion(y_a, m_a, order)
    contact = coefficient * (spurion @ moment + moment @ spurion.T)
    return delta_y, contact, schur_variation(y_a, m_a, delta_y)


def conditional_rational(y_a: np.ndarray, m_a: np.ndarray, y_d: np.ndarray,
                         m_d: np.ndarray, coefficient: complex = 3) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Rational sum of the assumed local matrices, not resummed loop matching."""
    delta_y = np.zeros_like(y_a)
    contact = np.zeros((y_a.shape[0], y_a.shape[0]), dtype=complex)
    for a, mass in enumerate(m_a):
        kernel = tree_resolvent(y_d, m_d, mass**2)
        column = y_a[:, a]
        delta_y[:, a] = -coefficient * mass * kernel @ column.conj()
        outer = np.outer(column, column.conj())
        contact += coefficient * (outer @ kernel + kernel @ outer.T)
    return delta_y, contact, schur_variation(y_a, m_a, delta_y)
