#!/usr/bin/env python3
"""Exact polynomial Ward-kernel classification for the declared odd-chain space.

No on-shell, integration-by-parts, or equations-of-motion quotient is imposed.
The coefficient field is constant complex/rational scalars. A separate note
must justify the declared Lorentz/locality premises; this is a linear proof
of the Ward kernel inside that space, not a complete EFT operator basis.
"""
from fractions import Fraction as F
import json

LABELS = (
    'gamma*p2', 'gamma*r2', 'gamma*k2', 'gamma*(p.r)',
    'gamma*(p.k)', 'gamma*(r.k)',
    'p^nu*slashp', 'p^nu*slashr', 'p^nu*slashk',
    'r^nu*slashp', 'r^nu*slashr', 'r^nu*slashk',
    'k^nu*slashp', 'k^nu*slashr', 'k^nu*slashk',
    'Gamma(p,nu,r)', 'Gamma(k,nu,p)', 'Gamma(k,nu,r)',
    'm2*gamma',
)


def rref(matrix):
    a = [[F(x) for x in row] for row in matrix]
    pivots, i = [], 0
    for j in range(len(a[0])):
        pivot = next((k for k in range(i, len(a)) if a[k][j]), None)
        if pivot is None:
            continue
        a[i], a[pivot] = a[pivot], a[i]
        d = a[i][j]
        a[i] = [v/d for v in a[i]]
        for k in range(len(a)):
            if k != i and a[k][j]:
                f = a[k][j]
                a[k] = [v-f*w for v, w in zip(a[k], a[i])]
        pivots.append(j)
        i += 1
        if i == len(a):
            break
    return a, pivots


def column(terms):
    return [F(terms.get(label, 0)) for label in LABELS]


def apply(a, v):
    return [sum((x*y for x, y in zip(row, v)), F()) for row in a]


def check(ok, label):
    if not ok:
        raise RuntimeError(label)


def main():
    # The targets are the independent slash(p/r/k) scalar monomials and
    # Gamma(p,k,r). Contracting k into the public numerator basis yields:
    rows = []
    for idx in (0, 1, 3, 6, 7, 9, 10, 12, 13, 15, 18):
        row = [F(0)]*len(LABELS)
        row[idx] = 1
        rows.append(row)
    for left, right in ((2, 14), (4, 8), (5, 11)):
        row = [F(0)]*len(LABELS)
        row[left] = row[right] = 1
        rows.append(row)
    reduced, pivots = rref(rows)
    check(len(pivots) == 14, 'Ward map has exact rank fourteen including mass spurion')
    check(len(LABELS)-len(pivots) == 5, 'five-dimensional transverse space')
    basis = [column(t) for t in (
        {'p^nu*slashk':1, 'gamma*(p.k)':-1},
        {'r^nu*slashk':1, 'gamma*(r.k)':-1},
        {'k^nu*slashk':1, 'gamma*k2':-1},
        {'Gamma(k,nu,p)':1},
        {'Gamma(k,nu,r)':1},
    )]
    check(all(not any(apply(rows, v)) for v in basis), 'five proposed tensors are transverse')
    check(len(rref(list(map(list, zip(*basis))))[1]) == 5, 'five tensors are independent')
    check(all(v[-1] == 0 for v in basis), 'no independent local m2 transverse tensor')
    # Exhaust the nullspace from row reduction rather than assuming these five
    # basis vectors span it. Each free variable generates an independent vector.
    for j in set(range(len(LABELS)))-set(pivots):
        v = [F(0)]*len(LABELS)
        v[j] = 1
        for i, pivot in enumerate(pivots):
            v[pivot] = -reduced[i][j]
        check(not any(apply(rows, v)), 'constructed nullspace vector')
    # D=1/2 T_r+1/3 T_k+1/2 Gamma(k,nu,r).
    coords = (F(0), F(1,2), F(1,3), F(0), F(1,2))
    d = [sum((c*v[i] for c,v in zip(coords,basis)),F()) for i in range(len(LABELS))]
    check(not any(apply(rows,d)), 'computed loop remainder is in Ward kernel')
    check(d == column({'r^nu*slashk':F(1,2),'gamma*(r.k)':-F(1,2),
                      'k^nu*slashk':F(1,3),'gamma*k2':-F(1,3),
                      'Gamma(k,nu,r)':F(1,2)}), 'loop coefficients in complete declared kernel')
    check(sum(bool(c) for c in coords) == 3, 'three nonzero coordinates; compact two-term expression is not a full basis')
    bad = list(d)
    bad[-1] = 1
    check(any(apply(rows,bad)), 'mass-spurion mutation violates Ward')
    omitted = list(d)
    omitted[17] = 0
    check(not any(apply(rows,omitted)) and omitted != d, 'omitted grade-three term evades Ward')
    bad = list(d)
    bad[8] += 1
    check(any(apply(rows,bad)), 'nontransverse extra momentum tensor detected')
    print(json.dumps({'status':'PASS','declared_coefficients':19,'ward_rank':14,
          'transverse_dimension':5,'basis':['T_p','T_r','T_k','Gamma(k,nu,p)','Gamma(k,nu,r)'],
          'loop_coordinates':[str(c) for c in coords],'mass_spurion_in_kernel':False,
          'negative_controls':3,
          'scope':'local degree-two odd open chain with p,r,k; no EOM/IBP quotient or full EFT basis claim'},indent=2))


if __name__ == '__main__':
    main()
