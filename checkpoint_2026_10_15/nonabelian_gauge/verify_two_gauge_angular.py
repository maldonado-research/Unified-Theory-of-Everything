#!/usr/bin/env python3
"""Independent exact large-q and color-index extraction of two-gauge-leg matter graphs.

Only Python's standard library is used. No producer or repository module is
imported. Clifford blades implement Cl(1,3) over Q; UV extraction uses the
homogeneous radial logarithm and exact d=4 angular moments, with an implicit
infrared separator. This never identifies an unseparated scaleless integral
with a UV pole. All returned poles strip i/(16*pi**2*epsilon).
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import prod
import argparse
import json
import random

ETA = (1, -1, -1, -1)
Z = (0, 0, 0, 0)
E = tuple(tuple(int(i == j) for i in range(4)) for j in range(4))


def va(*vectors):
    return tuple(sum(v[i] for v in vectors) for i in range(4))


def vs(a, b):
    return tuple(a[i] - b[i] for i in range(4))


def vm(c, a):
    return tuple(c * x for x in a)


def dot(a, b):
    return sum(ETA[i] * a[i] * b[i] for i in range(4))


@lru_cache(None)
def blade(a, b):
    """Ordered gamma^0...gamma^3 blades; return product sign and mask."""
    inversions = sum((b & ((1 << i) - 1)).bit_count()
                     for i in range(4) if a & (1 << i))
    sign = (-1) ** inversions
    sign *= prod(ETA[i] for i in range(4) if a & b & (1 << i))
    return sign, a ^ b


def clean(a):
    return {k: v for k, v in a.items() if v}


def add(*items):
    out = {}
    for item in items:
        for k, v in item.items():
            out[k] = out.get(k, 0) + v
    return clean(out)


def scale(c, a):
    return clean({k: c * v for k, v in a.items()})


def mul(a, b):
    out = {}
    for ia, ca in a.items():
        for ib, cb in b.items():
            sign, k = blade(ia, ib)
            out[k] = out.get(k, 0) + sign * ca * cb
    return clean(out)


def chain(*items):
    out = {0: 1}
    for item in items:
        out = mul(out, item)
    return out


def slash(v):
    return clean({1 << i: ETA[i] * v[i] for i in range(4)})


def comm(a, b):
    return add(mul(a, b), scale(-1, mul(b, a)))


def pc(a):
    return {(Z, mask): value for mask, value in a.items()}


def qdot(v):
    return clean({(E[i], 0): ETA[i] * v[i] for i in range(4)})


QS = {(E[i], 1 << i): ETA[i] for i in range(4)}
Q2 = {(vm(2, E[i]), 0): ETA[i] for i in range(4)}


def pmul(a, b):
    out = {}
    for (ea, ia), ca in a.items():
        for (eb, ib), cb in b.items():
            sign, mask = blade(ia, ib)
            key = (va(ea, eb), mask)
            out[key] = out.get(key, 0) + sign * ca * cb
    return clean(out)


def pchain(*items):
    out = {(Z, 0): 1}
    for item in items:
        out = pmul(out, item)
    return out


@lru_cache(None)
def moment(exp):
    """Normalized Minkowski angular moment from Wick pairings at d=4."""
    if any(x % 2 for x in exp):
        return F(0)
    n = sum(exp) // 2
    pairings = prod(prod(range(1, x, 2)) for x in exp)
    signs = prod(ETA[i] ** (exp[i] // 2) for i in range(4))
    return F(signs * pairings, prod(4 + 2 * j for j in range(n)))


def angular(poly):
    out = {}
    for (exp, mask), value in poly.items():
        out[mask] = out.get(mask, 0) + value * moment(exp)
    return clean(out)


def jmul(a, b):
    """Top two homogeneous numerator degrees, descending in q."""
    return pmul(a[0], b[0]), add(pmul(a[0], b[1]), pmul(a[1], b[0]))


def jchain(*items):
    out = ({(Z, 0): 1}, {})
    for item in items:
        out = jmul(out, item)
    return out


def seed_jet(v):
    # Full seed is q^2 - 2 q.v + v^2 + 2 t^2. The constant is two
    # degrees below its top, and hence below the logarithmic UV threshold.
    return Q2, scale(-2, qdot(v))


def fprop_jet(shift):
    return QS, pc(scale(-1, slash(shift)))


def scalar_vertex_jet(u, eps):
    return scale(-2, qdot(eps)), pc({0: dot(u, eps)})


def vertex_jet(eps):
    return pc(slash(eps)), {}


def pole(jet, shifts, denominators):
    # Box: N4/(q^2)^4 + 2(q.sum shifts)N5/(q^2)^5.
    # Triangle: N2/(q^2)^3 + 2(q.sum shifts)N3/(q^2)^4.
    expected = 2 * denominators - 3
    if any(sum(exp) != expected for exp, _ in jet[0]):
        raise RuntimeError("unexpected leading numerator degree")
    if any(sum(exp) != expected - 1 for exp, _ in jet[1]):
        raise RuntimeError("unexpected subleading numerator degree")
    return add(angular(jet[1]),
               angular(pmul(jet[0], scale(2, qdot(va(*shifts))))))


NAMES = ("FF12", "FF21", "SS12", "SS21", "FS12", "FS21", "SG")
WEIGHTS = (1, 1, 1, 1, -1, -1, -2)


@lru_cache(maxsize=4096)
def diagrams(p, h, k1, k2, e1, e2):
    K = va(k1, k2)
    t = vs(p, K)
    out = []
    # Index j labels the B emitted first, adjacent to the Yukawa vertex.
    for kj, ki, ej, ei in ((k1, k2, e1, e2), (k2, k1, e2, e1)):
        nj = jchain(seed_jet(va(h, K)), fprop_jet(K), vertex_jet(ei),
                    fprop_jet(kj), vertex_jet(ej), fprop_jet(Z))
        out.append(pole(nj, (Z, kj, K, p), 4))
    for kj, ki, ej, ei in ((k1, k2, e1, e2), (k2, k1, e2, e1)):
        u, w = vs(vm(2, p), kj), vs(vs(vm(2, p), vm(2, kj)), ki)
        nj = jchain(seed_jet(h), fprop_jet(Z),
                    scalar_vertex_jet(u, ej), scalar_vertex_jet(w, ei))
        out.append(pole(nj, (Z, p, vs(p, kj), t), 4))
    for ki, kj, ei, ej in ((k1, k2, e1, e2), (k2, k1, e2, e1)):
        nj = jchain(seed_jet(va(h, ki)), fprop_jet(ki), vertex_jet(ei),
                    fprop_jet(Z), scalar_vertex_jet(vs(vm(2, p), kj), ej))
        out.append(pole(nj, (Z, ki, p, vs(p, kj)), 4))
    # SG includes eps1.eps2 but not its separate phase/identical-field -2.
    nj = jchain(seed_jet(h), fprop_jet(Z))
    out.append(scale(dot(e1, e2), pole(nj, (Z, p, t), 3)))
    return out


def total(parts, weights=WEIGHTS):
    return add(*(scale(w, part) for w, part in zip(weights, parts)))


def closed_diagrams(p, h, k1, k2, e1, e2):
    """Individual analytic angular results, separately from target action."""
    out = []
    K = va(k1, k2)
    for b, a, eb, ea in ((k1, k2, e1, e2), (k2, k1, e2, e1)):
        d = vs(va(b, p), h)
        out.append(add(scale(F(1, 2), chain(slash(K), slash(ea), slash(eb))),
                       scale(F(1, 2), chain(slash(eb), slash(b), slash(ea))),
                       scale(-F(1, 6), add(chain(slash(d), slash(ea), slash(eb)),
                                           chain(slash(eb), slash(d), slash(ea)),
                                           chain(slash(ea), slash(eb), slash(d))))))
    for b, a, eb, ea in ((k1, k2, e1, e2), (k2, k1, e2, e1)):
        u = vs(vm(2, p), b)
        w = vs(vs(vm(2, p), vm(2, b)), a)
        d = vs(vs(vs(vm(3, p), vm(2, b)), a), h)
        out.append(add(scale(-F(1, 2) * dot(u, eb), slash(ea)),
                       scale(-F(1, 2) * dot(w, ea), slash(eb)),
                       scale(F(1, 3), add(scale(dot(ea, eb), slash(d)),
                                           scale(dot(d, ea), slash(eb)),
                                           scale(dot(d, eb), slash(ea))))))
    for a, b, ea, eb in ((k1, k2, e1, e2), (k2, k1, e2, e1)):
        u, d = vs(vm(2, p), b), va(h, vm(-2, p), b)
        out.append(add(scale(-F(1, 2) * dot(u, eb), slash(ea)),
                       scale(F(1, 2), chain(slash(a), slash(ea), slash(eb))),
                       scale(F(1, 3), add(scale(dot(ea, eb), slash(d)),
                                           scale(dot(ea, d), slash(eb)),
                                           scale(-2 * dot(eb, d), slash(ea))))))
    out.append(scale(F(1, 2) * dot(e1, e2), slash(va(vm(2, p), vm(-1, K), vm(-1, h)))))
    return out


def action(p, k1, k2, e1, e2):
    return add(scale(dot(e1, e2), slash(p)),
               scale(-F(1, 4), add(mul(comm(slash(k1), slash(e1)), slash(e2)),
                                   mul(comm(slash(k2), slash(e2)), slash(e1)))))


def one_b(p, h, k, eps, m2):
    """Previously audited one-B pole; used only for longitudinal controls."""
    r = vs(vs(p, h), k)
    gj = add(scale(-F(1, 3) * dot(va(vm(2, r), k), eps)
                   + F(1, 6) * dot(va(vm(2, h), k), eps), slash(p)),
             scale((dot(p, p) / 6 - m2 / 2), slash(eps)))
    dt = add(scale(F(1, 4), mul(comm(slash(k), slash(eps)), slash(r))),
             scale(F(1, 3), add(scale(dot(k, eps), slash(k)),
                                 scale(-dot(k, k), slash(eps)))))
    return add(gj, dt)


# Gaussian rationals and two-dimensional weak matrices, entirely exact.
def ga(r=0, i=0):
    return (F(r), F(i))


def gadd(*xs):
    return (sum(x[0] for x in xs), sum(x[1] for x in xs))


def gm(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def gs(c, a):
    return (c*a[0], c*a[1])


ZERO_G = ga()
I2 = ((ga(1), ga()), (ga(), ga(1)))
T = (
    ((ga(),ga(F(1,2))), (ga(F(1,2)),ga())),
    ((ga(),ga(0,-F(1,2))), (ga(0,F(1,2)),ga())),
    ((ga(F(1,2)),ga()), (ga(),ga(-F(1,2)))),
)


def ma(*xs):
    return tuple(tuple(gadd(*(x[i][j] for x in xs)) for j in range(2)) for i in range(2))


def ms(c, x):
    return tuple(tuple(gs(c,x[i][j]) for j in range(2)) for i in range(2))


def mm(a, b):
    return tuple(tuple(gadd(*(gm(a[i][k], b[k][j]) for k in range(2))) for j in range(2)) for i in range(2))


def mt(a):
    return tuple(tuple(a[j][i] for j in range(2)) for i in range(2))


def contract_seed(cchain, schain, crossed=True):
    """Literal K_{ab;cd} C_{be} S_{df} delta_{ef}, not a product shortcut."""
    out = []
    for a in range(2):
        row = []
        for c in range(2):
            s = ga()
            for b,d,e,f in product(range(2), repeat=4):
                pairing = (a == d and b == c) if crossed else (a == c and b == d)
                if pairing and e == f:
                    s = gadd(s,gm(cchain[b][e], schain[d][f]))
            row.append(s)
        out.append(tuple(row))
    return tuple(out)


def colors(m1, m2, wrong_transpose=False):
    c1 = ms(-1, m1 if wrong_transpose else mt(m1))
    c2 = ms(-1, m2 if wrong_transpose else mt(m2))
    anti = ma(mm(m1,m2),mm(m2,m1))
    return (
        contract_seed(mm(c2,c1),I2), contract_seed(mm(c1,c2),I2),
        contract_seed(I2,mm(m2,m1)), contract_seed(I2,mm(m1,m2)),
        contract_seed(c1,m2), contract_seed(c2,m1),
        ms(-1,contract_seed(I2,anti)),
    )


def tensor(color, cliff):
    return clean({(a,b,z,k): color[a][b][z]*v for a,b,z in product(range(2),repeat=3)
                  for k,v in cliff.items()})


def colored(parts, cs):
    return add(*(tensor(c, part) for c, part in zip(cs,parts)))


def ordered(parts):
    ff12,ff21,ss12,ss21,fs12,fs21,sg = parts
    return (add(ff12,ss21,scale(-1,fs21),scale(-1,sg)),
            add(ff21,ss12,scale(-1,fs12),scale(-1,sg)))


def action_pieces(p,h,k1,k2,e1,e2):
    r = vs(vs(vs(p,h),k1),k2)
    sigma_nonlinear = scale(F(1,4), mul(comm(slash(e1),slash(e2)),slash(r)))
    v = va(vm(dot(va(vm(2,k1),k2),e2),e1),
           vm(-dot(va(k1,vm(2,k2)),e1),e2),
           vm(dot(e1,e2),vs(k2,k1)))
    # Split partial F_nonlin and i[A,F_linear] to independently mutate them.
    partial = scale(F(1,3),slash(va(vm(dot(va(k1,k2),e2),e1),
                                        vm(-dot(va(k1,k2),e1),e2))))
    connection = add(scale(F(1,3),slash(v)),scale(-1,partial))
    return sigma_nonlinear, partial, connection


def action_ordered(p,h,k1,k2,e1,e2,mutant=None):
    sigma, partial, connection = action_pieces(p,h,k1,k2,e1,e2)
    if mutant == 'omit_nonlinear_curvature':
        sigma,partial = {},{}
    elif mutant == 'reverse_curvature_sign':
        sigma,partial = scale(-1,sigma),scale(-1,partial)
    elif mutant == 'omit_adjoint_connection':
        connection = {}
    extra = add(sigma,partial,connection)
    base = scale(F(1,2)*dot(e1,e2),slash(p))
    a = add(base,scale(-F(1,4),mul(comm(slash(k2),slash(e2)),slash(e1))),extra)
    b = add(base,scale(-F(1,4),mul(comm(slash(k1),slash(e1)),slash(e2))),scale(-1,extra))
    return a,b


def action_color(p,h,k1,k2,e1,e2,m1,m2,mutant=None):
    a,b = action_ordered(p,h,k1,k2,e1,e2,mutant)
    return add(tensor(mm(m1,m2),a),tensor(mm(m2,m1),b))


def ward(p,h,k1,k2,e2,m1,m2,m2scalar,omit_comm=False):
    first=tensor(mm(m2,m1),one_b(p,va(h,k1),k2,e2,m2scalar))
    second=tensor(mm(m1,m2),one_b(p,h,k2,e2,m2scalar))
    commterm={} if omit_comm else tensor(ma(mm(m1,m2),ms(-1,mm(m2,m1))),
                                        one_b(p,h,va(k1,k2),e2,m2scalar))
    return add(first,scale(-1,second),commterm)


CHECKS=0

def check(condition,label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


def encode(a):
    return {str(k):str(v) for k,v in sorted(a.items())}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--samples',type=int,default=16)
    parser.add_argument('--skip-basis',action='store_true')
    args=parser.parse_args()
    groups=[('W%dW%d'%(i+1,j+1),a,b) for i,a in enumerate(T) for j,b in enumerate(T)]
    y=ms(-F(1,2),I2)
    groups += [('BW%d'%(i+1),y,a) for i,a in enumerate(T)]
    groups += [('W%dB'%(i+1),a,y) for i,a in enumerate(T)]
    groups += [('BB',y,y)]
    for name,m1,m2 in groups:
        a,b=mm(m1,m2),mm(m2,m1)
        expected=(a,b,b,a,ms(-1,b),ms(-1,a),ms(-1,ma(a,b)))
        actual=colors(m1,m2)
        for n,c,d in zip(NAMES,actual,expected):
            check(c==d,name+' literal index contraction '+n)
        for cross in (False,True):
            cchain=mm(ms(-1,mt(m2)),ms(-1,mt(m1)))
            contracted=contract_seed(cchain,I2,cross)
            target=a if cross else ms(sum(a[i][i][0] for i in range(2)),I2)
            check(contracted==target,name+' FF seed pairing index contraction')
    rng=random.Random(0x57574220261015)
    mutations={n:0 for n in ('omit_nonlinear_curvature','reverse_curvature_sign',
                             'omit_adjoint_connection','omit_lower_Ward_kernel',
                             'wrong_fermion_transpose','reverse_FF_color_order',
                             'omit_seagull','wrong_mixed_charge_sign')}
    witnesses=[]
    noncomm=0
    for _ in range(args.samples):
        p,h,k1,k2,e1,e2=[tuple(F(rng.randint(-7,7),rng.randint(1,4)) for _ in range(4)) for _ in range(6)]
        parts=diagrams(p,h,k1,k2,e1,e2)
        check(parts==closed_diagrams(p,h,k1,k2,e1,e2),'all seven individual angular reductions')
        ords=ordered(parts)
        check(ords==action_ordered(p,h,k1,k2,e1,e2),'universal free-color action coefficients')
        for name,m1,m2 in groups:
            raw=colored(parts,colors(m1,m2))
            target=action_color(p,h,k1,k2,e1,e2,m1,m2)
            check(raw==target,name+' full matter integral/action equality')
            swapped=colored(diagrams(p,h,k2,k1,e2,e1),colors(m2,m1))
            check(raw==swapped,name+' Bose symmetry')
            for ki,kj,ei,m_i,m_j in ((k1,k2,e2,m1,m2),(k2,k1,e1,m2,m1)):
                contracted=colored(diagrams(p,h,ki,kj,ki,ei),colors(m_i,m_j))
                for mass in (F(-5,3),F(0),F(11,7)):
                    check(contracted==ward(p,h,ki,kj,ei,m_i,m_j,mass),name+' non-Abelian Ward recursion')
            is_noncomm = mm(m1,m2) != mm(m2,m1)
            if is_noncomm:
                noncomm+=1
                for mut in ('omit_nonlinear_curvature','reverse_curvature_sign','omit_adjoint_connection'):
                    mutations[mut] += raw != action_color(p,h,k1,k2,e1,e2,m1,m2,mut)
                ct=colored(diagrams(p,h,k1,k2,k1,e2),colors(m1,m2))
                mutations['omit_lower_Ward_kernel'] += ct != ward(p,h,k1,k2,e2,m1,m2,F(5,7),True)
                cs=list(colors(m1,m2));cs[0],cs[1]=cs[1],cs[0]
                mutations['reverse_FF_color_order'] += colored(parts,cs)!=target
            mutations['wrong_fermion_transpose'] += colored(parts,colors(m1,m2,True))!=target
            cs=list(colors(m1,m2));cs[-1]=ms(0,cs[-1])
            mutations['omit_seagull'] += colored(parts,cs)!=target
            cs=list(colors(m1,m2));cs[4]=ms(-1,cs[4]);cs[5]=ms(-1,cs[5])
            mutations['wrong_mixed_charge_sign'] += colored(parts,cs)!=target
            if len(witnesses)<2 and name in ('W1W2','BW2'):
                witnesses.append({'sector':name,'inputs':{n:list(map(str,v)) for n,v in zip(('p','h','k1','k2','eps1','eps2'),(p,h,k1,k2,e1,e2))},
                                  'ordered_coefficients':[encode(x) for x in ords],
                                  'full_kernel':encode(raw)})
    basis=0
    if not args.skip_basis:
        for slot,component,a,b in product(range(4),repeat=4):
            momenta=[Z,Z,Z,Z];momenta[slot]=E[component]
            p,h,k1,k2=momenta
            parts=diagrams(p,h,k1,k2,E[a],E[b])
            for n,raw,red in zip(NAMES,parts,closed_diagrams(p,h,k1,k2,E[a],E[b])):
                check(raw==red,n+' exhaustive angular/closed coefficient')
            check(ordered(parts)==action_ordered(p,h,k1,k2,E[a],E[b]),'exhaustive free-color coefficient')
            basis+=1
    for key,count in mutations.items():
        check(count>0,key+' has nonzero exact generic witness')
    for key in ('omit_nonlinear_curvature','reverse_curvature_sign','omit_adjoint_connection','omit_lower_Ward_kernel','reverse_FF_color_order'):
        check(mutations[key]==noncomm,key+' fails at every noncommuting generic case')
    print(json.dumps({'status':'PASS','checks':CHECKS,'generic_samples':args.samples,
                      'gauge_generator_pairs':len(groups),'noncommuting_generic_cases':noncomm,
                      'literal_color_contraction_checks':16*9,'complete_polynomial_points':basis,
                      'mutations_detected':mutations,'standalone':True,'repository_imports':False,
                      'route':'homogeneous large-q angular extraction plus literal weak-index sums',
                      'scope':'restricted seed/conjugate-Yukawa one-loop matter 1PI UV pole only',
                      'prior_angular_core':'adapted from independent verify_two_b_angular.py at 5c77bc5',
                      'witnesses':witnesses},indent=2))

if __name__=='__main__':
    main()
