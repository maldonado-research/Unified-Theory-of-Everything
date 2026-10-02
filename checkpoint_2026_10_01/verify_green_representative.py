#!/usr/bin/env python3
"""Exact algebra checks for GREEN_REPRESENTATIVE.md; no external packages."""
from dataclasses import dataclass
from fractions import Fraction as F
import json
import random

COUNT = 0


def check(condition, label):
    global COUNT
    COUNT += 1
    if not condition:
        raise RuntimeError(label)


@dataclass(frozen=True)
class G:
    """Gaussian rational: no floating point in Clifford or trace checks."""
    re: F = F(0)
    im: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, "re", F(self.re))
        object.__setattr__(self, "im", F(self.im))

    @staticmethod
    def coerce(x):
        return x if isinstance(x, G) else G(x)

    def __add__(self, x):
        x = self.coerce(x)
        return G(self.re + x.re, self.im + x.im)

    __radd__ = __add__

    def __neg__(self):
        return G(-self.re, -self.im)

    def __sub__(self, x):
        return self + -self.coerce(x)

    def __rsub__(self, x):
        return self.coerce(x) + -self

    def __mul__(self, x):
        x = self.coerce(x)
        return G(self.re*x.re-self.im*x.im,
                 self.re*x.im+self.im*x.re)

    __rmul__ = __mul__

    def conj(self):
        return G(self.re, -self.im)


I = G(0, 1)
Z = G()


def mat(rows):
    return [[G.coerce(x) for x in row] for row in rows]


ONE = mat([[int(a == b) for b in range(4)] for a in range(4)])
ZERO = mat([[0]*4 for _ in range(4)])


def add(a, b):
    return [[a[i][j]+b[i][j] for j in range(4)] for i in range(4)]


def scale(c, a):
    return [[c*x for x in row] for row in a]


def mul(a, b):
    return [[sum((a[i][k]*b[k][j] for k in range(4)), Z)
             for j in range(4)] for i in range(4)]


def trace(a):
    return sum((a[i][i] for i in range(4)), Z)


def adjoint(a):
    return [[a[j][i].conj() for j in range(4)] for i in range(4)]


g = [
    mat([[1,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,-1]]),
    mat([[0,0,0,1],[0,0,1,0],[0,-1,0,0],[-1,0,0,0]]),
    mat([[0,0,0,-I],[0,0,I,0],[0,I,0,0],[-I,0,0,0]]),
    mat([[0,0,1,0],[0,0,0,-1],[-1,0,0,0],[0,1,0,0]]),
]
metric = [1,-1,-1,-1]
g5 = scale(I, mul(mul(mul(g[0], g[1]), g[2]), g[3]))
PL = scale(F(1,2), add(ONE, scale(-1,g5)))
PR = scale(F(1,2), add(ONE,g5))


def dot(a,b):
    return sum(metric[j]*a[j]*b[j] for j in range(4))


def slash(v):
    answer = ZERO
    for mu in range(4):
        answer = add(answer, scale(metric[mu]*v[mu], g[mu]))
    return answer


def dirac_adjoint(a):
    return mul(mul(g[0], adjoint(a)), g[0])


def routed(p,h,m2):
    """Original independently expressed p,h tensor, before using r=p-h."""
    p2,h2,ph = dot(p,p),dot(h,h),dot(p,h)
    return [p[j]*(p2+m2-F(2,3)*ph+F(1,2)*h2)
            +h[j]*(F(1,6)*p2-F(1,2)*m2) for j in range(4)]


def representative(p,r,h,m2, row_sign=1, scalar_sign=1):
    """Evaluate six ordered Fourier monomials from their derivative action."""
    p2,r2,h2 = dot(p,p),dot(r,r),dot(h,h)
    coeff = (F(5,6),F(1,3),F(1,6),-F(1,6),F(1,2),F(1,2))
    operators = ([p2*x for x in p], [r2*x for x in p],
                 [scalar_sign*h2*x for x in p],
                 [row_sign*p2*x for x in r], [m2*x for x in p],
                 [row_sign*m2*x for x in r])
    return [sum(coeff[a]*operators[a][j] for a in range(6))
            for j in range(4)]


# Every ordinary derivative acts on its declared Fourier wave.
for q in [F(-3,7),F(0),F(2,5)]:
    check(I*(-I*q) == G(q), "incoming i partial sign")
    check((-I)*(I*q) == G(q), "outgoing -i partial sign")
    check(-(I*q)*(I*q) == G(q*q), "outgoing -partial squared")

rng = random.Random(20261002)
for case in range(256):
    r = [F(rng.randrange(-15,16),rng.randrange(1,8)) for _ in range(4)]
    h = [F(rng.randrange(-15,16),rng.randrange(1,8)) for _ in range(4)]
    p = [r[j]+h[j] for j in range(4)]
    m2 = F(rng.randrange(-10,11),rng.randrange(1,8))
    check(routed(p,h,m2) == representative(p,r,h,m2),
          "routed/representative equality %d" % case)
    reverse = representative([-x for x in p],[-x for x in r],
                             [-x for x in h],m2)
    check(reverse == [-x for x in routed(p,h,m2)],
          "simultaneous momentum reversal %d" % case)

p,h,m2 = [F(3),F(0),F(0),F(0)],[F(1),F(1),F(0),F(0)],F(0)
r = [p[j]-h[j] for j in range(4)]
check(routed(p,h,m2) == [F(45,2),F(3,2),F(0),F(0)],
      "public nonproportional counterexample")
check(representative(p,r,h,m2,row_sign=-1) != routed(p,h,m2),
      "mutation: incorrect row Fourier sign must fail")
# A non-null Higgs wave tests the sign of its Laplacian independently.
h_massive = [F(2),F(0),F(0),F(0)]
r_massive = [p[j]-h_massive[j] for j in range(4)]
check(representative(p,r_massive,h_massive,F(1),scalar_sign=-1)
      != routed(p,h_massive,F(1)),
      "mutation: incorrect scalar Laplacian sign must fail")

# Apply free external equations only after matching the full polynomial.
# bar u_L slash r = 0; h^2=m2 gives 5 p^2/6 + 2 m2/3.
for p2 in map(F, [-7,0,2,19]):
    for m2 in map(F, [-2,0,3]):
        projected = F(5,6)*p2+F(1,6)*m2+F(1,2)*m2
        check(projected == F(5,6)*p2+F(2,3)*m2,
              "massive scalar on-shell projection")
        if m2:
            check(projected != F(5,6)*p2+F(1,2)*m2,
                  "mutation: dropped scalar mass descendant must fail")
check(F(5,6)/F(3,2) == F(5,9), "massless light projection ratio")

# Symbolic sterile source tower: keys (field, derivative degree, M degree).
def deriv(poly):
    out = {}
    def term(key, value):
        out[key] = out.get(key,F(0))+value
    for (field,a,m), coefficient in poly.items():
        if field == "NL" and a == 0:
            term(("NR",0,m+1),coefficient)
            term(("SR",0,m),coefficient)
        elif field == "NR" and a == 0:
            term(("NL",0,m+1),coefficient)
            term(("SL",0,m),coefficient)
        else:
            term((field,a+1,m),coefficient)
    return out

lhs2 = deriv(deriv({("NL",0,0):F(1)}))
rhs2 = {("NL",0,2):F(1),("SL",0,1):F(1),("SR",1,0):F(1)}
lhs3 = deriv(lhs2)
rhs3 = {("NR",0,3):F(1),("SR",0,2):F(1),
        ("SL",1,1):F(1),("SR",2,0):F(1)}
check(lhs2 == rhs2, "sourced A squared")
check(lhs3 == rhs3, "sourced A cubed")
for removed in [("SR",0,2),("SL",1,1),("SR",2,0)]:
    mutation = dict(rhs3)
    del mutation[removed]
    check(lhs3 != mutation, "mutation: dropped sterile source contact")

# Exact Clifford check and curvature sign from the coefficient of each D_mu D_nu.
for mu in range(4):
    for nu in range(4):
        product = mul(g[mu],g[nu])
        reverse = mul(g[nu],g[mu])
        metric_term = scale(metric[mu] if mu == nu else 0,ONE)
        check(add(product,reverse) == scale(2,metric_term),
              "Clifford relation")
        sigma = scale(I*F(1,2),add(product,scale(-1,reverse)))
        # -D^2 - (1/2) sigma F, with F=-i[D,D].
        correct = add(scale(-1,metric_term),scale(I,sigma))
        check(scale(-1,product) == correct, "curvature sign")
        if mu != nu:
            wrong = add(scale(-1,metric_term),scale(-I,sigma))
            check(scale(-1,product) != wrong,
                  "mutation: reversed curvature sign must fail")
check(mul(PL,PL) == PL and mul(PR,PR) == PR, "chiral projectors")
for mu in range(4):
    check(mul(PR,mul(g[mu],PL)) == mul(g[mu],PL),
          "mixed bilinear chirality")

# Gauge-completion ambiguity: transverse vertex and nonzero on-shell trace.
p = [F(2),F(0),F(0),F(0)]
r = [F(1),F(0),F(0),F(1)]
h = k = [F(1,2),F(0),F(0),-F(1,2)]
eps = [F(0),F(1),F(0),F(0)]
check(p == [r[j]+h[j]+k[j] for j in range(4)], "four-point conservation")
check(dot(p,p) == 4 and dot(r,r) == dot(h,h) == dot(k,k) == 0,
      "four-point on-shell conditions")
check(dot(eps,k) == 0 and dot(eps,eps) == -1, "physical gauge polarization")
def gauge_kernel(p,k,eps):
    return mul(add(scale(dot(p,eps),slash(k)),
                   scale(-dot(p,k),slash(eps))),PL)
check(gauge_kernel(p,k,k) == ZERO, "Ward identity")
check(gauge_kernel(p,[F(0)]*4,eps) == ZERO, "zero gauge kernel")
vertex = gauge_kernel(p,k,eps)
density_l = mul(PL,slash(r))
density_n = add(slash(p),scale(F(2),ONE))
spin_sum = trace(mul(mul(mul(density_l,vertex),density_n),
                     dirac_adjoint(vertex)))
check(spin_sum == G(4), "nonzero on-shell ambiguity spin sum")

print(json.dumps({
    "status":"PASS",
    "checks":COUNT,
    "arithmetic":"exact rational and Gaussian-rational; standard library",
    "seed":20261002,
    "representative_coefficients":["5/6","1/3","1/6","-1/6","1/2","1/2"],
    "scalar_mass_descendant":"2/3",
    "gauge_ambiguity_spin_sum":"4",
    "scope":"local representative, sourced identities, nonunique completion; no beta function"
},indent=2))
