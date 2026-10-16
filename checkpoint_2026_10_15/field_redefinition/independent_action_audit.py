#!/usr/bin/env python3
"""Independent local audit: no imports from the research repository.

Checks action-variation, Grassmann, chirality, and Fourier-IBP lemmas needed
to turn the already-derived equation residual into a field redefinition.
This is not a new loop calculation or a complete EFT canonicalization.
"""
import json
import sympy as s

I = s.I
R = s.Rational
checks = []
controls = []

def check(name, condition):
    if not condition:
        raise RuntimeError(name)
    checks.append(name)

def control(name, condition):
    if condition:
        raise RuntimeError("negative control was not rejected: " + name)
    controls.append(name)

def zero_mat(x):
    return all(s.expand(z) == 0 for z in x)

class G:
    """Exact finite exterior algebra; monomials use ordered generator ids."""
    def __init__(self, terms=None):
        self.t = {m: s.expand(c) for m, c in (terms or {}).items() if c != 0}
        self.t = {m: c for m, c in self.t.items() if c != 0}
    @staticmethod
    def coerce(x):
        return x if isinstance(x, G) else G({(): s.sympify(x)})
    def __add__(self, other):
        d = self.t.copy()
        for m, c in G.coerce(other).t.items():
            d[m] = d.get(m, 0) + c
        return G(d)
    __radd__ = __add__
    def __neg__(self):
        return G({m: -c for m, c in self.t.items()})
    def __sub__(self, other):
        return self + -G.coerce(other)
    def __rsub__(self, other):
        return G.coerce(other) - self
    def __mul__(self, other):
        d = {}
        for a, x in self.t.items():
            for b, y in G.coerce(other).t.items():
                if set(a) & set(b):
                    continue
                inversions = sum(j < i for i in a for j in b)
                m = tuple(sorted(a + b))
                d[m] = d.get(m, 0) + (-1)**inversions * x * y
        return G(d)
    __rmul__ = __mul__
    def dagger(self):
        # Hermitian conjugation reverses every product of odd fields.
        return G({m: (-1)**(len(m)*(len(m)-1)//2)*s.conjugate(c)
                  for m, c in self.t.items()})
    def derivative(self, j, side):
        d = {}
        for m, c in self.t.items():
            if j in m:
                k = m.index(j)
                parity = k if side == "left" else len(m)-1-k
                d[m[:k] + m[k+1:]] = (-1)**parity*c
        return G(d)
    def __eq__(self, other):
        return not (self - other).t

def theta(j):
    return G({(j,): 1})

def gm(A, v):
    return [sum((A[i,j]*v[j] for j in range(A.cols)), G())
            for i in range(A.rows)]

def row(v, A):
    return [sum((v[i]*A[i,j] for i in range(A.rows)), G())
            for j in range(A.cols)]

def dot(a, b):
    return sum((x*y for x,y in zip(a,b)), G())

pauli = [s.Matrix([[0,1],[1,0]]), s.Matrix([[0,-I],[I,0]]),
         s.diag(1,-1)]
z2 = s.zeros(2)
g0 = s.BlockMatrix([[z2,s.eye(2)],[s.eye(2),z2]]).as_explicit()
gam = [g0] + [s.BlockMatrix([[z2,p],[-p,z2]]).as_explicit() for p in pauli]
g5 = I*gam[0]*gam[1]*gam[2]*gam[3]
PL, PR = (s.eye(4)-g5)/2, (s.eye(4)+g5)/2
C = I*gam[2]*gam[0]
B = C*g0.T

def slash(p):
    return p[0]*gam[0] - sum((p[i]*gam[i] for i in range(1,4)), s.zeros(4))

def adj(v):
    return row([x.dagger() for x in v],g0)

def cc(v):
    return gm(B,[x.dagger() for x in v])

def veq(a,b):
    return all(x == y for x,y in zip(a,b))

check("charge conjugation matrix", zero_mat(C.T+C) and zero_mat(B*B.conjugate()-s.eye(4)))
for mu in range(4):
    check("Clifford chirality %d"%mu, zero_mat(PL*gam[mu]-gam[mu]*PR))
    check("charge conjugation gamma %d"%mu,
          zero_mat(C*gam[mu].T*C.inv()+gam[mu]))

# A fully generic complex Grassmann spinor W, and independent Majorana E.
W = [theta(2*j)+I*theta(2*j+1) for j in range(4)]
n = [theta(8)+I*theta(9), theta(10)+I*theta(11)]
E = [n[1].dagger(),-n[0].dagger(),n[0],n[1]]
check("Majorana E reality",veq(E,cc(E)))
check("Majorana scalar flip", dot(adj(W),E) == dot(adj(E),cc(W)))
check("Hermitian scalar conjugate", dot(adj(W),E).dagger() == dot(adj(E),W))
D = dot(adj(W),E)+dot(adj(W),E).dagger()
DeltaN = [-(w+c) for w,c in zip(W,cc(W))]
check("Majorana shift preserves reality",veq(DeltaN,cc(DeltaN)))
check("Majorana cancellation",dot(adj(E),DeltaN)+D == 0)
control("spurious Majorana half",dot(adj(E),[x*R(1,2) for x in DeltaN])+D == 0)
control("omit charge-conjugate shift",dot(adj(E),[-x for x in W])+D == 0)
control("charge-conjugate sign reversed",veq([w-c for w,c in zip(W,cc(W))],cc([w-c for w,c in zip(W,cc(W))])))

# Left and right derivatives differ for odd variables. A wrong ordering
# produces a measurable sign defect, rather than being papered over.
x, xb, eta = theta(12), theta(13), theta(14)
bil = xb*x
actual = xb*eta
check("right Grassmann derivative",bil.derivative(12,"right")*eta == actual)
check("left Grassmann derivative",eta*bil.derivative(12,"left") == actual)
control("left derivative used as right",bil.derivative(12,"left")*eta == actual)

# Check the 1/2 in a generic Majorana quadratic form by differentiating
# the action itself, not by postulating its variation.
A = s.Matrix([[0,2+I,3,-5],[-2-I,0,7,11],[-3,-7,0,13],[5,-11,-13,0]])
t = [theta(15+j) for j in range(4)]
d = [theta(19+j) for j in range(4)]
S = dot(t,gm(A,t))*R(1,2)
var_from_action = sum((d[j]*S.derivative(15+j,"left") for j in range(4)),G())
var_expected = dot(d,gm(A,t))
check("literal Majorana half cancels",var_from_action == var_expected)
control("retain half in first variation",var_from_action == var_expected*R(1,2))

# Literal scalar source and constrained tilde map.
phi, phic, m2, lam, left_source, right_source = s.symbols("phi phic m2 lam left_source right_source")
L_alg = -m2*phic*phi-lam*(phic*phi)**2-left_source*phi-phic*right_source
check("scalar Euler source sign",s.expand(s.diff(L_alg,phic)+m2*phi+2*lam*phic*phi**2+right_source)==0)
control("wrong scalar source sign",s.expand(s.diff(L_alg,phic)+m2*phi+2*lam*phic*phi**2-right_source)==0)
eps = s.Matrix([[0,1],[-1,0]])
tt = s.Matrix([[2+3*I,5-7*I]])
dphi_bar = -tt
dH = eps*tt.T
check("tilde scalar constrained conjugation",zero_mat(-(dH.T)*eps-dphi_bar))
control("independent H and tildeH shift sign",zero_mat(-((-dH).T)*eps-dphi_bar))

# Fourier IBP checks with generic off-shell momenta. The row residual
# carries +ir.x, phi +ih.x, and N -ip.x, with p=r+h.
momenta = [([2,-1,3,1],[1,2,-2,4]),([3,4,1,-2],[-1,2,3,5]),
           ([5,2,-3,1],[2,-4,1,-1]),([1,3,2,6],[4,-1,5,2])]
for trial,(r,h) in enumerate(momenta):
    p=[r[j]+h[j] for j in range(4)]
    rS,hS,pS = slash(r),slash(h),slash(p)
    literal = rS*pS/3-pS*pS/6
    covector = -hS*pS/3+pS*pS/6
    check("lepton IBP trial %d"%trial,zero_mat(literal-covector))
    control("omit derivative on phi %d"%trial,zero_mat(literal-pS*pS/6))
    control("wrong first IBP sign %d"%trial,zero_mat(literal-(hS*pS/3+pS*pS/6)))
    check("lepton left chirality %d"%trial,zero_mat(PR*covector*PL))
    # On a neutral outgoing row Q, K0Q=Q slash p, so two integrations
    # by parts preserve rather than reverse the sign of A^2.
    check("sterile second IBP %d"%trial,zero_mat(pS*pS-(p[0]**2-sum(x*x for x in p[1:]))*s.eye(4)))
    control("wrong second IBP sign %d"%trial,zero_mat(-pS*pS-pS*pS))

# Noncommuting flavor data makes right-acting mass order testable.
F=s.Matrix([[1+I,2-I],[3+2*I,-1+I]])
M=s.diag(2,5)
check("ordered mass identity",F*M*M == (F*M)*M)
control("commuted flavor mass",F*M*M == M*M*F)

# Gauge reality: an allowed shift is proportional to current+h.c.
z=s.Rational(2,3)+I*s.Rational(5,7)
realshift=z+s.conjugate(z)
check("gauge shift real",s.conjugate(realshift)==realshift)
control("gauge shift lacks h.c.",s.conjugate(z)==z)

result={"status":"pass","independent_checks":len(checks),
        "negative_controls_rejected":len(controls),"checks":checks,
        "negative_controls":controls,
        "scope":"Independent action-variation/IBP/Grassmann lemmas only; no loop recomputation or full EFT canonicalization."}
print(json.dumps(result,indent=2))
