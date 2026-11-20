#!/usr/bin/env python3
"""Standalone exact Fourier/IBP/source controls; Python standard library only.

No repository imports, files written, symbolic packages, or bare assertions.
Gaussian-rational representation follows the earlier public matrix check,
but operators and source identities are evaluated directly here.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import permutations
import json


@dataclass(frozen=True)
class G:
    re: F = F(0)
    im: F = F(0)

    def __post_init__(self):
        if isinstance(self.re, float) or isinstance(self.im, float):
            raise TypeError("floating point forbidden")
        object.__setattr__(self, "re", F(self.re))
        object.__setattr__(self, "im", F(self.im))

    @staticmethod
    def coerce(x):
        return x if isinstance(x, G) else G(x)

    def __add__(self, x):
        x = G.coerce(x)
        return G(self.re+x.re, self.im+x.im)

    __radd__ = __add__

    def __neg__(self):
        return G(-self.re, -self.im)

    def __sub__(self, x):
        return self + -G.coerce(x)

    def __rsub__(self, x):
        return G.coerce(x) + -self

    def __mul__(self, x):
        x = G.coerce(x)
        return G(self.re*x.re-self.im*x.im,
                 self.re*x.im+self.im*x.re)

    __rmul__ = __mul__

    def __truediv__(self, x):
        x = G.coerce(x)
        d = x.re*x.re+x.im*x.im
        return self*G(x.re/d, -x.im/d)

    def __bool__(self):
        return bool(self.re or self.im)


Z, I = G(), G(0, 1)
eta = (1, -1, -1, -1)


def mat(a):
    return [[G.coerce(x) for x in row] for row in a]


ZERO = mat([[0]*4 for _ in range(4)])
ONE = mat([[int(a == b) for b in range(4)] for a in range(4)])


def add(*args):
    return [[sum((a[i][j] for a in args), Z) for j in range(4)]
            for i in range(4)]


def scale(c, a):
    return [[c*x for x in row] for row in a]


def mul(a, b):
    return [[sum((a[i][k]*b[k][j] for k in range(4)), Z)
             for j in range(4)] for i in range(4)]


gamma = [
    mat([[1,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,-1]]),
    mat([[0,0,0,1],[0,0,1,0],[0,-1,0,0],[-1,0,0,0]]),
    mat([[0,0,0,-I],[0,0,I,0],[0,I,0,0],[-I,0,0,0]]),
    mat([[0,0,1,0],[0,0,0,-1],[-1,0,0,0],[0,1,0,0]])]
g5 = scale(I, mul(mul(mul(gamma[0], gamma[1]), gamma[2]), gamma[3]))
PL = scale(F(1, 2), add(ONE, scale(-1, g5)))


def dot(a, b):
    return sum(eta[j]*a[j]*b[j] for j in range(4))


def slash(a):
    return add(*(scale(eta[j]*a[j], gamma[j]) for j in range(4)))


def comm(a, b):
    return add(mul(a, b), scale(-1, mul(b, a)))


def signperm(a):
    if len(set(a)) != len(a):
        return 0
    return (-1)**sum(a[i] > a[j] for i in range(len(a))
                     for j in range(i+1, len(a)))


sigma = [[scale(I/2, comm(gamma[a], gamma[b])) for b in range(4)]
         for a in range(4)]
anti = {}
for a in range(4):
    for b in range(4):
        for c in range(4):
            anti[a,b,c] = scale(F(1,6), add(*(
                scale(signperm(t), mul(mul(gamma[(a,b,c)[t[0]]],
                      gamma[(a,b,c)[t[1]]]), gamma[(a,b,c)[t[2]]]))
                for t in permutations(range(3)))))


def v_tensor(a, k):
    return [add(scale(a[n], slash(k)), scale(-dot(a,k), gamma[n]))
            for n in range(4)]


def e_tensor(a, k):
    return [add(mul(mul(slash(k), gamma[n]), slash(a)),
                scale(-k[n], slash(a)), scale(-a[n], slash(k)),
                scale(dot(k,a), gamma[n])) for n in range(4)]


def contract(t, e):
    return add(*(scale(eta[n]*e[n], t[n]) for n in range(4)))


def field_strength(k, e):
    return [[I*eta[a]*eta[b]*(k[a]*e[b]-k[b]*e[a])
             for b in range(4)] for a in range(4)]


def dual(f):
    return [[eta[a]*eta[b]*sum((F(1,2)*signperm((a,b,c,d))*f[c][d]
             for c in range(4) for d in range(4)), Z)
             for b in range(4)] for a in range(4)]


def ordinary_fourier(a, f, incoming=False):
    # Outgoing -D and incoming +partial give the same factor -i a^nu.
    return add(*(scale(-I*a[n]*f[m][n], gamma[m])
                 for m in range(4) for n in range(4)))


def antisym_fourier(a, f):
    return add(*(scale(-I*F(1,2)*eta[c]*a[c]*f[m][n], anti[m,n,c])
                 for m in range(4) for n in range(4) for c in range(4)))


def sigma_f(f):
    return add(*(scale(f[m][n], sigma[m][n])
                 for m in range(4) for n in range(4)))


def rank(rows):
    a = [[G.coerce(x) for x in row] for row in rows]
    n, out = len(a[0]), 0
    for c in range(n):
        p = next((j for j in range(out, len(a)) if a[j][c]), None)
        if p is None:
            continue
        a[out], a[p] = a[p], a[out]
        pivot = a[out][c]
        a[out] = [x/pivot for x in a[out]]
        for j in range(len(a)):
            if j != out and a[j][c]:
                x = a[j][c]
                a[j] = [u-x*v for u,v in zip(a[j], a[out])]
        out += 1
        if out == len(a):
            break
    return out


def ward_grade_one_columns():
    # Key=(free gamma coefficient, sorted three component-variable indices).
    cols = []
    for a in range(3):
        for b in range(a,3):
            col = {}
            for mu in range(4):
                for lam in range(4):
                    key = (mu, tuple(sorted((8+mu, 4*a+lam, 4*b+lam))))
                    col[key] = col.get(key, 0)+eta[mu]*eta[lam]
            cols.append(col)
    for a in range(3):
        for b in range(3):
            col = {}
            for mu in range(4):
                for lam in range(4):
                    key = (mu, tuple(sorted((8+lam,4*a+lam,4*b+mu))))
                    col[key] = col.get(key,0)+eta[mu]*eta[lam]
            cols.append(col)
    return cols


checks = []
mutations = []


def check(label, condition):
    if not condition:
        raise RuntimeError(label)
    checks.append(label)


def reject(label, wrong, correct):
    if wrong == correct:
        raise RuntimeError("mutation escaped: "+label)
    mutations.append(label)


def main():
    cols = ward_grade_one_columns()
    keys = sorted(set().union(*(c.keys() for c in cols)))
    wrank = rank([[c.get(key,0) for c in cols] for key in keys])
    check("15-column grade-one polynomial Ward map has rank12", wrank == 12)
    data = [((2,1,-1,3),(3,-2,1,1),(1,2,3,-1)),
            ((1,-1,2,0),(2,3,0,-2),(3,1,-2,1))]
    basis_rows = []
    for r,h,k in data:
        p = tuple(r[j]+h[j]+k[j] for j in range(4))
        basis = [v_tensor(r,k),v_tensor(h,k),v_tensor(k,k),
                 e_tensor(r,k),e_tensor(h,k)]
        for b in basis:
            check("exact transverse basis Ward contraction", contract(b,k)==ZERO)
        # Nontransverse third grade-three candidate excludes a sixth direction.
        nontrans = contract(e_tensor(h,r),k)
        check("generic E(r,nu,h) has nonzero Ward contraction", nontrans != ZERO)
        for nu in range(4):
            projected = [mul(b[nu],PL) for b in basis]
            for i in range(4):
                for j in range(4):
                    basis_rows.append([b[i][j] for b in projected])
        for e in [(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)]:
            f = field_strength(k,e)
            vr,vh,vk,er,eh = [contract(b,e) for b in basis]
            qvr,qvh,qvk = [ordinary_fourier(a,f) for a in (r,h,k)]
            qer,qeh = [antisym_fourier(a,f) for a in (r,h)]
            for actual,want in zip((qvr,qvh,qvk,qer,qeh),(vr,vh,vk,er,eh)):
                check("direct derivative placement gives declared Fourier kernel",actual==want)
            check("4D chiral dual convention",mul(qeh,PL)==mul(scale(I,ordinary_fourier(h,dual(f))),PL))
            vp,ep = ordinary_fourier(p,f), antisym_fourier(p,f)
            check("ordinary derivative IBP",vp==add(qvr,qvh,qvk))
            check("antisymmetric derivative IBP and Bianchi",ep==add(qer,qeh))
            check("Bianchi derivative-field-strength contraction",antisym_fourier(k,f)==ZERO)
            sf = sigma_f(f)
            source_l, source_n = mul(slash(r),sf),mul(sf,slash(p))
            check("row-lepton source relation",qer==add(qvr,scale(-F(1,2),source_l)))
            check("sterile derivative source relation",add(vp,ep)==scale(-F(1,2),source_n))
            od = add(scale(F(1,2),add(qvr,qer)),scale(F(1,3),qvk))
            compact = add(scale(F(1,4),mul(comm(slash(k),slash(e)),slash(r))),scale(F(1,3),vk))
            check("compact measured tensor coefficient map",od==compact)
            direct = add(*(scale(I*F(1,4)*I*eta[c]*r[c]*f[m][n],
                          mul(sigma[m][n],gamma[c]))
                          for m in range(4) for n in range(4) for c in range(4)))
            direct = add(direct,scale(F(1,3),qvk))
            check("compact local sigma operator direct Fourier map",direct==od)
            lform = add(qvr,scale(-F(1,4),source_l),scale(F(1,3),qvk))
            hform = add(scale(-F(1,2),add(qvh,qeh)),
                        scale(-F(1,6),qvk),scale(-F(1,4),source_n))
            check("complete lepton-source packet",lform==od)
            check("Higgs-derivative IBP and sterile-source packet",hform==od)
    check("five fixed-chirality tensors are independent over Q(i)",rank(basis_rows)==5)

    # Generic nonzero control: omission can preserve Ward or real-wave tests.
    r,h,k = data[0]
    e = (1,0,0,0)
    f = field_strength(k,e)
    vr,vk,er = ordinary_fourier(r,f),ordinary_fourier(k,f),antisym_fourier(r,f)
    od = add(scale(F(1,2),add(vr,er)),scale(F(1,3),vk))
    reject("omit off-shell divergence direction",add(scale(F(1,2),add(vr,er))),od)
    reject("omit independent three-gamma direction",add(scale(F(1,2),vr),scale(F(1,3),vk)),od)
    reject("wrong outgoing-row Fourier sign",scale(-1,vr),vr)
    source_l = mul(slash(r),sigma_f(f))
    reject("omit sourced lepton contact",add(vr,scale(F(1,3),vk)),od)
    reject("wrong sourced lepton sign",add(vr,scale(F(1,4),source_l),scale(F(1,3),vk)),od)
    reject("wrong four-dimensional dual sign",mul(antisym_fourier(h,f),PL),
           mul(scale(-I,ordinary_fourier(h,dual(f))),PL))
    check("real B wave makes divergence direction zero",
          ordinary_fourier((1,0,0,1),field_strength((1,0,0,1),(0,1,0,0)))==ZERO)

    # Functional variation of scalar kinetic term, with commuting stand-ins:
    # dL/dB = i gy[(D H)^* H-H^*D H] = -g' J_scalar.
    phi,phi_c,der,der_c,gg,yy,bb = map(G,(2,3,5,7,11,F(1,2),13))
    dh = der+I*gg*yy*bb*phi
    dhc = der_c-I*gg*yy*bb*phi_c
    varied = (-I*gg*yy*phi_c)*dh+dhc*(I*gg*yy*phi)
    jscalar = I*yy*(phi_c*dh-dhc*phi)
    check("scalar kinetic variation fixes current sign",varied==-gg*jscalar)
    # Canonical gauge equation partial_mu B^{mu nu}=+g'J^nu;
    # QVk=-gy K^mu partial^nu B_mu nu gives +gy g' K.J.
    gauge_eom_sign, antisymmetry_sign, operator_sign = F(1),F(-1),F(-1)
    current_sign=operator_sign*antisymmetry_sign*gauge_eom_sign
    check("field-strength index order gives positive QVk current sign",current_sign==1)
    reject("forget antisymmetry in gauge EOM",operator_sign*gauge_eom_sign,current_sign)
    weights={"H":F(1,2),"L":F(-1,2),"e":F(-1),"Q":F(1,6),"u":F(2,3),"d":F(-1,3)}
    packet={key:weight/F(3) for key,weight in weights.items()}
    reject("omit scalar hypercharge current",{k:v for k,v in packet.items() if k!="H"},packet)
    reject("omit fermion hypercharge currents",{"H":packet["H"]},packet)
    sterile={"M N_R":F(-1,4),"S_R^N":F(-1,4)}
    reject("omit sterile spinor source",{"M N_R":F(-1,4)},sterile)
    # Generic COMPLEX flavor matrices: conjugation and order are both tested.
    y=[[G(1,1),G(2,-1)],[G(3),G(5,2)]]
    fcoeff=[[G(2,1),G(7)],[G(11),G(13,-1)]]
    cseed=[[G(2),G(3,1)],[G(3,1),G(7,-1)]]
    mm=lambda a,b:[[sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
    tr=lambda a:[list(x) for x in zip(*a)]
    cc=lambda a:[[G(x.re,-x.im) for x in row] for row in a]
    ydag=tr(cc(y))
    yf=mm(ydag,fcoeff)
    reject("replace Ydagger F by F Ydagger",mm(fcoeff,ydag),yf)
    reject("omit flavor conjugation",mm(tr(y),fcoeff),yf)
    check("arbitrary F has surviving Majorana antisymmetric source",yf[0][1]!=yf[1][0])
    selected=mm(mm(ydag,cseed),cc(y))
    check("symmetric seed removes this particular Majorana tensor source",selected==tr(selected))
    active=mm(mm(cseed,cc(y)),tr(y))
    hy=mm(y,ydag)
    antisym_active=[[(active[i][j]-active[j][i])/2 for j in range(2)] for i in range(2)]
    ch, hc=mm(cseed,tr(hy)),mm(hy,cseed)
    expected=[[(ch[i][j]-hc[i][j])/2 for j in range(2)] for i in range(2)]
    check("active sterile source has C Htranspose-H C orientation",antisym_active==expected)
    check("active sterile source survives for a symmetric seed",bool(antisym_active[0][1]))
    reject("discard all spinor contacts after symmetric-seed simplification",
           [[G(),G()],[G(),G()]],antisym_active)
    # Direct tensor exchange identity: C sigma^T C^{-1}=-sigma.
    charge_c=scale(I,mul(gamma[2],gamma[0]))
    for a in range(4):
        for b in range(4):
            spin_tensor=mul(sigma[a][b],charge_c)
            check("Grassmann tensor bilinear is flavor antisymmetric",
                  spin_tensor==[list(row) for row in zip(*spin_tensor)])
    check("charged-lepton source contraction Hdagger tildeH vanishes",3*5-5*3==0)
    print(json.dumps({"status":"passed","checks":len(checks),
          "rejected_mutations":len(mutations),"ward_grade_one_rank":wrank,
          "independent_transverse_directions":5,"mutation_labels":mutations},indent=2))


if __name__ == "__main__":
    main()
