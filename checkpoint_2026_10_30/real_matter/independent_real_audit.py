#!/usr/bin/env python3
"""Independent first-jet/exterior audit of the complete selected matter map.

No project imports. Scalar identities are polynomial identities in arbitrary
real first jets. Odd fields and their derivative jets are independent exterior
generators. This checks classical first-order transport, not a new loop.
"""
from collections import Counter
from pathlib import Path
import json
import sympy as s

I = s.I
positive = Counter()
negative = Counter()


def expand(value):
    return s.expand(value)


class Exterior:
    mates = {}
    next_generator = 0

    def __init__(self, terms=None):
        self.terms = {}
        for key, value in (terms or {}).items():
            value = expand(value)
            if value != 0:
                self.terms[key] = value

    @classmethod
    def scalar(cls, value):
        return value if isinstance(value, cls) else cls({(): s.sympify(value)})

    @classmethod
    def pair(cls):
        left = cls.next_generator
        right = left + 1
        cls.next_generator += 2
        cls.mates[left] = right
        cls.mates[right] = left
        return cls({(left,): 1}), cls({(right,): 1})

    def __add__(self, other):
        terms = dict(self.terms)
        for key, value in self.scalar(other).terms.items():
            terms[key] = terms.get(key, 0) + value
        return Exterior(terms)

    __radd__ = __add__

    def __neg__(self):
        return Exterior({key: -value for key, value in self.terms.items()})

    def __sub__(self, other):
        return self + -self.scalar(other)

    def __rsub__(self, other):
        return self.scalar(other) + -self

    def __mul__(self, other):
        terms = {}
        for left, value in self.terms.items():
            for right, coefficient in self.scalar(other).terms.items():
                if set(left).intersection(right):
                    continue
                sign = (-1) ** sum(a > b for a in left for b in right)
                key = tuple(sorted(left + right))
                terms[key] = terms.get(key, 0) + sign * value * coefficient
        return Exterior(terms)

    __rmul__ = __mul__

    def __truediv__(self, other):
        return self * (s.Integer(1)/s.sympify(other))

    def dagger(self):
        result = Exterior()
        for key, value in self.terms.items():
            term = self.scalar(s.conjugate(value))
            for generator in reversed(key):
                term = term * Exterior({(self.mates[generator],): 1})
            result += term
        return result


def is_zero(value):
    if isinstance(value, Exterior):
        return not value.terms
    if isinstance(value, s.MatrixBase):
        return all(expand(entry) == 0 for entry in value)
    return expand(value) == 0


def check(name, difference):
    if not is_zero(difference):
        raise RuntimeError(name + ' failed')
    positive[name] += 1


def reject(name, difference):
    if is_zero(difference):
        raise RuntimeError(name + ' mutation escaped')
    negative[name] += 1


def total(terms):
    return sum(terms, Exterior())


def main():
    # All scalar identities are established as exact real-component polynomials.
    x, y, u, v, a, b, c, d, g, t = s.symbols('x y u v a b c d g t', real=True)
    ph = s.Matrix([x + I*y, u + I*v])
    dp = s.Matrix([a + I*b, c + I*d])
    rho = (ph.H*ph)[0]
    dr = (dp.H*ph + ph.H*dp)[0]
    current = I*((ph.H*dp)[0] - (dp.H*ph)[0])
    normjet = (dp.H*dp)[0]
    z = (ph.H*dp)[0]
    raw_scalar = (g + I*t)*(z*z-rho*normjet) + (g-I*t)*(s.conjugate(z)**2-rho*normjet)
    decomposed = g*(dr**2-current**2)/2 - 2*g*rho*normjet + t*dr*current
    check('general scalar pole decomposition', raw_scalar-decomposed)
    dph = rho*(g+I*t)*ph
    ddp = (g+I*t)*(dr*ph+rho*dp)
    kinetic_change = (ddp.H*dp+dp.H*ddp)[0]
    check('general scalar kinetic product rule', kinetic_change-(2*g*rho*normjet+g*dr**2-t*dr*current))
    residual = g*(3*dr**2-current**2)/2
    check('general transported scalar derivative packet', raw_scalar+kinetic_change-residual)
    check('general scalar radial norm change', (dph.H*ph+ph.H*dph)[0]-2*g*rho**2)
    eps = s.Matrix([[0, 1], [-1, 0]])
    higgs = -eps*ph.conjugate()
    dhiggs = -eps*dph.conjugate()
    check('inverse tilde Higgs convention', eps*higgs.conjugate()-ph)
    check('full Higgs radial and opposite phase', dhiggs-rho*(g-I*t)*higgs)
    reject('radial scalar coefficient halved', raw_scalar+kinetic_change/2-residual)
    reject('drop scalar norm-gradient residual', residual+g*current**2/2)
    reject('reverse scalar current residual', residual-g*(3*dr**2+current**2)/2)
    reject('give H the same phase as phi', dhiggs-rho*(g+I*t)*higgs)
    # A concrete nonzero-jet witness is recorded in the receipt.
    witness = {x:1,y:2,u:2,v:-1,a:2,b:-1,c:-1,d:2,g:5,t:5}
    rho_value = expand(rho.subs(witness))
    dr_value = expand(dr.subs(witness))
    current_value = expand(current.subs(witness))
    reject('zero norm-gradient witness', dr_value)
    reject('zero scalar-current witness', current_value)
    reject('zero radial residual witness', residual.subs(witness))
    # A polynomial potential is varied directly rather than using a kinetic EOM.
    q, e, m2, lam, w6, b3 = s.symbols('q e m2 lam w6 b3', real=True)
    potential = m2*q+lam*q**2+w6*q**3
    direct = s.diff(potential.subs(q,q+e*2*g*q*q),e).subs(e,0)
    check('general potential chain rule through rho cubed', direct-2*g*q*q*s.diff(potential,q))
    net_potential = 2*b3*q*q - 2*g*q*q*(m2+2*lam*q)
    check('raw potential pole plus declared reference variation', net_potential-(2*(b3-g*m2)*q*q-4*g*lam*q**3))
    reject('discard induced sextic scalar interaction', net_potential-2*(b3-g*m2)*q*q)

    # Independent complex Hermitian flavor matrices; their commutator is nonzero.
    GG = s.Matrix([[2,1+2*I],[1-2*I,3]])
    AA = s.Matrix([[1,2-I],[2+I,4]])
    gv = s.trace(GG)
    tv = s.trace(AA)
    F = GG/2+I*AA
    cv = gv+I*tv
    B = F+s.conjugate(cv)*s.eye(2)
    D = B.H
    check('full lepton map conjugacy', F.H-GG/2+I*AA)
    check('trace invariant of selected Hermitian parts', s.trace(GG+I*AA)-cv)
    reject('commuting-flavor witness', GG*AA-AA*GG)
    # Use simple complex scalar components, retaining arbitrary rho and dr in jets.
    phv = s.Matrix([1+2*I,2-I])
    hv = -eps*phv.conjugate()
    rr, dd = s.symbols('rho drho', real=True)
    dphv = rr*cv*phv
    dhv = -eps*dphv.conjugate()
    L={}; BL={}; DL={}; DBL={}; eta={}; beta={}
    for f in range(2):
        for weak in range(2):
            for spin in range(2):
                L[f,weak,spin], BL[f,weak,spin] = Exterior.pair()
                DL[f,weak,spin], DBL[f,weak,spin] = Exterior.pair()
                eta[f,weak,spin], beta[f,weak,spin] = Exterior.pair()
    def left_shift(fields, matrix, factor):
        return {(f,w,z):total(factor*matrix[f,h]*fields[h,w,z] for h in range(2)) for f in range(2) for w in range(2) for z in range(2)}
    def row_shift(fields, matrix, factor):
        return {(f,w,z):total(factor*fields[h,w,z]*matrix[h,f] for h in range(2)) for f in range(2) for w in range(2) for z in range(2)}
    dL=left_shift(L,F,rr); dBL=row_shift(BL,F.H,rr)
    ddL={key:value+left_shift(L,F,dd)[key] for key,value in left_shift(DL,F,rr).items()}
    ddBL={key:value+row_shift(BL,F.H,dd)[key] for key,value in row_shift(DBL,F.H,rr).items()}
    sigma=[s.eye(2),s.Matrix([[0,-1],[-1,0]]),s.Matrix([[0,I],[-I,0]]),s.diag(-1,1)]
    for mu, sig in enumerate(sigma):
        actual=Exterior(); target=Exterior(); wronghalf=Exterior()
        for f in range(2):
            for w in range(2):
                for z1 in range(2):
                    for z2 in range(2):
                        actual += I*sig[z1,z2]*(dBL[f,w,z1]*DL[f,w,z2]+BL[f,w,z1]*ddL[f,w,z2]-ddBL[f,w,z1]*L[f,w,z2]-DBL[f,w,z1]*dL[f,w,z2])/2
                        for h in range(2):
                            radial = I*rr*GG[f,h]*sig[z1,z2]*(BL[f,w,z1]*DL[h,w,z2]-DBL[f,w,z1]*L[h,w,z2])/2
                            phase = -dd*AA[f,h]*sig[z1,z2]*BL[f,w,z1]*L[h,w,z2]
                            target += radial+phase
                            wronghalf += radial/2+phase
        check('literal symmetric lepton kinetic variation '+str(mu),actual-target)
        reject('half lepton normalization '+str(mu),actual-wronghalf)
        check('kinetic variation is Hermitian '+str(mu),actual-actual.dagger())
    # Literal interaction monomials, each including independent conjugate fields.
    N={}; BN={}; E={}; BE={}; Q={}; BQ={}; U={}; BU={}; V={}; BV={}
    for f in range(2):
        for z0 in range(2):
            N[f,z0], BN[f,z0] = Exterior.pair()
            E[f,z0], BE[f,z0] = Exterior.pair()
    for w in range(2):
        for z0 in range(2): Q[w,z0], BQ[w,z0] = Exterior.pair()
    for z0 in range(2):
        U[z0], BU[z0] = Exterior.pair()
        V[z0], BV[z0] = Exterior.pair()
    Y=s.Matrix([[1+I,2-I],[3-2*I,-1+I]])
    Ye=s.Matrix([[2-I,1+I],[-1+2*I,3]])
    def yuk(row, col, scalar, matrix):
        return total(row[f,w,z0]*matrix[f,h]*scalar[w]*col[h,z0] for f in range(2) for h in range(2) for w in range(2) for z0 in range(2))
    nu_change = -yuk(dBL,N,phv,Y)-yuk(BL,N,dphv,Y)
    nu_target = -rr*yuk(BL,N,phv,(F.H+cv*s.eye(2))*Y)
    e_change = -yuk(dBL,E,hv,Ye)-yuk(BL,E,dhv,Ye)
    e_target = -rr*yuk(BL,E,hv,(F.H+s.conjugate(cv)*s.eye(2))*Ye)
    check('neutrino Yukawa direct variation',nu_change-nu_target)
    check('charged-lepton Yukawa direct variation',e_change-e_target)
    check('neutrino Yukawa with Hermitian conjugate',nu_change+nu_change.dagger()-nu_target-nu_target.dagger())
    check('charged-lepton Yukawa with Hermitian conjugate',e_change+e_change.dagger()-e_target-e_target.dagger())
    reject('put lepton flavor matrix after Yukawa',nu_change+rr*yuk(BL,N,phv,Y*(F.H+cv*s.eye(2))))
    reject('drop charged-lepton source',e_change)
    def quark(scalar,col):
        return total(BQ[w,z0]*scalar[w]*col[z0] for w in range(2) for z0 in range(2))
    up_change=-(1+2*I)*quark(dphv,U)
    down_change=-(2-I)*quark(dhv,V)
    check('up-quark Yukawa direct variation',up_change+rr*cv*(1+2*I)*quark(phv,U))
    check('down-quark Yukawa direct variation',down_change+rr*s.conjugate(cv)*(2-I)*quark(hv,V))
    reject('discard radial up-quark variation',up_change+I*rr*tv*(1+2*I)*quark(phv,U))
    reject('discard radial down-quark variation',down_change-I*rr*tv*(2-I)*quark(hv,V))
    probe=total(beta[f,w,z0]*dL[f,w,z0]+dBL[f,w,z0]*eta[f,w,z0] for f in range(2) for w in range(2) for z0 in range(2))
    probe_target=rr*total(beta[f,w,z0]*F[f,h]*L[h,w,z0]+BL[f,w,z0]*F.H[f,h]*eta[h,w,z0] for f in range(2) for h in range(2) for w in range(2) for z0 in range(2))
    check('lepton independent probe direct variation',probe-probe_target)
    reject('drop independent lepton probe variation',probe)
    jrow=s.Matrix([[2-I,3+I]])
    scalar_probe=(jrow*dphv)[0]+s.conjugate((jrow*dphv)[0])
    scalar_probe_target=rr*(cv*(jrow*phv)[0]+s.conjugate(cv)*s.conjugate((jrow*phv)[0]))
    check('scalar independent probe direct variation',scalar_probe-scalar_probe_target)
    reject('discard scalar radial probe term',scalar_probe-I*rr*tv*((jrow*phv)[0]-s.conjugate((jrow*phv)[0])))
    # Direct exterior variation of the actual barX C5 Xc/2 monomial.
    C5=s.Matrix([[1+I,2-I],[2-I,-1+2*I]])
    xb={(f,z0):total(BL[f,w,z0]*phv[w] for w in range(2)) for f in range(2) for z0 in range(2)}
    dxb={(f,z0):total(dBL[f,w,z0]*phv[w]+BL[f,w,z0]*dphv[w] for w in range(2)) for f in range(2) for z0 in range(2)}
    xc={(f,z0):total(-eps[z0,z1]*xb[f,z1] for z1 in range(2)) for f in range(2) for z0 in range(2)}
    dxc={(f,z0):total(-eps[z0,z1]*dxb[f,z1] for z1 in range(2)) for f in range(2) for z0 in range(2)}
    change=total(C5[f,h]*(dxb[f,z0]*xc[h,z0]+xb[f,z0]*dxc[h,z0])/2 for f in range(2) for h in range(2) for z0 in range(2))
    def insertion(matrix):
        return rr*total(matrix[f,h]*xb[f,z0]*xc[h,z0]/2 for f in range(2) for h in range(2) for z0 in range(2))
    matrix=D*C5+C5*D.T
    target=insertion(matrix)
    check('conditional C5 literal exterior variation',change-target)
    check('conditional C5 variation with Hermitian conjugate',change+change.dagger()-target-target.dagger())
    check('conditional C5 coefficient remains symmetric',matrix-matrix.T)
    reject('replace charge-conjugate transpose by adjoint',change-insertion(D*C5+C5*D.H))
    reject('omit scalar map from conditional C5 insertion',change-insertion(F.H*C5+C5*F.conjugate()))
    reject('omit lower-order insertion half',change-2*target)
    # A map of order kappa acting on an already one-loop action is order kappa^2.
    loop, shift=s.symbols('loop shift')
    check('one-loop four-lepton packet unaffected at first order',s.diff(e*loop*(1+e*shift)**4,e).subs(e,0)-loop)
    result={
        'status':'PASS',
        'input_commit':'5094584ebf402736b0bc5b597f6769f640178fd1',
        'scope':'Independent classical first-order transport of selected real and imaginary matter-pole pieces, declared interactions/probes and conditional C5 insertion; no new loop, quantum map or total physical coefficient.',
        'positive_checks':sum(positive.values()),
        'negative_checks':sum(negative.values()),
        'positive_categories':dict(positive),
        'negative_categories':dict(negative),
        'scalar_witness':{'rho':str(rho_value),'norm_gradient':str(dr_value),'scalar_current':str(current_value),'g':5,'t':5,'residual':str(expand(residual.subs(witness)))},
        'repository_imports':False,
        'assertions_required':False,
    }
    Path(__file__).with_name('INDEPENDENT_REAL_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
