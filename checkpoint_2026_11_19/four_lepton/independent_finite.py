#!/usr/bin/env python3
"""Angular audit of the zero-momentum-subtracted primitive four-lepton graph."""
import json
from collections import Counter
from pathlib import Path

from mpmath import mp


mp.dps = 48
checks = Counter()
adverse = Counter()


def check(category, condition):
    if not condition:
        raise ValueError(f"failed {category}")
    checks[category] += 1


def close(category, left, right, tolerance="1e-35"):
    check(category, abs(left - right) <= mp.mpf(tolerance) * (1 + abs(right)))


def reject(category, survives):
    if survives:
        raise ValueError(f"mutation escaped {category}")
    adverse[category] += 1


def angular_difference(t, b, s):
    """S³ average of 1/((k+p)²+b) - 1/(k²+b)."""
    if s == 0 or (b == 0 and t >= s):
        return mp.zero
    baseline = t + b
    if b == 0:
        return (t - s) / (s * t)
    offset = baseline - s
    root = mp.sqrt(offset**2 + 4*b*s)
    den = baseline * (baseline + s + root)
    if offset < 0:
        return (offset - root) / den
    return -4*b*s / (den * (root + offset))


def radial(a, b, s):
    if s == 0:
        return mp.zero
    def integrand(t):
        if t == 0:
            return mp.zero
        return (t / (t + a))**2 * angular_difference(t, b, s)
    if b == 0:
        return mp.quad(integrand, [0, s])
    return mp.quad(integrand, sorted({0, a, b, s, mp.inf}))


def parameter(a, b, s, interior=(mp.mpf("0.5"),)):
    """Bubble identity k²/(k²+a)²=(1+a d/da)/(k²+a)."""
    def integrand(x):
        d = x*a + (1-x)*b
        if d == 0:
            return -mp.log1p(s/a) + s / (s+a)
        ratio = x*(1-x)*s / d
        return -mp.log1p(ratio) + (x*a / d) * ratio / (1+ratio)
    return mp.quad(integrand, [0, *interior, 1])


def massless(a, s):
    return 2 - (1 + 2*a/s) * mp.log1p(s/a)


def equal_mass(a, s):
    r = mp.sqrt(1 + 4*a/s)
    return mp.mpf("2.5") - (2*r + 2*a/(s*r))*mp.atanh(1/r)


def check_integrals():
    for t, b, s in [(2, 3, 7), (5, 0, 2), (1, 0, 4), (4, 0, 4)]:
        t, b, s = map(mp.mpf, (t,b,s))
        # Half-angle form handles the b=0, t=s angular endpoint.
        def angular(theta):
            c = mp.cos(theta/2)
            if t == s and b == 0:
                return mp.sin(theta/2)**2 / t
            den = (mp.sqrt(t)-mp.sqrt(s))**2 + b + 4*mp.sqrt(t*s)*c*c
            return mp.sin(theta)**2 / den
        direct = 2 * mp.quad(angular, [0, mp.pi/2, mp.pi]) / mp.pi
        close("direct angular quadrature", direct, 1/(t+b) + angular_difference(t,b,s))

    for values in [(2,3,7), (3,0,2), (3,3,4), (7,"0.001",9), (1,13,"0.1")]:
        a,b,s = map(mp.mpf, values)
        value = radial(a,b,s)
        close("radial versus independent bubble parameter", value, parameter(a,b,s))
        check("spacelike sign", value < 0)
        close("common squared-scale invariance", radial(11*a,11*b,11*s), value)
        check("zero momentum subtraction", radial(a,b,0) == 0)
    for s in map(mp.mpf, ["0.001","2","100"]):
        a = mp.mpf(3)
        close("massless closed form", radial(a,0,s), massless(a,s))
        close("equal-mass closed form", radial(a,a,s), equal_mass(a,s))
    s = mp.mpf("1e-5")
    for a,b in [(mp.mpf(2),mp.mpf(3)), (mp.mpf(3),mp.mpf(3))]:
        slope = -b * mp.quad(lambda t: t**2 / ((t+a)**2*(t+b)**3), [0,a,mp.inf])
        if a == b:
            analytic = -1/(12*a)
        else:
            analytic = b*(5*a+b)/(2*(a-b)**3) + a*b*(a+2*b)*mp.log(b/a)/(a-b)**4
        close("massive low-momentum slope", slope, analytic)
        close("low-momentum finite difference", radial(a,b,s)/s, analytic, "1e-6")
    close("equal-mass cubic after zero quadratic", (radial(1,1,s)+s/12)/s**3, mp.mpf(1)/840, "1e-8")
    close("massless leading quadratic", radial(1,0,s)/s**2, -mp.mpf(1)/6, "2e-6")
    close("continuity from positive mass", radial(2,mp.mpf("1e-12"),3), massless(mp.mpf(2),3), "1e-10")


def mix(phi, c, s):
    """Unitary with a complex phase and a real rotation."""
    return mp.matrix([[phi*c,phi*s],[-s,c]])


def spectral_matrix(h, distinct, values):
    """Lagrange projectors using only H and its distinct eigenvalues."""
    identity = mp.eye(h.rows)
    total = mp.zeros(h.rows)
    for eigenvalue in distinct:
        projector = identity
        for other in distinct:
            if other != eigenvalue:
                projector = projector * (h - other*identity) / (eigenvalue-other)
        total += projector * values[eigenvalue]
    return total


def check_masses():
    i = mp.j
    w = mix((3+4*i)/5, mp.mpf(5)/13, mp.mpf(12)/13)
    u = mix(i, mp.mpf(3)/5, mp.mpf(4)/5).T
    v = mix((5+12*i)/13, mp.mpf(8)/17, mp.mpf(15)/17)
    y = mp.matrix([[1+2*i,4-i],[-2+i,i]])
    a, s = mp.mpf(3), mp.mpf(2)
    weights = {b: radial(a,b,s) for b in [0,1,4]}
    for masses in [(1,2),(0,2),(2,2),(0,0)]:
        d = mp.diag(masses)
        m = w*d*w.T
        h = m*m.H
        distinct = sorted({int(mass**2) for mass in masses})
        numerator = spectral_matrix(h, distinct, weights)*m
        takagi = w*mp.diag([mass*weights[mass**2] for mass in masses])*w.T
        close("spectral functional calculus versus Takagi", mp.norm(numerator-takagi), 0)
        close("sterile tensor symmetry", mp.norm(numerator.T-numerator), 0)
        z = mp.mpf("0.7")
        direct = (z*mp.eye(2)+h)**-1 * m
        diagonal = w*mp.diag([mass/(z+mass**2) for mass in masses])*w.T
        close("full anomalous resolvent", mp.norm(direct-diagonal), 0)
        mp_rot = u.H*m*u.conjugate()
        yp = v.H*y*u.conjugate()
        tensor = y.conjugate()*numerator*y.H
        primed = yp.conjugate()*spectral_matrix(mp_rot*mp_rot.H, distinct, weights)*mp_rot*yp.H
        close("passive active and sterile covariance", mp.norm(primed-v.T*tensor*v), 0)
        if masses == (0,2):
            reject("mass adjoint orientation", mp.norm(numerator-spectral_matrix(h,distinct,weights)*m.H) < mp.mpf("1e-15"))
            reject("noncommuting denominator order", mp.norm(numerator-m*spectral_matrix(h,distinct,weights)) < mp.mpf("1e-15"))
        if masses == (1,2):
            reject("single shared finite kernel", mp.norm(numerator-weights[1]*m) < mp.mpf("1e-15"))
        if masses == (0,0):
            check("massless full triangle", mp.norm(tensor) == 0)


def check_distinct_thresholds():
    # Jacobi localization at the unique interior zero of the parameter
    # denominator. Scaling the threshold separation regulates the real integral.
    a = mp.mpf(1)
    for m in map(mp.mpf, [1,2]):
        h = mp.sqrt(a) + m
        x0 = m / h
        u0 = x0*(1-x0)
        delta = mp.mpf("1e-16")
        width = mp.sqrt(delta*u0) / h
        locations = (x0-8*width, x0, x0+8*width)
        near_threshold = parameter(a,m*m,-h*h+delta,locations)
        close("threshold singular coefficient", mp.sqrt(delta)*near_threshold,
              -mp.pi*a*x0/(h*mp.sqrt(u0)), "1e-7")
    # Active row Y*=(1,i/2), masses (1,4): the zero local mass
    # moment is exact, but the two finite functions have distinct slopes.
    check("null UV flavor source", 1 + (mp.j/2)**2*4 == 0)
    finite = radial(1,1,mp.mpf(1)) - radial(1,16,mp.mpf(1))
    check("nonzero finite kernel for null UV example", finite < -mp.mpf("0.01"))


def mirror_slope(a, b):
    if b == a:
        return 1/(2*a)
    if b == 0:
        return 1/a
    return (a-b + b*mp.log(b/a))/(a-b)**2


def spin_product(barred, unbarred, dagger=False):
    """Exterior components of P at two fixed weak indices."""
    terms = {}
    epsilon = [(0,1,1), (1,0,-1)]
    for r,t,e in epsilon:
        for v,w,f in epsilon:
            fields = [(1,barred[0],0,r), (1,barred[1],1,t),
                      (0,unbarred[0],0,v), (0,unbarred[1],1,w)]
            if dagger:
                fields = [(1-kind,flavor,weak,spin) for kind,flavor,weak,spin in reversed(fields)]
            swaps = sum(a > b for j,a in enumerate(fields) for b in fields[j+1:])
            terms[tuple(sorted(fields))] = e*f*(-1)**swaps
    return terms


def check_mirror():
    for a,b in [(1,1),(2,0),(3,7),(2,16)]:
        a,b = mp.mpf(a),mp.mpf(b)
        def integrand(t):
            if b == 0:
                return 1/(t+a)**2
            return t / ((t+a)**2*(t+b))
        direct = mp.quad(integrand, sorted({mp.zero,a,b,mp.inf}))
        close("mirror angular radial integral", direct, mirror_slope(a,b))
        check("positive mirror slope", direct > 0)
    k,p = mp.mpf(2),mp.mpf(3)
    for direction in [1,-1]:
        averaged = (2/mp.pi)*mp.quad(lambda th: mp.sin(th)**2 *
                                      (p**2+direction*2*k*p*mp.cos(th)), [0,mp.pi])
        close("two seed assignments have same mirror numerator", averaged,p**2)
    adjoint = spin_product((0,1),(2,3),dagger=True)
    exchanged = spin_product((2,3),(0,1))
    check("Grassmann adjoint exchanges the two flavor pairs", adjoint == exchanged)
    reject("minus under Hermitian exchange", adjoint == {k:-v for k,v in exchanged.items()})
    b = mp.mpf(16)
    f16 = -b * mp.quad(lambda t: t**2/((t+1)**2*(t+b)**3), [0,1,b,mp.inf])
    total_slope = -mp.mpf(1)/12 - f16 + mirror_slope(mp.mpf(1),1) - mirror_slope(mp.mpf(1),b)
    check("null UV witness finite slope including mirror", total_slope > mp.mpf("0.25"))


def main():
    check_integrals()
    check_masses()
    check_distinct_thresholds()
    check_mirror()
    # One Wick crossing. Phases use eta=(+---), with the final common
    # Wick rotation i stripped when translating the amplitude to the action.
    assignment = (-1j/2) * (-1j)**2 * 1j**2 * (-1j) * (-1) / 2
    check("ordered Wick factor", assignment == 1/4)
    check("paired P-action factor", 2*assignment == 1/2)
    reject("extraneous Majorana half", 2*assignment/2 == 1/2)
    reject("wrong Wick sign", -2*assignment == 1/2)
    reject("omitted crossed assignment", assignment == 1/2)
    reject("finite difference sign reversed", parameter(mp.mpf(2),3,7) > 0)
    result = {
        "status": "PASS",
        "method": "S3 angular integration versus subtracted bubble parameterization",
        "normalization": "F=16*pi^2*(I(Q^2)-I(0))",
        "benchmarks": [
            {"a": a, "b": b, "s": s,
             "F": mp.nstr(radial(*map(mp.mpf,(a,b,s))), 26)}
            for a,b,s in [(2,3,7),(3,0,2),(3,3,4)]
        ],
        "positive_checks": sum(checks.values()),
        "positive_categories": dict(checks),
        "negative_checks": sum(adverse.values()),
        "negative_categories": dict(adverse),
        "repository_imports": False,
        "assertions_required": False,
    }
    data = json.dumps(result, indent=2) + "\n"
    Path(__file__).with_name("INDEPENDENT_FINITE_RESULTS.json").write_text(data, encoding="utf-8")
    print(data, end="")


if __name__ == "__main__":
    main()
