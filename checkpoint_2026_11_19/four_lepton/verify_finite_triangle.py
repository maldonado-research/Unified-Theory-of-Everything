"""Verify the subtracted Euclidean four-lepton graph in its ordered channel."""

import json
from pathlib import Path

import mpmath as mp
import sympy as sp


positive = []
adverse = []


def check(label, condition):
    if not condition:
        raise ValueError(label)
    positive.append(label)


def reject(label, wrong_matches):
    if wrong_matches:
        raise ValueError("Undetected mutation: " + label)
    adverse.append(label)


def zero(expr):
    return sp.cancel(expr) == 0


def close(left, right, tol=mp.mpf("1e-35")):
    return abs(left - right) <= tol * max(1, abs(left), abs(right))


def parameter_form(a, b, z):
    def integrand(x):
        # At b=0 the endpoint x=0 is a removable limit.
        w = z * (1 - x) / a if b == 0 else x * (1 - x) * z / (a*x + b*(1-x))
        fraction = 1 if b == 0 else a*x/(a*x + b*(1-x))
        return -mp.log1p(w) + fraction*w/(1+w)

    return mp.quad(integrand, [0, 1])


def radial_form(a, b, z):
    # Here u=k^2. At b=0 the angular difference is supported on u<z.
    if z == 0:
        return mp.mpf(0)
    if b == 0:
        return mp.quad(lambda u: (u*u/z-u)/(u+a)**2, [0, z])

    def integrand(u):
        base = u+b
        # Rationalize the difference of the two angular averages at large u.
        root = mp.sqrt((base-z)**2+4*b*z)
        difference = -4*b*z/((base-z+root)*base*(base+z+root))
        return u**2/(u+a)**2 * difference

    return mp.quad(integrand, [0, z, mp.inf])


def slope_by_ratio(a, b):
    r = b/a
    if not r:
        return mp.mpf(0)
    delta = r-1
    if abs(delta) < mp.mpf("0.25"):
        series = mp.nsum(lambda n: (-delta)**n*(n+1)/((n+3)*(n+4)), [0, mp.inf])
        return -r*series/a
    return -r*(r**2+4*r-5-(4*r+2)*mp.log(r))/(2*a*delta**4)



def mirror_slope(a, b):
    r = b/a
    if r == 0:
        return 1/a
    delta = r-1
    if abs(delta) < mp.mpf("0.25"):
        return mp.nsum(lambda n: (-delta)**n/((n+1)*(n+2)), [0, mp.inf])/a
    return (1-r+r*mp.log(r))/(a*delta**2)


def algebra():
    a, b, x, z, u, r = sp.symbols("a b x z u r", positive=True)
    c = x*(1-x)
    origin = x*a + (1-x)*b
    shifted = origin+c*z
    log_ratio = sp.log(shifted/origin)
    diff = -log_ratio + a*x*c*z/(origin*shifted)
    check("double scalar propagator from the mass derivative",
          zero((1/(u+a) + a*sp.diff(1/(u+a), a))-u/(u+a)**2))
    check("subtracted bubble and derivative",
          zero(diff - (-log_ratio + a*sp.diff(-log_ratio, a))))
    shifted_numerator = 2*x*(-log_ratio + (1-x)**2*z/(2*shifted))
    check("direct shifted triangle agrees by an endpoint derivative",
          zero(shifted_numerator-diff-sp.diff(c*log_ratio, x)))
    check("total derivative has vanishing massive boundary",
          (c*log_ratio).subs(x, 0) == (c*log_ratio).subs(x, 1) == 0)
    check("dimensionless pole residue from two scalar lines",
          sp.integrate(2*x, (x, 0, 1)) == 1)
    check("Wick measure and signs of numerator and three denominators",
          sp.I*(-1)/(-1)**3 == sp.I)
    sign_kernel = -c*((1-x)*b+c*z)/shifted**2
    check("spacelike derivative has definite sign", zero(sp.diff(diff, z)-sign_kernel))
    check("degenerate-mass slope", sp.integrate(sign_kernel.subs({b: a, z: 0}), (x, 0, 1)) == -1/(12*a))
    # An explicit primitive for the dimensionless slope after h=1+(r-1)(1-x).
    h = sp.symbols("h", positive=True)
    primitive = -h**2/2 + (r+2)*h -(2*r+1)*sp.log(h) - r/h
    check("mass-ratio slope primitive", zero(sp.diff(primitive, h)-(h-1)**2*(r-h)/h**2))
    check("mass-ratio logarithm coefficients",
          zero(2*(primitive.subs(h, r)-primitive.subs(h, 1))-(r**2+4*r-5-(4*r+2)*sp.log(r))))
    massless_function = 2-(1+2/z)*sp.log(1+z)
    check("massless kernel starts at fourth momentum order",
          sp.series(massless_function, z, 0, 3).removeO().expand() == -z**2/6)
    dot = sp.symbols("k_dot_p", real=True)
    check("opposite mirror routings average to a pure z numerator difference",
          sp.expand(((u+2*dot+z)+(u-2*dot+z))/2-u-z) == 0)
    mirror_primitive = (origin-b*sp.log(origin))/(a-b)**2
    check("mirror primitive", zero(sp.diff(mirror_primitive, x)-x/origin))
    mirror_integral = mirror_primitive.subs(x, 1)-mirror_primitive.subs(x, 0)
    mirror_exact = (a-b-b*sp.log(a/b))/(a-b)**2
    check("mirror analytic slope at unequal positive masses", zero(sp.expand_log(mirror_integral-mirror_exact, force=True)))
    paired_density = x**2/origin + c**2*z/(origin*shifted) + c*a*x/shifted**2
    check("pointwise positive decomposition of the sum of scalar slopes",
          zero(x/origin+sign_kernel-paired_density))
    check("mirror equal-mass slope", sp.integrate((x/origin).subs(b, a), (x, 0, 1)) == 1/(2*a))
    check("two real ordered orientations at equal masses",
          sp.integrate((x/origin+sign_kernel).subs({b: a, z: 0}), (x, 0, 1)) == 5/(12*a))
    sample = {x: sp.Rational(1, 3), a: 2, b: 3, z: 5}
    reject("drop the doubled scalar propagator", zero((diff+log_ratio).subs(sample)))
    reject("omit Euclidean numerator shift in the direct triangle",
           zero((2*x*(-log_ratio)-diff-sp.diff(c*log_ratio, x)).subs(sample)))


def flavor():
    # Fully complex unitary matrix for M=U diag(m_i) U^T, and an independent
    # active transformation. Rational weights stand for arbitrary real kernels.
    I = sp.I
    U = sp.Matrix([[3, 4*I], [4*I, 3]])/5
    V = sp.diag(I, 1)*sp.Matrix([[3, -4], [4, 3]])/5
    D = sp.diag(0, 3)
    Y = sp.Matrix([[1+I, 2-I], [3*I, 2+I]])
    M = U*D*U.T
    check("complex Takagi change with a singular symmetric mass",
          U.H*U == sp.eye(2) and M == M.T and M.rank() == 1)
    t = sp.Rational(7, 3)
    # Direct comparison to the resolvent of the squared mass.
    lhs = (t*sp.eye(2)+M*M.H).inv()*M
    rhs = U*D*(t*sp.eye(2)+D**2).inv()*U.T
    check("exact complex anomalous propagator order", lhs.equals(rhs))
    reject("adjoint in the final Takagi slot", lhs.equals(U*D*(t*sp.eye(2)+D**2).inv()*U.H))
    reject("commute mass through its left squared-mass resolvent",
           lhs.equals(M*(t*sp.eye(2)+M*M.H).inv()))
    weights = sp.diag(sp.Rational(-1, 8), sp.Rational(-2, 5))
    tensor = Y.conjugate()*U*D*weights*U.T*Y.H
    Yp = V.H*Y*U.conjugate()
    check("four-lepton tensor is symmetric", tensor.equals(tensor.T))
    check("mass-basis rotation and active covariance",
          (Yp.conjugate()*D*weights*Yp.H).equals(V.T*tensor*V))
    other_weights = sp.diag(11, weights[1, 1])
    check("zero-mode weight can vary without affecting the triangle",
          (Y.conjugate()*U*D*other_weights*U.T*Y.H).equals(tensor))
    reject("substitute a Yukawa transpose for its adjoint",
           (Y.conjugate()*U*D*weights*U.T*Y.T).equals(tensor))
    # At fixed labels the mirror conjugates and exchanges the pair slots.
    seed = sp.Matrix([1+I, 2-I, 3+I])
    def pairs(mat):
        return sp.Matrix([mat[0, 0], sp.sqrt(2)*mat[0, 1], mat[1, 1]])
    direct = pairs(tensor)
    mirror = pairs(Y.conjugate()*U*D*sp.diag(2, 3)*U.T*Y.H)
    soft_barred = seed*direct.T + mirror.conjugate()*seed.H
    soft_unbarred = seed*mirror.T + direct.conjugate()*seed.H
    check("adjoint exchanges the two asymmetric momentum assignments",
          soft_barred.H.equals(soft_unbarred))
    reject("impose Hermiticity at one fixed asymmetric momentum assignment",
           soft_barred.H.equals(soft_barred))


def numerical():
    mp.mp.dps = 55
    cases = [("massless", 2, 0, 3), ("equal", 2, 2, 3),
             ("heavy", 2, 50, 7), ("light", 4, "0.04", 11),
             ("nearly equal", 3, "3.000003", "0.5")]
    values = {}
    for name, a, b, z in cases:
        a, b, z = map(mp.mpf, (a, b, z))
        value = parameter_form(a, b, z)
        radial_value = radial_form(a, b, z)
        check(name + " direct four-dimensional radial integral", close(value, radial_value))
        check(name + " negative subtracted kernel", value < 0)
        check(name + " strict decrease at positive momentum", parameter_form(a, b, 2*z) < value)
        check(name + " scale homogeneity", close(parameter_form(3*a, 3*b, 3*z), value))
        slope = slope_by_ratio(a, b)
        if b:
            slope_integral = mp.quad(lambda x: -b*x*(1-x)**2/(a*x+b*(1-x))**2, [0, 1])
            check(name + " analytic slope against exact-integrand quadrature", close(slope, slope_integral))
        else:
            s = z/a
            check("massless elementary expression", close(value, 2-(1+2/s)*mp.log1p(s)))
        mirror = mirror_slope(a, b)
        # Integral of the mirror numerator difference z over the three
        # unshifted propagators, with the independent 4D radial measure.
        if b:
            radial_mirror = mp.quad(lambda u: u/((u+a)**2*(u+b)), [0, a, mp.inf])
        else:
            radial_mirror = mp.quad(lambda u: 1/(u+a)**2, [0, a, mp.inf])
        check(name + " mirror slope against the original radial integral", close(mirror, radial_mirror))
        check(name + " positive subtracted sum of both scalar kernels", z*mirror+value > 0)
        reject(name + " double the direct graph at fixed external momenta", close(z*mirror, value))
        values[name] = {"a": str(a), "b": str(b), "Q_squared": str(z),
                        "F": mp.nstr(value, 24), "F_prime_0": mp.nstr(slope, 24),
                        "H": mp.nstr(z*mirror, 24), "H_prime_0": mp.nstr(mirror, 24)}
    reject("flip the relative Wick sign", close(-parameter_form(2, 2, 3), radial_form(2, 2, 3)))
    reject("use a mass-independent weight for split eigenmasses",
           close(parameter_form(2, 2, 3), parameter_form(2, 9, 3)))
    check("origin subtracts the same integrand", parameter_form(2, 2, 0) == 0)
    return values


if __name__ == "__main__":
    algebra()
    flavor()
    report = {"status": "PASS", "kinematics": "ordered external (0,0,+p,-p), spacelike",
              "normalization": "F = 16 pi^2 [I_E(Q^2) - I_E(0)]"}
    report["witnesses"] = numerical()
    report.update(positive_checks=positive, rejected_mutations=adverse)
    Path(__file__).with_name("FINITE_TRIANGLE_RESULTS.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Selected triangle: {len(positive)} checks passed, {len(adverse)} mutations rejected")
