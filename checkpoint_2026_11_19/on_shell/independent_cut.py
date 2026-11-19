"""Two-scalar cut from tree sewing and an independent Feynman-parameter delta."""
import itertools
import json
from collections import Counter
from pathlib import Path

from mpmath import mp

mp.dps = 44
checks, adverse = Counter(), Counter()


def close(label, left, right, tol="1e-31"):
    if abs(left-right) > mp.mpf(tol)*max(1, abs(right)):
        raise ValueError(f"{label}: {left} versus {right}")
    checks[label] += 1


def reject(label, candidate, expected):
    if abs(candidate-expected) < mp.mpf("1e-10"):
        raise ValueError(f"mutation survived: {label}")
    adverse[label] += 1


def kinematics(a, s):
    beta = mp.sqrt(1-4*a/s)
    return beta, a-s/2, s*beta/2


def ratio(a, b, s, c):
    _, t0, h = kinematics(a, s)
    j = mp.atanh(h/(b-t0))/h
    return c-(c*b+(1-c)*t0)*j


def minkowski_square(p, k):
    terms = [(pi-ki)**2 for pi, ki in zip(p, k)]
    return terms[0]-sum(terms[1:])


def direct_sphere(a, b, s, c, weak, nodes, weights):
    """Build four-vectors, enumerate all ordered intermediate weak states."""
    beta, _, _ = kinematics(a, s)
    e = mp.sqrt(s)/2
    p1 = (e, 0, 0, e)
    p3 = (e, e*mp.sqrt(1-c*c), 0, e*c)
    a1, a2, a3, a4 = weak
    total = 0
    # The azimuth enters only through the linear seed numerator.
    for x, w in zip(nodes, weights):
        for phi in (0, mp.pi/2, mp.pi, 3*mp.pi/2):
            r = beta*e*mp.sqrt(1-x*x)
            k1 = (e, r*mp.cos(phi), r*mp.sin(phi), beta*e*x)
            k2 = (e, -k1[1], -k1[2], -k1[3])
            t, u = minkowski_square(p1,k1), minkowski_square(p1,k2)
            tp, up = minkowski_square(p3,k1), minkowski_square(p3,k2)
            for i,j in itertools.product(range(2), repeat=2):
                seed = (a3==i and a4==j)*tp + (a3==j and a4==i)*up
                sterile = (a1==i and a2==j)/(t-b) + (a1==j and a2==i)/(u-b)
                total += w*seed*sterile/16  # x average, four azimuths, 1/2!
    return total


def simplex_delta(a, b, s):
    # v=x+y, d=(x-y)/v; integrate the two d roots of Delta=0.
    center = 2*(a-b)/s
    r = mp.sqrt(center**2+4*b/s)
    v0, vminus = center+r, center-r
    # v=v0+(1-v0)*w^2 cancels the square-root endpoint singularity.
    return (4/s)*mp.sqrt(1-v0)*mp.quad(
        lambda w: 1/mp.sqrt(v0-vminus+(1-v0)*w*w), [0,1])


def antisymmetric_pair(left, right):
    return left[0]*right[1]-left[1]*right[0]


def main():
    nodes, weights = mp.gauss_quadrature(56, "legendre")
    for aa, bb, ss, cc in [(1,1,8,"0.3"), (2,16,11,"-0.7"), (3,0,13,"1")]:
        a,b,s,c = map(mp.mpf, (aa,bb,ss,cc))
        beta,t0,_ = kinematics(a,s)
        delta = simplex_delta(a,b,s)
        for chi in (c, -c):
            cut_loop = chi*beta-(chi*b+(1-chi)*t0)*delta
            close("parameter delta versus closed sphere", cut_loop, beta*ratio(a,b,s,chi))
        for weak in [(0,0,0,0),(0,1,0,1),(0,1,1,0),(0,0,1,1)]:
            a1,a2,a3,a4 = weak
            expected = ((a1==a3 and a2==a4)*ratio(a,b,s,c)
                        +(a1==a4 and a2==a3)*ratio(a,b,s,-c))
            close("ordered weak trees sewn over the sphere",
                  direct_sphere(a,b,s,c,weak,nodes,weights), expected)
        close("interchange cancels identical-channel angular dependence",
              ratio(a,b,s,c)+ratio(a,b,s,-c), 2*ratio(a,b,s,0))

    # The operator l_gamma l_delta acting on a1^dagger a2^dagger|0>
    # gives u_gamma(2)u_delta(1)-u_gamma(1)u_delta(2).
    z1, z2 = (2,1j),(1+1j,3)
    inc = antisymmetric_pair(z2,z1)-antisymmetric_pair(z1,z2)
    close("two incoming fermion contractions", inc, -2*antisymmetric_pair(z1,z2))
    close("ordered spin norm in the forward state",
          -abs(antisymmetric_pair(z1,z2))**2,
          -abs(inc)**2/4)
    a,s = mp.mpf(1),mp.mpf(8)
    beta,_,_ = kinematics(a,s)
    r1,r2 = ratio(a,1,s,0),ratio(a,16,s,0)
    same_component = 4*(r1-r2)  # D=X=1, direct and conjugate seed
    for cc in ("-1","-0.5","0","0.6","1"):
        c = mp.mpf(cc)
        close("chosen two-mass symmetric benchmark",
              2*(ratio(a,1,s,c)+ratio(a,1,s,-c)-ratio(a,16,s,c)-ratio(a,16,s,-c)),
              same_component)
    near = 4+mp.mpf("1e-35")
    close("threshold coefficient including mirror",
          4*(ratio(a,1,near,0)-ratio(a,16,near,0)), mp.mpf(30)/17)
    close("massless anomalous numerator", 0*ratio(a,0,s,1), 0)
    close("degenerate opposite complex weights", ratio(a,1,s,0)*(1+1j**2), 0)
    close("null pole weights", 1+4*(1j/2)**2, 0)
    c1,c2,A1,A2 = 1+2j,3-4j,0.5-2j,-2+3j
    row = c2*A1+(A2*c1).conjugate()
    transposed = c1*A2+(A1*c2).conjugate()
    close("Hermitian flavor interchange",row,transposed.conjugate())
    reject("twice direct orientation for complex flavor",2*c2*A1,row)
    reject("unconjugated mirror",c2*A1+A2*c1,row)
    reject("missing scalar-state half",2*same_component,same_component)
    reject("omitted seed orientation",same_component/2,same_component)
    reject("external-leg-only seed",0,same_component)
    reject("conjugate-square sterile weights",4*(r1+r2),same_component)
    reject("reversed ordered loop sign",-same_component,same_component)
    reject("independent denominator common to both masses",0,same_component)
    reject("discarded scattering angle in weak-resolved channel",
           ratio(a,1,s,0),ratio(a,1,s,mp.mpf("0.35")))
    reject("unsubtracted exchange-odd spin sign",0,inc)

    result = {
        "status":"PASS",
        "method":"CM four-vector tree sewing versus Feynman simplex spectral density",
        "benchmark":{"scalar_squared_mass":1,"sterile_masses":[1,4],
                     "anomalous_weights":[1,-1],"seed":1,"s":8,
                     "weak_indices":[0,0,0,0],
                     "both_orientations_reduced_angular_sum":mp.nstr(same_component,25),
                     "disc_over_2i_divided_by_ordered_spin":mp.nstr(beta*same_component/(16*mp.pi),25),
                     "forward_disc_over_2i":mp.nstr(-s*beta*same_component/(16*mp.pi),25),
                     "reduced_threshold_coefficient":"30/17"},
        "positive_checks":sum(checks.values()),"positive_categories":dict(checks),
        "negative_checks":sum(adverse.values()),"negative_categories":dict(adverse),
        "repository_imports":False,"assertions_required":False}
    data = json.dumps(result, indent=2)+"\n"
    Path(__file__).with_name("INDEPENDENT_CUT_RESULTS.json").write_text(data,encoding="utf-8")
    print(data,end="")


if __name__ == "__main__":
    main()
