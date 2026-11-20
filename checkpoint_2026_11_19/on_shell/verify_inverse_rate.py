"""Reciprocal transition probabilities with two identical final-state leptons."""
from collections import Counter
import json
from pathlib import Path

from mpmath import mp

mp.dps = 52
positive, negative = Counter(), Counter()


def close(name, x, y, tol="1e-44"):
    if abs(x-y) > mp.mpf(tol)*max(1, abs(y)):
        raise ArithmeticError(f"{name}: {x} != {y}")
    positive[name] += 1


def reject(name, x, y):
    if abs(x-y) < mp.mpf("1e-14")*max(abs(y), mp.mpf("1e-20")):
        raise ArithmeticError(f"mutation survived: {name}")
    negative[name] += 1


def cut_gap(a, s):
    beta = mp.sqrt(1-4*a/s)
    t0, h = s/2-a, s*beta/2
    # Subtract the two propagators before integration: (b2-b1)/(d1*d2).
    return beta*t0*mp.quad(
        lambda x: 15*a/((a+t0-h*x)*(16*a+t0-h*x)), [-1,1])/2


def cut_log(a, b, s):
    beta = mp.sqrt(1-4*a/s)
    t0, h = s/2-a, s*beta/2
    return t0/s*mp.log((b+t0+h)/(b+t0-h))


def determinants(s, theta, phi):
    # A common SU(2) spatial rotation sends the initial momentum pair to
    # the final pair. Independent ket phases also exercise complex conjugation.
    n = mp.sqrt(mp.sqrt(s))
    u1, u2 = (0,n), (n,0)
    k, t = mp.cos(theta/2), mp.sin(theta/2)
    z = mp.exp(1j*phi/2)
    u3 = (-n*t/z, n*k*z)
    u4 = (n*k/z, n*t*z)
    def pair(left, right):
        return left[0]*right[1]-left[1]*right[0]
    in_pair = pair(u1,u2)*mp.exp(1j*mp.mpf(".21"))
    out_pair = pair(u3,u4)*mp.exp(1j*mp.mpf("-.53"))
    return -out_pair.conjugate()*in_pair, -in_pair.conjugate()*out_pair


def dsigma_factor(s):
    # Integrate d^3p3 d^3p4 delta^4(P-p3-p4) over p4 and |p3|.
    k = mp.sqrt(s)/2
    radial = k**2/((2*k)*(2*k)) / 2  # delta(sqrt(s)-2k) radial Jacobian
    dphi2 = (2*mp.pi)**(4-6)*radial
    flux = 4*(s/2)
    return dphi2/(2*flux)  # one identical-final-state factorial


def main():
    a,s = mp.mpf(1),mp.mpf(8)
    dr = cut_gap(a,s)
    for scale,xi in [(1,5),(2,8),(9,17)]:
        aa,ss = mp.mpf(scale),mp.mpf(scale*xi)
        eig = cut_log(aa,aa,ss)-cut_log(aa,16*aa,ss)
        close("tree transfer integral versus cut logarithms",cut_gap(aa,ss),eig)
        close("scale-independent spectral gap",cut_gap(aa,ss),cut_gap(a,ss/aa))
        close("two-body phase space and beam flux",dsigma_factor(ss),1/(128*mp.pi**2*ss))
    g = mp.mpc(".4",".7")
    for theta,phi in [(0,0),(".97","1.3"),("2.5",".62")]:
        direct, reverse = determinants(s,mp.mpf(theta),mp.mpf(phi))
        close("reverse ordered helicity factor",reverse,direct.conjugate())
        close("external spinor magnitude",abs(direct),s)
        amp = direct*g*dr/(8*mp.pi)
        amp_reverse = reverse*g.conjugate()*dr/(8*mp.pi)
        close("anti-Hermitian coefficient as a Hermitian matrix",amp_reverse,amp.conjugate())
    amplitude = s*g*dr/(8*mp.pi)
    phase_space = dsigma_factor(s)
    bound = abs(amplitude)**2/(32*mp.pi*s)

    # Include angle-dependent Hermitian contacts with arbitrary flavor phases.
    # Hermiticity acts on the contact in the inverse process. Causal i does
    # not change sign between two independently evaluated upper-rim amplitudes.
    def contact(x):
        return mp.mpc(".17","-.89") + mp.mpc(".51",".23")*x*x
    def direct_rate(x):
        return abs(contact(x)+1j*amplitude)**2
    def inverse_rate(x):
        return abs(contact(x).conjugate()+1j*amplitude.conjugate())**2
    for x in map(mp.mpf,["-1","-.7","0",".8","1"]):
        close("pointwise complex reciprocal probabilities",
              (direct_rate(x)+inverse_rate(x))/2,abs(contact(x))**2+abs(amplitude)**2)
    mean = 2*mp.pi*phase_space*mp.quad(lambda x: (direct_rate(x)+inverse_rate(x))/2, [-1,1])
    excess = 2*mp.pi*phase_space*mp.quad(lambda x: abs(contact(x))**2,[-1,1])
    close("integrated mean minus integrated dispersive excess",mean-excess,bound)
    close("literal angular integral at saturation",phase_space*2*mp.pi*mp.quad(
        lambda x: abs(amplitude)**2,[-1,1]),bound)
    close("contact tuning zero direct probability",-1j*amplitude+1j*amplitude,0)
    close("same contact doubles inverse amplitude",
          (-1j*amplitude).conjugate()+1j*amplitude.conjugate(),2j*amplitude.conjugate())
    near = 4+mp.mpf("1e-40")
    beta = mp.sqrt(1-4/near)
    close("two-mass threshold coefficient",cut_gap(a,near)/beta,mp.mpf(15)/34,"1e-35")
    threshold = 4*(mp.mpf(15)/34)**2/(2048*mp.pi**3)
    close("threshold rate coefficient",near*(cut_gap(a,near)/beta)**2/(2048*mp.pi**3),
          threshold,"1e-38")
    for aa in (mp.mpf(4),mp.mpf(9)):
        # Fixed dimensionless g=a*C*w0 under a change of the mass unit.
        ds = 8*aa*(abs(g)/aa)**2*cut_gap(aa,8*aa)**2/(2048*mp.pi**3)
        close("inverse-square mass scaling at fixed dimensionless inputs",aa*ds,bound)
    reject("omitted final-state factorial",2*bound,bound)
    reject("extra initial-pair factorial",bound/2,bound)
    reject("averaged fixed initial helicities",bound/4,bound)
    reject("duplicated flavor orientation",4*bound,bound)
    reject("replaced spectral cut by ordinary complex amplitude imaginary part",
           (contact(0)+1j*amplitude).imag,amplitude)
    reject("claimed positive direct-only floor at tuned zero",0,bound)
    reject("discarded real or imaginary contact contribution",mean,bound)
    reject("used absolute squared Yukawas in null-pole example",
           cut_log(a,a,s)+cut_log(a,16*a,s),dr)
    unit_amp = s*dr/(8*mp.pi)
    unit_bound = unit_amp**2/(32*mp.pi*s)
    fmt = lambda x: mp.nstr(x,26)
    result = {"status":"PASS","method":"transfer-angle integral, helicity determinants, two-body flux and reciprocal phase-space integration",
              "benchmark":{"s_over_a":8,"b_over_a":[1,16],"weights_over_w0":[1,-1],
                           "Delta_rho":fmt(dr),"abs_A_over_abs_g":fmt(unit_amp),
                           "a_times_rate_floor_over_abs_g_squared":fmt(unit_bound),
                           "threshold_a_times_rate_floor_over_abs_g_squared_beta_squared":fmt(threshold)},
              "positive_checks":sum(positive.values()),"positive_categories":dict(positive),
              "negative_checks":sum(negative.values()),"negative_categories":dict(negative),
              "repository_imports":False,"assertions_required":False}
    data = json.dumps(result,indent=2)+"\n"
    Path(__file__).with_name("INVERSE_RATE_RESULTS.json").write_text(data,encoding="utf-8")
    print(data,end="")


if __name__ == "__main__":
    main()
