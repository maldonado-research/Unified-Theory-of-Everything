"""A charged-fermion phase-space test on the gas-loaded HDblast wall."""
import argparse
import gzip
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq


HERE = Path(__file__).resolve().parent
COLD_BUDGET = 0.1
DIRAC_SPECIES = 3


def spectrum(ell, shift):
    """Spatial Dirac eigenvalues in a single orbital multiplet; each repeated ell times."""
    magnetic_twice = np.arange(2-ell,ell,2)
    # Radicand in a positive form for either sign of the gauge amplitude.
    root = np.sqrt((ell-abs(shift))**2 + 2*abs(shift)*(ell-np.sign(shift)*magnetic_twice))
    return np.r_[ell+0.5-shift, ell+0.5+shift, 0.5+root, 0.5-root]


def direct_spectrum(ell, shift):
    """Independent spin and left-invariant orbital matrices in the uncoupled basis."""
    j = (ell-1)/2
    m = np.arange(ell)-j
    plus = np.diag(np.sqrt((j-m[:-1])*(j+m[:-1]+1)), -1)
    jx, jy = (plus+plus.T)/2,(plus-plus.T)/(2j)
    sx = np.array([[0,1],[1,0]])
    sy = np.array([[0,-1j],[1j,0]])
    sz = np.diag([1,-1])
    op = np.kron(sx,2*jx) + np.kron(sy,2*jy) + np.kron(sz,np.diag(2*m-shift))
    return np.linalg.eigvalsh(op+1.5*np.eye(2*ell))


def ceilings(radius, mass, field, response, delta=COLD_BUDGET):
    """Fill complete/partial levels in increasing kinetic energy; allow mixtures."""
    shift = mass*response*radius**2*field/2
    # Include every level with kinetic energy <= 0.6 times rest energy.
    cutoff = 0.6
    momentum = radius*mass*math.sqrt(cutoff*(2+cutoff))
    last_ell = math.ceil(abs(shift) + momentum + 1)
    levels, degeneracies = [], []
    for ell in range(1,last_ell+1):
        p_over_m = spectrum(ell,shift)/(radius*mass)
        kinetic = p_over_m**2/(np.hypot(1,p_over_m)+1)
        keep = kinetic <= cutoff
        levels.extend(kinetic[keep])
        degeneracies.extend([2*DIRAC_SPECIES*ell]*np.count_nonzero(keep))
    if not levels:
        return 0., 0., 0., shift
    idx = np.argsort(levels)
    energy = np.array(levels)[idx]
    capacity = np.array(degeneracies)[idx]
    balance = np.cumsum((energy-delta)*capacity)
    expensive = np.flatnonzero(balance>0)
    if not len(expensive):
        raise ArithmeticError("increase energy cutoff to enclose the first over-budget level")
    i = expensive[0]
    spare = -balance[i-1] if i else 0.
    filled = capacity[:i].sum() + spare/(energy[i]-delta)
    # Linear-programming dual: fill every level below this marginal energy.
    dual = np.maximum(energy[i]-energy,0.) @ capacity / (energy[i]-delta)
    if not math.isclose(filled,dual,rel_tol=2e-12,abs_tol=1e-11):
        raise ArithmeticError("primal/dual fermion count mismatch")
    return float(filled),float(energy[i]),float(dual),shift


def free_count(momentum):
    shells = max(0,math.floor(momentum - 1.5)+1)
    return 4*shells*(shells+1)*(shells+2)


def perturbative_demand(inputs, K):
    # Divide out mu from the minimum comoving count resulting from local alpha < 1.
    return 3*math.pi*inputs["initial_radius"]**3*inputs["initial_cold_rest_density"]*inputs["Q"]**2/(2*K)


def mass_bound(radius, field, K, inputs):
    """Conservative continuous-mass obstruction from free singular-value counting."""
    c = radius*math.sqrt(2*COLD_BUDGET*(2+2*COLD_BUDGET)) + inputs["Q"]*radius**2*abs(field)/2
    demand = perturbative_demand(inputs,K)
    # If c*mu<1.5 no low-energy states exist. Above it this capacity/mass grows.
    def difference(mu):
        return 8*(c*mu+0.5)**3/mu-demand
    min_mass=1.5/c
    if difference(min_mass)>=0:
        raise ArithmeticError("benchmark parameter outside quoted root bracket")
    return brentq(difference,min_mass,1e9)


def sqrt_scaled_floor(x, scale):
    """Floor of scale*sqrt(x) for a nonnegative Fraction."""
    return math.isqrt(x.numerator*scale**2//x.denominator)


def lattice_count(shift, momentum, scale):
    """Integer spectral count at rational conformal |s| and |lambda| cutoff.

    Arguments shift and momentum are already scaled integers. Multiplications
    on orbital labels use int64 only after bounding the largest intermediate.
    """
    if momentum<0:
        return 0
    if scale%2 or shift<=0:
        raise ValueError("positive shift and even grid denominator required")
    half = scale//2
    hi = (shift+momentum+half)//scale+2
    lo = max(1,(shift-momentum-half)//scale-1)
    if (hi*scale+shift+half)**2 >= 2**62 or 12*hi**3 >= 2**62:
        raise OverflowError("revise integer grid for this parameter range")
    ell = np.arange(lo,hi+1,dtype=np.int64)
    edges = ((abs(scale*ell+half-shift)<=momentum).astype(np.int64)
             + (abs(scale*ell+half+shift)<=momentum).astype(np.int64))
    branches=np.zeros_like(ell)
    base=(ell*scale-shift)**2
    denominator=4*shift*scale
    for lower,upper in ((0,momentum-half),(max(0,half-momentum),half+momentum)):
        if upper<lower:
            continue
        jmin = np.maximum(1,-((base-lower**2)//denominator))
        jmax = np.minimum(ell-1,(upper**2-base)//denominator)
        branches += np.maximum(0,jmax-jmin+1)
    return int((2*DIRAC_SPECIES*ell*(edges+branches)).sum())


def mass_window(inputs):
    """Whole-mass lower certificate and an explicit admissible state at mu=2000."""
    F=Fraction
    pt=inputs["switch"]
    a=F(str(pt["a"]))
    d=F(str(inputs["Q"]))*a*a*abs(F(str(pt["B"])))/2
    a_hi, d_lo, d_hi = F(9554,1000), F(4396,100), F(4397,100)
    if not (a<a_hi and d_lo<d<d_hi):
        raise ArithmeticError("switch parameter bounds")
    # Bound A=N_req/(mu*S0/S0min) with rational brackets for pi.
    pref=3*inputs["initial_radius"]**3*F(str(inputs["initial_cold_rest_density"]))*F(str(inputs["Q"]))**2/(2*F(str(pt["K"])))
    A_lo, A_hi=2855000000,2856000000
    if not (A_lo<F(333,106)*pref<F(355,113)*pref<A_hi):
        raise ArithmeticError("local charge ceiling")
    delta=F(1,10)
    # First eliminate all mu <= 52 by the singular-value displacement bound.
    c=a_hi*F(664,1000)+d_hi  # sqrt(.44) < .664
    if not (8*(c*52+F(1,2))**3 < A_lo*52):
        raise ArithmeticError("starting mass exclusion")
    z=F(163,1000)
    bins=4096
    for mu_floor,exclude_through,inv_alpha in ((52,1900,1),(1900,1970,1),(1970,6240,10)):
        cubic=quadratic=F(0)
        for i in range(1,bins+1):
            y=z*i/bins
            c=a_hi*F(sqrt_scaled_floor(y*(2+y),10**6)+1,10**6)
            u=c+F(1,2*mu_floor)
            if u>=d_lo:
                raise ArithmeticError("bulk support crosses other endpoint")
            cubic+=4*u**3
            quadratic+=6*(d_hi+c)*(2*c+F(1,mu_floor)) + 3*(1+u/d_lo)*u**2
        # Direct rational comparison at the upper endpoint; the RHS increases in mu.
        capacity=(cubic*exclude_through**2+quadratic*exclude_through)*z/bins
        if not capacity < inv_alpha*A_lo*(z-delta):
            raise ArithmeticError("continuous mass exclusion")
    if not (cubic*2000**2+quadratic*2000)*z/bins < F(103,100)*A_lo*(z-delta):
        raise ArithmeticError("coupling bound for the mass-2000 trial")
    # Construct an explicit trial: fill all states below z, uniformly thin
    # if necessary, and enclose their mean energy by level counting.
    def trial(mu,alpha_inv):
        scale=1000
        stop=F(13,80)  # individual kinetic energy threshold .1625
        s_grid=math.floor(scale*d*mu)
        def bracket(y):
            w=sqrt_scaled_floor(a*a*mu**2*y*(2+y),scale)
            return lattice_count(s_grid,w-1,scale),lattice_count(s_grid,w+2,scale)
        low,high=bracket(stop)
        nbin=384
        integ_low=stop*sum(bracket(stop*i/nbin)[0] for i in range(nbin))/nbin
        mean_upper=stop-integ_low/high
        if not (low>alpha_inv*A_hi*mu and mean_upper<delta):
            raise ArithmeticError("feasible ensemble certificate")
        return dict(mass_L=mu,S0_over_local_minimum=float(alpha_inv),
                    individual_kinetic_cut=float(stop),count_lower=low,count_upper=high,
                    mean_kinetic_over_rest_upper=math.ceil(mean_upper*10**7)/10**7)
    permissive=trial(2000,F(101,100))
    weak=trial(6500,F(101,10))
    # Cross-check integer counting against the explicit spectrum at small inputs.
    scale=1000
    for ell_shift in (1250,43965,175858):
        mom=15001
        direct=sum(6*ell*np.count_nonzero(abs(spectrum(ell,ell_shift/scale))<=mom/scale)
                   for ell in range(1,math.ceil((ell_shift+mom)/scale)+3))
        if lattice_count(ell_shift,mom,scale)!=direct:
            raise ArithmeticError("rank count cross-check")
    # Only the exceptional edge |s|=ell+1/2 can have zero kinetic energy.
    # Bulk R²=(ell-s)²+4sj for j>=1 has minimum >=4 for ell>=2.
    if not 6*d/A_lo < F(925,10**10):
        raise ArithmeticError("zero-kinetic state capacity")
    return dict(strictly_zero_kinetic_fraction_upper=9.25e-8,
                excluded_mass_L_at_or_below=1970,admissible_trial=permissive,
                at_mass_2000_inverse_alpha_strictly_below=1.03,
                local_alpha_at_most_point_one=dict(excluded_mass_L_at_or_below=6240,
                                                    admissible_trial=weak),
                lower_certificate_quadrature_bins=bins,trial_certificate_quadrature_bins=384)


def checked_inputs(source_trace, source_gas):
    p = json.loads((HERE/"HD_INPUTS.json").read_text())
    geom = gzip.decompress((HERE/"trace_geometry.csv.gz").read_bytes())
    if hashlib.sha256(geom).hexdigest()!=p["geometry_csv_uncompressed_sha256"]:
        raise ArithmeticError("geometry digest")
    if source_trace:
        original = source_trace.read_bytes()
        trimmed = b"\n".join(b",".join(line.split(b",",3)[:3]) for line in gzip.decompress(original).splitlines())+b"\n"
        if hashlib.sha256(original).hexdigest()!=p["background_source_gzip_sha256"] or trimmed!=geom:
            raise ArithmeticError("geometry differs from original archive")
    if source_gas:
        original = source_gas.read_bytes()
        j=json.loads(original)
        case=p["source_case"]
        run=j["runs"][case]
        params=j["parameters"][case]
        start=dict(zip(j["trace_columns"],run["trace"][0]))
        if (start["t"] != 0 or start["a"] != p["initial_radius"] or start["K"] != p["initial_K"]
                or start["rest"] != p["initial_cold_rest_density"] or params["D"] != start["rest"]
                or params["Q"] != p["Q"] or params["activation"] != p["switch"]["t"]):
            raise ArithmeticError("initial gas normalization or switch time")
        point=next(dict(zip(j["trace_columns"],r)) for r in run["trace"] if r[0]==p["switch"]["t"])
        if hashlib.sha256(original).hexdigest()!=p["charged_json_sha256"] or any(point[k]!=v for k,v in p["switch"].items()):
            raise ArithmeticError("gas source")
        if point["u"] != 0 or point["theta"] != 0:
            raise ArithmeticError("mobility switch preparation")
        if any(run[k]!=v for section in ("magnetic_zero","heating_interval") for k,v in p[section].items()):
            raise ArithmeticError("gas event extracts")
    return p


def verify(args):
    p = checked_inputs(args.source_trace,args.source_gas)
    point=p["switch"]
    for ell in (1,2,3,11,43,45):
        for shift in (0., 1.25, -43.7):
            if not np.allclose(np.sort(spectrum(ell,shift)),direct_spectrum(ell,shift),rtol=1e-12,atol=2e-13):
                raise ArithmeticError("Dirac spectrum disagrees with orbital-spin matrix")
    for k in (1.49, 1.5, 3., 15.5):
        count=sum(2*DIRAC_SPECIES*ell*np.count_nonzero(abs(spectrum(ell,0))<=k)
                  for ell in range(1,math.ceil(k+3)))
        if count!=free_count(k):
            raise ArithmeticError("round three-sphere state count")
    h=p["heating_interval"]
    r=h["maximum_radius"]/h["first_ninety_intrinsic_radius"]
    bound=(math.sqrt(r)-1)**2/(2*(1+r))
    sites=[("initial_if_already_charged",p["initial_radius"],0.,p["initial_K"]),
           ("reference_magnetic_zero",p["magnetic_zero"]["postexit_magnetic_zero_radius"],0.,point["K"]),
           ("actual_mobility_switch",point["a"],point["B"],point["K"])]
    results={}
    for name,a,B,K in sites:
        need=perturbative_demand(p,K)
        rows={}
        for mu in (1.,4.):
            n,edge,dual,shift=ceilings(a,mu,B,p["Q"])
            markov=2*free_count(a*mu*math.sqrt(.44)+abs(shift))
            if not (0<=n<=markov+1e-8):
                raise ArithmeticError("exclusion-principle cross-check")
            rows[str(int(mu))]=dict(gauge_shift=round(shift,8),maximum_cold_count=round(n,6),
                marginal_kinetic_per_rest=round(edge,9),required_count_strictly_above=math.floor(need*mu),
                cold_fraction_upper=float(f"{n/(need*mu):.9g}"),conservative_count_ceiling=markov)
        results[name]=dict(three_Dirac_species=rows,necessary_mass_L_above=float(f"{mass_bound(a,B,K,p):.10g}"))
    out=dict(assumptions="three free Dirac species (particles and antiparticles), common mass and response, instantaneous positive energies, mean kinetic/rest <= 0.1; local alpha < 1 for demand bound",
             neutral_single_mode_late_segment_probability_bound=float(f"{bound:.12g}"),sites=results, switch_mass_window=mass_window(p))
    return out


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-trace",type=Path)
    parser.add_argument("--source-gas",type=Path)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    result=json.dumps(verify(args),indent=2,allow_nan=False)+"\n"
    if args.output:
        args.output.write_text(result,encoding="utf-8")
    else:
        print(result,end="")
