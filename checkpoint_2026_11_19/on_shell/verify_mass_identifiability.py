"""Check what the two-scalar cut resolves for finite sterile mass spectra."""
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


def check(ok, name):
    if not ok:
        raise ArithmeticError(name)


def angular_average(s, a, weights, cosine):
    """Direct weak routing by angular quadrature, with spin/phase factors stripped."""
    beta = mp.sqrt(1 - 4*a/s)
    n = s/2-a
    h = s*beta/2
    # After the azimuth integral, the projection of the outgoing direction
    # onto the cut momentum is cosine*x. Scalar masses are identical.
    def integrand(x):
        return sum(w*(n - h*cosine*x)/(n+b-h*x) for b, w in weights)
    return mp.quad(integrand, [-1, 0, 1])/2


def main():
    mp.mp.dps = 55
    y, r, x = sp.symbols("y r x", positive=True)
    ratio = (r+y)/(r+1/y)
    check(sp.cancel(ratio*ratio.subs(r,1/r)-y**2)==0, 'reciprocal log pair')
    log_derivative = 1/y + 1/(y+r)-1/(y+1/r)
    check(sp.cancel(sp.diff(ratio,y)/ratio-log_derivative)==0, 'logarithmic derivative')
    # The simple poles of a weighted sum carry the entries of this linear map.
    spectrum = [sp.Rational(1,4), sp.Integer(1), sp.Integer(4), sp.Integer(9)]
    residual = []
    for p in [sp.Integer(0)] + sorted(set(spectrum+[1/b for b in spectrum])):
        residual.append([sp.limit((y+p)*log_derivative.subs(r,b), y, -p) for b in spectrum])
    even_map = sp.Matrix(residual)
    weights = sp.Matrix([1,-2,1,0])
    check(even_map.rank()==3 and even_map*weights==sp.zeros(even_map.rows,1),
          'null space of equal-component channel')
    heavy_columns = [i for i,rval in enumerate(spectrum) if rval>=1]
    check(even_map[:,heavy_columns].rank()==len(heavy_columns), 'single-sided positive mass spectrum')
    weighted_map = even_map*sp.diag(*spectrum)
    check(even_map.col_join(weighted_map).rank()==len(spectrum), 'angular resolution')
    check(weighted_map*weights!=sp.zeros(even_map.rows,1), 'odd projection sees blind vector')
    # The proposed spectrum is realizable by complex Yukawas in one incoming
    # flavor, with the derivative seed confined to a different outgoing flavor.
    masses = sp.diag(sp.Rational(1,2),1,2)
    row = sp.Matrix([[sp.sqrt(2), sp.I*sp.sqrt(2), 1/sp.sqrt(2)]])
    actual = sp.Matrix([sp.expand(masses[i,i]*row[i]**2) for i in range(3)])
    check(actual==weights[:3,:], 'coherent anomalous mass weights')
    check((row*masses*row.H)[0]==4, 'modulus-square mutation erases cancellation')
    # The coefficient of beta^3*cos(theta) after including phase space.
    coeff = -sp.Rational(4,3)*sum(weights[i]/(1+spectrum[i])**2 for i in range(3))
    check(coeff==-sp.Rational(6,25), 'odd threshold term')
    family = [(mp.mpf(1)/4,1),(mp.mpf(1),-2),(mp.mpf(4),1)]
    sample_s = mp.mpf(9)
    cosine = mp.mpf('0.5')
    beta = mp.sqrt(1-4/sample_s)
    fwd = angular_average(sample_s, 1, family, cosine)
    back = angular_average(sample_s, 1, family, -cosine)
    check(abs(fwd+back)<mp.mpf('1e-48') and abs(fwd)>mp.mpf('1e-3'),
          'component assignment and opposite angles')
    for s in map(mp.mpf, ('4.0001', '5', '21', '175')):
        d = angular_average(s, 1, family, cosine)
        e = angular_average(s, 1, family, -cosine)
        check(abs(d+e)<mp.mpf('1e-46') and abs(d-e)>mp.mpf('1e-8'),
              'direct solid-angle quadrature away from threshold')
    small_beta = mp.mpf('1e-5')
    s = 4/(1-small_beta**2)
    reduced_odd = angular_average(s, 1, family, cosine)
    limit = reduced_odd/(small_beta**2*cosine)
    check(abs(limit+mp.mpf(6)/25)<mp.mpf('1e-9'), 'threshold limit from angular integral')
    two_mass = [(mp.mpf(1),1),(mp.mpf(16),-1)]
    check(abs(angular_average(s,1,two_mass,0)-mp.mpf(15)/34)<mp.mpf('1e-9'),
          'pole cancellation survives as nonzero two-mass cut')
    result = {
        'scope':'isolated flavor channel with seed on outgoing pair; common scalar squared mass a=1; ordered amplitude discontinuity',
        'reciprocal_example': {'sterile_masses':['1/2','1','2'],'pole_weights':[1,-2,1],
                               'sum_same_weak_components':'0 for all s>4',
                               's':9,'cos_theta':'1/2',
                               'beta_times_direct_weak_kernel':mp.nstr(beta*fwd,20),
                               'beta_times_exchanged_weak_kernel':mp.nstr(beta*back,20),
                               'direct_kernel_including_beta_over_cos_theta_beta_cubed_at_threshold':'-6/25'},
        'symbolic_residue_map_example':{'squared_masses':['1/4','1','4','9'],
                                       'rank_even':even_map.rank(),'rank_even_and_odd':even_map.col_join(weighted_map).rank()},
    }
    Path(__file__).with_name('MASS_IDENTIFIABILITY_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: reciprocal identity, residues, angular resolution, coherent weights and threshold limit')


if __name__=='__main__':
    main()
