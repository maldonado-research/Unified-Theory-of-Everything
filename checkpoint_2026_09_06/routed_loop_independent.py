"""Independent exact Feynman-parameter and routing check.

This script tests the algebra in routed_loop_independent.md. It does not
claim a complete EFT basis or EOM reduction. Tensor integrals use the
original Zhang loop denominator k^2[(k+p)^2-m_H^2].
"""
from fractions import Fraction as F
from pathlib import Path
import json


def require(ok, message):
    if not ok:
        raise ArithmeticError(message)


def integrate(poly):
    return sum((c/F(n+1) for n,c in poly.items()), F(0))


i1_p = integrate({1: -F(1)})
i2_pp = integrate({2: F(1)})
i2_g_m2 = integrate({1: F(1,2)})
i2_g_p2 = integrate({1: -F(1,2), 2: F(1,2)})
i3_ppp = integrate({3: -F(1)})
i3_gp_m2 = integrate({2: -F(1,2)})
i3_gp_p2 = integrate({2: F(1,2), 3: -F(1,2)})
require(4*i2_g_m2 == 1, 'I2 trace must reproduce massive tadpole')
require(4*i2_g_p2 + i2_pp == 0, 'I2 p^2 trace must cancel')
require(6*i3_gp_m2 == -1, 'I3 trace has -mH^2 p sign')
require(6*i3_gp_p2 + i3_ppp == 0, 'I3 p^3 trace must cancel')

eps=((0,1),(-1,0))
for a in range(2):
    for c in range(2):
        crossed=sum(eps[a][d]*eps[b][c]*eps[b][d] for b in range(2) for d in range(2))
        direct=eps[a][c]*sum(eps[b][d]**2 for b in range(2) for d in range(2))
        require(crossed == eps[a][c], 'Crossed SU2 weight is one')
        require(direct == 2*eps[a][c], 'Direct SU2 weight is two')

# Original k convention: 2 p^2 I_mu + trace(I3)_mu + 2 h^nu I2_mu,nu + h^2 I1_mu.
original = [2*i1_p, 2*i2_pp, i1_p, 2*i2_g_p2, 6*i3_gp_m2, 2*i2_g_m2]
expected=[F(-1),F(2,3),F(-1,2),F(-1,6),F(-1),F(1,2)]
require(original == expected, 'Independent weighted tensor assembly')
root_path=Path(__file__).parent/'routed_loop/ROUTED_POLE_CHECKS.json'
root=json.loads(root_path.read_text())
require(list(map(F,root['coefficients'])) == [-x for x in original], 'ell=-k must reverse the complete numerator sign')

out={
    'status':'PASS_INDEPENDENT_EXACT_FEYNMAN_PARAMETER_ROUTING_AND_MASS_CHECK',
    'original_denominator':'k^2[(k+p)^2-mH^2]',
    'removed_factor':'i/(16*pi^2*epsilon), d=4-2epsilon',
    'basis':root['basis'],
    'original_Zhang_k_coefficients':list(map(str,original)),
    'root_ell_equals_minus_k_coefficients':list(map(str,[-x for x in original])),
    'tensor_coefficients':{
        'I1_p':str(i1_p),'I2_pp':str(i2_pp),'I2_g_mH2':str(i2_g_m2),
        'I2_g_p2':str(i2_g_p2),'I3_ppp':str(i3_ppp),
        'I3_symgp_mH2':str(i3_gp_m2),'I3_symgp_p2':str(i3_gp_p2)
    },
    'checks':['exact parameter moments','two tensor traces','SU2 contractions','weight-two direct plus weight-one crossed routing','signed comparison with independently produced parent result'],
    'limitation':'Declared one-seed cubic UV polynomial; not full Green-basis matching, sourced EOM, a physical RGE coefficient, or a finite threshold.'
}
Path(__file__).with_name('ROUTED_LOOP_INDEPENDENT_CHECK.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
