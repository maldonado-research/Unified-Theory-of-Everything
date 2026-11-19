"""Profile a single incident flavor over four energies and arbitrary complex nuisances."""
import json
from pathlib import Path
from runpy import run_path

from mpmath import mp


def upper_gap(xi):
    """Dimensionless difference of two triangles from the physical spectral density."""
    beta = mp.sqrt(1-4/xi)
    k = mp.mpf(15)/17
    value = 2*mp.atanh(k*beta)

    def regular(u):
        # atanh subtraction identity removes the apparent pole at u=beta.
        if u == beta:
            return 2*k/(1-(k*beta)**2)
        argument = k*(u-beta)/(1-k*k*u*beta)
        return 4*u*mp.atanh(argument)/((u-beta)*(u+beta))

    principal = mp.quad(regular, [0, beta, 1])
    return (principal+value*(mp.log((1-beta**2)/beta**2)+mp.j*mp.pi))/xi


def close(name, actual, expected, tolerance='1e-42'):
    if abs(actual-expected) > mp.mpf(tolerance)*max(1, abs(expected)):
        raise ArithmeticError(f'{name}: {actual} versus {expected}')


def fmt(value):
    return mp.nstr(value, 25)


def complex_parts(z):
    return [fmt(z.real), fmt(z.imag)]


def calculate():
    mp.dps = 64
    points = list(map(mp.mpf, (5, 6, 9, 11)))
    weights = [mp.sqrt(x/4) for x in points]
    shapes = [upper_gap(x) for x in points]
    triangle = run_path(str(Path(__file__).with_name('verify_lower_operator.py')))['triangle']
    deviations = []
    for x, f in zip(points, shapes):
        direct = triangle(x, 1, 1)-triangle(x, 1, 16)
        close('dispersion against two simplex integrals', f, direct)
        deviations.append(abs(f-direct))

    # Column order f, 1, xi; the physical nuisance coefficients are -g5, d0, d1.
    design = mp.matrix([[w*f, w, w*x] for x, f, w in zip(points, shapes, weights)])
    target = mp.matrix([w*(x/2-1)*f for x, f, w in zip(points, shapes, weights)])
    enlarged = mp.matrix([list(design[i,:])+[target[i]] for i in range(4)])
    signed_minors = []
    for i in range(4):
        minor = mp.matrix([list(design[j,:]) for j in range(4) if j != i])
        signed_minors.append((-1)**(i+3)*mp.det(minor))
    norm = mp.norm(mp.matrix(signed_minors))
    q = mp.matrix([m.conjugate()/norm for m in signed_minors])
    close('left null vector', mp.norm(q.H*design), 0)
    close('null normalization', mp.norm(q), 1)
    determinant = mp.det(enlarged)
    close('oriented cofactor expansion', determinant, norm*(q.H*target)[0])

    solution, qr_residual = mp.qr_solve(design, target)
    residual = target-design*solution
    profile = abs((q.H*target)[0])**2
    close('QR profile', qr_residual**2, profile)
    close('normal equations', mp.norm(design.H*residual), 0)
    close('residual direction', mp.norm(residual-q*(q.H*target)[0]), 0)
    gram_nuisance = mp.det(design.H*design)
    close('Cauchy-Binet', gram_nuisance, norm**2)
    gram_all = mp.det(enlarged.H*enlarged)
    close('Gram ratio', gram_all/gram_nuisance, profile)
    scale = 2048*mp.pi**5
    close('single-direction phase space', (mp.mpf(2)/(16*mp.pi**2))**2/(32*mp.pi), 1/scale)
    # The unprofiled rate includes a genuinely complex offset in every column.
    g7 = mp.mpc('0.7','-0.4')
    offsets = mp.matrix([mp.mpc('0.5','1.1'), mp.mpc('-0.3','0.7'), mp.mpc('0.2','-0.8')])
    test_residual = g7*target+design*offsets
    close('complex Pythagorean identity', mp.norm(test_residual)**2,
          abs(g7)**2*profile+mp.norm(design*(offsets+g7*solution))**2)
    sig_x = mp.svd(design, compute_uv=False)
    sig_a = mp.svd(enlarged, compute_uv=False)
    if not all(abs(m) > mp.mpf('0.05') for m in signed_minors) or profile <= 0:
        raise ArithmeticError('degenerate four-energy scan')
    return {
        'status': 'PASS',
        'xi': [int(x) for x in points],
        'exposure_per_energy': '1/4',
        'sterile_masses_over_sqrt_a': [1, 4], 'sterile_weights_over_w0': [1,-1],
        'basis': ['sqrt(xi/4)*f', 'sqrt(xi/4)', 'sqrt(xi/4)*xi'],
        'target': 'sqrt(xi/4)*(xi/2-1)*f',
        'complex_numbers': '[real, imaginary]',
        'f': [complex_parts(z) for z in shapes],
        'max_simplex_dispersion_mismatch': mp.nstr(max(deviations), 5),
        'signed_nuisance_minors': [complex_parts(z) for z in signed_minors],
        'det_enlarged': complex_parts(determinant),
        'det_gram_nuisance': fmt(gram_nuisance.real), 'det_gram_enlarged': fmt(gram_all.real),
        'singular_values_design': [fmt(s) for s in sig_x],
        'singular_values_enlarged': [fmt(s) for s in sig_a],
        'condition_design_2norm': fmt(sig_x[0]/sig_x[2]),
        'condition_enlarged_2norm': fmt(sig_a[0]/sig_a[3]),
        'q_left_null': [complex_parts(z) for z in q],
        'minimizing_coefficients_g5_d0_d1_per_g7': [
            complex_parts(solution[0]), complex_parts(-solution[1]), complex_parts(-solution[2])],
        'weighted_squared_residual': fmt(profile),
        'a_times_mean_direct_cross_section_floor_per_abs_g7_squared': fmt(profile/scale),
    }


if __name__ == '__main__':
    data = json.dumps(calculate(), indent=2)+'\n'
    Path(__file__).with_name('SINGLE_DIRECTION_RESULTS.json').write_text(data, encoding='utf-8')
    print(data, end='')
