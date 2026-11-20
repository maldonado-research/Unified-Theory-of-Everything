"""Check the positive-action C5 and C7 triangles in the charged pair channel."""
import json
from pathlib import Path

from mpmath import mp


def close(label, actual, expected):
    if abs(actual-expected) > mp.mpf('1e-34') * max(1, abs(expected)):
        raise ArithmeticError(f'{label}: {actual} != {expected}')


def triangle(s, a, b):
    """Upper lip, with y integrated exactly on the Feynman simplex."""
    knots = [mp.mpf(0), mp.mpf(1)]
    if s >= 4*a:
        beta = mp.sqrt(1-4*a/s)
        # 2a/[s(1+beta)] avoids cancellation at high energy.
        lower, upper = 2*a/(s*(1+beta)), (1+beta)/2
        knots = sorted(set([*knots, lower, upper]))

    def at_x(x):
        start = b+(a-b)*x
        slope = b-a+s*x
        end = s*(x-lower)*(x-upper) if s >= 4*a else a-s*x*(1-x)
        if not end:
            return 0  # Boundary values at integrable logarithmic singularities.
        if not slope:
            return (1-x)/start
        if abs(slope*(1-x)) < start/2:
            return -mp.log1p(-slope*(1-x)/start)/slope
        return (mp.log(start/abs(end)) + (mp.j*mp.pi if end < 0 else 0))/slope

    return mp.quad(at_x, knots)


def transfer_log(t, a, b):
    """s times the spectral density of T at s=4a/t, t in (0,1]."""
    beta = mp.sqrt(1-t)
    h_plus_v = a*(1+beta)**2/t
    h_minus_v = a*a/h_plus_v
    return mp.log1p((4*a*beta/t)/(b+h_minus_v))


def dispersion(s, a, b):
    """Unsubtracted dispersion integral from the angular logarithm."""
    edge = 4*a
    if s < edge:
        return mp.quad(lambda t: transfer_log(t,a,b)/(edge-s*t), [0, 1])
    if s == edge:
        # t=1-u²; retain the beta=u coordinate when u is extremely small.
        def endpoint(u):
            if u == 0:
                return 2/(a+b)
            if u == 1:
                return 0
            t = (1-u)*(1+u)
            h_minus_v = a*(1-u)/(1+u)
            return 2*mp.log1p(4*a*u/(t*(b+h_minus_v)))/(edge*u)
        return mp.quad(endpoint, [0,1])
    pivot = edge/s
    at_pole = transfer_log(pivot,a,b)

    def regular(t):
        if abs(t-pivot) < mp.mpf('1e-26'):
            return -mp.diff(lambda z: transfer_log(z,a,b), pivot)
        return (transfer_log(t,a,b)-at_pole)/(pivot-t)

    principal = mp.quad(regular, [0, pivot, 1])+at_pole*mp.log(pivot/(1-pivot))
    return (principal+mp.j*mp.pi*at_pole)/s


def transfer_gap(s,a):
    """Analytic difference of logs with positive terms throughout the cut."""
    h = s/2-a
    v = mp.sqrt(s*(s-4*a))/2
    return mp.log1p(2*v*15*a/((a+a*a/(h+v))*(16*a+h+v)))/s


def sew(s,a,weights,c5,c7,spin):
    """Compute A from the two ordered trees, with the full scalar phase space."""
    beta = mp.sqrt(1-4*a/s)
    v = s*beta/2
    t0 = a-s/2
    # B12 and B34 are taken so W=-B34* B12=spin. Their arbitrary split
    # cancels from the product, including when an external pair is rephased.
    b12 = mp.sqrt(s)
    out_bar = -spin/b12
    def product(x):
        # Choose an outgoing direction perpendicular to the incoming beam.
        # cos(phi)=1 suffices: the two contact routings cancel its azimuth.
        y = mp.sqrt(1-x*x)
        t3, u3 = t0+v*y, t0-v*y
        electron = -out_bar*((c5+c7*t3)+(c5+c7*u3))
        muon = b12*sum(w*(1/(t0+v*x-b)+1/(t0-v*x-b)) for b,w in weights)
        return electron*muon
    # 1/2 from unitarity, 1/2! from identical scalars; dPhi2=beta dOmega/(32pi²).
    return beta*(2*mp.pi)/(4*32*mp.pi**2)*mp.quad(product, [-1, 1])


def fmt(x):
    return mp.nstr(x, 25)


def complex_record(z):
    return {'real': fmt(mp.re(z)), 'imag': fmt(mp.im(z))}


def main():
    mp.dps = 54
    a = mp.mpf(1)
    rows = []
    for b_int,s_int in [(1,-3),(16,0),(16,3),(1,4),(16,4),(1,8),(16,8),(16,19)]:
        b,s = mp.mpf(b_int),mp.mpf(s_int)
        direct = triangle(s,a,b)
        close('unsubtracted dispersion of T',direct,dispersion(s,a,b))
        if b == a and s >= 4*a:
            angle = mp.pi/2+mp.j*mp.acosh(mp.sqrt(s/(4*a)))
            close('equal masses analytic continuation', direct, 2*angle**2/s)
        rho = transfer_log(4*a/s,a,b)/s if s > 4*a else 0
        close('spectral density from simplex', mp.im(direct)/mp.pi, rho)
        rows.append({'a':1,'b':b_int,'s':s_int,'K5':complex_record(-direct),'rho5':fmt(-rho)})

    s = 8*a
    h = s/2-a
    dt = triangle(s,a,a)-triangle(s,a,16*a)
    db = mp.quad(lambda x: mp.log(x+16*(1-x)),[0,1])
    close('finite bubble difference',db,mp.mpf(16)*mp.log(16)/15-1)
    gap = transfer_gap(s,a)
    close('benchmark cut log difference', mp.im(dt)/mp.pi, gap)
    c5,c7,w0 = mp.mpc('.8','-.4'),mp.mpc('.1','.3'),mp.mpc('1.4','.25')
    spin = s*mp.exp(mp.j*mp.mpf('.23'))
    weights = [(a,w0),(16*a,-w0)]
    def upper(g5,g7,w,W):
        return W*w*(g7*db+(h*g7-g5)*dt)/(8*mp.pi**2)
    forward = upper(c5,c7,w0,spin)
    inverse = upper(c5.conjugate(),c7.conjugate(),w0.conjugate(),spin.conjugate())
    absorptive = (forward-inverse.conjugate())/(2*mp.j)
    close('two full upper rims',absorptive,spin*w0*(h*c7-c5)*gap/(8*mp.pi))
    close('two tree unitarity sewing',absorptive,sew(s,a,weights,c5,c7,spin))
    close('pure constant-contact sewing',sew(s,a,weights,1,0,spin),-spin*w0*gap/(8*mp.pi))
    close('pure derivative sewing',sew(s,a,weights,0,1,spin),spin*w0*h*gap/(8*mp.pi))
    close('crossing cancels scalar tree',sew(s,a,weights,h*c7,c7,spin),0)
    close('crossing leaves local bubble',upper(h*c7,c7,w0,spin),spin*w0*c7*db/(8*mp.pi**2))
    if abs(forward.imag-absorptive) < mp.mpf('.001'):
        raise ArithmeticError('complex-phase reciprocal test lost sensitivity')

    for scale in [mp.mpf('.25'),mp.mpf(9)]:
        close('mass-dimension minus two', scale*triangle(8*scale,scale,16*scale),triangle(s,a,16*a))
    for xi in ['4.00000000000000000000001','8','37','1e30']:
        energy = mp.mpf(xi)*a
        gap_x = transfer_gap(energy,a)
        close('positive stable difference',gap_x,
              (transfer_log(4*a/energy,a,a)-transfer_log(4*a/energy,a,16*a))/energy)
        if gap_x <= 0:
            raise ArithmeticError('sterile transfer gap must be strictly positive')
    receipt = {
        'status':'PASS',
        'normalization':'M = W/(8 pi^2) sum_i w_i (C7 B_i + (h C7 - C5) T_i), h=s/2-a',
        'methods':['upper-lip Feynman simplex','unsubtracted dispersion from transfer logarithm',
                   'physical product of trees including identical-intermediate-state factor',
                   'inverse upper-rim conjugation and isolated zero of the tree'],
        'individual_lower_triangles': rows,
        'two_mass_benchmark':{'a':1,'s':8,'b':[1,16],'relative_weights':[1,-1],
                              'Delta_K5':complex_record(-dt),'Delta_B':fmt(db),
                              'Delta_K7_unsubtracted':complex_record(db+h*dt),
                              'Delta_P':fmt(gap),'rho7_gap':fmt(h*gap),
                              'tuned_C5_over_C7':3},
    }
    path = Path(__file__).with_name('LOWER_OPERATOR_RESULTS.json')
    path.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print('PASS positive C5/C7 normalization, dispersion, cut, reciprocal upper rims, cancellation')


if __name__ == '__main__':
    main()
