"""Rational enclosures for a four-setting, single-direction cross-section bound."""
from dataclasses import dataclass
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path
import json


@dataclass(frozen=True)
class Interval:
    lo: Q
    hi: Q

    def __add__(self, other):
        other = enclosure(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-enclosure(other))

    def __mul__(self, other):
        other = enclosure(other)
        products = [a*b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return Interval(min(products), max(products))

    __rmul__ = __mul__

    def square(self):
        top = max(self.lo**2, self.hi**2)
        low = 0 if self.lo <= 0 <= self.hi else min(self.lo**2, self.hi**2)
        return Interval(Q(low), top)


def enclosure(value):
    if isinstance(value, Interval):
        return value
    return Interval(Q(value), Q(value))


def outward(value, digits=26):
    scale = 10**digits
    return Interval(Q(value.lo*scale//1,scale), Q(-((-value.hi*scale)//1),scale))


def sqrt_bounds(value, digits=28):
    scale = 10**digits
    low = isqrt(value*scale**2//1)
    return Interval(Q(low,scale), Q(low+1,scale))


def atan_unit_fraction(q, n=32):
    partial = sum(Q((-1)**j, (2*j+1)*q**(2*j+1)) for j in range(n))
    next_sum = partial + Q((-1)**n,(2*n+1)*q**(2*n+1))
    return Interval(min(partial,next_sum),max(partial,next_sum))


def atanh_over_argument(z, n):
    """sum z^r/(2r+1), including a positive geometric tail bound."""
    partial = sum(z**r/(2*r+1) for r in range(n))
    return Interval(partial, partial+z**n/((2*n+1)*(1-z)))


def reciprocal(x):
    if x.lo <= 0:
        raise ArithmeticError('Positive divisor required')
    return Interval(1/x.hi, 1/x.lo)


def triangle_gap(xi, pi, count=180):
    k = Q(15,17)
    inv_x = reciprocal(xi)
    b = enclosure(1) - 4*inv_x
    if not 0 < b.lo <= b.hi < 1:
        raise ArithmeticError('Settings must lie above two-scalar threshold')
    j = enclosure(1) - b*Interval(atanh_over_argument(b.lo,count).lo,
                                  atanh_over_argument(b.hi,count).hi)
    j = outward(j)
    coefficient = k
    real = enclosure(0)
    for r in range(count):
        real += (4*coefficient/(2*r+1))*j
        j = outward(Q(1,2*r+3) + b*j)
        coefficient *= k*k
    tail = 4*coefficient/((2*count+1)*(1-k*k)*(1-b.hi))
    real = (real + Interval(-tail,tail))*inv_x
    beta = Interval(sqrt_bounds(b.lo).lo, sqrt_bounds(b.hi).hi)
    series = Interval(atanh_over_argument(k*k*b.lo,count).lo,
                      atanh_over_argument(k*k*b.hi,count).hi)
    imag = 2*k*inv_x*pi*beta*series
    return (outward(real),outward(imag))


def add(z, w):
    return z[0]+w[0], z[1]+w[1]


def scale(c, z):
    return c*z[0], c*z[1]


def multiply(z, w):
    return z[0]*w[0]-z[1]*w[1], z[1]*w[0]+z[0]*w[1]


def modulus_squared(z):
    return z[0].square()+z[1].square()


def rational_strings(interval, digits):
    ends = outward(interval,digits)
    return [str(ends.lo),str(ends.hi)]


def coefficient(energies, pi):
    values = [triangle_gap(x, pi) for x in energies]
    # Constant and linear nuisance columns allow fixed coordinate shifts.
    shift = (enclosure(Q(2,5)), enclosure(Q(11,20)))
    centered = [add(z,scale(-1,shift)) for z in values]
    top = (enclosure(0),enclosure(0))
    norm = enclosure(0)
    for omitted in range(4):
        i,j,k = [row for row in range(4) if row!=omitted]
        # Determinant of the remaining three rows with columns [f,1,xi].
        l = add(scale(energies[k]-energies[j], add(centered[i],scale(-1,centered[j]))),
                scale(energies[j]-energies[i], add(centered[k],scale(-1,centered[j]))))
        l = scale((-1)**omitted, l)
        energy = energies[omitted]
        if modulus_squared(l).lo <= 0:
            raise ArithmeticError('A three-setting nuisance minor could vanish')
        top = add(top, multiply(l, scale(energy*Q(1,2)-Q(7,2),centered[omitted])))
        norm += 4*reciprocal(energy)*modulus_squared(l)
    numerator = modulus_squared(top)
    if norm.lo <= 0 or numerator.lo <= 0:
        raise ArithmeticError('Four-energy annihilator or projection unresolved')
    bound = Interval(numerator.lo/(norm.hi*2048*pi.hi**5),
                     numerator.hi/(norm.lo*2048*pi.lo**5))
    return bound, values


def main():
    # Machin's identity with explicit alternating series brackets.
    pi = outward(16*atan_unit_fraction(5)-4*atan_unit_fraction(239))
    points = [5,6,9,11]
    bound, values = coefficient([enclosure(x) for x in points],pi)
    lower,upper=Q('1.33764501629e-9'),Q('1.33764501631e-9')
    if not lower < bound.lo <= bound.hi < upper:
        raise ArithmeticError('Certified interval exceeds published bracket')
    halfwidth = Q(1,1000)
    windows = [Interval(x-halfwidth, x+halfwidth) for x in points]
    uniform_bound, _ = coefficient(windows,pi)
    if uniform_bound.lo <= Q('1e-10'):
        raise ArithmeticError('Unresolved finite-width result: '+str(float(uniform_bound.lo)))
    result = {
        'method':'exact Fraction arithmetic, explicit two-sided tails and rational square-root brackets',
        'settings':points,'exposures':['1/4']*4,
        'principal_value_terms':180,
        'coefficient_bracket_strict':['1.33764501629e-9','1.33764501631e-9'],
        'outward_rational_coefficient':rational_strings(bound,23),
        'four_windows_halfwidth_in_xi':'1/1000',
        'uniform_coefficient_lower_bound':str(outward(uniform_bound,18).lo),
        'triangle_gap_rectangles':[
            {'xi':x, 'Re':rational_strings(f[0],19), 'Im':rational_strings(f[1],19)}
            for x,f in zip(points,values,strict=True)],
    }
    Path(__file__).with_name('SCAN_CERTIFICATE_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: rational bounds on four complex triangles and the single-direction profile coefficient')


if __name__=='__main__':
    main()
