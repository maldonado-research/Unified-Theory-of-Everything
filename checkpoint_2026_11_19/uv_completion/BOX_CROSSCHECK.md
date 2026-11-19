# Direct Feynman-parameter check of the absolute box rate

For the four-singlet benchmark with $M/m=10$ and $y=g=1$, direct integration gives the leading polarized cross sections $m^2\sigma=4.009254\times10^{-11}$ at $s=6m^2$ and $7.174007\times10^{-11}$ at $s=10m^2$. These include the dispersive part of both exchanged boxes and all their angular dependence. The imaginary part, computed from the same Feynman integral, independently reproduces the signed scalar-pair cut. [Numerical receipt](BOX_DIRECT.json) · [On-shell sewing](CUT.json)

## The box integral

Use units $m=1$ and the external phase $B_{34}^*B_{12}=s$. Denote the external scattering angle by $z=\cos\theta$, with $t=-s(1-z)/2$ and $u=-s(1+z)/2$. For muon squared mass and residue $(b,w_b)=(1,1),(16,-1)$ and electron squared mass and residue $(c,w_c)=(100,10),(400,-40)$, the amplitude is
$$
\mathcal M(s,z)=\frac{s}{16\pi^2}\sum_{b,c} w_b w_c
  \big[D(s,t;b,c)+D(s,u;b,c)\big],\qquad
D=\int_{x_i\ge0}\! d^4x\,
 \frac{\delta(1-\sum_i x_i)}
 {(x_1+x_2+b x_3+c x_4-s x_1 x_2-t x_3 x_4-i0)^2}.
$$
The two scalar lines have parameters $x_1,x_2$ and the muon and electron singlets have $x_3,x_4$. The equal-chirality propagators furnish the two masses in $w_b w_c$. Exchanging the identical external leptons interchanges $t$ and $u$. The exchange sign combines with the spinor sign to give the displayed sum. The normalized four-propagator loop and the three-dimensional simplex supply $1/(16\pi^2)$ in the amplitude. Power counting gives a UV-convergent scalar integral for each individual mass pair. All internal masses are positive.

Set $x=x_2$, $y=x_3$, $L=1-x$, and $x_4=L-y-x_1$. At fixed $x,y$ the denominator before squaring is linear in $x_1$. Its real parts at $x_1=0$ and $x_1=L-y$ are
$$
 A(y)=x+cL+(b-c-tL)y+t y^2,\qquad d+k y=1-sxL+(b-1+sx)y .
$$
Integration of its inverse square gives $(L-y)/[(A(y)-i0)(d+k y-i0)]$. The first endpoint is positive throughout the simplex for physical $t\le0$. For $0<x<1$, $k=b-1+sx>0$. Writing $p=-d/k$ and $H(y)=(L-y)/A(y)$, the remaining $y$ integral near a pole is
$$
 \frac1k\left[
 \int_0^L\frac{H(y)-H(p)}{y-p}\,dy+
 H(p)\left(\log\left|\frac{L-p}{p}\right|+
 i\pi\,\mathbf1_{0<p<L}\right)\right].
$$
The checker forms the divided difference algebraically. For a distant pole outside the interval it integrates the smooth original expression. The physical pole lies in the integration interval for $(1-\beta)/2<x<(1+\beta)/2$, where $\beta=\sqrt{1-4/s}$. Splitting the $x$ interval at those two logarithmic endpoints and using a squared-sine substitution on each panel regularizes the quadrature. This calculation uses the $-i0$ prescription directly and requires no subtraction constant or dispersive input. [Executable](check_box_direct.py)

## Normalization and numerical convergence

For fixed incoming helicities the flux and identical-electron phase space give
$$
 m^2\sigma=\frac{1}{32\pi s}\int_0^1 |\mathcal M(s,z)|^2\,dz
$$
in the dimensionless variables above. The integrand is even in $z$ because the exchanged boxes are summed. The separate $J=0$ lower bound uses the square of the mean imaginary part. Integrating the square of the pointwise imaginary part keeps its higher partial waves.

| $s/m^2$ | $\mathcal M(s,0)$ | $\langle\mathrm{Im}\,\mathcal M\rangle$ | $m^2\sigma$ from the imaginary part | Full $m^2\sigma$ |
|---:|---:|---:|---:|---:|
| 6 | $-1.4143754\,10^{-4}-6.5067129\,10^{-5}i$ | $-6.5042088\,10^{-5}$ | $7.0135500\,10^{-12}$ | $4.0092537\,10^{-11}$ |
| 10 | $-1.9338113\,10^{-4}-1.8757071\,10^{-4}i$ | $-1.8715511\,10^{-4}$ | $3.4842175\,10^{-11}$ | $7.1740071\,10^{-11}$ |

The Gauss–Legendre refinements $(N_{x,y},N_z)=(64,6),(128,10),(224,16)$ give the following total rates. The last two refinements differ by less than $4\times10^{-9}$ in relative terms. These are numerical convergence comparisons. The receipt also records the right-angle and mean amplitudes on each grid.

| $s/m^2$ | $(64,6)$ | $(128,10)$ | $(224,16)$ |
|---:|---:|---:|---:|
| 6 | $4.00925569\,10^{-11}$ | $4.00925371\,10^{-11}$ | $4.00925370\,10^{-11}$ |
| 10 | $7.17400981\,10^{-11}$ | $7.17400710\,10^{-11}$ | $7.17400709\,10^{-11}$ |

Run `python -I -B -O check_box_direct.py` in this directory. The checker verifies convergence, compares the integrated imaginary part with the analytic scalar-pair cut, and checks both integrated bounds. To write the same receipt to a chosen file, append `--output PATH`. The reported absolute rates use the specified renormalizable interactions at their leading nonzero amplitude order. [Checker](check_box_direct.py) · [Receipt](BOX_DIRECT.json)
