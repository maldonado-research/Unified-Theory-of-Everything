# External QCDLoop check of the scalar continuum

For the scalar propagator continuum insertion alone, the normalized rate slopes at $s=8m^2$ are $2.78740068115\times10^{-12}$ for $|g|^2$ and $4.84448606787\times10^{-14}$ for $|y|^2$. The electron block contributes through interference with the real part of the leading amplitude. Its first continuum threshold is $s=121m^2$, above this scan. [Numerical receipt](QCDLOOP_CONTINUUM.json) · [Standalone calculation](qcdloop_continuum.cpp)

## Normalization and calculation

Set $m=1$, $M=10$, and strip the common external coupling and spinor phase as in the leading box benchmark, retaining the spinor modulus $s$. Write $t=-s(1-z)/2$, $u=-s(1+z)/2$. The external library computes the finite scalar box $D(U,b,1,d;s,t)$ with squared masses $(U,b,1,d)$ and invariants $(0,0,0,0,s,t)$ on the upper rim. The driver checks both Laurent poles at each point. Define

$$
A_U=\frac{s}{16\pi^2}\sum_{b,d} w_b w_d
  [D(U,b,1,d;s,t)+D(U,b,1,d;s,u)].
$$

Here $(b,w_b)=(1,1),(16,-1)$ and $(d,w_d)=(100,10),(400,-40)$. The leading amplitude is $A_1$. The coupling factors in the continuum have positive weights because they come from a scalar self-energy cut. With $f_B(U)=\Theta(U-B)(U-B)^2/[U(U-1)^2]$, the coefficients of the two extra coupling moduli are

$$
\Delta_g=\frac{2}{16\pi^2}\int_1^\infty dU\, A_U[f_1(U)+f_{16}(U)/4],\qquad
\Delta_y=\frac{2}{16\pi^2}\int_1^\infty dU\, A_U[f_{100}(U)+2f_{400}(U)].
$$

The factor two counts the two scalar lines. Interchanging their Feynman parameters leaves each box unchanged because the four external masses vanish. For $U>1$, $f_1(U)=1/U$. For $f=g,y$ the reported coefficient $S_f$ in the normalized rate is $2\int_0^1 dz\,\mathrm{Re}(A_1^*\Delta_f)/(32\pi s)$. This flux includes the identical final electrons and the specified incident polarization. The absolute rate contains the common factor $|g|^4|y|^4$. Other two-loop diagrams would be needed for a prediction of the complete Yukawa correction. [Leading normalization and QCDLoop routing](../uv_completion/QCDLOOP.md) · [QCDLoop publication](https://doi.org/10.1016/j.cpc.2016.07.033)

| $s/m^2$ | $\mathrm{Re}\,\Delta_g(s,0)$ | $\mathrm{Im}\,\Delta_g(s,0)$ | $\mathrm{Re}\,\Delta_y(s,0)$ | $S_g$ | $S_y$ |
|---:|---:|---:|---:|---:|---:|
| 6 | $-3.866832770\times10^{-6}$ | $-3.997025444\times10^{-7}$ | $-8.244942824\times10^{-8}$ | $1.895420940\times10^{-12}$ | $3.859105087\times10^{-14}$ |
| 8 | $-5.608333492\times10^{-6}$ | $-1.264685994\times10^{-6}$ | $-1.133695856\times10^{-7}$ | $2.787400681\times10^{-12}$ | $4.844486068\times10^{-14}$ |
| 10 | $-7.283110684\times10^{-6}$ | $-2.462922767\times10^{-6}$ | $-1.459595699\times10^{-7}$ | $3.698825419\times10^{-12}$ | $5.583439398\times10^{-14}$ |

## Integration and ultraviolet check

The driver splits at $U=1,(\sqrt{s}-1)^2,16,100,400$. Only the first interval can contribute an imaginary part. On each finite interval $[l,h]$, the map $U=l(h/l)^{\sin^2(\pi x/2)}$ resolves the square-root cusp at the physical endpoint from both sides. The map $U=400/x^2$ integrates the entire tail on $0<x<1$. The signed flavor sum and exchanged boxes are evaluated in GNU binary128 before conversion to binary64 for the mass and angular integrals. The mass and angular Gauss–Legendre orders are $(10,4),(18,6),(28,10)$. The last two grids differ by less than $2.1\times10^{-14}$ relative for the reported nonzero central components, angular means and rate coefficients. The raw imaginary parts for the closed electron block are below $10^{-36}$ in magnitude. Grid differences measure numerical convergence.

The large-$U$ convergence follows directly from the positive Feynman denominator. Fix the sum of the three light-line parameters $x_s+x_\mu+x_e=1$ and let the heavy parameter $v$ range from zero to infinity. With $L=x_s+b x_\mu+d x_e$ and $C=L-t x_\mu x_e$, the scalar box is

$$
D=\int_{x_\mu+x_e\le1}dx_\mu\,dx_e\int_0^\infty
\frac{dv}{[Uv^2+(U+L-sx_s)v+C]^2}.
$$

The triangle domain has nonnegative coordinates. For $U>s$, $C\ge1$ and $U+L-sx_s\ge U+1-s$. Integrating the upper bound obtained by dropping $Uv^2$ gives $0<D\le 1/[2(U+1-s)]$. Every active spectral weight obeys $0\le f_B\le 1/U$. Each flavor term in the mass integral therefore has an absolutely integrable tail of order $U^{-2}$.

To check cancellation in the signed sum at large $U$, the driver also computes $UD$ by putting $V=Uv$, $q=1/U$, $H=1+q(L-sx_s)$ and $K=H^2-4qC$. At the three tested large masses $K>0$. The $V$ integral is

$$
\frac{H/C-2q\,\log[(H+\sqrt K)^2/(4qC)]/\sqrt K}{K}.
$$

A 64-point Gauss–Legendre rule per axis integrates the remaining triangle with $x_e=h^2$, $x_\mu=(1-x_e)k$. At $s=8$, $z=0$ and $U=4\times10^4,10^7,2\times10^8$, the relative difference of the signed sum from QCDLoop is below $9\times10^{-14}$. The highest test mass exceeds the largest quadrature node. [Check implementation](qcdloop_continuum.cpp) · [Raw comparisons](QCDLOOP_CONTINUUM.json)

## Reproduction

Use x86-64 Linux with Bash, Git, GNU C++ and libquadmath. The recorded run used GCC 14.2.0. Start in a scratch directory containing the driver. The upstream checkout and compiler are needed to regenerate the external receipt. [Pinned QCDLoop 2.1.0 source](https://github.com/scarrazza/qcdloop/tree/52d94ce9e628c1b1f9784f66a6ef9de285359ae3)

```bash
git clone --depth 1 --branch 2.1.0 https://github.com/scarrazza/qcdloop.git
test "$(git -C qcdloop rev-parse HEAD)" = 52d94ce9e628c1b1f9784f66a6ef9de285359ae3 || exit 1
mkdir -p build/qcdloop
sed 's/@VERSION@/2.1.0/g' qcdloop/src/qcdloop/config.h.in > build/qcdloop/config.h
g++ -std=gnu++11 -fext-numeric-literals -O2 -Ibuild -Iqcdloop/src \
  qcdloop_continuum.cpp qcdloop/src/{box,cache,topology,tools,types}.cc \
  -lquadmath -o build/qcdloop_continuum
./build/qcdloop_continuum --output QCDLOOP_CONTINUUM.json
```
