# How a one-percent change in the electron couplings changes the rate

Fix the four singlet masses, the muon Yukawa row and the first electron Yukawa coupling in the [four-singlet benchmark](../uv_completion/00_READ_FIRST.md). Vary the second **squared** electron Yukawa through $Y_{e2}^2=-2y_e^2(1+\delta)$. At the same leading four-Yukawa amplitude order, the direct cross-section ratio $\sigma(10m^2)/\sigma(6m^2)$ ranges from 1.48822 to 2.53590 for a complex $\delta$ of magnitude up to 0.01. The tuned value is 1.78936. This interval describes a choice of model inputs. An experimental uncertainty would require measured masses and couplings.

## Rate as a quadratic function of complex detuning

Let $\xi=s/m^2$, choose $M=10m$, and strip $\gamma=y_e^2(y_\mu^*)^2$ from the amplitude for $\mu^-_L\mu^-_L\to e^-_Le^-_L$. The electron and muon mass residues in scalar-mass units are $(10,-40[1+\delta])$ and $(1,-1)$. Write $A_0(z)$ for the tuned direct amplitude at scattering cosine $z$ and $A_1(z)$ for its contribution from the second electron singlet at the tuned residues. The exact leading amplitude divided by $\gamma$ is $A_0+\delta A_1$. Both functions include the two scalar pairings. The matching coefficients in the original fields become $C_5=-\delta y_e^2/M$ and $C_7=(3-\delta)y_e^2/(4M^3)$. [Tree matching](../uv_completion/MATCHING.md)

With the physical upper-rim propagator prescription, define the angular inner product by
$$
 \langle f,g\rangle_\xi=\frac{1}{32\pi\xi}\int_0^1 f(z)^*g(z)\,dz.
$$
The scalar functions are even in $z$ after removing the external spinor phase. The prefactor includes the identical-electron factor at fixed incident helicities. Put $a_\xi=\langle A_0,A_0\rangle_\xi$, $b_\xi=\langle A_0,A_1\rangle_\xi$ and $c_\xi=\langle A_1,A_1\rangle_\xi$. Then
$$
 Q_\xi(\delta)\equiv\frac{m^2\sigma_{\mu\mu\to ee}}{|\gamma|^2}
 =a_\xi+2\operatorname{Re}(\delta b_\xi)+|\delta|^2c_\xi .
$$

| $\xi$ | $a_\xi$ | $\operatorname{Re} b_\xi$ | $\operatorname{Im} b_\xi$ | $c_\xi$ |
|---:|---:|---:|---:|---:|
| 6  | $4.00925370\times10^{-11}$ | $1.34882615\times10^{-9}$ | $6.41879498\times10^{-10}$ | $5.56548417\times10^{-8}$ |
| 8  | $5.61821543\times10^{-11}$ | $1.47545212\times10^{-9}$ | $8.15740365\times10^{-10}$ | $5.05926136\times10^{-8}$ |
| 10 | $7.17400709\times10^{-11}$ | $1.58692497\times10^{-9}$ | $9.10283235\times10^{-10}$ | $4.66542038\times10^{-8}$ |

At zero detuning the direct log-rate derivatives with respect to $\operatorname{Re}\delta$ are $2\operatorname{Re}b_\xi/a_\xi$. They are 67.28564698, 52.52387139 and 44.24096474 in ascending energy order. The reciprocal mean has the same real derivative and zero imaginary-detuning derivative at the origin.

For real detuning the rates in the two directions coincide. The entries below include all scattering angles and both the real and imaginary parts of the loop.

| $\delta$ | $10^{11}Q_6$ | $10^{11}Q_{10}$ | $Q_{10}/Q_6$ |
|---:|---:|---:|---:|
| $-0.01$ | 1.86815 | 4.46670 | 2.39097 |
| $0$ | 4.00925 | 7.17401 | 1.78936 |
| $+0.01$ | 7.26345 | 10.81440 | 1.48888 |

In the inverse reaction the Yukawa factors are conjugated and each scalar integral keeps the upper-rim prescription. Its scaled rate is $Q_\xi(\delta^*)$. For a complex perturbation the difference between the direct and inverse rates is $-4\operatorname{Im}\delta\operatorname{Im} b_\xi$ in these units. At $\delta=+0.01i$, the difference divided by the sum is $-0.28117$ at $\xi=6$ and $-0.23828$ at $\xi=10$. This exclusive directional asymmetry requires both the relative Yukawa phase and the relative absorptive phases of the two mass contributions. The input perturbation introduces that Yukawa phase. A real scalar quartic with the tuned relative Yukawas supplies no such phase.

## Range over a complex disk

Take $\rho=0.01$. The denominator is positive throughout the disk since $\sqrt{a_6}>\rho\sqrt{c_6}$ and the angular norm obeys the triangle inequality. On its boundary set $L=a_6+\rho^2c_6$ and $H=a_{10}+\rho^2c_{10}$. A trial real ratio $r$ is reached if $H-rL+2\operatorname{Re}[\delta(b_{10}-rb_6)]=0$. Its two boundary extrema therefore solve
$$
 (H-rL)^2=4\rho^2 |b_{10}-rb_6|^2 .
$$
At both roots the difference $Q_{10}-rQ_6$ has negative curvature $c_{10}-rc_6$ as a function of the real and imaginary parts of $\delta$. Its minimum over the disk lies on the boundary, establishing the lower endpoint. At the upper root its unconstrained maximum is at radius $|b_{10}-rb_6|/|c_{10}-rc_6|=0.0208395>\rho$, establishing the upper endpoint on the same boundary. Complex perturbations $0.00994327-0.00106365i$ and $-0.00931247+0.00364390i$ attain the lower and upper ratios, respectively, up to rounding. Continuity supplies the intermediate values.

The numerical calculation integrates the two electron-mass contributions separately using the fixed-transfer dispersion formula and compares the entire complex Gram matrix with direct Feynman-simplex quadrature. A second dispersion grid checks refinement. The tuned rates at the two previously published energies reproduce the QCDLoop values. The comparisons apply to the declared model at leading four-Yukawa amplitude order. Computing a radiative correction to the same model requires the amplitudes and renormalization conditions at the relevant loop order. Substituting both a loop-induced low-momentum contact and the same loop graph would count that effect twice. [Dispersion integral](../uv_completion/DISPERSION.md) · [Direct box integration](../uv_completion/BOX_CROSSCHECK.md) · [Checker](check_detuning.py) · [Results](DETUNING.json)
