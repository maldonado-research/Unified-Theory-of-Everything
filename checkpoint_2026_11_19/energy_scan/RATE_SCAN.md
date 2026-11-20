# A reciprocal-rate bound across collision energies

In the flavor-isolated channel $\mu^-_L\mu^-_L\leftrightarrow e^-_Le^-_L$, an independent complex Weinberg coefficient $C_{5,ee}$ can cancel the scalar-pair cut at any chosen energy. Across a scan, the same coefficient must cancel a term linear in $s$ at every energy. For nonzero $C_{7,ee}w_0$, two distinct energies above threshold already suffice to leave a strictly positive, calculable contribution to the exposure-averaged reciprocal rate.

## The spectral weight and the observable

Use massless external leptons, common scalar squared mass $a=m_\phi^2$, sterile masses $(m_\phi,4m_\phi)$ and weights $w_j=m_j(Y^*_{\mu j})^2=(w_0,-w_0)_j$ in the positive-mass basis. Only the muon Yukawa row and the two electron-pair seeds enter. Use the seed sign convention in which their common on-shell factor is $(s/2-a)C_{7,ee}-C_{5,ee}$. Define

$$
\xi=s/a,\quad n=\xi/2-1,\quad \beta=\sqrt{1-4/\xi},\quad
h=\xi\beta/2,\quad g_7=aC_{7,ee}w_0,\quad g_5=C_{5,ee}w_0 .
$$

Both $g_i$ are dimensionless. In terms of the transfer-angle integral $J_b$ with squared sterile mass $b$, the dimensionless spectral shape is

$$
\begin{aligned}
j(\xi)&=a\beta(J_a-J_{16a})\\
 &=\frac{\beta}{2}\int_{-1}^1
 \frac{15\,dx}{(1+n-hx)(16+n-hx)}\\
 &=\frac{2}{\xi}\left[
 \operatorname{artanh}\frac{h}{1+n}
 -\operatorname{artanh}\frac{h}{16+n}\right].
\end{aligned}
$$

Since $1+n-h=\xi(1-\beta)/2>0$, the integral proves $j>0$ at every finite $\xi>4$. For $g_5=0$ the derivative-seed spectral gap is $n j$, which reproduces the preceding same-component cut. With the ordered external spinor factor of magnitude $s$, the anti-Hermitian amplitude has magnitude $|A_{fi}|=\xi j |n g_7-g_5|/(8\pi)$. [Transfer integral](../on_shell/INDEPENDENT_CUT.md) · [Channel and normalization](../on_shell/INVERSE_RATE_BOUND.md)

Let $\bar\sigma=(\sigma_{\mu\mu\to ee}+\sigma_{ee\to\mu\mu})/2$ for fixed incident helicities and identical outgoing leptons. Hermitian analyticity gives $M_{fi}=H_{fi}+iA_{fi}$ and $M_{if}^*=H_{fi}-iA_{fi}$ with the full two-particle states reversed. Thus the mean of squared amplitudes at each angle is $|H_{fi}|^2+|A_{fi}|^2$. The angle-independent cut contributes $|A_{fi}|^2/(32\pi s)$ after full-sphere integration including the identical-final-state factor. For a normalized nonnegative common exposure $dP(\xi)$ in both directions, the result is

$$
 a\langle\bar\sigma\rangle_P\ \geq\
 Q(g_5)=\int k(\xi)|g_7n(\xi)-g_5|^2\,dP(\xi),
 \qquad k(\xi)=\frac{\xi j(\xi)^2}{2048\pi^3}.
$$

The observable averages cross sections incoherently across runs. Hermitian local four-lepton amplitudes contribute through the nonnegative integrated $|H_{fi}|^2$ term.

## Minimization over a complex Weinberg coefficient

Set $S_r=\int k n^r\,dP$ for $r=0,1,2$, assuming finite moments and $S_0>0$. For the full complex plane of $g_5$, completion of the square yields

$$
 Q(g_5)=S_0\left|g_5-g_7\frac{S_1}{S_0}\right|^2
   + |g_7|^2 R_P,\qquad R_P=S_2-\frac{S_1^2}{S_0}.
$$

The unique minimizer of the cut contribution is $g_5=g_7 S_1/S_0$. The minimizing ratio is real because $n$ is real and the exposure is positive. A relative imaginary part increases the rate contribution. Taking $d\mu=k\,dP$ gives a direct positivity test,

$$
 R_P=\int (n-S_1/S_0)^2\,d\mu
 =\frac{1}{2S_0}\iint [n(\xi)-n(\xi')]^2\,d\mu(\xi)d\mu(\xi').
$$

Thus $a\langle\bar\sigma\rangle_P\ge |g_7|^2 R_P>0$ whenever $g_7\ne0$ and the effective measure samples more than one energy above threshold. At a single energy $\xi_*$ the choice $C_{5,ee}=a n(\xi_*)C_{7,ee}$ nulls the cut there and gives $R_P=0$. The dispersive rate may still be nonzero. For a nontrivial scan the optimized residual is sharp for the cut contribution. Saturation of the bound on the entire reciprocal mean also requires $H_{fi}=0$ almost everywhere in angle and energy. Optimizing the cut alone makes no claim of such simultaneous saturation.

## Equal-exposure benchmarks

| Normalized exposure in $\xi$ | Optimal $g_5/g_7$ | $R_P$ in $a\langle\bar\sigma\rangle_P\geq |g_7|^2R_P$ |
|---|---:|---:|
| Point masses $1/2$ at 6 and 10 | 3.14037330148265348 | $3.78076496757922739\times10^{-6}$ |
| Uniform density $1/4$ on $[5,7]\cup[9,11]$ | 3.17083201734730432 | $3.87456007022534449\times10^{-6}$ |

For the point scan, $n_1=2$, $n_2=4$. Writing $v_i=k(\xi_i)/2$ gives $R_P=v_1v_2(n_2-n_1)^2/(v_1+v_2)$. Each finite bin carries half of the total exposure, uniform in squared collision energy. The finite-bin calculation integrates the actual kernel and the variation of $n$ within each bin. A scan uniform in $\sqrt{s}$ would use a different measure.

The checker integrates the analytic logarithms by tanh–sinh quadrature and evaluates the centered square independently using Gauss–Legendre quadrature over both transfer angle and scan energy. The two rate-floor values agree between methods to absolute differences below $4\times10^{-38}$. It also tests the complex square completion and recovers the derivative-only benchmark at $\xi=8$. Run `python verify_rate_scan.py` in this directory to reproduce the receipt. [Executable](verify_rate_scan.py) · [Numerical moments and residuals](RATE_SCAN_RESULTS.json)

## Scope of the rate statement

This calculation uses the square of the leading flavor-changing amplitude, linear in either seed and quadratic in the Yukawa couplings. The two coefficients are fixed at a common reference scale for all scan energies. Allowing a separate arbitrary $g_5$ at every energy removes the cross-energy constraint. Corrections from running enter at their respective perturbative orders. In a general EFT expansion the square of a dimension-seven contribution can occur at the same inverse-scale order as interference between dimension-five and dimension-nine amplitudes. Applying the bound to an observable requires control of any additional terms retained at the claimed order.

In the specified sector a scalar decay to a massless lepton and the lighter sterile state has zero phase space, and the heavier sterile state lies above that decay threshold. The exchanged sterile propagator is spacelike on the scalar cut. These facts permit the fixed-order, stable-scalar cut used here. A model with an open scalar decay channel or additional coherent cuts requires the discontinuity for that spectrum, including physical decay products where appropriate. Additional independent on-shell operator structures can also change the projection and its rank. The result constrains the reciprocal mean of the restricted amplitudes and sets no separate floor on either direction. [Channel and threshold scope](../on_shell/CHANNEL_SCOPE.md)
