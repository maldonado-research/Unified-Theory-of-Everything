# How the scalar continuum enters the earlier four-rate test

A finite change in the relative electron Yukawa coupling cancels in a particular combination of three direct rates and one reverse rate at leading order. The scalar propagator insertion calculated here survives that cancellation. Its contribution to the combination has opposite signs for the muon and electron continua. A residual in this test can therefore receive Yukawa corrections even when the scalar quartic vanishes. [Leading consistency relation](../radiative_stability/CLOSURE.md)

Use $s/m^2=6,8,10$, in that order, and keep the scalar pole mass fixed. The leading direct rate after removing the common coupling strength is $Q_i=a_i+2\operatorname{Re}(b_i\delta)+c_i|\delta|^2$, where $\delta$ changes the second electron singlet's squared Yukawa coupling. Define $d_i=2\operatorname{Re}b_i/a_i$, $e_i=c_i/a_i$, and $f_i=-2\operatorname{Im}b_i/a_i$. The normalized direct rates can be expressed as $r_i=U+d_i V+e_i W+f_i T$. One inverse rate at the first setting determines $T$ from the signed difference. At leading order $U=N$, $V=N\operatorname{Re}\delta$, $W=N|\delta|^2$, $T=N\operatorname{Im}\delta$ for any common normalization $N>0$. Hence $\mathcal I=W/U-(V/U)^2-(T/U)^2=0$ throughout the single-detuning model.

With the tuned complex phases, both the continuum insertion and the explicit renormalization-scale logarithm give equal reciprocal rates at the indicated order. Introduce a bookkeeping parameter $h$ multiplying one such calculated contribution, and its fractional shifts $\ell_i$ through $r_i=N(1+h\ell_i)$. Then $T=0$ and the remaining coefficients are determined by the matrix $D$ with rows $(1,d_i,e_i)$. At $h=0$, $V=W=0$, so only $W$ contributes to the derivative of $\mathcal I$,
$$
\left.\frac{d\mathcal I}{dh}\right|_0=(D^{-1}\ell)_3
=\frac{(d_8-d_6)(\ell_{10}-\ell_6)-(d_{10}-d_6)(\ell_8-\ell_6)}
{(d_8-d_6)(e_{10}-e_6)-(d_{10}-d_6)(e_8-e_6)} .
$$
If a correction multiplies all rates by the same factor, the numerator vanishes. The rate ratio $R=\sigma_{10}/\sigma_6$ has derivative $R_0(\ell_{10}-\ell_6)$ with $R_0=1.7893622161$.

| Perturbation with its own coefficient | Shift in $R$ per unit coefficient | Shift in $\mathcal I$ per unit coefficient |
|---|---:|---:|
| Finite scalar insertion proportional to $\lvert y_\mu\rvert^2$ | $+0.00766304$ | $+2.70327\times10^{-5}$ |
| Finite scalar insertion proportional to $\lvert y_e\rvert^2$ | $-0.000329712$ | $-1.18345\times10^{-6}$ |
| Previously calculated scalar quartic proportional to $\lambda$ | $+0.725044$ | $+5.21219\times10^{-5}$ |

The finite rows are coefficients from the specified diagrams in the pole-normalized scalar convention. At nonzero Yukawa couplings the same coupling powers occur in additional diagrams. Individual rows cannot be used as physical derivatives from varying a Yukawa coupling across the entire theory. The values quantify how accurately those other contributions must be controlled before attributing a measured residual to $\lambda$. Both the direct simplex calculation and the separately evaluated QCDLoop spectrum give the finite projections in the table. [Finite-amplitude results](CONTINUUM.json) · [External integration](QCDLOOP_CONTINUUM.json)

For the separate **explicit logarithm** $L=\ln(\mu/\mu_{\mathrm{ref}})$, at fixed renormalized parameters and to this order the coefficients of $L\lvert y_\mu\rvert^2$ and $L\lvert y_e\rvert^2$ in $\mathcal I$ are $+2.43232\times10^{-6}$ and $-1.04009\times10^{-5}$. These use the running of the entire leading amplitude. The implicit running of the physical leading rates provides the opposite scale dependence when compared to a fixed reference design. The scale logarithms constrain the complete two-loop answer independently of the finite scalar insertion. [Renormalization conventions](RG_LOG.md)

The verifier checks Cramer's rule against a numerical derivative of the nonlinear rate combination, tests the linear combination of log rates from the previous checkpoint, and recovers its quartic reference. [Program](check_projection.py) · [Numerical results](PROJECTION.json)
