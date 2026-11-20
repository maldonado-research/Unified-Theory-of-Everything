# Internal sterile self-energy in the two-scalar cut

With all four Majorana masses and Yukawas in full-theory $\overline{\rm MS}$ at $\mu=m$, a sterile self-energy insertion contributes $-2.615099745\,10^{-8}|y_e|^2-3.041359246\,10^{-8}|g|^2$ to the reduced $J=0$ two-scalar cut at $s=8m^2$. The reduction divides by $\gamma=y_e^2(g^*)^2$ and uses $B_{34}^*B_{12}=s$. The scalar pole mass $m$ and the external kinematics are fixed. This result concerns the internal sterile lines of the two production amplitudes. [Numerical receipt](STERILE.json) · [Leading cut and partial-wave convention](../../checkpoint_2026_11_19/uv_completion/CUT.md)

## Renormalized two-point function

Write $p^2=-q<0$ and $H=Y^\dagger Y$. For the action $-\bar L Y\phi N_R+\mathrm{h.c.}$, the subtracted correction to the inverse Majorana propagator with its overall $i$ stripped is
$$
 \delta K(p)=\frac{f(q;\mu)}{16\pi^2}
 \not p\,[H P_R+H^T P_L],\qquad
 f(q;\mu)=-2\int_0^1 dx\,(1-x)
 \log\frac{(1-x)(m^2+xq)}{\mu^2}.
$$
The active line is massless. Combining its denominator with the scalar denominator shifts its numerator to $(1-x)\not p$ under the integral. The factor two sums the two weak components of that internal doublet. The charge-conjugate Majorana chain supplies the transposed matrix multiplying the opposite projector. The massless numerator and the complex-scalar contractions give zero proper chirality-flipping mass correction. The sterile kinetic pole cancels against $\delta Z_R$, and the $\overline{\rm MS}$ mass-parameter counterterm cancels the field-counterterm contribution to the mass coefficient in the action. Both counterterms therefore leave the finite kinetic integral above. [Counterterm conventions and proper two-point poles](../../checkpoint_2026_10_16/yukawa_counterterms/YUKAWA_TWO_POINT.md)

For $r=q/m^2$, integration gives
$$
 f(q;\mu)=\log\frac{\mu^2}{m^2}+2+\frac1r-
                \frac{(1+r)^2}{r^2}\log(1+r).
$$
At $q=0$ its continuous value is $\log(\mu^2/m^2)+1/2$. The integral is real throughout the physical scalar-pair cut since every tree exchange is spacelike. At fixed spacelike transfer the massless-scalar limit gives $2-\ln(q/\mu^2)$. Buchmüller and Plümacher use $S^{-1}=\not p-M-\Sigma$ and obtain the same inverse-kernel coefficient from their massless bubble. [Buchmüller and Plümacher, Eqs. (3)–(7)](https://arxiv.org/abs/hep-ph/9710460)

## Inserting the Majorana propagator

In the positive Takagi basis, $S_i=(\not p+M_i)/(p^2-M_i^2)$. Stripping propagator factors of $i$, expansion of $(K_0+\delta K)^{-1}$ gives $\delta S=-S\delta K S$. With outgoing electron endpoints the equal-chirality block of $S_i\not p[H_{ij}P_R+H_{ji}P_L]S_j$ has numerator
$$
 p^2\,(M_jH_{ij}+M_iH_{ji})P_R.
$$
In particular an ordered pair $i,j$ can mix different Takagi mass eigenstates. For the reciprocal incoming tree one conjugates the complete outgoing flavor contraction. This is also the conjugate kinematic kernel at these spacelike transfers.

In each disjoint flavor block write $Y_{\ell i}=y_\ell c_i$, where $y_\mu=g$, $c_e=(1,i\sqrt2)$ and $c_\mu=(1,-i/2)$. With the tree written as $-B^* (y_\ell)^2 F_\ell(q)$ per outgoing scalar assignment, define
$$
 F_\ell(q)=\sum_i\frac{M_i c_i^2}{M_i^2+q},\qquad
 L_\ell(q)=\sum_i\frac{|c_i|^2}{M_i^2+q}.
$$
The sums $F_\ell$ are real in this texture. Contracting the two matrix orderings with $c_ic_j$ gives twice $F_\ell L_\ell$ because $H_{ij}=|y_\ell|^2 c_i^* c_j$. This product includes both diagonal and off-diagonal contributions. Combining the propagator-expansion sign, $p^2=-q$, and the tree convention gives the replacement
$$
 F_\ell(q)\ \longrightarrow\ F_\ell(q)+|y_\ell|^2 D_\ell(q;\mu),
 \qquad D_\ell(q;\mu)=-\frac{2q}{16\pi^2}\,F_\ell(q)L_\ell(q) f(q;\mu).
$$
Thus at zero transfer the entire correction from this kinetic insertion vanishes, even for a block with nonzero inverse-mass moment. The factor in front of $F_\ell$ is dimensionless. [Exchange signs and mass texture](../../checkpoint_2026_11_19/uv_completion/MATCHING.md)

## Angular projection and scale check

Let $\beta=\sqrt{1-4m^2/s}$, $q(x)=s/2-m^2-s\beta x/2$, and $\langle G\rangle=\int_{-1}^{1}G(q(x))\,dx/2$. The two identical-scalar assignments give equal means on each side. The $J=0$ average over the outgoing lepton direction factorizes the two production means. The ordered two-body measure, the identical-pair factor $1/2!$ and the discontinuity factor $1/2$ yield
$$
 \frac{A_{0,N}^{(1)}}{\gamma}=|y_e|^2 C_e+|g|^2 C_\mu,\qquad
 C_e=\frac{s\beta}{8\pi}\langle D_e\rangle\langle F_\mu\rangle,
 \qquad C_\mu=\frac{s\beta}{8\pi}\langle F_e\rangle\langle D_\mu\rangle .
$$

| $s/m^2$ | $C_e$ at $\mu=m$ | $C_\mu$ at $\mu=m$ | $\partial_{\ln\mu}C_e$ | $\partial_{\ln\mu}C_\mu$ |
|---:|---:|---:|---:|---:|
| 6 | $-2.775394959\,10^{-9}$ | $+3.445874617\,10^{-8}$ | $+6.029910349\,10^{-8}$ | $+9.468543378\,10^{-7}$ |
| 8 | $-2.615099745\,10^{-8}$ | $-3.041359246\,10^{-8}$ | $+1.776719649\,10^{-7}$ | $+1.930750408\,10^{-6}$ |
| 10 | $-8.040621480\,10^{-8}$ | $-1.649178225\,10^{-7}$ | $+3.549311204\,10^{-7}$ | $+3.008703911\,10^{-6}$ |

For the explicit slopes all renormalized inputs are fixed and $\partial_{\ln\mu} f=2$. The signs can be verified from the running of the tree. Assign to the sterile-field part of the one-loop flow $D_N Y=YH$ and $D_N M=H^T M+MH$, with $16\pi^2$ removed from both beta functions. In a general sterile basis the outgoing tree per assignment, with its common spinor sign stripped, is $K=Y(q+M^\dagger M)^{-1}M^\dagger Y^T$. Differentiating the inverse matrix and using the rank-one block gives
$$
 (D_N K)_{\ell\ell}=4q\,(y_\ell)^2|y_\ell|^2 F_\ell L_\ell .
$$
This cancels the explicit log derivative of the insertion after the loop factor is restored. The other terms in the Yukawa beta function belong to external active and scalar field renormalization. [Full-theory beta functions](../../checkpoint_2026_11_19/radiative_stability/RUNNING.md)

Run `python -I -B -O check_sterile.py --output STERILE.json` with mpmath. The checker integrates the original bubble, compares 40- and 64-node angular rules, and checks the tree mean against its elementary log integral. It independently perturbs the full chiral inverse kernel in a complex nondiagonal mass basis. At three transfers per flavor it checks the insertion against that inverse, including a light-row phase, and differentiates the mass and Yukawa matrices under $D_N$ to test the scale cancellation. The checks cover this single insertion class. [Checker](check_sterile.py)
