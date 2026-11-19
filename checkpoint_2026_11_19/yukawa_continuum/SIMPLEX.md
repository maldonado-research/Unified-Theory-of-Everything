# The finite two-loop scalar propagator insertion

Fix the scalar pole at $a=m^2$ and give the charged scalar unit pole residue. Evaluate the original four-lepton box with either scalar propagator replaced by the continuum from a massless lepton and a Majorana fermion. This determines its entire complex contribution to $\mu^-_L\mu^-_L\to e^-_Le^-_L$ at the order of six Yukawa vertices. The numerical result here is a contribution in the stated subtraction prescription. In the original model additional diagrams and the corresponding counterterms are required at that order. [Scalar normalization and spectral density](SPECTRUM.md) · [Original four-singlet action](../uv_completion/MATCHING.md)

## From the dressed propagator to a two-dimensional integral

Use units $m=1$ throughout this computation, with $b=(1,16)$ labeling the muon-block squared masses, $d=(100,400)$ the electron-block squared masses, and $B$ the squared mass in the scalar self-energy. The mass residues are $w_b=(1,-1)$ and $v_d=(10,-40)$. The weights of the **squared magnitudes** of the Yukawas inside the scalar self-energy are $1,1/4$ in the muon block and $1,2$ in the electron block, after their base Yukawas have been removed. For a unit squared coupling the continuum weight in $dU$ is $R_B(U)/(16\pi^2)$ with
$$
 R_B(U)=\frac{(U-B)^2}{U(U-1)^2}\,\mathbf1_{U\ge B}.
$$
For $B=1$ the endpoint factors cancel, leaving $1/U$ above the endpoint. There is no continuum weight from a state below its own threshold.

For one scalar site, one mass assignment and one external scalar pairing, the original scalar box has parameters $(x_1,x_2,x_3,x_4)$ for internal squared masses $(U,1,b,d)$. They are nonnegative and sum to one. Its denominator before squaring is $Ux_1+x_2+b x_3+d x_4-s x_1x_2-t x_3x_4-i0$. Set $x=x_2$, $y=x_3$, $L=1-x$, $z=L-y$, and integrate $0\le x_1\le z$. The denominator is affine in $x_1$. Its real endpoint values are
$$
 A_d=x+b y+(d-ty)z,\qquad
 z(U-P),\qquad P=sx-\frac{x+b y}{z}.
$$
Therefore the integral over $x_1$ is $1/[A_d(U-P-i0)]$ on the integration interior. This identity first holds at $s<0$ and extends by analytic continuation. For $t\le0$, $A_d$ is positive there. Folding in $R_B$ gives the dimensionless box integral
$$
 I_{Bbd}(s,t)=\int_0^1 dx\int_0^{1-x}dy\,\frac{J_B(P)}{A_d},\qquad
 J_B(P)=\int_B^\infty\frac{R_B(U)\,dU}{U-P-i0}.
$$
The imaginary part of $J_B$ is $+\pi R_B(P)$ inside the cut. For $B=1$ the full result is $J_1(P)=-\ln(1-P-i0)/P$, with $J_1(0)=1$. All other $B$ values are above the accessible range of $P$ for $s=6,8,10$. In fact, maximizing $P$ first over $y$ and then over $x$ gives $(\sqrt{s}-1)^2$ at $y=0$. A useful real form follows from substituting $v=B/U$,
$$
 J_B(P)=\frac1B\int_0^1\frac{(1-v)^2\,dv}{(1-v/B)^2(1-Pv/B)} .
$$
The checker evaluates this smooth integral or its partial fractions for large negative $P$. For the lightest state it splits the remaining simplex at $P=1$ and resolves the logarithmic boundary. Before quadrature, it combines the two signed heavy residues exactly,
$$
 \frac{10}{A_{100}}-\frac{40}{A_{400}}=
 \frac{-30(x+b y-tyz)}{A_{100} A_{400}} .
$$
This keeps a cancellation in the electron source explicit throughout the integral.

## Amplitude, rate and independent checks

Let $c$ denote the scattering cosine, $t=-s(1-c)/2$, $u=-s(1+c)/2$, $\gamma=y_e^2(y_\mu^*)^2$, and use the original phase $B_{34}^*B_{12}=s$. Define
$$
 D_B(s,c)=\frac{2s}{(16\pi^2)^2}\sum_{b,d} w_b v_d
      [I_{Bbd}(s,t)+I_{Bbd}(s,u)] .
$$
The factor two is exact for these external states. Exchanging $x_1$ and $x_2$ in the original simplex exchanges the two scalar squared masses while preserving both transfer terms. The separate pairings in $t,u$ already account for identical external electrons. The scalar-line amplitude is
$$
 \frac{\Delta\mathcal M_{\rm scalar}}{\gamma}
  =|y_\mu|^2 D_g+|y_e|^2 D_y,\qquad
 D_g=D_1+D_{16}/4,\quad D_y=D_{100}+2D_{400}.
$$
Keeping the interference with $A_0=\mathcal M_{\rm leading}/\gamma$, the coefficients of $|y_\mu|^2,|y_e|^2$ in the scaled rate $r=m^2\sigma/|\gamma|^2$ are $2\langle\operatorname{Re}(A_0^*D_{g,y})\rangle_c/(32\pi s)$ in these mass units. The angular average includes all partial waves, and the flux includes the identical-final-electron factor. The table gives these coefficients for the scalar-line insertion alone.

| $s/m^2$ | Leading $r_0$ | Coefficient of $\lvert y_\mu\rvert^2$ | Coefficient of $\lvert y_e\rvert^2$ |
|---:|---:|---:|---:|
| 6 | $4.00925370\times10^{-11}$ | $1.89542094\times10^{-12}$ | $3.85910509\times10^{-14}$ |
| 8 | $5.61821543\times10^{-11}$ | $2.78740068\times10^{-12}$ | $4.84448607\times10^{-14}$ |
| 10 | $7.17400709\times10^{-11}$ | $3.69882542\times10^{-12}$ | $5.58343940\times10^{-14}$ |

The corresponding contribution to $\Delta\ln[\sigma(10m^2)/\sigma(6m^2)]$ is $+0.0042825552\,|y_\mu|^2-0.0001842621\,|y_e|^2$. This small difference comes from ratios of energy-dependent corrections to the leading rates. The electron continuum is kinematically closed here, yet its dispersive part contributes to the interference. Other Yukawa corrections can modify both the size and sign of the total rate shift.

An independent evaluation keeps the continuum mass as an integration variable and calls the published QCDLoop massive box in binary128. It checks the large-mass tail separately. The Python checker eliminates the continuum mass first and compares the central amplitude, angular mean and integrated interference with the stored external data. Three-particle phase space gives an independent check on the imaginary part. [Python program](check_continuum.py) · [Numerical receipt](CONTINUUM.json) · [QCDLoop implementation and grid refinement](EXTERNAL.md) · [Three-body cut](THREE_BODY.md)
