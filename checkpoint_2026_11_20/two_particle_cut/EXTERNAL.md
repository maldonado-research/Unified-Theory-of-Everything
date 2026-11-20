# External pole residues and scale cancellation on the two-scalar cut

At $\mu=m$, external-state residues multiply the leading two-scalar cut by $0.229679017971\,E-0.0170128006694\,G$, with $E=|y_e|^2$ and $G=|g|^2$. This factor is common to all energies and external angles. It drops out of the ratio of two cut amplitudes at first order. The sterile self-energy gives an energy-dependent explicit scale derivative. Its sum with the residues reproduces the independent four-lepton RG calculation at all three tested energies. [Numerical receipt](EXTERNAL.json) · [Four-lepton RG calculation](../../checkpoint_2026_11_19/yukawa_continuum/RG_LOG.md)

## Residues of the renormalized fields

Take the unbroken gaugeless Yukawa model with massless lepton doublets, $a=m^2$ the scalar pole mass squared and $b_i=M_i^2\geq a$ the squared positive Takagi masses at the reference scale. Use full-theory $\overline{\rm MS}$ for the Yukawa and sterile-mass parameters. Denote $H=Y^\dagger Y$, $P=YY^\dagger$, $T=\operatorname{tr}H$. The pole residues of the $\overline{\rm MS}$ fields are $1+r/(16\pi^2)$ at first order. An external particle with unit normalization supplies the square root of this pole residue. In each production amplitude the pair of leptons of flavor $f$ and the pair of charged scalars therefore multiply the tree by $1+(r_{Lf}+r_\phi)/(16\pi^2)$. [Original action and two-point conventions](../../checkpoint_2026_10_16/yukawa_counterterms/YUKAWA_TWO_POINT.md)

For an external lepton the inverse kinetic term receives a positive open-chain bubble. The numerator $P_R(\not k+M_i)P_L$ reduces to $\not kP_L$. Combining denominators with weight $x$ on the sterile line shifts this to $(1-x)\not pP_L$ under integration. The finite bubble is minus the logarithm of its squared-mass polynomial. Inverting the resulting massless kinetic term gives
$$
 r_{Lf}=\sum_i |Y_{fi}|^2\int_0^1 dx\,(1-x)
       \ln\frac{x b_i+(1-x)a}{\mu^2}.
$$
The benchmark has diagonal $P$. For overlapping rows the weight becomes $Y_{fi}Y_{gi}^*$ in the residue matrix. Each fixed external lepton component couples to one scalar component in this loop.

For the scalar, let $z=p^2$. Its inverse propagator has the sign convention $z-m_R^2+\Gamma_\phi(z)/(16\pi^2)$. The closed fermion contraction contributes a minus sign. Its numerator trace is $2k\cdot(k-p)=k^2+(k-p)^2-z$. Cancelling one denominator at a time, the massless tadpole vanishes and the subtracted inverse kernel is
$$
\begin{aligned}
 \Gamma_\phi(z)&=\sum_i H_{ii}\big[(z-b_i)B(z;b_i)-A(b_i)\big],\\
 A(b)&=b[1-\ln(b/\mu^2)],\qquad
 B(z;b)=-\int_0^1dx\,\ln\frac{x[b-(1-x)z]-i0}{\mu^2}.
\end{aligned}
$$
One fixed external scalar component selects one lepton component. Its conjugate Yukawa endpoints have fixed orientations and there is one closed contraction. Enforcing the pole condition subtracts $\Gamma_\phi(a)$ from the inverse kernel. The pole residue follows by inversion and differentiation,
$$
 r_\phi=-\Gamma'_\phi(a)=-\sum_i H_{ii}
       \left[B(a;b_i)+(a-b_i)B'(a;b_i)\right].
$$
With $r=b/a>1$, the per-sterile contribution divided by $H_{ii}$ is
$$
 \ln(b/\mu^2)-1-r-(r^2-1)\ln(1-1/r).
$$
At $b=a$ this expression tends to $\ln(a/\mu^2)-2$. The logarithmic branch point in $B'$ is suppressed by its factor $a-b$. The on-shell two-point derivative exists and is real at the benchmark threshold.

The following coefficients include $1/(16\pi^2)$ and are evaluated at $\mu=m$.

| Factor | $E$ | $G$ |
|---|---:|---:|
| Scalar pole residue | $+0.0954830001545$ | $-0.00913487789313$ |
| Electron pole residue | $+0.0387130176622$ | $0$ |
| Muon pole residue | $0$ | $+0.00125695511685$ |
| External factors on $ee\to\phi\phi$ | $+0.134196017817$ | $-0.00913487789313$ |
| External factors on $\mu\mu\to\phi\phi$ | $+0.0954830001545$ | $-0.00787792277628$ |

## Fixed-parameter logarithm and running of the leading cut

Differentiating the finite two-point integrals directly gives $\partial_{\ln\mu}r_{Lf}=-P_{ff}$ and $\partial_{\ln\mu}r_\phi=-2T$. The sterile inverse kinetic bubble has coefficient $f(q;\mu)=-2\int_0^1(1-x)\ln[(1-x)(a+xq)/\mu^2]\,dx$ for spacelike transfer $p^2=-q$. Its scale derivative is $2$. With the loop factor stripped, the insertion on the incoming tree per transfer has logarithmic derivative
$$
 [\partial_{\ln\mu}K_{\rm insertion}]_{ij}
 =-\frac{2q\,[M_iH_{ij}+H_{ji}M_j]}{(q+b_i)(q+b_j)}.
$$
Here $K$ is the sterile kernel between the incoming Yukawa rows $Y^*$ and $Y^\dagger$. The sign comes from the first-order expansion of the inverse propagator, $\delta S=-S\,\delta S^{-1}S$, with the same propagator sign stripped as the tree. The opposite-chirality proper vertex is UV convergent. The target-orientation triangle is forbidden by the complex scalar and massless lepton contractions. Those vertices therefore supply no explicit $\mu$ derivative. A one-loop proper four-point graph, when allowed, also has negative superficial degree and no proper UV subdivergence. [Sterile insertion and its sign](STERILE.md) · [Proper vertex contraction count](../../checkpoint_2026_10_16/yukawa_counterterms/YUKAWA_VERTEX.md)

One can independently vary the tree without reassigning any graph's sign. In fixed field coordinates write $F=Y^* K Y^\dagger$, $K=M(q+M^\dagger M)^{-1}$. At the diagonal positive reference mass and for general symmetric $\delta M$,
$$
 \delta K_{ij}=\frac{q\,\delta M_{ij}-M_iM_j\,\delta M_{ij}^*}
 {(q+b_i)(q+b_j)} .
$$
Use the one-loop flow with $16\pi^2$ stripped, $DY=\frac12PY+TY+YH$ and $DM=H^T M+MH$. The first two Yukawa terms give $\frac12(P^*F+FP)+2TF$, which cancels the explicit two-lepton and two-scalar residue derivatives as matrices. The two $YH$ terms at the endpoints give $H^TK+KH$ inside the sterile kernel. Adding the displayed $\delta K$ with $\delta M=DM$ cancels the mass-cubic terms and leaves $+2q(M_iH_{ij}+H_{ji}M_j)/[(q+b_i)(q+b_j)]$. This cancels the explicit kinetic-bubble derivative. Complex conjugation gives the outgoing flavor contraction. These identities preserve Majorana mixing and the phases of the squared Yukawa couplings. [Mass-matrix and Yukawa flow](../../checkpoint_2026_11_19/radiative_stability/RUNNING.md)

The identities hold at each transfer. They commute with adding the two boson assignments, multiplying the two trees with the outgoing one conjugated, and integrating the on-shell sphere. The measure depends only on $s$ and the fixed scalar pole mass. Consequently $\partial_{\ln\mu}A_2^{(1)}=-D A_2^{(0)}/(16\pi^2)$ for the Hermitian cut. The three-particle cut first uses two tree amplitudes at this order and has zero explicit scale derivative. All quoted cut coefficients divide by the *reference* coupling $\gamma=y_e^2(g^*)^2$ in the phase convention $B_{34}^*B_{12}=s$. The reference normalization is held fixed when differentiating the full tree product. [Sewing normalization](../../checkpoint_2026_11_19/uv_completion/CUT.md) · [Three-particle cut](../../checkpoint_2026_11_19/yukawa_continuum/FULL_THREE_BODY.md)

| $s/m^2$ | $A_{2,0}^{(0)}/\gamma$ | Explicit log per $E$, $J=0$ | Explicit log per $G$, $J=0$ | Explicit log per $E$, $c=0$ | Explicit log per $G$, $c=0$ |
|---:|---:|---:|---:|---:|---:|
| 6 | $-6.504208766\,10^{-5}$ | $6.238556580\,10^{-6}$ | $3.521128286\,10^{-6}$ | $6.240216451\,10^{-6}$ | $3.521982615\,10^{-6}$ |
| 8 | $-1.248688279\,10^{-4}$ | $1.203878831\,10^{-5}$ | $6.872882221\,10^{-6}$ | $1.204872527\,10^{-5}$ | $6.878462856\,10^{-6}$ |
| 10 | $-1.871551124\,10^{-4}$ | $1.813253504\,10^{-5}$ | $1.041603888\,10^{-5}$ | $1.816074926\,10^{-5}$ | $1.043271192\,10^{-5}$ |

For the quotient of the $J=0$ cut amplitudes $Q=A_{2,0}(10m^2)/A_{2,0}(6m^2)$, $Q^{(0)}=2.87744626858$. The explicit derivative from the two-particle diagrams is $-0.00278932760\,E-0.00436919293\,G$. Differentiating the leading quotient using the complete running matrices gives the opposite derivative. This quotient concerns angular means of discontinuities.

Run **python -I -B -O check_external.py --output EXTERNAL.json** in this directory. The check differentiates the scalar inverse kernel with a high-precision one-sided difference at the pole, independently integrates the on-shell lepton residue, and refines the angular quadrature. It compares the central-angle leading cut and both summed logarithms with the archived imaginary parts of an independently calculated upper-rim box and its RG derivative. Away from the benchmark it checks all complex-matrix entries with two overlapping lepton rows. At the benchmark it evolves the nondiagonal mass and complex Yukawa matrices directly, recomputes the projected leading cut and its quotient at four displacements, and compares against the explicit slopes for two separate nonzero strength pairs. [Checker](check_external.py) · [Receipt](EXTERNAL.json)
