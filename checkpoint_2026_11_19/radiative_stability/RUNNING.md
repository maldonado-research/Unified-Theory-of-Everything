# One-loop running of the two mass moments

In the four-singlet benchmark, the electron inverse-mass moment and the muon positive-mass moment both remain zero under one-loop pure-neutrino-Yukawa running with all four singlets active. The mass ratios and the Yukawa ratios each change. A rotation that keeps the Majorana mass matrix diagonal makes their compensation explicit. With successive tree-level matching, the dimension-five electron cancellation also survives the pure-Yukawa logarithm between the two electron masses. A scalar quartic produces a different result in that interval.

## With all singlets active

Use the original unbroken-theory, barred-composite convention. The symmetric, invertible Majorana mass matrix is $M_N$, and $Y$ has light flavor rows. Write $D=16\pi^2\,d/d\ln\mu$, $H=Y^\dagger Y$, $P=YY^\dagger$ and $T=\operatorname{tr}H$. In fixed field coordinates the one-loop minimal-subtraction equations restricted to neutrino Yukawas are
$$
 D Y=\tfrac32 Y H+TY,\qquad
 D M_N=H^T M_N+M_N H .
$$
These equations apply to the renormalizable theory above the largest mass. The earlier counterterm calculation derived both from the action. Antusch et al. give the equivalent equations with conjugated Yukawa and mass conventions, and Zhang gives the one-loop poles in the present convention. [Action and texture](../uv_completion/MATCHING.md) · [Counterterm derivation](../../checkpoint_2026_10_16/yukawa_counterterms/YUKAWA_VERTEX.md) · [Source convention conversion](../../checkpoint_2026_10_16/yukawa_counterterms/LITERATURE_YUKAWA_POLES.md) · [Antusch et al., Appendix D.1](https://arxiv.org/abs/hep-ph/0501272) · [Zhang, Appendix A](https://arxiv.org/html/2504.00792v2#A1)

The tensors in the outgoing inverse-mass and incoming positive-mass projections are
$$
 K=Y M_N^{-1}Y^T,\qquad R=Y^* M_N Y^\dagger .
$$
In positive Takagi coordinates $K_{ee}=\sum_iY_{ei}^2/M_i$ and $R_{\mu\mu}=\sum_iM_i(Y^*_{\mu i})^2$. These quantities have mass dimensions $-1$ and $+1$, respectively. Under $Y\mapsto U^\dagger Y V$, $M_N\mapsto V^T M_N V$ with constant unitary $U,V$, they transform as $K\mapsto U^\dagger K U^*$ and $R\mapsto U^T R U$. Thus both are independent of sterile coordinates, while the subscripts refer to the specified external light flavors.

Use $D M_N^{-1}=-M_N^{-1}(D M_N)M_N^{-1}$ in $D K$. Each external Yukawa contributes $3/2$ to its side of the $P$ action. Differentiating the inverse mass contributes $-1$ on each side. For $D R$, differentiating the numerator mass supplies $+1$ on each side. This gives the two tensor equations
$$
 \boxed{D K=2TK+\tfrac12(PK+KP^T)},\qquad
 \boxed{D R=2TR+\tfrac52(P^*R+RP)}.
$$
The conjugate of $R$ gives the positive moment for the opposite external orientation. The equations are homogeneous on whole tensors. For a single flavor entry, other entries can contribute if $P$ mixes the selected flavor. For example, with $M_N=\operatorname{diag}(1,4)$ in common mass units, choose electron and muon rows $(1,2i)$ and $(1,0)$. Then $K_{ee}=0$ but $D K_{ee}=1$. Substituting $(1,i/2)$ as the electron row gives $R_{ee}=0$ and $D R_{ee}=5$ in the corresponding units.

For the benchmark use $(M,2M,m,4m)$ and the two disjoint rows $(y,i\sqrt{2}\,y,0,0)$ and $(0,0,g,-ig/2)$. Initially $P=\operatorname{diag}(3|y|^2,5|g|^2/4)$ and $T=P_{ee}+P_{\mu\mu}$. The differential identity $D P=3P^2+2TP$ preserves diagonal $P$. In particular the two entries satisfy
$$
 D K_{ee}=(2T+P_{ee}) K_{ee},\qquad
 D R_{\mu\mu}=(2T+5P_{\mu\mu}) R_{\mu\mu}.
$$
Both initial conditions are zero, so both entries vanish throughout the common full-theory flow to this accuracy. The other two diagonal moments initially equal $K_{\mu\mu}=15g^2/(16m)$ and $R_{ee}=-3M(y^*)^2$. The surviving finite-momentum amplitudes of the benchmark are compatible with these componentwise zeros. [Tree kernels and cut weights](../uv_completion/MATCHING.md)

## Following the mass eigenstates

The initial mass beta function has nonzero off-diagonal entries within both pairs,
$$
 (D M_N)_{12}=-i\sqrt{2}\,M|y|^2,\qquad
 (D M_N)_{34}=\tfrac32 i\,m|g|^2 .
$$
For positive distinct masses, let $A$ be the anti-Hermitian generator of the infinitesimal Takagi rotation, so the Yukawa derivative in diagonal coordinates is $D Y+YA$ at this point. Writing $B=D M_N$ in the initial basis gives
$$
 B_{ij}+M_i A_{ij}-M_j A_{ij}^*=0,\qquad
 A_{ij}=-\frac{\operatorname{Re}B_{ij}}{M_i-M_j}
          -i\frac{\operatorname{Im}B_{ij}}{M_i+M_j}\quad (i\ne j).
$$
Here $B_{ii}=2M_i H_{ii}$ is real, so one can set $A_{ii}=0$. The diagonal derivatives give $D\ln(M_2/M_1)=2(H_{22}-H_{11})$ within a pair with local labels $1,2$.

| Pair | $D\ln(M_2/M_1)$ at the benchmark | $D\ln\lvert Y_2/Y_1\rvert^2$ in the running Takagi basis |
|---|---:|---:|
| Electron | $2\lvert y\rvert^2$ | $2\lvert y\rvert^2$ |
| Muon | $-3\lvert g\rvert^2/2$ | $3\lvert g\rvert^2/2$ |

For example, the electron column ratio tracks the mass ratio because its squared complex value is $Y_2^2/Y_1^2=-M_2/M_1$. In the muon block $(Y_2^*/Y_1^*)^2=-M_1/M_2$. Both relative phases stay fixed in this texture, up to discrete Takagi column signs. Evolving only the diagonal masses and using the fixed-coordinate Yukawa beta would miss the $Y A$ term.

## Between the electron thresholds

When the heavier electron state is removed, the zero-momentum dimension-five tensor in the intermediate theory is
$$
 K_{\rm tot}=C+K_a,\qquad K_a=Y_a M_a^{-1}Y_a^T .
$$
Here $C$ is the local barred-composite coefficient and $a$ denotes active singlets. Tree matching in the mass basis transfers the removed column's contribution from $K_a$ to $C$. At the upper electron threshold the two electron terms are individually nonzero, $C_{ee}=-y^2/M$ and $(K_a)_{ee}=y^2/M$, for the benchmark imposed there. The symbol $M$ denotes the lighter electron mass at that reference scale. Using running parameters evaluated at a different common scale affects the one-loop logarithm below only at higher order when the benchmark separations are used at leading order.

Zhang's intermediate-theory equations include a contribution $-3C Y_a^* M_a$ to $D Y_a$, in addition to the renormalizable terms. Write $P_a=Y_aY_a^\dagger$, $T_a=\operatorname{tr}P_a$ and $\mathcal L_{P_a}(X)=P_a X+X P_a^T$. Projecting the published one-loop equations onto the neutrino Yukawas and, optionally, the scalar quartic $V\supset\lambda(H^\dagger H)^2$ gives
$$
\begin{aligned}
 D C &= (2T_a+4\lambda)C+\tfrac72\mathcal L_{P_a}(C),\\
 D K_a &= 2T_a K_a+\tfrac12\mathcal L_{P_a}(K_a)-3\mathcal L_{P_a}(C),\\
 D K_{\rm tot} &= 2T_a K_{\rm tot}
                 +\tfrac12\mathcal L_{P_a}(K_{\rm tot})+4\lambda C .
\end{aligned}
$$
These are dimension-five equations to first order in the heavy-scale expansion for the removed modes. In particular, the $7/2$ term and the feedback through $D Y_a$ both enter the evolution of the sum. The published tensor equation applies directly in these conventions, independently of the paper's further determinant results. [Zhang, Eqs. (7)–(11)](https://arxiv.org/html/2405.18017v3#S3)

Setting aside terms proportional to $\lambda$, tree matching and the one-loop running of $K_{\rm tot}$ preserve its $ee$ zero across the electron thresholds at dimension five. The electron and muon supports stay disjoint because both $C$ and $P_a$ are diagonal and $M_a$ splits into sterile flavor blocks. While both muon singlets remain active and $C_{\mu\mu}=0$, the extra term in $D Y_a$ has zero muon row, so the muon positive moment still obeys its homogeneous equation. Removing a muon state changes the large-transfer expansion itself. Its local inverse-mass expansion applies below that mass and cannot be used to continue a two-pole large-transfer residue through the threshold.

If an independent quartic is included, running down from $2M$ to $M$ at fixed first loop order gives
$$
 \left.\Delta (K_{\rm tot})_{ee}\right|_{\lambda y^2,\ \mathrm{log}}
 =\frac{4\lambda y^2}{16\pi^2 M}\ln 2 .
$$
The sign follows from the negative upper-threshold $C_{ee}$ and the negative running interval. This coefficient is the dimension-five logarithm, with couplings evaluated at a common reference scale. In the renormalizable full theory at $\lambda=0$ with gauge and other Yukawa couplings switched off, $D\lambda=-2\operatorname{tr}[(Y^\dagger Y)^2]$. This follows from Eq. (8) with $C=0$ in the convention above. A quartic generated from a zero boundary value contributes through the displayed logarithm at one additional loop order. In an intermediate EFT, Eq. (8) also contains terms in $D\lambda$ involving $C$ and retained masses. Powers of light or retained masses over the removed mass belong to higher orders in the sequential expansion. Finite threshold constants require matching. The companion scalar-quartic calculation performs this matching for its specified interaction and reproduces the logarithm. A complete physical neutrino mass or all-orders symmetry would require further inputs. [Quartic matching and finite scalar-mass dependence](QUARTIC.md)

## Reproduction

Run `python -I -B -O check_running.py --output RUNNING.json` using the pinned repository requirements. The checker differentiates both moment definitions, changes both sterile and active coordinates, rediagonalizes the benchmark mass matrix at first order and separately combines the two intermediate-theory derivatives. All matrix checks use exact complex algebra. [Checker](check_running.py) · [Receipt](RUNNING.json)
