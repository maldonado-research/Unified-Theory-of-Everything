# The charged three-particle cut and its threshold coefficient

In the four-singlet benchmark the charged intermediate state $\phi^-\mu^-N_1$ opens at $s=4m^2$, the same energy as the scalar pair. Its complete contribution to the Hermitian cut at order $\gamma |g|^2$ starts as $-1.097860094\times10^{-7}(s/m^2-4)^{3/2}$ in the phase $B_{34}^* B_{12}=s$. Both the right-angle value and the angular mean have this leading coefficient. The coefficient follows from sewing the two scalar-decay trees. The other connected tree routings first enter at order $(s/m^2-4)^{5/2}$. The two-scalar cut also contributes to the full two-loop amplitude. [Benchmark action and mass residues](../uv_completion/MATCHING.md) · [Spin and unitarity conventions](../uv_completion/CUT.md)

## The three-particle measure

Let $a=m^2$, $B=M_i^2$, $U=(p_L+p_N)^2$ and $U_-=(\sqrt{s}-\sqrt a)^2$. For one fixed charged scalar component the spin sum of its chiral decay tree is $|Y_{\ell i}|^2\,2p_L\!\cdot p_N=|Y_{\ell i}|^2(U-B)$. The decay into a definite charged lepton counts a single orientation of the Majorana line. The integrated $L N$ phase space is $(U-B)/(8\pi U)$. Factorizing $d\Phi_3=dU\,d\Phi_2(P;k,Q)\,d\Phi_2(Q;p_L,p_N)/(2\pi)$ and including the scalar propagator on each side gives the weight per unit squared Yukawa coupling
$$
\rho_B(U)=\frac{(U-B)^2\,\theta(U-B)}{16\pi^2 U(U-a)^2}.
$$
For $B=a$ the ratio of the two squared factors is one on the continuum. This leaves a finite density at its lower endpoint. The pole-subtracted insertion keeps the scalar pole at its specified mass with unit residue. Its three-particle cut has the same weight as this direct physical sewing.

Write $h=(s-a-U)/2$, $v=\sqrt{h^2-aU}$ and $\beta_U=2v/s$. The two exchanged momentum transfers on a given side are $-h\pm vx$ even when the scalar masses differ. Energy conservation cancels their mass difference in these combinations. Here $x$ is the cosine between the incident axis and the pole scalar, $y$ the corresponding cosine for the outgoing axis, and $c$ the cosine between the two axes. In mass units $m=1$, the benchmark kernels and tree sums are
$$
 F_\mu(q)=\frac{15}{(1+q)(16+q)},\quad
 F_e(q)=\frac{-30q}{(100+q)(400+q)},\quad
 S_\ell(x)=F_\ell(h-vx)+F_\ell(h+vx).
$$
Restore units by expressing $s,a,U,B,q$ in $m^2$ and $F$ in $m^{-1}$. Define $A_B=(\mathcal M_{fi}-\mathcal M_{if}^*)/(2i\gamma |Y_{\ell i}|^2)$ for the scalar self-energy insertion from this column only. Its signed discontinuity is
$$
 A_B(s,c)=\frac{s}{16\pi}\int_B^{U_-}dU\,\rho_B(U)\,\beta_U
 \int\!\frac{d\Omega_k}{4\pi}\,S_\mu(x)S_e(y),\qquad U_->B,\quad \sqrt{s}>\sqrt a+\sqrt B .
$$
It vanishes below the physical threshold. The three final particles are distinguishable. The $1/2$ from unitarity, combined with the angle-integrated $\beta_U/(8\pi)$, gives $1/(16\pi)$ here. Relative to the identical scalar-pair cut this doubles the prefactor. In an insertion calculation the same factor counts the two scalar lines that can carry the self-energy. Both tree assignments have already entered each $S_\ell$. The spin factor $-W=s$ fixes the sign, and $h^2-v^2=aU>0$ gives $F_\mu>0$, $F_e<0$ throughout the physical angular integral. The reduced discontinuity is negative wherever the channel is open.

## Projection and numerical values

At fixed $U$ and scalar direction, averaging the outgoing axis makes $y$ uniform on $[-1,1]$. The remaining sphere integral averages the incoming tree. Thus the angular factor in the $J=0$ cut is $4\bar F_\mu\bar F_e$, where $\bar F_\ell=\sum_i w_{\ell i}J_{b_i}$ and
$$
 J_b=\frac{1}{2v}\log\frac{b+h+v}{b+h-v},\qquad J_b(v=0)=\frac{1}{b+h}.
$$
The two means share the integration variable $U$. At $c=0$ the azimuthal average of a single electron denominator is $[(b+h)^2-v^2(1-x^2)]^{-1/2}$. These expressions give the following values for **the scalar self-energy insertions**, with $M/m=10$ and $\gamma |g|^2$ stripped.

| $s/m^2$ | $A_{N_1}(s,0)$ | $\langle A_{N_1}\rangle_c$ |
|---:|---:|---:|
| 6 | $-3.997025444\times10^{-7}$ | $-3.996201790\times10^{-7}$ |
| 8 | $-1.264685994\times10^{-6}$ | $-1.263796084\times10^{-6}$ |
| 10 | $-2.462922767\times10^{-6}$ | $-2.459469805\times10^{-6}$ |

The open column is $N_1$ at each energy. The thresholds in $s/m^2$ are $4,25,121,441$ for $N_1,N_2,E_1,E_2$ respectively. Their squared coupling weights are $|g|^2,|g|^2/4,|y|^2,2|y|^2$. The heavy columns can contribute to the dispersive amplitude below their cut thresholds.

A second integration takes the pole scalar's recoil momentum $p$ as the independent variable. It uses $d^3k/[(2\pi)^3 2E_\phi]$, the integrated decay phase space, the spin trace and both tree transfers computed from the center-of-mass energies and momenta. In this form the radial measure for the cut, including the external spin and the two propagators, is $s p^2 dp\,[(U-B)/(U-a)]^2/(64\pi^3 E_\phi U)$. Numerical quadrature over the production sphere and external angle reproduces both columns above. The program also tests the nondegenerate $N_2$ cut above threshold. [Executable](check_threebody.py) · [Numerical receipt](THREE_BODY.json)

## Threshold laws and the other tree routings

Set $m=1$, $r=\sqrt B$, $s_0=(1+r)^2$, $\delta=s-s_0$, $t=U-B$ and $\kappa=r/(1+r)$. The integration interval has length $\kappa\delta+O(\delta^2)$. To leading order $s\beta_U=2\sqrt r\sqrt{\delta-t/\kappa}$ and both scalar directions yield the common tree values $2F_\ell(r)$. One may take $t=(U_--B)z$, $0\le z\le 1$, to expand uniformly under the integral. The exact phase-space square root factors as $\sqrt{(U_--U)((\sqrt{s}+1)^2-U)}$.

For $B=1$, the constant endpoint density gives
$$
 A_1(s,c)=\frac{F_\mu(1)F_e(1)}{96\pi^3}\,\delta^{3/2}+O(\delta^{5/2}),\qquad
 F_\mu(1)=\frac{15}{34},\quad F_e(1)=-\frac{30}{101\cdot401}.
$$
For $B>1$, each isolated insertion density vanishes quadratically. Integrating $z^2\sqrt{1-z}$ gives $16/105$. Its leading term per unit column strength is
$$
 A_B(s,c)=\frac{F_\mu(r)F_e(r)\sqrt r\,\kappa^3}{32\pi^3 B(B-1)^2}\,
 \frac{16}{105}\,\delta^{7/2}+O(\delta^{9/2}).
$$
These higher thresholds are fixed-order cuts with free sterile lines. Heavy sterile decay and finite-width effects require their own treatment. Lower channels are already open at these energies.

To identify which coefficient survives in the complete $\phi^-\mu^-N_1$ state sum, label the three Yukawa vertices in each tree by the active particle attached to them. The two incoming leptons and the outgoing muon exhaust the active slots. Electric charge attaches the outgoing scalar to either incoming vertex. When the external $N_1$ attaches to the outgoing muon, the other two sterile ends join with a mass numerator. The two choices of scalar assignment are exactly the decay trees above. In an electron-initiated tree the flavor-diagonal blocks force this attachment of $N_1$. For a muon-initiated tree the external $N_1$ can instead attach to either incoming muon. Connectivity then fixes the outgoing scalar to the other incoming vertex. These two extra routings have an internal scalar of transfer $(p_i-k_N)^2$ and an internal sterile of transfer $(p_j-k_\phi)^2$ with a momentum numerator.

As $s\downarrow4m^2$, the two massive outgoing particles approach rest, so both extra transfers approach $-m^2$ uniformly in angle. In the rest frame of the $L N_1$ pair the lepton has energy $\omega=(U-a)/(2\sqrt U)$. Its external spinor is $O(\sqrt\omega)$, and the boost to the overall center of mass stays bounded. The decay tree is $O(\omega^{-1/2})$ from its scalar denominator. Each extra tree is $O(\sqrt\omega)$ because both propagators remain uniformly bounded. The interference between an extra muon tree and the electron decay tree is therefore bounded even as $U\downarrow a$. Its phase-space integral carries $\int_0^{O(\delta)}dt\,t\sqrt{U_--a-t}=O(\delta^{5/2})$. This proves the coefficient of $\delta^{3/2}$ above for the **complete contribution of this physical three-particle state**. The tabulated finite-energy values apply to the scalar insertion alone. The [complete three-particle sewing](FULL_THREE_BODY.md) gives the finite-energy exchange contributions. A full four-lepton rate additionally requires the other Yukawa diagrams, their renormalization and interference with the leading complex amplitude.
