# Scalar-quartic correction to the electron vertex

For the electron singlets with masses $M,2M$ and Yukawa couplings $y,i\sqrt{2}\,y$, the one-loop scalar quartic generates a nonzero zero-momentum coefficient even though the two tree terms cancel. With $V=aH^\dagger H+\lambda(H^\dagger H)^2$ and $M=10\sqrt a$, its value is $0.0168504162035541\,\lambda y^2/M$. Unequal singlet masses give the two tree terms different loop weights. The leading dimension-five result in the heavy-mass expansion is $2\ln(4)\,\lambda y^2/(16\pi^2 M)$. [Numerical receipt](QUARTIC.json)

## Vertex and counting

Work at zero background in the unbroken gaugeless model, with positive renormalized $a=m^2$. The selected coupling order is $\lambda y^2$. The muon columns have no electron Yukawa and do not enter this order. Use the preceding checkpoint's metric $(+---)$, $S=1+iT$ and literal action term $+\tfrac12\bar X_e C X_e^c+\mathrm{h.c.}$ with $X_e=\widetilde H^\dagger L_e$. Its outgoing-electron projection for equal charged weak components is $-2B_{34}^* C$, where $B_{34}$ is the ordered antisymmetric spinor contraction. Define $C^{(0)}$ by projecting the amputated open spinor tensor and then setting all external momenta to zero. This avoids treating a massless zero-momentum spinor as an external state. Adjoint fields give the conjugate coupling $y^{*2}$. [Tree action and spin conventions](../uv_completion/MATCHING.md)

The charged projection selects a vertex $-4i\lambda$ from $-\lambda(h^\dagger h)^2$. The doublet vertex is $-2i\lambda(\delta_{ac}\delta_{bd}+\delta_{ad}\delta_{bc})$ for two $H$ indices $a,b$ and two conjugate indices $c,d$. Both charge routings end on the same charged component. To count directly in the full theory, label the two Yukawa vertices temporarily. Two assignments of the external fermions cancel the $1/2!$ from expanding the Yukawa interactions. Their Wick exchange sign cancels the sign from the antisymmetric spin tensor. At tree level there are two assignments of the external scalar legs. With a quartic vertex there are two assignments of the internal scalar lines to its conjugate slots and two assignments of its remaining slots to the external scalars. Dividing four by the two already present in the normalized tree gives a relative factor two. Equivalently the two-scalar bubble joining an already symmetrized tree to the quartic has symmetry factor $1/2$. The singlet forms an open fermion line.

For each mass the analytic tree kernel at spacelike transfer is $Y_{ei}^2 M_i/(M_i^2+k_E^2)$. Multiplying its Minkowski vertex by $(-4i\lambda)/2$, the two scalar propagator numerators $i^2$ and the Wick rotation $i$ gives a relative $-2\lambda$. Thus in the stated action normalization the triangle contributes
$$
\delta C^{(0)}=-2\lambda\sum_i Y_{ei}^2
\int\frac{d^4k_E}{(2\pi)^4}\,
\frac{M_i}{(k_E^2+a)^2(k_E^2+M_i^2)} .
$$
The remaining quartic contractions either put a scalar mass tadpole on an external leg or give a disconnected vacuum contribution. Amputation with the full propagator removes the external insertions. Equivalently, a scalar on-shell mass counterterm cancels the momentum-independent one-loop quartic self-energy. Its derivative vanishes, and neither the external lepton wavefunction nor the singlet masses or Yukawa vertices receive a one-loop pure-$\lambda$ renormalization. Replacing $a$ by its one-loop mass shift in the displayed triangle first changes it at the next loop order.

## Exact scalar mass and sign

Set $r_i=a/M_i^2$. Combining the two light denominators with the heavy one gives
$$
16\pi^2 M_i\int\frac{d^4k_E}{(2\pi)^4}
\frac{M_i}{(k_E^2+a)^2(k_E^2+M_i^2)}
=f(r_i)=\int_0^1\frac{x\,dx}{1-x+r_i x}
=\frac{r_i-1-\ln r_i}{(1-r_i)^2},\qquad f(1)=\tfrac12 .
$$
For a mass ratio $R=M_2/M_1>1$ write $M_1=M$, $r=a/M^2$, $K=R^2$ and impose $Y_{e2}^2=-R y^2$ so the tree inverse-mass moment vanishes. Then
$$
\delta C^{(0)}=\frac{2\lambda y^2}{16\pi^2 M}\,[f(r/K)-f(r)] .
$$
As a separate calculation, set $t=k_E^2/M^2$ directly in the pair's subtracted radial integral. The bracket equals
$$
\int_0^\infty\frac{(K-1)t^2\,dt}{(t+r)^2(1+t)(K+t)} .
$$
The integrand is positive when $K>1$ and $r>0$. Its $r=0$ limit integrates to $\ln K$ and gives a strict upper bound at positive $r$. For equal singlet masses the tuned squares are opposite and the entire integrand vanishes at every momentum. Equality of a scalar and a singlet mass is regular and gives the single-mass value $f(1)=1/2$. For $r\ll1$, $f(r)=\ln(1/r)-1+r[2\ln(1/r)-1]+O(r^2\ln(1/r))$. Individual terms have a light-mass logarithm while their tuned difference tends to $\ln K$. At $R=2$, $M/m=10$, the difference is $1.33045553538228$, or $95.9721\%$ of $\ln4$. [Analytic and quadrature checks](check_quartic.py) · [Results](QUARTIC.json)

## Dimension-five matching and running

At fixed positive $a$, the exact zero-momentum vertex includes fluctuations of the propagating light scalar and powers of $a/M_i^2$. Define the leading dimension-five Wilson coefficient by matching in the original fields, expanding to leading inverse mass and subtracting the light EFT amplitude in $\overline{\mathrm{MS}}$. For $c_i=Y_{ei}^2/M_i$, the finite single-pole separation is
$$
\begin{aligned}
 \delta C^{(0)}_{i,\mathrm{leading}}
 &=\frac{2\lambda c_i}{16\pi^2}\,[1-\ln(M_i^2/a)],\\
 \delta C_{5,i}^{\mathrm{hard}}(\mu)
 &=\frac{2\lambda c_i}{16\pi^2}\,[1+\ln(\mu^2/M_i^2)],\\
 \delta C^{\mathrm{light}}_{i}(\mu)
 &=-\frac{2\lambda c_i}{16\pi^2}\ln(\mu^2/a).
\end{aligned}
$$
The light term follows from the renormalized zero-momentum bubble with a local tree insertion $c_i$. The last two lines sum to the first. For the tuned pair $\sum c_i=0$, the constants and common scale logs cancel. The surviving hard difference is $2\lambda y^2\ln K/(16\pi^2M)$. These signs and the single-pole finite term agree with the published type-I hard matching at the selected coupling order. [Zhang & Zhou, Eqs. (47), (48), (88)](https://arxiv.org/pdf/2107.12133v3)

The coefficient below one threshold obeys $16\pi^2\,d C_5/d\ln\mu=4\lambda C_5$ in this quartic convention. Sequential matching first adds $c_2=-y^2/M$ at $M_2$, so running down to $M_1$ gives a positive increment $4\lambda(y^2/M)\ln R/(16\pi^2)$. Adding $c_1$ cancels the tree term, and the two one-loop constants at their own thresholds also cancel. This reproduces the summed hard result. The older sequential-threshold calculation uses $-\lambda_{\mathrm A}(H^\dagger H)^2/4$ in its Lagrangian, giving $\lambda_{\mathrm A}=4\lambda$ here. Both the effective bubble and its finite singlet counterpart were previously identified. [Antusch et al., Fig. 3 and Eq. (D.33)](https://arxiv.org/pdf/hep-ph/0501272v3)

The correction diagnoses the electron tree tuning for the chosen renormalized parameters and $\lambda\ne0$. Its exact positive-mass value is an off-shell four-field coefficient. A light-neutrino pole mass or broken-phase rate requires additional physical inputs. The calculation selects one independent coupling monomial.

## Normalization check in the four-lepton channel

There is a further check on the quartic factor using the polarized box process. Factor out the common couplings $\gamma=y^2(g^*)^2$ and choose the spin phase $B_{34}^*B_{12}=s$. For each flavor let $S_\alpha(s)=\sum_i w_{\alpha i}T_{M_i^2}(s)$, with real residues $w_{\alpha i}$ after stripping its coupling and
$$
 T_b(s)=\int_{x,y\ge0,\ x+y\le1}\frac{dx\,dy}{a(x+y)+b(1-x-y)-sxy-i0}.
$$
The charged corner normalized to its ordered spin and two scalar assignments receives $-2\lambda S_\alpha(s)/(16\pi^2)$. At zero total momentum this triangle agrees with the integral above after restoring the mass numerator. At order $\lambda y^2(g^*)^2$ all four Yukawa active legs are external. Connecting all four scalar slots to the quartic produces two triangles, one per flavor. A single self-contraction at the quartic inserts a momentum-independent tadpole in either scalar line of a box, canceled in the scalar pole-mass convention. Two self-contractions give a disconnected vacuum graph. In the pole-mass convention the quartic term is therefore
$$
 \delta\mathcal M(s)=-\frac{4\lambda s\gamma}{(16\pi^2)^2}\,S_e(s)S_\mu(s).
$$
The central quartic supplies four scalar contractions. Both triangle integrals contain one ordered internal scalar routing. The product is independent of the external scattering angle.

To check the sign, write $\bar F_\alpha=\sum_i w_{\alpha i}\langle[b+s/2-a-(s\beta/2)x]^{-1}\rangle_x$ for $b=M_i^2$, uniform $x\in[-1,1]$, and $\beta=\sqrt{1-4a/s}$. The upper-rim two-scalar cut gives $\operatorname{Im}S_\alpha=\pi\beta\bar F_\alpha$. Replace the normalized tree in turn on each side of the physical two-scalar cut by the one-loop charged corner. The imaginary parts of those corners cancel in the interference between the two cut orderings. For real reduced residues, the result is
$$
\operatorname{Im}(\delta\mathcal M/\gamma)=
-\frac{2\lambda s\beta}{8\pi(16\pi^2)}
\,[\bar F_\mu\operatorname{Re}S_e+\bar F_e\operatorname{Re}S_\mu].
$$
It equals the discontinuity divided by $2i$ of the displayed product. The ordered phase space, identical intermediate scalar factor and tree signs are those of the preceding cut calculation. [Four-singlet cut normalization](../uv_completion/CUT.md)
