# Tree matching for four Majorana singlets

Two unequal positive-mass singlets coupled to the electron can cancel the tree-level Weinberg coefficient while producing a nonzero two-derivative amplitude. A second pair coupled to the muon realizes the opposite *mass-weighted* residues of the earlier cut. The two cancellations involve different moments of the masses.

## Model and external states

Use the unbroken gaugeless theory with massless external left-handed leptons, $\phi=\widetilde H$, $a=m_\phi^2=m^2>0$, and $X_\alpha=\phi^\dagger L_\alpha$. Choose $M\gg m$, with $m,4m$ retained in propagators. For the singlets in the column order $(E_1,E_2,N_1,N_2)$, the positive Takagi masses and Yukawa rows are
$$
 (M_i)=(M,2M,m,4m),\qquad
 Y=\begin{pmatrix}
 y & i\sqrt{2}\,y & 0 & 0\\
 0 & 0 & g & -ig/2
 \end{pmatrix}.
$$
Here $y,g$ may be complex. The sign of an individual Takagi column can be reversed. Canonical light kinetic terms, a scalar potential with positive quadratic mass, and the following singlet terms specify a renormalizable model. Select the term of order $y^2(g^*)^2$ in the direct charged process $\mu_L^-\mu_L^-\to e_L^-e_L^-$. Couplings in other sectors are switched off for the calculation.
$$
 \mathcal L_N=\tfrac12\sum_i\bar N_i(i\not\partial-M_i)N_i
       -\big(\bar X_\alpha Y_{\alpha i}N_{Ri}+\mathrm{h.c.}\big),
 \qquad N_i=N_i^c .
$$
The gauge-singlet composite permits an ordinary derivative acting on the whole $X^c$. The literal EFT action in these original fields uses
$$
 \tfrac12\bar X_e\,[C_5+C_7(-\Box)+C_9(-\Box)^2+\cdots]X_e^c+\mathrm{h.c.}
$$
The signs refer to $\eta=(+---)$ and $S=1+iT$. The ordered spin determinant is $B_{ij}=u_i^T\epsilon u_j$, $\epsilon_{12}=1$. [Literal singlet and seed action](../../checkpoint_2026_10_15/operator_map/ACTION_SIGN.md) · [Dimension-five sign and charged tree](../energy_scan/LOWER_OPERATOR.md)

## Matching a single transfer

Write $z=q^2$, where $q=p-k$ is the difference of a lepton momentum and a scalar momentum. The incoming same-charge pair has a sterile tree $+B_{12}\,M_i(Y^*_{\mu i})^2/(z-M_i^2)$ per scalar assignment, including the open-fermion ordering sign. Reversing this tree for the outgoing electron pair on a spacelike transfer gives $+B_{34}^*\,M_iY_{ei}^2/(z-M_i^2)$. The local electron tree per assignment is $-B_{34}^*(C_5+C_7z+C_9z^2+\cdots)$. Thus its exact reduced kernel is
$$
 H_e(z)=\sum_{i\in E}\frac{M_iY_{ei}^2}{M_i^2-z}
       =\sum_{n\ge0} z^n\sum_{i\in E}\frac{Y_{ei}^2}{M_i^{2n+1}},
 \qquad |z|<M^2 \quad\text{for the series.}
$$
Equivalently, eliminating the heavy fields gives a minus one-half times the source contracted with the inverse free Majorana kernel. The equal-chirality part of $(\not q-M_i)^{-1}$ has numerator $M_i$ and denominator $q^2-M_i^2$. The two source signs multiply to plus. The resulting one-half times $H_e$ multiplies the barred composite bilinear, with the conjugate coefficient in its adjoint. The opposite-chirality propagator also generates lepton-number-conserving interactions starting at dimension six. They have different external field content from this two-electron tree. The half in the action cancels the two equivalent external fermion contractions. Both scalar assignments survive with the same sign, so the normalized electron factor is their **mean** if the overall factor $-2B_{34}^*$ is stripped. [Ordered fermion contractions and tree normalization](../on_shell/INDEPENDENT_CUT.md)

For the specified electron columns this gives
$$
 H_e(z)=y^2\!\left(\frac{M}{M^2-z}-\frac{4M}{4M^2-z}\right)
       =\frac{C_7 z}{(1-z/M^2)(1-z/(4M^2))}.
$$
In particular,
$$
 C_5=0,\qquad C_7=\frac{3y^2}{4M^3},\qquad C_9=\frac{15y^2}{16M^5},
 \qquad C_{5+2n}=\frac{y^2}{M^{2n+1}}(1-4^{-n}).
$$
The literal $C_7$ has dimension $-3$. The muon residues are $m_i(Y^*_{\mu i})^2=(w_0,-w_0)$, $w_0=m(g^*)^2$. Their sum vanishes, giving the zero derivative-triangle UV residue of the EFT. The retained muon exchange has a nonzero inverse mass moment $\sum (Y^*_{\mu i})^2/m_i=15(g^*)^2/(16m)$ and must keep both denominators at scan energies. In the inverse process the coupling product and spin factors are conjugated. Causal amplitudes in both directions use the upper rim. [Two-mass cut and phase conventions](../on_shell/ON_SHELL_TRIANGLE.md)

## What the on-shell slope measures

For equal scalar and charged weak components put $t=(p_3-k_1)^2$, $u=(p_3-k_2)^2$, so $t+u=2a-s$. The outgoing UV tree equals $-B_{34}^*[H_e(t)+H_e(u)]$. Through dimension seven its mean after stripping $-2B_{34}^*$ is $F(s)=C_5+C_7(a-s/2)$. The coefficient of $s$ is **$-3y^2/(8M^3)$**. The apparent intercept $aC_7$ at $s=0$ coexists with a vanishing single-transfer constant $C_5$. In this model, eliminating the actual heavy fields fixes the representative with the derivative on the whole composite. Measuring this charged on-shell slope alone fixes the sum of contributions in that projection. For example a separately matched term $D\,k_1\cdot k_2$ in the normalized tree would give a slope $-(C_7-D)/2$. [General on-shell operator projection](../energy_scan/CONTACT_SCOPE.md)

With $x=-z/M^2\ge0$, the exact per-transfer multiplier relative to dimension seven is
$$
 f(x)=\frac4{4+5x+x^2},\quad
 1-f(x)=\frac{x(5+x)}{4+5x+x^2}\le \frac54 x,\quad
 f(x)-(1-5x/4)=\frac{x^2(21+5x)}{4(4+5x+x^2)}\le\frac{21}{16}x^2 .
$$
The second identity isolates the remainder after dimension nine. The dimension-nine contribution to the *mean* is $C_9(t^2+u^2)/2=C_9[(s/2-a)^2+(s\beta/2)^2\cos^2\chi]$, with $\beta=\sqrt{1-4a/s}$ and $\chi$ the angle between $p_3$ and $k_1$ in the center-of-mass frame. Thus higher orders introduce angular dependence. For all cut directions $x\le x_{\max}=(s/2-a+s\beta/2)/M^2$. The Taylor series converges uniformly for $x_{\max}<1$ and a small $x_{\max}$ controls the omission. The rational inequalities themselves apply to every spacelike transfer. Both electron routes carry a common phase and the same spacelike sign, so their mean at dimension seven is multiplied by a weighted average of $f(-t/M^2)$ and $f(-u/M^2)$ with positive weights proportional to $-t,-u$. This supplies a uniform tree error estimate before angular sewing.

## Stability of the two cancellations

At fixed masses, changing $Y_{e2}^2=-2y^2(1+\delta)$ gives $C_5=-\delta y^2/M$ and $C_7=(3-\delta)y^2/(4M^3)$. A small fractional mismatch can dominate the derivative contact at $|z|\ll M^2$. For an electron pair with generic $y\ne0$, unequal masses and a surviving $\Delta L_e=2$ derivative amplitude, the specified field content and couplings supply no conserved electron lepton number that enforces this tuning. An internal charge of $X_e$ assigns the same charge to the Weinberg bilinear and its two-derivative version. Unequal thresholds admit different radiative corrections. A claim of protection or a loop-level prediction for $C_5$ would require additional matching or a specified symmetry. For comparison, two equal masses with opposite squared Yukawas would cancel $H_e(z)$ identically, including its slope. The muon condition independently fixes a mass-weighted moment. Perturbing the relative Yukawa generically restores that moment as well.

The construction establishes a renormalizable example of the stated lepton-scalar matching in the unbroken theory. It provides no measured neutrino masses, absolute broken-phase rate, ultraviolet theory of gravity, or test that selects the broader hypothesis over other singlet models.

Run `python -I -B -O check_matching.py` in this directory for the JSON checks. To write the same JSON to a chosen path use `--output PATH`. [Matching checker](check_matching.py) · [Receipt](MATCHING.json)
