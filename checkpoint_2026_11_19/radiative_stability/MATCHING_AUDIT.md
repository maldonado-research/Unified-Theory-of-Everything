# The one-loop electron Weinberg coefficient

In the unbroken gaugeless model, imposing the two electron Yukawas and masses at one common $\overline{\rm MS}$ scale gives zero one-loop pure-Yukawa shift of the electron Weinberg coefficient at order $1/M$. This includes both $y^2|y|^2$ and $y^2|g|^2$. Adding $V_4=\lambda(H^\dagger H)^2$ gives a positive shift relative to $\lambda y^2/M$,

$$
\Delta C^{ee}_{5,\lambda}=\frac{2\lambda}{16\pi^2}\frac{y^2}{M}\ln4 .
$$

Here $y$ is the electron–singlet Yukawa, $g$ is the muon–singlet Yukawa, and $C_5$ multiplies $\bar X_e X_e^c/2$ plus its adjoint. The complex phase is that of $y^2$. The scalar has positive pole mass squared $a=m^2$ and the light leptons are massless. The coefficient is evaluated to first order in the loop expansion and to leading order in light scales over $M$. [Benchmark action and sign](../uv_completion/MATCHING.md) · [Independent executable](check_audit.py) · [Receipt](AUDIT.json)

## Common-scale renormalized matching

Let $v_i$ be a column of $Y$, $M_i>0$ its Takagi mass, $C_i=v_i v_i^T/M_i$ and $C=\sum_i C_i$. In the simultaneous type-I seesaw matching, write $L_i=\ln(\mu^2/M_i^2)$ after modified minimal subtraction. Zhang and Zhou give the gaugeless hard vertex and field normalizations through this order as

$$
\begin{aligned}
16\pi^2\,\delta C^{\rm G}&=2\lambda\sum_i(1+L_i)C_i,\\
16\pi^2\,\delta Z_H&=\tfrac12\sum_i(1+2L_i) v_i^\dagger v_i,\\
16\pi^2\,\delta Z_L&=\tfrac14\sum_i(3+2L_i)v_i v_i^\dagger,\\
C_{5,\rm eff}&=C+\delta C^{\rm G}-\delta Z_H C
              -\tfrac12(\delta Z_L C+C\delta Z_L^T).
\end{aligned}
$$

The finite constants in these expressions matter for generic textures. Canonical normalization of the two scalar legs and two lepton legs supplies the last three field terms. The field redefinition from eliminating the dimension-six scalar derivative operator contributes to $C_5$ proportional to a light mass squared times $C/M^2$. The authors explicitly discard this higher-order term in the displayed dimension-five matching. Equations here follow the v3 PDF, where the potential, tree coefficient, one-loop fields and logarithms, hard vertex, EOM term, field replacements and final coefficient occur at Eqs. (34), (38), (46)–(48), (85), (86) and (88), respectively. [Zhang & Zhou, arXiv:2107.12133v3](https://arxiv.org/pdf/2107.12133v3)

For the electron singlets $C_{1,ee}=y^2/M$ and $C_{2,ee}=-y^2/M$. Every $C_{e\alpha}$ is zero in the sum. This removes every field term in the electron diagonal entry, even when the muon inverse mass moment is nonzero. The hard quartic term is proportional to $L_1-L_2=\ln4$. The two additive finite constants cancel in the texture, and the matching scale drops out at this order. All parameters specifying the texture are renormalized at the same scale. A mass definition or a finite redefinition of the Yukawas changes the input relation that expresses this tuning.

The published SMEFT result integrates out all contributing singlets. At the electron threshold of this benchmark the $m,4m$ muon singlets remain dynamical. The hard quartic triangle contains only one electron singlet. The graph argument below treats the muon masses exactly. Purely light external-leg corrections are common to the full and partially matched theories, and a scalar field factor induced by the muon sector multiplies the zero electron tree coefficient. These facts justify the displayed electron entry in the partially matched EFT.

There is a separate check of the logarithm between the two electron thresholds. At order $\lambda y^2$, the already integrated $2M$ column gives $C_{ee}=-y^2/M$ just below its threshold. Its quartic running is $16\pi^2 dC_{ee}/d\ln\mu=4\lambda C_{ee}$. Running to $M$ adds $+4\lambda(y^2/M)\ln2/(16\pi^2)$. The separate quartic matching constants at their own thresholds sum to zero. For sequential calculations involving Yukawas, the intermediate theory has additional terms in the running of both the surviving Yukawa and the Wilson coefficient. Zhang derives their cancellation in the evolution of the combined tree-plus-local coefficient. [Zhang, arXiv:2405.18017v3, Eqs. (8)–(11)](https://arxiv.org/pdf/2405.18017v3)

## Why the pure-Yukawa vertex vanishes

Use two-component spinors and choose the orientation with external $LLHH$. In a four-Yukawa graph there are three $LHN$ vertices and one conjugate vertex. The two external leptons and two external scalars all attach to the first three. Label by $h,l,n$ the vertices joined to the scalar, lepton and sterile lines of the conjugate vertex. The two other sterile legs must be paired. This accounts for all contractions up to external permutations and conjugation.

* If $n=h=l$ the graph is disconnected.
* If $n=h\ne l$ or $n=l\ne h$ the single light line running from the conjugate vertex to the remaining endpoint is a bridge. Cutting it disconnects the graph. Such terms are removed by light-particle amputation or handled by field and parameter matching.
* In every remaining graph, vertex $n$ holds both an external lepton and an external scalar. Its sole internal line is a mixed-chirality sterile propagator with numerator $q\cdot\sigma$, where $q$ is their total incoming momentum. It vanishes when the external momenta are set to zero. At positive scalar and sterile masses the remaining loop has at most one massless lepton line and is infrared regular there.

The 27 assignments give 3 disconnected graphs, 12 connected light-reducible graphs and 12 one-light-particle-irreducible graphs. This proves the absence of a pure-Yukawa zero-momentum bare vertex before expanding any massive propagator. The $\overline{\rm MS}$ renormalized coefficient additionally requires the canonical factors above. As a check on the UV counterterms, above the thresholds the purely singlet-Yukawa RGEs imply $16\pi^2\,dC/d\ln\mu=2\operatorname{tr}(P)C+(P C+C P^T)/2$ with $P=YY^\dagger$. Its $ee$ entry vanishes for this texture. Higher loops, dimension-seven matching and finite-transfer on-shell matrix elements involve additional structures. The existing on-shell tree projection contains the constant $a C_7$. [Antusch et al., Table 1 and Appendix D.1](https://arxiv.org/pdf/hep-ph/0501272v3) · [Tree-level on-shell projection](../uv_completion/MATCHING.md)

## Direct zero-momentum quartic loop at finite scalar mass

For the tuned electron pair, the two zero-momentum Euclidean triangle integrands may be combined before taking $a\to0$. Strip the Yukawa $y^2$ and let $u=k_E^2$. After the angular integral the positive-simplex triangle sum becomes

$$
S_e(0)=\int_0^\infty du\,\frac{u}{(u+a)^2}
  \left(\frac{M}{u+M^2}-\frac{4M}{u+4M^2}\right)
 =-3M\int_0^\infty du\,
  \frac{[u/(u+a)]^2}{(u+M^2)(u+4M^2)}.
$$

The integrand has definite negative sign. For every $a\ge0$ the absolute value of the last integrand is bounded by its integrable $a=0$ value. Therefore $\lim_{a\to0} M S_e(0)=-\ln4$. The quartic correction to the zero-momentum coefficient in the original fields is $-2\lambda y^2 S_e(0)/(16\pi^2)$. The charged-scalar rule gives $-4i\lambda$ and the tree normalization of the local two-lepton, two-scalar vertex supplies the factor of two in the denominator. This yields the common-scale result directly.

At positive $a$ the closed form is $S_e(0)=M[T(M^2,a)-4T(4M^2,a)]$, with
$$
T(B,a)=\frac{B\ln(B/a)-B+a}{(B-a)^2},\qquad T(a,a)=\frac{1}{2a}.
$$
The numerical integral gives $M S_e(0)=-1.33045553538$ at $M/m=10$, compared to the leading term $-1.38629436112$. The finite-$a$ result is a specified off-shell full-theory vertex. Its power corrections include soft contributions. Matching higher-dimensional Wilson coefficients requires subtracting the corresponding EFT amplitude and choosing an operator basis.

## Completeness of the single-quartic four-lepton coefficient

For $\mu^-_L\mu^-_L\to e^-_Le^-_L$ at order $\lambda\gamma$, $\gamma=y^2(g^*)^2$, the four active-lepton legs of the four Yukawa vertices are external. The sterile pairings are fixed by flavor. The four charged scalar ends on the Yukawa vertices and four scalar ends on the quartic vertex exhaust the remaining possibilities. Classify them by the number of quartic self-contractions.

| Self-contractions | Charged Wick contractions with labeled slots | Effect |
|---|---:|---|
| 0 | 4 | Each Yukawa scalar joins the central quartic. The two independent loop momenta give a product of triangles. |
| 1 | 16 | A scalar tadpole corrects one propagator in either one-loop box. |
| 2 | 4 | A disconnected quartic vacuum bubble cancels in normalized amplitudes. |

The neutral component of the doublet also contributes to the scalar tadpole. The tadpole is independent of momentum, so fixing the charged scalar pole mass $a$ cancels it with the mass counterterm, and its derivative gives zero residue correction. At this coupling order there is no lepton or sterile self-energy proportional to $\lambda$ alone and no quartic-only one-loop renormalization of a Yukawa vertex in the zero-background theory. A calculation holding a scalar running mass fixed instead must keep the mass-conversion term acting on the one-loop boxes.

Choose external spinor phases $B_{34}^* B_{12}=s$ as in the previous box calculation. Write $S_\ell(s)=\sum_i w_{\ell i}T_i(s)$, using the real mass-weighted residues after removing $\gamma$ and positive-simplex upper-rim triangles
$$
T_i(s)=\int_{x,y\ge0,\ x+y\le1}
 \frac{dx\,dy}{a(x+y)+M_i^2(1-x-y)-sxy-i0}.
$$
A box carries $i^4$ from its propagators and $i(-1)^4/(16\pi^2)$ from one Wick-rotated integral. The six propagators of the double triangle give $i^6$, the quartic gives $-4i\lambda$, and each of its two integrations gives $i(-1)^3/(16\pi^2)$. The ordered spinor factors and four Yukawa vertices agree with the box. Dividing by the common $i$ in $i\mathcal M$ gives
$$
\Delta\mathcal M_\lambda(s)
 =-\frac{4\lambda s\gamma}{(16\pi^2)^2}S_e(s)S_\mu(s).
$$

For an $s$-channel unitarity cut, the muon and electron Yukawa vertices must be on their respective external sides. Assigning the central quartic to either side gives the two scalar-pair cuts. After removing the common weak phase, the upper-rim product has imaginary part $\mathrm{Re}\,S_e\,\mathrm{Im}\,S_\mu+\mathrm{Im}\,S_e\,\mathrm{Re}\,S_\mu$. This graph count fixes the coefficient linear in $\lambda$ at the stated Yukawa degree and pole-mass convention. Other coupling orders and the rest of the two-loop amplitude require their own calculation. [One-loop box normalization](../uv_completion/BOX_CROSSCHECK.md) · [Two-body cut convention](../uv_completion/CUT.md)

Run `python -I -B -O check_audit.py` in this directory. The checker enumerates both graph sets, evaluates the finite common-scale field factors for a complex texture and a detuned control, and integrates the signed Euclidean combination independently of its analytic primitive. The optional `--output PATH` writes the same JSON.
