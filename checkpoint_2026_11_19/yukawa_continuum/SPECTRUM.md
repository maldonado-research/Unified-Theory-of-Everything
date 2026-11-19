# The charged-scalar pole and its Yukawa continuum

In the unbroken gaugeless four-singlet model, each charged scalar line in the leading four-lepton box can propagate a massless charged lepton and a Majorana singlet. For one coupling $Y$, let $a=m^2>0$ be the physical squared scalar mass and $b=M_N^2\ge a$ the singlet squared mass. With scalar pole residue fixed to one, the positive continuum density per unit invariant mass squared $v$ is
$$
 \rho_b(v)=\frac{|Y|^2}{16\pi^2}\,
  \frac{(v-b)^2}{v(v-a)^2}\,\theta(v-b).
$$
For the column with $b=a$, the two powers cancel exactly and $\rho_a(v)=|Y|^2\theta(v-a)/(16\pi^2v)$. The two scalar lines in the box give identical continuum insertions for every mass assignment. [Original action and flavor texture](../uv_completion/MATCHING.md) · [Leading box](../uv_completion/DISPERSION.md)

## Spin sum, weak component and sign

Use $\eta=(+---)$, $S=1+iT$ and the interaction $-\bar L^A Y\phi^A N_R+\mathrm{h.c.}$, with $\phi=\widetilde H$. For a charged scalar of virtual mass $\sqrt v$, label the outgoing massless lepton momentum $p$ and the Majorana momentum $k$. Cutting the charged scalar two-point function gives
$$
 \sum_{\rm spins}|\mathcal A|^2
 =|Y|^2\operatorname{tr}[\not p P_R(\not k+\sqrt b)P_L]
 =2|Y|^2p\cdot k=|Y|^2(v-b).
$$
The scalar has one spin state. The chiral projector selects one massless helicity and the Majorana spin sum includes both physical spin states. The lepton and singlet are distinguishable. The adjoint Yukawa vertex gives the decay of the opposite scalar charge. In the $H$ variables, the two antisymmetric weak tensors contract to $\sum_C\epsilon_{AC}\epsilon^*_{BC}=\delta_{AB}$. For fixed charged $A=B$ only the charged lepton runs. A trace over both scalar components or an extra sum over conjugate decays would double the answer. In the original flavor basis new thresholds add incoherently in this two-point function, with weights $|Y_{\alpha i}|^2$.

Write the causal propagator as $iD(z)$ with $D(z)=[z-a+\Pi_R(z)]^{-1}$, where $i\Pi$ is the amputated loop. The integrated two-body phase space is $(v-b)/(8\pi v)$ above threshold. The two-point optical theorem gives
$$
 \operatorname{Im}\Pi(v+i0)
 =\tfrac12\int d\Phi_2\sum|\mathcal A|^2
 =\frac{|Y|^2}{16\pi}\frac{(v-b)^2}{v}.
$$
It also gives the positive virtual decay width $\operatorname{Im}\Pi(v+i0)/\sqrt v$, in agreement with the published charged-scalar width for this Yukawa interaction. The sign in the inverse propagator now fixes $\rho_b=-\operatorname{Im}\delta D/\pi$ away from the pole. [Charged-scalar decay, Eq. (5.6)](https://arxiv.org/pdf/2207.07382v2)

## Pole subtraction and a spacelike loop check

The two local counterterms impose $\Pi_R(a)=\Pi_R'(a)=0$. For $b=a$ take both limits from below. The dispersion integral and the pole-normalized propagator through order $|Y|^2$ are
$$
 \Pi_R(z)=(z-a)^2\int_b^\infty\frac{\rho_b(v)\,dv}{v-z},
 \qquad D(z)=\frac1{z-a}+\int_b^\infty\frac{\rho_b(v)\,dv}{z-v}
 +O(|Y|^4).
$$
For $z=-Q^2<0$, $-D(-Q^2)$ increases by $\int\rho_b(v)dv/(Q^2+v)>0$. At $b=a$ this increment equals $|Y|^2\ln(1+Q^2/a)/(16\pi^2 Q^2)$. The endpoint continuum gives a logarithm in $D$ as $z\to a$. The loop itself vanishes as $(z-a)^2\ln(1-z/a)$, so the coefficient of the simple pole stays one. The perturbative spectral measure has a delta function of unit weight at $a$ and a continuum starting there. If $b<a$, the pole is inside an open decay channel and these stable-pole formulas do not apply.

An independent spacelike calculation starts with the closed chiral loop. For massive denominator $d_1=k^2-b$ and massless denominator $d_2=(k-q)^2$, its trace is $2k\cdot(k-q)=d_1+d_2+b-z$ with $z=q^2$. The fermion loop sign gives $(z-b)B_0(z;0,b)-A_0(b)$ in units of $|Y|^2/(16\pi^2)$, with the scaleless tadpole zero in dimensional regularization. Its UV residue is $z-2b$ in $d=4-2\epsilon$, in agreement with the earlier independent two-point calculation. Up to an additive constant, $B_0=-\int_0^1 dx\,\ln(1-xz/b-i0)$. Integrating $x$ gives the one-loop inverse kernel up to an affine polynomial,
$$
 F_b(z)=-\frac{(z-b)^2}{z}\ln(1-z/b),\qquad
 \frac{16\pi^2}{|Y|^2}\Pi_R(z)=F_b(z)-F_b(a)-(z-a)F_b'(a).
$$
The apparent singularity at $z=0$ is removable, with $F_b(0)=b$. At $b=a$ the endpoint values of this representative satisfy $F_a(a)=F_a'(a)=0$. The upper rim of the same logarithm has imaginary part $-\pi$ for $z>b$. These signs reproduce the positive cut. [Literal UV polynomial and weak multiplicities](../../checkpoint_2026_10_16/yukawa_counterterms/YUKAWA_TWO_POINT.md)

At $z=-2a$ the positive continuum integral, multiplied by $16\pi^2a/|Y|^2$, takes the following values. The checker evaluates the cut integral, the subtracted closed form and the subtracted original bubble parameter integral separately. [Checker](check_spectrum.py) · [Numerical receipt](SPECTRUM.json)

| $b/a$ | Continuum integral |
|---:|---:|
| 1 | 0.549306144334055 |
| 16 | 0.0208573166820733 |
| 100 | 0.00333343300850122 |
| 400 | 0.000833334894539603 |

## Insertion into the four-lepton box

For two scalar squared masses $a_1,a_2$ and two singlet squared masses $d_i,e_j$, the ordered scalar box is
$$
 D_4(s,t;a_1,a_2,d_i,e_j)
 =\int\frac{[dx]}
 {[a_1x_1+a_2x_2+d_ix_3+e_jx_4-sx_1x_2-tx_3x_4-i0]^2},
 \quad [dx]=\prod_{k=1}^4 dx_k\,\delta(1-\textstyle\sum x_k),\quad x_k\ge0 .
$$
The self-energy is diagonal in the scalar indices and carries no external spin tensor. The equal-chirality singlet numerators and the external factor $B_{34}^*B_{12}=s$ are those of the leading amplitude. In that amplitude make the replacement, individually for transfers $t$ and $u$,
$$
 D_4(s,t;a,a,d_i,e_j)\ \longmapsto\
 2\sum_{\alpha,n}\int_{b_n}^\infty dv\,\rho_{b_n}(v;Y_{\alpha n})\,
 D_4(s,t;v,a,d_i,e_j)
$$
to obtain the correction alone. Expanding the product of the two scalar propagators places one continuum on either line. Interchanging $x_1,x_2$ preserves both kinematic products and proves equality of these insertions before the external exchange. The factor two counts scalar sites. The leading $t,u$ sum still counts the two scalar pairings. The signed outer mass residues in the benchmark are $w_i=(m,-m)$ and $v_j=(M,-4M)$ after stripping $\gamma=y^2(g^*)^2$. The positive inner weights, in order of increasing singlet mass at $M=10m$, are $|g|^2(1,1/4)$ and $|y|^2(1,2)$. The mass factors or relative phases of the outer residues must not be copied into the inner weights. [Original simplex and identical-lepton exchange](../uv_completion/BOX_CROSSCHECK.md)

At large $v$, $\rho_b(v)\sim |Y|^2/(16\pi^2 v)$. Its integral alone diverges logarithmically, as expected for a pole-normalized field at this order. Each propagator insertion converges because of the additional $1/(v-z)$. At fixed $s,t\le0$ and all other masses squared $\ge a$, the denominator before squaring in $D_4$ is at least $a+(v-a)x_1$. Integrating the other parameters over their triangle gives $D_4\le\tfrac12\int_0^1 dx_1/[a+(v-a)x_1]^2=1/(2av)$. This proves absolute convergence for the Euclidean box convolution before the finite flavor sum. Timelike values are the causal continuation. The mass-$\sqrt v$ line together with the pole-mass line has its ordinary threshold at $(\sqrt v+\sqrt a)^2$, so the continuum contribution starts at $(\sqrt b+\sqrt a)^2$.

In the original model, the monomials $\gamma|g|^2$ and $\gamma|y|^2$ also receive diagrams with active or singlet fermion self-energies, Yukawa-vertex or other connected contractions and the associated parameter and external-state renormalization where applicable. A correction from the scalar two-point function alone therefore cannot determine either full two-loop coefficient. Finite changes in normalization of the scalar field or of the input Yukawas redistribute terms proportional to the leading box.

One can isolate the insertion by a formal flavor extension. Add a massless lepton doublet $L_f$ and a distinct Majorana $N_f$ coupled only to one another by $h\bar L_f\phi N_{Rf}$ and its adjoint, with diagonal, unmixed mass terms and $b_f\ge a$. At order $\gamma|h|^2$ the only two spectator vertices have to contract both fermion slots with each other. In a connected four-lepton graph the two remaining scalar ends lie on one of the original scalar lines. Fix the scalar pole mass and residue and hold the active amputated Yukawa couplings in those field coordinates fixed. The spectator loop’s contribution to the scalar mass and field subtractions is included in $\Pi_R$ for the scalar line. The insertion gives the entire explicit spectator contribution at this monomial. Reexpressing the result in another Yukawa input scheme can add a finite multiple of the leading box. To count the contractions explicitly, pair three charged scalar slots with three conjugate slots. Of the six permutations, two close the spectator lines into a disconnected vacuum graph. Each of the two original box pairings occurs twice among the four connected assignments. The extra row and column define a diagnostic model. Attaching a new lepton to an existing singlet would allow a spectator correction to that singlet propagator as well.
