# The two-scalar cut from physical tree amplitudes

For massless external leptons and a scalar doublet of common squared mass
$a>0$ treated with fixed-order free propagators, the retained-Majorana
triangle has an on-shell two-scalar
discontinuity for $s>4a$. It can survive an exact cancellation of the
primitive four-lepton UV pole. The calculation here sews two tree amplitudes.
A separate Feynman-parameter computation checks the sign and normalization
in the ordered action convention. [October four-lepton pole](../../checkpoint_2026_10_23/seed_yukawa/FOUR_LEPTON_POLE.md)

## External states and the two trees

Use $S=1+iT$, strip the momentum-conserving delta function from $T$,
and order the incoming state as $a_1^\dagger a_2^\dagger|0\rangle$, with
the outgoing ket ordered likewise. Let $p_1,p_2$ be incoming and
$p_3,p_4$ outgoing future-directed momenta. Write
$B_{ij}=u_i^T\epsilon u_j$, where $\epsilon_{12}=1$. Spinors have
$|B_{ij}|^2=2p_i\cdot p_j$. The external spin factor of October's
ordered polynomial $P=(\bar l\epsilon\bar l^T)(l^T\epsilon l)$ is
$W=-B_{34}^* B_{12}$ per unsigned permutation. The minus follows because
$l_\gamma l_\delta$ on the ordered incoming state gives the wavefunction
of particle 2 in the first slot minus that of particle 1 in the first
slot. The analogous barred contractions on the outgoing bra have the
opposite order. Antisymmetry of the epsilon matrix absorbs the permutation
of spinors within each bilinear. Exchanging two complete external labels
changes the sign of the physical matrix element. In the forward elastic
matrix element $W=-s$.

Write $a_i$ for external weak components, $I,J$ for intermediate weak
components, $M=U\operatorname{diag}(m_j)U^T$, and
$A^{(j)}_{rs}=m_j(Y^*U)_{rj}(Y^*U)_{sj}$, with flavor labels $r,s$
and $b_j=m_j^2$. The charge of the two scalar states is the one connected
to the chosen two-lepton state by the displayed vertices. The charge
conjugate process is obtained by conjugating the external states as well.
Let $k_1,k_2$ be outgoing from the initial pair into the intermediate
state and define $t_i=(p_i-k_1)^2, u_i=(p_i-k_2)^2$.
For fixed components the displayed seed tree and sterile tree are

$$
\begin{aligned}
{\cal M}^{C}_{34,IJ}
 &= -B_{34}^* C_{34}
    (\delta_{a_3 I}\delta_{a_4 J}t_3+
     \delta_{a_3 J}\delta_{a_4 I}u_3),\\
{\cal M}^{Y}_{IJ,12}
 &= B_{12}\sum_j A^{(j)}_{12}
  \left(\frac{\delta_{a_1 I}\delta_{a_2 J}}{t_1-b_j}
       +\frac{\delta_{a_1 J}\delta_{a_2 I}}{u_1-b_j}\right).
\end{aligned}
$$

Numerical subscripts on $A,C$ in these formulas indicate the corresponding
external flavor labels. The seed action has a Weyl minus from charge
conjugation and a half for the symmetric seed pair. External contractions
absorb the half. The derivative acts on a single *whole* composite. The
momentum in it is $p_i-k_j$, and the other composite carries its negative.
In particular $p_i^2=0$ does not set the seed to zero. The sterile tree
has an anomalous numerator $m_j$, and its propagator denominator depends
on which scalar goes to which lepton. The same-sign sum of the two weak
assignments follows from factoring out the antisymmetric external spinor
product. In the October ordering the two Yukawa vertices, anomalous
propagator and Wick crossing give the ordered two-unbarred-field action
coefficient $-A^{(j)}/[2(t-b_j)]$ before external permutations.
The minus from the incoming state then gives the sign displayed above.

## Angular kernel and Hermitian cut

Let $P=p_1+p_2, s=P^2, \beta=\sqrt{1-4a/s}$,
$T_0=a-s/2, h=s\beta/2$. In the center of mass define
$\boldsymbol n=\hat{\boldsymbol k}_1, x=\boldsymbol n\cdot\hat{\boldsymbol p}_1, y=\boldsymbol n\cdot\hat{\boldsymbol p}_3, c=\hat{\boldsymbol p}_1\cdot\hat{\boldsymbol p}_3$.
Then $t_1=T_0+hx, u_1=T_0-hx, t_3=T_0+hy, u_3=T_0-hy$.
Since $T_0<-h$, neither exchanged propagator has a physical pole on
this cut, even for $b=0$. The Lorentz-invariant measure for two ordered
momenta is $d\Phi_2=\beta d\Omega/(32\pi^2)$.
Unitarity gives the difference $T-T^\dagger$ divided by $2i$ as
one half of the sewn product integrated over intermediate states.
Completeness of the scalar Fock states also gives $1/2!$ when summing
over both ordered weak labels and the full sphere.

Set $D=\delta_{a_3a_1}\delta_{a_4a_2}$ and
$X=\delta_{a_3a_2}\delta_{a_4a_1}$. At fixed $\boldsymbol n$, define
the seed-on-outgoing scalar kernel including the intermediate-state half,

$$
\begin{aligned}
G_b(x,y)=
&\frac D2\left[\frac{T_0+hy}{T_0+hx-b}
                +\frac{T_0-hy}{T_0-hx-b}\right]\\
+&\frac X2\left[\frac{T_0-hy}{T_0+hx-b}
                +\frac{T_0+hy}{T_0-hx-b}\right].
\end{aligned}
$$

The differential cut is $W\,d\Phi_2/2$ times
$\sum_j[C_{34}A^{(j)}_{12}G_{b_j}(x,y) +(A^{(j)}_{34})^* C_{12}^*G_{b_j}(y,x)]$.
The second term sews the reversed tree orientations. The rotation
exchanging the two external axes gives them equal integrated angular
kernels. For $b\geq0$ set

$$
J_b=\frac1{2h}\log\frac{b-T_0+h}{b-T_0-h},\qquad
K_b(c)=c-[cb+(1-c)T_0]J_b .
$$

Indeed the azimuthal average of $y$ at fixed $x$ is $cx$. The
averages of $1/(T_0+hx-b)$ and $hx/(T_0+hx-b)$ are respectively
$-J_b$ and $1-(b-T_0)J_b$. With lower and upper rims related by
Hermitian analyticity, the projected cut reads

$$
\frac{\operatorname{Disc}_s {\cal M}_{34,12}}{2i}
=\frac{\beta W}{16\pi}\sum_j
  [C_{34}A^{(j)}_{12}+(A^{(j)}_{34})^* C_{12}^*]
  [D K_{b_j}(c)+X K_{b_j}(-c)] . \tag{1}
$$

For different external states the left side is
$[ {\cal M}_{34,12}(s+i0)-{\cal M}_{12,34}(s+i0)^*]/(2i)$
with the incoming and outgoing states interchanged in the second term.
For a forward diagonal matrix element it is the imaginary part. A complex
coefficient in a single ordered graph cannot in general be replaced by
twice its real part. There is no weak trace or additional Majorana symmetry
factor. Differentiating the October contact $+C_{34}A_{12}P/2$ gives
$W C_{34}A_{12}(D+X)$ before its conjugate and loop normalization.
Thus a unit scalar integral multiplies the same external tensor in both
calculations.

## An independent sign and normalization check on the loop

For the direct assignment use $D_1=k^2-a, D_2=(P-k)^2-a, D_N=(p_1-k)^2-b$, with the usual $+i0$.
The transverse part of $p_3$ integrates to zero by reflection about
the $P,p_1$ plane. In the triangle numerator this gives the replacement

$$
(p_3-k)^2\longrightarrow c(D_N+b)
 +(1-c)[(D_1+D_2)/2+T_0].
$$

The terms cancelling $D_1,D_2$ leave scalar–sterile bubbles at external
invariant zero. The two-scalar bubble, in units $i/(16\pi^2)$, has
absorptive part $+\pi\beta$, obtained from the interval on which
$a-sx(1-x)<0$ in its negative logarithm. The scalar triangle in
the same units is $-\int dx\,dy/(\Delta-i0)$ on the unit simplex,
where $\Delta=a(x+y)+b(1-x-y)-sxy$.
Its absorptive part is $-\pi\int dx\,dy\,\delta(\Delta)$.
To evaluate this integral set $v=x+y, d=(x-y)/v$ and integrate over
both roots of the quadratic in $d$. With
$v_\pm=2(a-b)/s\pm\sqrt{4(a-b)^2/s^2+4b/s}$,

$$
\int dx\,dy\,\delta(\Delta)=\frac{2}{s}\int_{v_+}^1
 \frac{dv}{\sqrt{(v-v_+)(v-v_-)}}=\beta J_b.
$$

Combining the bubble and triangle gives $+\pi\beta K_b(c)$.
After restoring $1/(16\pi^2)$, this agrees with (1).
The crossed assignment corresponds to $c\to-c$. The executable builds
all CM four-vectors, sums the intermediate weak labels explicitly, and
compares the sphere quadrature to the angular formula. It also integrates
the parameter-space delta independently. [Standalone checker](independent_cut.py)
and [receipt](INDEPENDENT_CUT_RESULTS.json)

## Threshold and pole-zero benchmark

As $s\downarrow4a$, $J_b\to1/(a+b)$, and
$K_b(c)\to a/(a+b)$, independently of angle. The reduced cut therefore
opens as $\beta$ if its mass and flavor sum survives. Tensor components
can cancel the leading coefficient. The angular part
$1-(b-T_0)J_b=-h^2/[3(b-T_0)^2]+O(h^4/(b-T_0)^4)$
contributes at order $\beta^3$ to (1).

For identical weak components $D=X=1$, the sum is $-2T_0 J_b$.
Take $a=1$, one active flavor, $C=1$, $m_1=1, m_2=4$,
$(Y^*U)=(1,i/2)$. The pole weights are $1,-1$.
After removing $\beta W/(16\pi)$, the full two-orientation cut
tends to $4(1/2-1/17)=30/17$ at threshold. It stays strictly positive
for each finite $s>4$, since the angular average
$J_b=\langle(b-T_0-hx)^{-1}\rangle$ strictly decreases with $b$.
At $s=8$ the reduced coefficient is $3.1030440581306561$.
The forward spin factor is $-8$, giving the contribution
$-0.3492151494111737$ to (1) in the chosen units.

Only the term linear in the seed and quadratic in the Yukawa couplings
is retained. The signed interference may be negative, while the complete
forward optical sum contains the squared modulus of the sum of amplitudes.
The benchmark has no open two-body scalar decay into a lepton and sterile
fermion through the specified Yukawa vertices. More generally, if a
coupled sterile mode has $m_j<\sqrt a$, that scalar decay is kinematically
allowed. Equation (1) then describes the perturbative cut of the free
propagators at this order. A physical rate needs the appropriate
asymptotic states, including decay products for unstable scalars, and the
relevant interactions. Local polynomial counterterms have zero
discontinuity across this cut.
