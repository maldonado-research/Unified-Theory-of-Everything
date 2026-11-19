# A two-scalar cut in a physical four-lepton channel

In the unbroken gaugeless theory, consider $\mu^-_L(p_1)\mu^-_L(p_2)\to e^-_L(p_3)e^-_L(p_4)$ with massless external leptons and both scalars of physical squared mass $a>0$. For the primitive triangle containing one derivative seed and two Yukawa vertices, opposite UV weights at two unequal positive sterile masses leave a discontinuity at every $s>4a$. For scalar mass one, sterile masses 1 and 4, and weights $+1,-1$, the reduced discontinuity divided by $2\pi i$ is **0.5485458739562 at $s=8$**. Both modes have the same two-scalar threshold. Their different weights above threshold come from the *uncut* sterile propagator.

## External states and normalization

Use metric $(+---)$, $S=T\exp(i\int\mathcal L_{\rm int})$, $\epsilon_{12}=+1$, $\phi=\widetilde H$, and the literal action
$$
 \mathcal L_{\rm seed}=\tfrac12 C_{\alpha\beta}
 (\bar L_\alpha^a\phi^a)(-\Box)(\phi^b L_\beta^{c,b})+\mathrm{h.c.},
 \qquad \mathcal L_Y=-\bar L_\alpha^aY_{\alpha i}\phi^aN_{Ri}+\mathrm{h.c.}
$$
Work in a positive Takagi mass basis. Select $C_{ee}$ as the only nonzero seed entry and the $\mu$ row as the only nonzero Yukawa row. Both external pairs use the charged lower weak component. The displayed orientation can connect these states. The conjugate seed annihilates two electrons and therefore has zero matrix element in this channel. For each sterile mode put $b=m_i^2$ and $w_i=m_i(Y^*_{\mu i})^2$. The square comes from two conjugate Yukawas and the anomalous Majorana propagator. A massless sterile mode has zero weight.

Use $|1,2\rangle=c_1^\dagger c_2^\dagger|0\rangle$, $|3,4\rangle=c_3^\dagger c_4^\dagger|0\rangle$ and commuting negative-helicity wavefunctions. Factor the ordered tensor
$$
 S_{34;21}=(\bar u_3\epsilon\bar u_4^T)(u_2^T\epsilon u_1).
$$
Annihilating the ket fixes the incoming order $2,1$. With $u^\dagger u=2E$ for each leg, $|S_{34;21}|=s$ for $s>0$. For example, choose $\mathbf p_1$ along $+z$, $\mathbf p_3$ at polar angle $\theta$ in the $xz$ plane, and unit wavefunctions $(0,1),(1,0),(-\sin(\theta/2),\cos(\theta/2)),(\cos(\theta/2),\sin(\theta/2))$ in the order $1,2,3,4$. Restoring their $\sqrt{2E}$ factors gives $S_{34;21}=-s$ in this phase convention.

Write a Feynman amplitude insertion as $i\mathcal M$, with the overall momentum-conserving delta stripped. The primitive contributes
$$
 \mathcal M(s)=\frac{2 C_{ee}S_{34;21}}{16\pi^2}
                    \sum_i w_i K(s;a,m_i^2).                                          \tag{1}
$$
Here $K$ has UV pole $+1/\epsilon$ in $d=4-2\epsilon$. At fixed labeled momenta the coefficient of each scalar assignment, before the loop integral and external contractions, is
$$
  (-i/2)(-i)^2 i^2 (-i)(-1)/2! = +1/4 .
$$
These factors include the seed spin-conjugation sign, both signed Yukawas, two scalar propagators, the stripped sterile propagator and the crossing of one open lepton by its Wick contraction. Each of the two scalar Wick assignments admits two assignments of the incoming and two of the outgoing external sources. Reversing either pair changes both the commuting spin determinant and the fermion permutation sign, so all eight terms have the same sign in $S_{34;21}$. Four give each of the two loop numerators below. This yields $2=8/4$ multiplying the *mean* of the two kernels in (1). No identical-particle phase-space factors enter the external amplitude. This normalization also follows by differentiating the ordered pole action. [Literal pole and Wick derivation](../../checkpoint_2026_10_23/seed_yukawa/FOUR_LEPTON_POLE.md)

For the earlier one-active-flavor benchmark with a real seed and real $w_i$, elastic same-weak scattering also receives its conjugate orientation. On shell, swapping the external pairs preserves $s$ and the reduced symmetrized kernel below. Both orientations thus give the same $K(s+i0)$, with total coefficient $4 C S_{34;21}/(16\pi^2)$ times $\sum_i w_i K$. Conjugation acts on couplings and external fields. Both causal loop integrals use $+i0$.

## Distinct scalar routings and reduction

Let $P=p_1+p_2=p_3+p_4$ and $s=P^2$. Route $k$ along the scalar ending at the $p_1$ Yukawa vertex. The other scalar has momentum $P-k$, and the sterile line has momentum $k-p_1$. Write
$$
 D_1=k^2-a,\qquad D_2=(P-k)^2-a,\qquad D_N=(k-p_1)^2-b .
$$
The whole neutral composite carries a charged external lepton momentum *and* a scalar momentum. When this scalar is paired with external $p_j$ at the seed, the numerator is $N_j=(k-p_j)^2$, $j=3,4$. The complementary composite has the opposite total momentum by conservation and gives the same square. The two integrated routings can differ at fixed angle. Their equal multiplicities in the equal-weak-component channel yield the exact identity
$$
 \tfrac12(N_3+N_4)=k^2-k\cdot P
              =\tfrac12(D_1+D_2)+a-s/2.                                              \tag{2}
$$
Keeping only derivatives on an external lepton would instead give $p_j^2=0$.

Normalize Minkowski scalar masters by extracting $i/(16\pi^2)$ from $\int d^d k/(2\pi)^d$, with the usual regulator scale and Feynman $+i0$. The two pinched bubble integrals in (2) have external invariant $p_1^2$ or $p_2^2$, both zero. In this normalization the three-denominator master is $-T(s)$, where $T$ below has a positive Euclidean parameter integrand. Thus
$$
 \begin{aligned}
 K(s)&=B_0(0;a,b)+(s/2-a)T(s),\\
 T(s)&=\int_{x,y\ge0,\ x+y\le1}
  \frac{dx\,dy}{d(x,y)-sxy-i0},\qquad
 d(x,y)=a(x+y)+b(1-x-y).
 \end{aligned}                                                                         \tag{3}
$$
The $x,y$ parameters multiply the two scalar denominators. In $\overline{\rm MS}$ conventions the bubble is $1/\bar\epsilon-\int_0^1 dx\,\log[(ax+b(1-x))/\mu^2]$ before subtraction, where $1/\bar\epsilon=1/\epsilon-\gamma_E+\log 4\pi$. Only this term has a UV pole. The finite reduced subtraction $\widehat K(s)=K(s)-K(0)$ compares form factors after factoring out the external spin tensor. Its subtraction point is the analytic limit of the form factor. At $s=0$ the physical spin tensor itself vanishes.

## Cut and real part

For $b>0$, $d(x,y)\ge4a xy$ throughout the simplex, since $4xy\le(x+y)^2\le x+y$. Equality at threshold occurs only at $x=y=1/2$, where the sterile Feynman parameter is zero. Write $\rho(s)=\operatorname{Disc}K(s)/(2\pi i)$, with upper lip minus lower lip. At $s>4a$, the imaginary part of $1/(d-sxy-i0)$ is $+\pi\delta(d-sxy)$. For fixed $x$, the denominator starts at $b+(a-b)x>0$ and ends at $a-sx(1-x)$ when $y=1-x$. Its zero lies on the integration segment precisely for $x_-\le x\le x_+$, where $x_\pm=(1\pm\beta)/2$ and $\beta=\sqrt{1-4a/s}$. Integrating the delta gives
$$
 \rho(s)=(s/2-a)\int_{x_-}^{x_+}\frac{dx}{b-a+sx}
     =\frac{s/2-a}{s}\log\frac{b-a+s(1+\beta)/2}{b-a+s(1-\beta)/2}.                 \tag{4}
$$
Below threshold $\rho=0$. Both logarithm arguments are positive on the physical cut.

An on-shell phase-space check uses $k^0=\sqrt{s}/2$ and $|\mathbf k|=\sqrt{s}\,\beta/2$. Put $h=s/2-a$, $v=s\beta/2$, and $u=\widehat{\mathbf k}\cdot\widehat{\mathbf p}_1$. With $\cos\chi=\widehat{\mathbf k}\cdot\widehat{\mathbf p}_3$, the two numerators are $-h\pm v\cos\chi$. The exchange denominator is $-b-h+vu<0$. The cut bubble has density $\beta$ in the same normalization, so the scalar cut gives
$$
 \rho=\beta\int_{-1}^1\frac{du}{2}\frac{h}{b+h-vu}.
$$
This integral equals (4). At $a=b=1$, $s=8$ and external $\cos\theta=0.6$, the separately integrated routed densities are $0.5564701068$ and $0.7655902738$. Their mean is $0.6610301903$. Replacing $N_j$ by the cut scalar mass $a$ would even reverse the cut sign in this projection.

The finite real and imaginary parts have a common once-subtracted representation,
$$
 \widehat K(s)=s\int_{4a}^\infty
              \frac{\rho(\sigma)\,d\sigma}{\sigma(\sigma-s-i0)}.                       \tag{5}
$$
To verify the subtraction without assuming its high-energy behavior, subtract (3) at $s=0$ inside the simplex. The integrand is $s(d/2-a xy)/[d(d-sxy-i0)]$ on the upper lip. Inserting $\sigma=d/(xy)$ in the delta-function form of (4) gives that same integrand in (5). Endpoint sets with $xy=0$ have zero measure. The dispersion integral converges because (4) grows at most logarithmically. Its above-threshold real part uses the Cauchy principal value.

At $s=4a$ the real part is finite and the discontinuity vanishes. In (3) the upper endpoint of the $y$ integral gives $a(2x-1)^2$. Splitting the remaining logarithmic integral at $x=1/2$ resolves the coincident roots. In (5) the substitution $4a/\sigma=1-u^2$ removes the inverse-square-root endpoint in the threshold integral. For $b=a$, direct integration gives $T(s)=2\arcsin^2\sqrt{s/(4a)}/s$ on $0<s\le4a$, so $\widehat K(4a)=1/2+\pi^2/8$. Both numerical integrals give $1.73370055013617$ at $a=b=1$ and $0.366535322282412$ at $a=1,b=16$.

## An exact UV cancellation

From the angular integral, $\rho>0$ and $\partial_b\rho<0$ for all $s>4a$. Near the common threshold,
$$
 \rho(s;a,b)=\beta\frac{a}{a+b}+O(\beta^3).
$$
Thus two nonzero, opposite weights at distinct masses cannot cancel this projected cut. Take $C_{ee}=1$ in the chosen units, $a=1$, sterile masses $(1,4)$ and conjugated muon-row Yukawas $(1,i/2)$. The weights are $(1,-1)$ and the primitive UV pole is zero. The difference $\Delta\rho=\rho(s;1,1)-\rho(s;1,16)$ is strictly positive above threshold and starts at $(15/34)\beta$. At $s=8$ the reduced results are
$$
 \Delta\rho=0.548545873956202579,\qquad
 \Delta\widehat K=1.189985107065383553+1.723307687777798713\,i.
$$
Here $\operatorname{Im}\Delta\widehat K=\pi\Delta\rho$. The mass-difference notation in these numbers does not include the coupling, ordered external spinor, or prefactor in (1). Three or more mass groups can have exact cut cancellations even with unequal masses. Reciprocal squared masses relative to $a$ give a precise example. [Mass-identifiability calculation](MASS_IDENTIFIABILITY.md)

A finite local four-lepton counterterm can shift the absolute form factor, and local derivative operators can shift polynomial terms in a larger EFT. Neither supplies a two-scalar discontinuity. This calculation fixes the discontinuity and subtracted real part of the specified primitive in the original action coordinates. A total scattering amplitude or rate would require the other contributions at the chosen EFT order. The scalar mass in (3)–(5) is an unbroken-phase input, so the location $4a$ has no electroweak-broken-phase interpretation here.

Run `python -I -B -O verify_on_shell.py` in this directory. The checker evaluates the upper-lip parameter integral and its principal-value dispersion integral separately, integrates the two routed phase-space cuts, and tests the literal normalization with explicit fermionic state operations. [Checker](verify_on_shell.py) · [Numerical receipt](ON_SHELL_TRIANGLE_RESULTS.json)
