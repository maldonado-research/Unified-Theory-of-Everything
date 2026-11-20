# Finite spacelike subtraction of the ordered four-lepton triangle

Consider the four-lepton graph with one derivative-seed insertion, two conjugate neutrino Yukawa vertices, and a retained Majorana propagator. Its ordered UV loop pole is $+C_{\alpha\beta}(Y^* M Y^\dagger)_{\gamma\delta}P_{\alpha\beta\gamma\delta}/(32\pi^2\epsilon)$, using the spinor order and field coordinates of the [October 23 derivation](../../checkpoint_2026_10_23/seed_yukawa/FOUR_LEPTON_POLE.md). That local pole cancels in the difference of the graph at two spacelike momentum configurations. The difference can be evaluated as a finite integral at the chosen scalar and sterile masses.

## Routing and normalization

Use a common positive scalar squared mass $a$. In this calculation it is an input to the propagator. Changing an auxiliary infrared mass would change the finite answer. Specify the four incoming momenta by vertex. Both barred leptons at the seed carry zero momentum. The unbarred leptons at the two Yukawa vertices carry $+p,-p$, with $z=Q^2=-p_M^2\geq0$. A scalar line entering the seed with momentum $k$ leaves it with the same momentum on the other scalar line. Momentum conservation then fixes the sterile line momentum to $k+p$, up to an equivalent reversal. The derivative on the entire second seed composite supplies $k_M^2$. Interchanging the two scalar contractions sends $p\to-p$ and preserves this scalar integral. External momenta label Grassmann sources and can be set off shell independently of the spinor basis.

For a Takagi mass $m\geq0$, set $b=m^2$ and strip the factor $m$ coming from the anomalous propagator. Define, using a common translation-invariant regulator for the two logarithmically divergent integrals,

$$
 I_E(z)=\int\frac{d^4 k}{(2\pi)^4}
  \frac{k^2}{(k^2+a)^2((k+p)^2+b)},\qquad
 F(a,b,z)=16\pi^2\,[I_E(z)-I_E(0)].
$$

Wick rotation of the Minkowski integral gives $+i I_E$. The Wick signs of the numerator and the three denominators cancel. The positive loop-action sign, seed half and two Wick assignments are the same as in the ordered pole above. The subtracted integrand falls at least as $k^{-5}$ at large Euclidean momentum. At $b=0$, the shifted massless singularity is locally integrable in four dimensions and the soft origin is integrable as well.

Write $M=U\,\mathrm{diag}(m_j)\,U^T$ with $U$ unitary in the prior left-sterile convention. The ordered difference multiplying $P_{\alpha\beta\gamma\delta}$ in the loop effective action is

$$
 \frac{C_{\alpha\beta}}{32\pi^2}
 \left[Y^* U\,\mathrm{diag}\{m_j F(a,m_j^2,z)\}\,U^T Y^\dagger\right]_{\gamma\delta}.
 \tag{1}
$$

Indeed the Euclidean mass factor before integration is
$(q^2+M M^\dagger)^{-1}M=U\,\mathrm{diag}\{m_j/(q^2+m_j^2)\}\,U^T$. Equal-mass sterile subspaces have equal weights, so (1) is independent of the Takagi choice within them. A massless subspace contributes zero. The conjugate action term conjugates the fields and the momentum assignment together. At fixed external barred and unbarred momenta a full Green function requires its other orientations and permutations on that same configuration. Doubling (1) at this asymmetric routing would count a different quantity.

## Integral and sign

The numerator reduces the triangle to a bubble and its mass derivative,

$$
 \frac{k^2}{(k^2+a)^2}=\left(1+a\,\partial_a\right)\frac{1}{k^2+a}.
$$

Let $x\in[0,1]$ be the Feynman weight of the scalar denominator, $D=ax+b(1-x)$, and $c=x(1-x)$. Completing the square gives $D+cz$ as the shifted mass squared. The regulated bubble difference in units of $1/(16\pi^2)$ is $-\int_0^1\log(1+cz/D)\,dx$. Acting with the indicated derivative gives

$$
 F(a,b,z)=\int_0^1 dx\left[-\log(1+cz/D)
                   +\frac{ax}{D}\frac{cz}{D+cz}\right]. \tag{2}
$$

At $b=0$, take the endpoint value by continuity. A direct calculation with the doubled scalar line gives the integrand $2x[-\log(1+cz/D)+(1-x)^2z/(2(D+cz))]$. It differs from the integrand of (2) by $\partial_x[c\log(1+cz/D)]$, with zero endpoint contribution.

The subtracted kernel is zero at the origin. Differentiation gives a sign test over the whole spacelike interval,

$$
 \partial_z F=-\int_0^1 dx\;
       \frac{c\,[b(1-x)+cz]}{(D+cz)^2}. \tag{3}
$$

It is strictly negative at $z>0$ for every $a>0, b\geq0$. For $b>0$ it is also negative at $z=0$. This sign belongs to the scalar kernel. The complex flavor tensor and additional graph contributions enter separately.

## Slope and mass limits

For $r=b/a>0$, the slope has the exact form

$$
 F'(0)=-\frac{r}{a}\int_0^1\frac{(1-y)y^2}{[1+(r-1)y]^2}\,dy
       =-\frac{r}{a}\frac{r^2+4r-5-(4r+2)\log r}{2(r-1)^4}.
$$

The integral is continuous at $r=1$, where $F'(0)=-1/(12a)$. For a stable numerical evaluation near equal masses use

$$
 F'(0)=-\frac{r}{a}\sum_{n=0}^\infty
      \frac{(-1)^n(n+1)(r-1)^n}{(n+3)(n+4)},\qquad |r-1|<1.
$$

In the massless sterile limit, (2) reduces to $F=2-(1+2/s)\log(1+s)$ with $s=z/a>0$. Its series begins at $-s^2/6$. The mass factor $m_j$ in (1) suppresses the entire corresponding eigenchannel, including at nonzero momentum.

For $a=b=2, z=3$ in common squared-mass units, (2) gives $-0.122398886130692358055401$, and $F'(0)=-1/24$ in reciprocal units. These values are a reproducible calculation point.

## Conjugate orientation at the same labelled momenta

Keep the barred external momenta at zero and the unbarred momenta at $+p,-p$. In the conjugate graph the unbarred leptons attach to the conjugate seed. The two Yukawa slots are now soft, so all three propagators carry the same loop momentum up to a reversal. The seed composite supplies $(k\pm p)^2$ in Euclidean units. After subtracting the zero-momentum graph, the term odd in $k$ integrates to zero and the numerator leaves $+z$. The Wick-rotation signs are the same as for (1). The normalized mirror kernel $H(a,b,z)$ is exactly linear in $z$ at this routing,

$$
 H=16\pi^2 z\int\frac{d^4k}{(2\pi)^4}\frac{1}{(k^2+a)^2(k^2+b)}
   =z h(a,b),\qquad h(a,b)=\int_0^1 \frac{x\,dx}{D}.
$$

For $b>0,\,b\ne a$ the slope is $h=(1-r+r\log r)/[a(1-r)^2]$. The continuous limits are $1/(2a)$ at $b=a$ and $1/a$ at $b=0$. Near $r=1$ use $a h=\sum_{n\geq0}(1-r)^n/[(n+1)(n+2)]$ for $|r-1|<1$.

Define $R[F]$ as the flavor matrix in brackets in (1), and $R[H]$ by substituting $H$ for $F$. In the ordered $\bar B_{\alpha\beta} B_{\gamma\delta}$ convention of the [spinor-pair derivation](../../checkpoint_2026_10_23/seed_yukawa/GREEN_MIXING_AND_SOURCES.md), the two orientations contribute the numerator

$$
 C_{\alpha\beta}R[F]_{\gamma\delta}
       +R[H]^*_{\alpha\beta}\,C^*_{\gamma\delta}
$$

over $32\pi^2$. The spin contraction is the antisymmetric two-component epsilon. Anticommuting spinors make each bilinear symmetric under exchange of its two **complete** field labels, including momentum and weak component. Adjointing one bilinear contributes a minus sign; adjointing their product cancels the two signs and interchanges the pair slots. Both scalar contractions have the same radial kernel under $p\to-p$; the soft pair contracted with a symmetric flavor matrix is weak-symmetric. No extra flavor multiplicity is introduced. Hermiticity interchanges the two momentum assignments, so the matrix at this single asymmetric assignment need not be Hermitian.

There is an exact sign check for the paired **scalar kernels**. With $\Delta=D+cz$, the integrand of $h+\partial_zF$ can be written as

$$
 \frac{x^2}{D}+\frac{c^2z}{D\Delta}+\frac{cax}{\Delta^2}.
$$

It is strictly positive in the interior of the parameter interval for $a>0$, $b\geq0$. Consequently $0<-F<H$ at $z>0$, and $F+H$ grows strictly from zero. For equal masses the slope of $F+H$ at the origin is $5/(12a)$. For real one-flavor inputs the same real flavor factor multiplies both kernels. Complex flavor tensors and additional contributions preclude reading this scalar inequality as a positivity bound on a physical amplitude. At $a=b=2,\,z=3$ the mirror gives $H=0.75$, while the direct kernel is negative. A zero sterile eigenvalue has a finite mirror kernel and still contributes zero after the mass factor.

## Four-dimensional check and scope

In the finite integral the four-dimensional angular average of the shifted propagator, with $u=k^2$, is

$$
 A(u,b,z)=\frac{2}{u+b+z+\sqrt{(u+b-z)^2+4bz}}.
$$
Thus $F=\int_0^\infty du\,u^2(u+a)^{-2}[A(u,b,z)-1/(u+b)]$. In the massless case the angular difference has support $u<z$. This verifies the scalar kernel independently of the Feynman parameters, including the near-degenerate case.

The [standalone checker](verify_finite_triangle.py) compares the parameter, radial and closed slope representations, verifies the matrix propagator order, and rejects sign, numerator and conjugation mutations. Its [receipt](FINITE_TRIANGLE_RESULTS.json) includes the mass and momentum values used in each comparison.

The subtraction removes this graph’s momentum-independent local ambiguity. The ordered spinor structure may be projected in four dimensions after the integral difference is finite. Independent higher-derivative contacts or changes of field coordinates can affect the same off-shell momentum dependence. A physical amplitude still needs the remaining diagrams, external states, renormalization and matching conditions. The calculation evaluates two ordered orientations of the selected primitive graph at the stated Euclidean kinematics.
