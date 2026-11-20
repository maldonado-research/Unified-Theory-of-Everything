# Four-dimensional angular audit of the finite four-lepton kernel

19 November 2026. This AI-assisted internal check starts from the ordered
Euclidean integrand of the retained-field triangle with one derivative seed
and two conjugate neutrino Yukawas. The barred seed leptons have zero
momentum. The unbarred Yukawa leptons carry opposite momenta with Euclidean
square $s=Q^2$. The scalar squared mass is $a>0$. In a Takagi mass
eigenstate write $b=m^2$ with $m\geq0$. All momenta and masses enter the
formulas below through these squared quantities except where indicated.

## Direct angular integration

Write $t=k^2$, and denote the subtracted ordered integral by
$F(a,b,s)/(16\pi^2)$, taking the value at momentum $p$ minus its value
at zero. The four-dimensional radial measure after averaging over the unit
three-sphere is $t\,dt/(16\pi^2)$. With polar angle measured from $p$,
the normalized angular density is $2\sin^2\theta\,d\theta/\pi$. Set
$A=t+s+b$, $B=2\sqrt{ts}$. Polynomial division in $\cos\theta$
and the elementary integral
$\int_0^\pi d\theta/(A+B\cos\theta)=\pi/\sqrt{A^2-B^2}$
give the angular average of the shifted propagator as
$2/(A+\sqrt{A^2-B^2})$. Continuity gives the endpoint value at $A=B>0$.
Consequently

```text
F = ∫₀∞ dt t²/(t+a)² [2/(t+s+b+√D) − 1/(t+b)],
D = (t+b−s)²+4bs.
```

The ordered difference is absolutely integrable in the UV even before averaging.
At large $k$ it falls as $O(|k|^{-5})$. For $b=0$ the isolated
massless-propagator singularities are locally integrable in four dimensions.
For $b,s>0$ the bracket has numerator
$t+b-s-\sqrt D<0$ when written over
$(t+b)(t+s+b+\sqrt D)$. At $b=0$ the shifted average is
$1/\max(t,s)$. These identities fix the strict negative sign of $F$
for every $s>0$, including $b=0$. A change of integration variable
shows that rescaling $a,b,s$ by the same positive constant preserves $F$.

For a second evaluation, split the numerator by
$k^2/(k^2+a)^2=(1+a\partial_a)/(k^2+a)$. The corresponding subtracted
bubble equals $-\int_0^1 dx\log(1+us/d)$ in the same normalization,
where $u=x(1-x)$ and $d=xa+(1-x)b$. Differentiation gives

```text
F = ∫₀¹ dx [−log(1+us/d) + ax us/{d(d+us)}].
```

The checker compares the two one-dimensional integrals with 48 decimal
digits of working precision. For example $F(2,3,7)=-0.2412060116194333$.
The angular identity is also tested by direct angular quadrature, including
the massless point $t=s$.

## Limits at small momentum and special masses

For $b>0$, differentiating the angular expression at $s=0$ gives

```text
F′(0) = −b ∫₀∞ dt t²/[(t+a)²(t+b)³]
      = b(5a+b)/[2(a−b)³] + ab(a+2b)log(b/a)/(a−b)⁴     (a≠b).
```

The slope is strictly negative. At equal masses the integral directly gives
$F'(0)=-1/(12a)$. At $b=a$ the parameter weight $x$ may be replaced
by $1/2$ against a symmetric function of $x(1-x)$. Performing the
integral with $r=\sqrt{1+4a/s}$ gives

```text
F(a,a,s) = 5/2 − [2r+2a/(sr)] atanh(1/r)
         = −s/(12a) + s³/(840a³) + O((s/a)⁴).
```

At $b=0$, the radial difference vanishes above $t=s$. Direct
integration over the remaining finite range gives

```text
F(a,0,s) = 2 − (1+2a/s) log(1+s/a)
         = −s²/(6a²) + s³/(6a³) + O((s/a)⁴).
```

This scalar kernel is continuous at $b=0$, as can also be seen by
dominated convergence in the parameter integral using $0\leq xa/d\leq1$
and $0\leq us/d\leq s/a$. The physical sterile numerator multiplies
it by $m$, so an exactly zero Takagi mode contributes zero to the graph.

## Literal action and arbitrary constant complex symmetric masses

With $M=W\,\mathrm{diag}(m_i)W^T$, form $H=MM^\dagger$. The full
anomalous Euclidean mass factor for a sterile momentum with square $z>0$
is $(z+H)^{-1}M$. Substituting the Takagi decomposition gives
$W\,\mathrm{diag}(m_i/(z+m_i^2))W^T$. Integration replaces this by
the symmetric matrix $K/(16\pi^2)$, where

```text
K = F(a,H,s) M = M F(a,M†M,s)
  = W diag[mᵢ F(a,mᵢ²,s)] Wᵀ,
R = Y* K Y†.
```

Here $F(a,H,s)$ is the spectral function of the Hermitian matrix $H$.
This ordering makes the formula unambiguous at a singular mass and inside
a degenerate eigenspace. Under passive rotations
$M'=U^\dagger M U^*,\ Y'=V^\dagger YU^*$, the transformations are
$K'=U^\dagger K U^*$, $R'=V^T R V$. The barred slots transform with
$C'=V^\dagger C V^*$. The checker constructs $F(a,H,s)$ independently
using the spectral projectors of $H$ and compares to the Takagi sum for
full rank, rank one, repeated positive masses, and zero.

The action phases follow from the literal one-assignment factors
$(-i/2)(-i)^2 i^2(-i)(-1)/2!=+1/4$. The last minus comes from the two
sterile endpoints straddling one external Grassmann lepton. The numerator
minus and the three denominator minuses cancel under Wick rotation, leaving the common factor $i$. The two scalar assignments are
equal after symmetric flavor sums and $p\mapsto-p$. Removing the common
amplitude phase gives the subtracted ordered action coefficient
$+C_{\alpha\beta}R_{\gamma\delta}/(32\pi^2)$ multiplying the ordered
spinor product $P_{\alpha\beta\gamma\delta}$ of the pole calculation,
for the displayed seed orientation. This sign refers to the $p$ minus
zero action, with the negative real kernel $F$ already in $R$.
[Literal action and ordered spinor product](../../checkpoint_2026_10_23/seed_yukawa/FOUR_LEPTON_POLE.md)

## Conjugate seed at the same labelled momenta

Keep the barred leptons soft and the unbarred leptons hard as above. The
conjugate seed now carries the two hard leptons, and each conjugate-graph
Yukawa vertex carries a soft barred lepton. Momentum conservation gives
the same loop momentum $k$ on both scalar lines and the sterile line. The
two seed assignments have derivatives $(k+p)^2$ and $(k-p)^2$ in Euclidean
signature. Subtract at $p=0$ before integrating. The resulting odd term
proportional to $k\cdot p$ integrates to zero for either routing. Wick
rotation of each quadratic numerator supplies the same minus as in the
displayed orientation, cancelling the common sign of the three denominators.

Denote the mirror kernel by $G(a,b,s)=sB(a,b)$. Including the angular
measure leaves

~~~text
B = ∫₀∞ dt t/[(t+a)²(t+b)] = ∫₀¹ dx x/[xa+(1−x)b]
  = [a−b + b log(b/a)]/(a−b)² (>0),     b>0 and a≠b,
B(a,a) = 1/(2a),      B(a,0) = 1/a.
~~~

The parameter equality follows by combining the scalar double pole with
the sterile simple pole and performing the elementary radial integral.
Hermitian conjugation of either two-spinor epsilon bilinear gives a minus
when reversing its odd fields. Their product gives a plus and interchanges
the two flavor pairs. Define $R_G$ by substituting $G$ for $F$ in $R$.
Then the contribution from both seed orientations to the coefficient of
the ordered $P_{\alpha\beta\gamma\delta}$ at these labelled momenta is

~~~text
[C_(αβ) R_(γδ) + (R_G)_(αβ)* C_(γδ)*]/(32π²).
~~~

The spin check compares explicitly sorted Grassmann monomials with
distinct barred and unbarred labels before and after the adjoint. The
radial integral checks $B$ independently of the closed expression. This
mirror term is a linear polynomial in $s$ on this slice and has no
singularity at a positive-mass threshold.

## Distinct mass groups

The separate proposal to distinguish positive masses by their thresholds
can be checked within the bubble representation. Write its rational term
as $ax/d-ax/(d+us)$. For $m>0$ put $h=\sqrt a+m$, $x_*=m/h$
and $u_*=x_*(1-x_*)$. At $s=-h^2+\delta$,

```text
d + us = h²(x−x*)² + δx(1−x).
```

The logarithm has an integrable limit. Rescaling the distance from $x_*$
by $\sqrt\delta$ in the pole term gives
$\lim_{\delta\downarrow0}\sqrt\delta F=-\pi ax_*/(h\sqrt{u_*})$.
The checker compares this coefficient against near-threshold integrals at
two positive masses. Larger positive masses give analytic functions at
this first threshold, since their parameter denominators remain strictly
positive. In a finite linear combination vanishing for all $s>0$,
analytic continuation and successive limits thus force every distinct
positive mass group's coefficient to vanish. This applies componentwise
to the ordered source, before adding its Hermitian conjugate.

A model with one active flavor with $Y^*=(1,i/2)$, $m=(1,4)$ and $a=1$
has zero local mass moment $1+(i/2)^2\,4=0$. Its finite scalar combination
is $F(1,1,s)-F(1,16,s)$. The first slope is $-1/12$, whereas the second
has absolute value below $1/32$, by bounding $t^2/(t+a)^2<1$ in the
radial slope integral. A real unit seed gives a nonzero ordered action
slope. Adding the mirror adds $1/2-(16\log16-15)/225$ to this difference.
The combined normalized slope is positive since $\log16<3$ and
$\partial_s F(1,16,0)<0$. Thus the sampled finite effect survives in both
orientations together for this choice.

## Interpretation and reproduction

The subtraction removes a momentum-independent local term in this ordered
graph for the chosen fields, action and routing. Finite local derivative
interactions and field redefinitions can alter the polynomial momentum
dependence in an off-shell effective action. Physical matching or a
scattering amplitude would require the other relevant diagrams, operators
and external-state treatment. The positive scalar parameter here enters
the finite function explicitly. Keeping it solely as an auxiliary IR mass
would describe a regulated kernel.

Run `python -I -B -O independent_finite.py` in this directory using the
root pinned dependencies. The program imports no repository code and
raises explicit errors on failed checks. The receipt records 72 checks
and eight mutation controls. These are analytic and numerical internal
checks of this one-loop contribution.
[Checker](independent_finite.py) · [Receipt](INDEPENDENT_FINITE_RESULTS.json)
