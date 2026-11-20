# Four-scalar pole from one derivative seed and two Yukawa vertices

23 October 2026 (Pacific). Input candidate:
`f1f9a892e15545c70e4fa05827ea7fefc5a438eb`.

This note evaluates an actual new primitive one-loop 1PI sector: one original
symmetric derivative seed and two conjugate neutrino-Yukawa vertices, with four
external scalar fields and zero external gauge field. The retained Majorana
field remains in the theory. This is an unreduced local effective-action pole,
not finite sterile-threshold matching, a physical Higgs or EFT beta function,
or an observable prediction.

## Result and coefficient conventions

Use the literal action of the earlier
[action-sign note](../../checkpoint_2026_10_15/operator_map/ACTION_SIGN.md):

```text
eta=(+---), d=4-2epsilon, kappa=1/(16 pi^2 epsilon),
phi=tilde H, X=phi^dagger L, P0=-Box,
L_seed=(1/2)bar X_alpha C7_alpha beta P0 X_beta^c+h.c.,
C7=C7^T,
L_Y=-bar L_alpha Y_alpha i phi N_Ri+h.c., N=N^c.
```

The sterile chiral mass action is
`-[bar N_L M N_R+h.c.]/2`, with constant complex symmetric `M`. Define

```text
rho=phi^dagger phi,
A=rho (partial_mu phi)^dagger partial^mu phi,
B=(phi^dagger partial_mu phi)(phi^dagger partial^mu phi),
B1=Tr(C7 Y* M Y^dagger),
B3=Tr(C7 Y* M M^dagger M Y^dagger).
```

The trace is over active flavor after the displayed matrix products. Neither
`B1` nor `B3` is assumed real or positive. The displayed seed orientation gives

```text
L_loop,div = kappa [ B3 rho^2 + B1(B-A) ] + h.c.,          (1)
L_ct       = -kappa [ B3 rho^2 + B1(B-A) ] + h.c.           (2)
```

Both the sign and the normalization include the literal seed half, the two
Yukawa expansion assignments, all active-fermion contractions and the Majorana
propagator normalization. The Hermitian conjugate is added once. In particular,
the potential coefficient in the full action is `2 kappa Re(B3)`, while the
complex derivative term is `B1 B+B1* B^dagger-2 Re(B1) A`. Taking a real part
before adding the conjugate source would lose the possible imaginary component
of `B1` multiplying `B-B^dagger`.

`B1` has mass dimension minus two and multiplies four-scalar dimension-six
terms. `B3` is dimensionless and multiplies a quartic-potential term. Thus a
higher-dimensional insertion generates a renormalizable-potential UV pole in
the presence of the retained mass scale. No mass inverse or sterile-field
elimination was used. There is no internal Higgs propagator in this triangle,
so its raw pole contains no `m_H^2`. Source-preserving EOM reductions can
subsequently rearrange the displayed terms; they have not been applied in (1).

These are zero-gauge local polynomials. Their global weak/unitary invariance
is explicit, but additional gauge amplitudes and the complete coefficient map
are outside this calculation.

## Fully labelled four-point polynomial

Treat `phi` and `phi^dagger` as independent Fourier sources. All momenta are
incoming with modes `exp(-ip.x)`:

```text
phi^a(p1), phi^b(p2), phi_dagger^c(k1), phi_dagger^d(k2),
p1+p2+k1+k2=0,
D=delta_ac delta_bd, X=delta_ad delta_bc,
K=k1+k2, w=k1-k2,
hD=p1+k1, hX=p1+k2,
JD=-K^2/2+hD^2-w.hD/2,
JX=-K^2/2+hX^2+w.hX/2.
```

For the displayed seed orientation, before adding its Hermitian conjugate,

```text
i Gamma_S,div = 2i kappa [B3(D+X)+B1(JD D+JX X)].          (3)
```

The conjugate action exchanges the scalar and conjugate-scalar sources and
conjugates the coefficients. Equation (3) is Bose symmetric within each pair
of equal fields; the full action obtained this way is Hermitian. A fixed set
of scalar indices is used throughout. No free weak-doublet trace supplies an
additional factor two.

For real diagonal retained masses, write

```text
w_i=(Y^dagger C7 Y*)_ii,
B1=sum_i w_i M_i, B3=sum_i w_i M_i^3.
```

The spinor/loop calculation can therefore be performed for one real positive
`M_i` and multiplied by `w_i M_i`, retaining its exact ordering when the
complex-mass dictionary is restored below.

## Literal Weyl/Wick normalization and closed-loop sign

Take chiral gamma matrices with
`C_spin=i gamma^2 gamma^0=diag(epsilon,-epsilon)`, where
`epsilon=[[0,1],[-1,0]]`. Write `L=(l,0)` and `N_L=(n,0)`. Then
`L^c=(0,-epsilon bar l^T)`. The relevant ordered terms are

```text
S = -(1/2) C7 phi phi bar l epsilon P0 bar l^T,
B_Y = -Y* phi_dagger l^T epsilon n,
<n n^T>(q) = -i M epsilon/(q^2-M^2+i0).
```

The last equation follows by projecting
`<N N^T>=i(slash q+M)C_spin^T/(q^2-M^2+i0)`; it is not an additional Majorana
vertex half. The two mixed active propagators in the ordered expectation
`<bar l(x)l(y)>` each have factor `-i`, with their respective Fourier momenta.

Expand `iS (iB_Y)^2/2!`. Before the spin contractions its coefficient is
`(-i/2)(-i)^2/2=i/4`. For the fermion order

```text
bar l_1, bar l_2, l_1, n_1, l_2, n_2,
```

the two connected active matchings have Wick signs plus and minus. The sterile
pair is always `n_1 n_2`. Their explicit spinor contractions are respectively

```text
 epsilon_ab epsilon_cd epsilon_ef epsilon_df
     A_ca B_eb = 2 a.b,
 epsilon_ab epsilon_cd epsilon_ef epsilon_df
     A_cb B_ea = -2 a.b,                                (4)
```

for `A=sigma.a`, `B=sigma.b`. Hence both Wick-weighted contributions have the
same sign. The three ordered contraction factors give `(-i)^3=i`; the first
active line has momentum `-q-k1`, while the second has `q-k2`. Factoring that
first numerator minus leaves

```text
(i/4) i (-1) [ (+1)2+(-1)(-2) ] = +1.                   (5)
```

This is the coefficient of the ordered scalar graph for `i Gamma`
below, before its ordinary loop integral. The latter supplies `i kappa`,
leaving the positive real action coefficient in (1). It includes the closed-fermion
sign through the actual Wick order and anomalous propagator. Appending a
second closed-loop sign or another Majorana half would change (5) incorrectly.
The later scalar functional derivatives generate the factor and weak
pairings in (3); they are not another unexplained graph multiplicity.

The executable verifies every Lorentz coefficient of (4) with explicit Pauli
matrices, both six-fermion matching signs, and the factors in (5).

## Actual Feynman-parameter UV integration

Route sterile momentum `q`. The active denominators are `(q+k1)^2` and
`(q-k2)^2`. For the direct scalar pairing, the derivative in the original seed
acts on a composite of momentum
`-q+k2+p2=-(q+k1+p1)`. Set

```text
r=k1, s=-k2, t=k1+p1=s-p2.
```

The normalized integral whose pole is needed is

```text
J(r,s,t;M)=integral d^d q/(2pi)^d
 ((q+r).(q+s)) (q+t)^2 /
 [(q^2-M^2+i0)((q+r)^2+i0)((q+s)^2+i0)].                 (6)
```

The crossed pairing uses `t=s-p1=-(p1+k2)` instead. The other seed placement
is related by the exchange of its identical fermions. The derivative acts on
the full scalar/lepton composite: replacing `t` by the bare active momentum
would omit the external scalar momentum and change the answer.

With Feynman parameters `a+b+c=1`, `a,b,c>=0`, define

```text
v=b r+c s, q=l-v,
Delta=a M^2-b r^2-c s^2+v^2,
R=r-v, S=s-v, T=t-v.
```

The denominator becomes `(l^2-Delta+i0)^3` with factor two. The even numerator
is

```text
l^4+l^2(T^2+R.S)+2(l.(R+S))(l.T)+(R.S)T^2.
```

The last term multiplies a UV-finite radial integral. Angular integration
uses `l_mu l_nu -> g_mu nu l^2/d`. Since only a simple UV pole is retained,
`2/d` may be replaced by `1/2` in its residue; its epsilon correction affects
finite terms, which are not calculated.

Gamma recursion gives, after factoring `i kappa`,

```text
I1(Delta)|UV=Delta, I2(Delta)|UV=1, I3(Delta)|UV=0,
integral l^2/(l^2-Delta)^3|UV=1,
integral l^4/(l^2-Delta)^3|UV=3Delta.
```

The last identity includes a tadpole term and two `Delta I2` terms. Dropping
the tadpole would lose a physical retained-mass contribution. Thus (6) has

```text
J|UV = i kappa 2 integral_simplex db dc
 [3Delta+T^2+R.S+(R+S).T/2]
     = i kappa [M^2-(r-s)^2/2+t^2-(r+s).t/2].             (7)
```

The program integrates the full abstract polynomial in the six independent
Gram entries `r^2,s^2,t^2,r.s,r.t,s.t`, using exact simplex moments
`integral b^j c^k=j!k!/(j+k+2)!`. It checks every coefficient. Together with
the degree-two UV bound, this proves the local polynomial rather than only
fitting selected momentum values.

For UV/IR separation take nonexceptional spacelike external invariants.
The massive/massless triangle is then IR safe; after Wick rotation `Delta`
is positive in the interior and its endpoint limits are integrable. The
retained `M` is a physical parameter, not an auxiliary IR regulator. No
scaleless integral is set to zero to infer a vanishing UV pole. When `M=0`,
the whole displayed sector vanishes because its anomalous numerator is
zero, not because a scaleless integral was discarded.

## Exact local-action reconstruction

At the scalar pairings `D,X`, the respective Fourier kernels of the two
derivative monomials are

| Monomial | D kernel | X kernel |
| --- | --- | --- |
| `A=rho partial phi_dagger.partial phi` | `-k1.p1-k2.p2` | `-k1.p2-k2.p1` |
| `B=(phi_dagger partial phi)^2` | `-2p1.p2` | `-2p1.p2` |
| `rho^2` | `2` | `2` |

Momentum conservation gives `B_D-A_D=2JD` and `B_X-A_X=2JX` identically. The
potential kernel is the literal fourth functional derivative of `rho^2`.
This proves that (7), including the derivative placement, matches exactly
(1) and (3). Twelve generic off-shell points with all sixteen weak-index
assignments separately check both Bose exchanges and the full Fourier
coefficient including the conjugate action. These are implementation
controls in addition to the complete polynomial identity.

As a constant-background sign/factor control, for one real active and sterile
flavor take the momentum-dependent chiral mass matrix

```text
M_R(q)=[[a,d],[d,M]], a=-c q^2, c=C7 phi^2, d=y phi.
```

The Majorana/Weyl Gaussian action has `-i integral log det(q^2-M_R^2)`.
The derivative of that determinant with respect to `a` at zero is `-2d^2M`.
At order `c d^2`, its logarithm therefore contributes
`2c d^2 M/(q^2-M^2)`. The tadpole pole gives
`+2 kappa C7 y^2 M^3 phi^4`, agreeing with (1) plus its conjugate for real
coefficients. This control fixes the potential normalization and sign; it
does not recompute the off-shell derivative polynomial. A constant scalar
background does not make the seed vanish inside the loop: its derivative
still acts on the internal lepton.

## General complex mass and flavor covariance

The holomorphic left-Weyl mass is `M*` in these conventions. Its anomalous
propagator has the ordered large-momentum expansion

```text
<n n^T>(q) = -i [M/q^2+M M^dagger M/q^4+...] epsilon.
```

The one-mass term controls the external-momentum part of (7); the next term
controls its local potential part. This gives precisely `B1,B3` above.
For a symmetric `M`, `M M^dagger M` is symmetric as required by the anomalous
fermion contraction. The identity
`(M M^dagger)M=M(M^dagger M)` retains the order without commuting `M` and
`M^dagger`.

Use the established passive coordinates

```text
L_old=V L_new, N_R,old=U* N_R,new,
Y'=V^dagger Y U*, C7'=V^dagger C7 V*, M'=U^dagger M U*.
```

The two active source matrices `Y* M Y^dagger` and
`Y* M M^dagger M Y^dagger` transform as `V^T(...)V`, so both traces are
invariant. The checker uses complex unitary rotations and noncommuting
symmetric complex masses, including singular and zero cases. Replacing `M`
by `M*`, or replacing the ordered cubic product by `M^dagger M M`, fails.

## Reproduction and limits

Run `python -I -B -O verify_four_scalar.py` in the root pinned symbolic
environment. The script imports no repository module and uses explicit
failures rather than removable assertions. Its
[receipt](FOUR_SCALAR_RESULTS.json) records **677 positive controls and
13 rejected mutations**. The count includes repeated weak-index and
momentum implementation comparisons, not 677 separate loop calculations.
The independent route in
[the adjacent audit](INDEPENDENT_FOUR_SCALAR.md) uses a different
large-momentum angular/radial integration; it is not imported by this
producer.

The calculation excludes the mixed lepton/scalar sector, other interactions,
other operator insertions, finite terms, on-shell amplitudes and the complete
source/EOM coefficient map. Neither the zero-gauge action nor its quartic
pole should be promoted to an observable or a complete physical Higgs/EFT
running prediction. The resulting counterterms are ingredients for that
remaining calculation.
