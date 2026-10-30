# Independent primitive four-lepton pole from the derivative seed

23 October 2026, Pacific. Input candidate
`f1f9a892e15545c70e4fa05827ea7fefc5a438eb`.

This audit was derived before reading the producer's coefficient. It uses a
large-momentum expansion with a common auxiliary infrared mass and a radial
Gamma-function residue, rather than the producer's Feynman-parameter
integration. The overall action factor is obtained independently from a
scalar Nambu determinant and checked against explicit Wick contractions.
The executable imports no producer or other repository code. These are
AI-assisted internal derivations and checks, not independent human peer review.

## Fixed action and result

Keep the [literal reference conventions](../../checkpoint_2026_10_15/operator_map/ACTION_SIGN.md),
with `eta=(+---)`, `d=4-2epsilon`, `S=T exp(i integral L_int)`,
`kappa=1/(16 pi^2 epsilon)` and `+i0` propagators. The interactions are

```text
L_seed=+(1/2) C_alpha beta (bar L_alpha^a phi^a)
                      (-Box)(phi^b L_beta^{c,b}) + h.c.,
L_Y= -bar L_alpha^a Y_alpha i phi^a N_Ri + h.c.,
L_M= -(bar N_L M N_R+h.c.)/2,
C=C^T, M=M^T, phi=tilde H, N_R=N_L^c.
```

`M` is constant and may be complex or singular. Other interactions and
higher-operator insertions are absent in the computed monomial. Exactly
one displayed seed and two conjugate Yukawa vertices occur. The internal
lines are two scalars and one anomalous sterile propagator; the four
external fermions are active leptons. The full scalar derivatives in the
seed are retained.

Use chiral two-component fields `l`, their row adjoints `bar l`, and
`epsilon=[[0,1],[-1,0]]`. Define the ordered local polynomial

```text
R_gamma delta=(Y* M Y^dagger)_gamma delta = R_delta gamma,

P_alpha beta gamma delta
  = sum_(a,b) (bar l_alpha^a epsilon bar l_beta^{b,T})
              (l_gamma^{a,T} epsilon l_delta^b).
```

Every repeated flavor index is summed; `a,b` are weak-doublet indices.
`Y*` means elementwise complex conjugation. The primitive loop action and
its cancelling counterterm are

```text
L_loop^UV = +kappa C_alpha beta R_gamma delta
                    P_alpha beta gamma delta/2 + h.c.,
L_ct      = -kappa C_alpha beta R_gamma delta
                    P_alpha beta gamma delta/2 + h.c.       (1)
```

There is no additional Majorana half or free weak-doublet multiplicity.
The Hermitian-conjugate seed gives the conjugate orientation, which must
be included once as indicated.

For direct comparison in four-component notation,

```text
P_alpha beta gamma delta
 = -sum_(a,b) (bar L_alpha^a L_beta^{c,b})
               (bar L_gamma^{c,a} L_delta^b).
```

In the four-dimensional current projection, with
`J_alpha gamma^mu=sum_a bar L_alpha^a gamma^mu L_gamma^a`, the actual
Grassmann Fierz relation is

```text
J_alpha gamma^mu J_beta delta,mu = -2 P_alpha beta gamma delta,
L_ct = +kappa C_alpha beta R_gamma delta
                  J_alpha gamma^mu J_beta delta,mu/4 + h.c. (2)
```

The minus in the Fierz relation includes the interchange of an unbarred
and a barred external fermion. Equation (2) is a four-dimensional chiral
projection of (1). No full dimensional/evanescent basis, associated finite
scheme conversion or two-loop use of that Fierz relation is supplied.

## Anomalous sterile numerator from the actual chiral inverse

In the chiral gamma basis,

```text
C_spin=i gamma^2 gamma^0=diag(epsilon,-epsilon),
L=(l,0), bar L=(0,bar l),
L^c=(0,-epsilon bar l^T),
bar L_alpha L_beta^c=-bar l_alpha epsilon bar l_beta^T.
```

Thus the displayed seed has Weyl coefficient
`-C(bar l epsilon bar l)phi phi/2`. Write `n=N_L`. The conjugate Yukawa
vertex is `-Y* phi^dagger n^T epsilon l`, equivalent to
`-Y* phi^dagger l^T epsilon n` for anticommuting fields. The left-holomorphic
sterile mass is `M*`, as fixed by the reference four-component action.

The inverse kernel is `slash q-M P_R-M^dagger P_L`. The same-chiral
Majorana contraction is obtained by taking the left-left part of its
inverse and multiplying by the actual `C_spin^T`. The leading slash term
vanishes between those endpoints. Its first mass term gives

```text
<n_i n_j>_UV = -i M_ij epsilon/q^2 + O(q^-4).            (3)
```

The checker constructs this term from
`slash q (M P_R+M^dagger P_L) slash q` before projecting; it does not
assign its complex mass numerator by analogy. For a real diagonal mass
the exact denominator in (3) is `q^2-M_i^2`. At zero momentum its value
is `+i epsilon/M_i`, agreeing with the two-variable Gaussian for the
mass action `-M_i* n^T epsilon n/2`. Neither that calibration nor the UV
computation requires an inverse for the general mass matrix.

Contracting the two Yukawa coefficients with (3) gives `R=Y* M Y^dagger`.
Replacing `M` by `M^dagger` is inconsistent with this orientation. Sterile
mass-squared insertions in the propagator decrease the UV degree of this
already logarithmic graph and do not contribute another local pole here.

## Independent homogeneous UV integration with IR separation

Let the three denominators have generic off-shell shifts and physical
masses. A routing of the differentiated neutral composite gives numerator
`(q+r)^2`, where `r` is a fixed external momentum. With loop momentum
scaled by `t`,

```text
(t q+r)^2 = t^2 q^2 + 2t q.r + r^2,
each quadratic denominator = t^2 q^2 + O(t).
```

The only potentially logarithmic term is therefore `q^2/(q^2)^3`.
The linear and constant numerator terms have UV degrees at most `-1`
and `-2`; denominator shifts and mass insertions lower them further.
The program constructs the actual polynomials for twelve non-null rational
UV directions and generic external shifts and checks their homogeneous
coefficients. This verifies the derivative implementation; the general
degree argument is the reason these tests apply to the whole declared
primitive graph. Applying derivatives only to the external lepton would
incorrectly remove the degree-two term and falsely report a UV zero.

Introduce an independent positive auxiliary mass squared `r_aux` in
`Q=q^2-r_aux+i0` before expanding. The primitive UV residue is

```text
I_UV = integral d^d q/(2pi)^d q^2/Q^3 |_UV = i kappa.    (4)
```

This is evaluated directly by the radial formula

```text
integral (q^2)^a/(q^2-r_aux+i0)^b
 = i(-1)^(a-b)/(4pi)^(d/2) r_aux^(d/2+a-b)
   Gamma(a+d/2) Gamma(b-a-d/2)/[Gamma(d/2) Gamma(b)].
```

At `a=1,b=3`, `Gamma(b-a-d/2)=Gamma(epsilon)` and the remaining ratio
approaches `Gamma(3)/[Gamma(2)Gamma(3)]=1`. The simple-pole coefficient
is independent of `r_aux`. The difference between the expanded and
original denominators is UV convergent by the same degree argument.
The auxiliary mass keeps the expanded vacuum integral infrared safe.
No scaleless integral is set to zero and used to infer absence of a UV
pole. The physical matrix `M` in (3) is not this auxiliary regulator.

## Independent action-factor calculation

Completing the sterile Gaussian with its Yukawa source produces, for
purposes of evaluating this one graph, the nonlocal source kernel

```text
- (1/2) R_gamma delta (phi^dagger l_gamma)^T epsilon
              (phi^dagger l_delta)/(q^2-M_i^2)
```

in a diagonal mass component. Its sign follows from (3) and from ordering
the two external odd sources. The complete mass-dependent kernel need
not be diagonalized to obtain its UV numerator. This intermediate
nonlocal Gaussian calculation is not a local EFT expansion or sterile
threshold matching; the result still concerns the original theory with
retained sterile fields.

Treat the external fermion bilinears as even background coefficients.
For the doubled complex-scalar Hessian, write the anomalous blocks as

```text
A_ab = -C_alpha beta (bar l_alpha^a epsilon bar l_beta^{b,T}) q^2,
B_ab = -R_gamma delta (l_gamma^{a,T} epsilon l_delta^b)/(q^2-M_i^2).
```

Their weak-index symmetry is checked explicitly using symmetric `C,R`
and the exterior algebra. The complex-scalar Gaussian contributes
`Gamma_scalar=i Tr log(D)/2`. The cross term in its expansion is

```text
-i Tr(K^-1 A K^-T B)/2,
```

because `-i Tr[(D0^-1 delta D)^2]/4` has two equal off-diagonal block
traces. Inserting (4) and the two negative anomalous blocks yields
`+kappa C R P/2`, the first line of (1).

A literal Wick calculation independently checks this normalization. For
one scalar assignment its factors are

```text
seed                  -i/2,
two conjugate Yukawas (-i)^2,
two scalar propagators i^2,
anomalous sterile     -i,
selected Wick sign    -1,
S B^2 expansion        3/3! = 1/2.
```

The two selected sterile fields straddle one external odd lepton, giving
the displayed Wick minus. Their epsilon contraction is evaluated
component by component. The product above is `+1/4` before the common
loop pole and coupling tensors. The two scalar Wick assignments are
distinct before symmetric flavor sums; each becomes the same ordered
polynomial afterward. They produce the factor two giving (1).

There is no additional closed-fermion-loop minus: the only internal
fermion propagator joins vertices with external fermion endpoints. The
triangle is not a closed fermion trace. Appending such a sign, another
Majorana half, omitting a scalar assignment, or canceling the seed half
twice fails the independent checks.

## Flavor, weak and identical-field controls

Use the passive coordinates already fixed for the previous complex-mass
calculation:

```text
L_old=V L_new, N_L,old=U N_L,new,
Y'=V^dagger Y U*, M'=U^dagger M U*,
C'=V^dagger C V*, R'=V^T R V.
```

The last relation follows directly from `R=Y* M Y^dagger`. The checker
substitutes the actual rotated exterior fields into (1), not just the
coefficient matrices, and verifies equality with the primed action. It
tests a general complex symmetric mass, a singular complex mass and zero
mass without a mass inverse. It also verifies the full weak-unitary
invariance, Hermitian conjugation of the ordered action and
counterterm cancellation. The sum over each weak index joins an external
pair; no unused weak trace supplies a factor two.

An identical-flavor limit does not force (1) to vanish. Even at one fixed
weak component the ordered monomial is
`(bar l epsilon bar l)(l epsilon l)=4 bar l_1 bar l_2 l_1 l_2`, which is
nonzero. The full one-active-flavor Hermitian action is checked explicitly.
The pole vanishes when `M=0` in this particular primitive sector, as its
anomalous numerator requires the mass. This is not a statement about all
loop sectors or all physical lepton-number-conserving operators.

## Reproduce and interpret

Run `python -I -B -O independent_four_lepton.py` with the root pinned
SymPy environment. The adjacent
[receipt](INDEPENDENT_FOUR_LEPTON_RESULTS.json) records **224 positive
checks and 25 adverse controls**. Of the positive checks, 128 compare
explicit Grassmann scalar/current components; they are not separate
loop integrations. Failures use explicit exceptions, so optimization
does not remove them.

This supports the primitive four-lepton UV action (1), its literal sign
and the declared four-dimensional current projection. Other diagrams,
counterterm insertions, source/field-coordinate transport and lower-order
operator effects can feed the same target operators. The result does not
give a total physical Wilson coefficient, complete anomalous dimension,
finite matching correction, scattering prediction, external novelty or
unification. A complete dimensional/evanescent prescription is still
needed before treating the projected current expression as part of a
physical coefficient vector.
