# Independent scalar-current invariant and source audit

16 October 2026. Independently derived from the scalar currents in the actual
[`ELECTROWEAK_REDUCTION.md`](../../checkpoint_2026_10_15/interacting_reduction/ELECTROWEAK_REDUCTION.md)
and [`CURRENT_DESCENDANTS.md`](../../checkpoint_2026_10_15/interacting_reduction/CURRENT_DESCENDANTS.md).
The independent derivation and executable were prepared before reading the
new scalar producer code; comparison with its final derivation finds no sign
or coefficient discrepancy. This is AI-assisted internal verification, not
external peer review or a further loop integration.

## Independent weak-invariant classification

Use the same `D=partial+i Acal`, `phi=tilde H`, `T=tau/2`,
`y_phi=y_L=-1/2`, and `F=C_seed Y_nu*`. Write

```
b_a^mu=F_(alpha i)bar L_alpha,a gamma^mu N_Li,
j^mu=b_a^mu phi_a,       rho=phi^dagger phi,
A=rho b_mu D^mu phi,
B=j_mu(phi^dagger D^mu phi),
D=j_mu((D^mu phi)^dagger phi).
```

Here `D` without an index in the third monomial is its label, not a new
derivative convention; it is producer `O_3`. This audit uses `A,B,D` to
keep its independently written formulas distinct. Their sum identity is
`B+D=j_mu partial^mu rho`.

For a fixed Lorentz component first regard `b` as an arbitrary even weak
row. A general contraction with derivative on `phi` has 16 coefficients
on `b_a phi_b phi_c^dagger (Dphi)_d`. Apply the three complexified weak
generators `E12,E21,diag(1,-1)` independently to every scalar and row slot.
Their stacked coefficient constraints have exact rank 14, giving two
invariants. `A` and `B` are both in this kernel.

With the derivative on `phi^dagger`, the two undifferentiated commuting
scalar columns are symmetric. The monomials
`b_a phi_b phi_c (Dphi)^dagger_d`, `b<=c`, have 12 coefficients.
The corresponding constraints have rank 11, giving one invariant, `D`.
The hypercharge sum vanishes in both modules. No derivative/EOM relation
has yet been applied in these algebraic ranks.

The checker additionally replaces `b` by its actual four-dimensional Weyl
bilinear and exhibits a diagonal three-component minor. Choose the odd
generators `bar L_(0,weak0,spin0)` and `N_(0,spin0)`. Scalar/jet components
that distinguish `A,B,D` are respectively

```
phi_1^dagger phi_1 (D_0 phi)_0,
phi_1^dagger phi_0 (D_0 phi)_1,
phi_0 phi_1 (D_0 phi)^dagger_1.
```

Each appears only in its designated monomial, with coefficient `F_00`.
The exact witness uses `F_00=1/2+2i`, giving a nonzero determinant
`F_00^3`. This is a component/exterior-algebra independence witness as well
as a weak-tensor nullity calculation. These ranks are restricted to this
scalar-jet placement sector, not full EFT EOM closure.

## Literal scalar coefficients and full integration by parts

The scalar pieces of the actual variational currents are

```
J_Y^mu=i y_phi[phi^dagger D^mu phi-(D^mu phi)^dagger phi],
J_I^mu=i phi^dagger T^I D^mu phi-i(D^mu phi)^dagger T^I phi.
```

Direct Pauli component expansion, with no weak or fermion equation, gives

```
j_I,mu J_I^mu=i A/2-i(B+D)/4,
P_s=-(g'^2 y_L j_mu J_Y^mu+g^2 j_I,mu J_I^mu)/6
   =i[-2g^2 A+(g^2-g'^2)B+(g^2+g'^2)D]/24.
```

The positive product `y_L y_phi=1/4` is essential for the hypercharge
sign. Let `d=(g^2+g'^2)/24`. Use the exact definitions
`A_N=i slash partial` and `K bar L=-i(D_mu bar L)gamma^mu`. Differentiating
the **whole** neutral current gives

```
partial_mu j^mu=b_mu D^mu phi
       +i[(K bar L)phi N_L]_F-i[bar L phi A_N N_L]_F.
```

The independent checker constructs ordinary derivatives from arbitrary
covariant jets and a completely generic connection matrix:
`partial phi=Dphi-i Acal phi`,
`partial phi^dagger=(Dphi)^dagger+i phi^dagger Acal`, and
`partial bar L=Dbar L+i bar L Acal`.
Their connection terms cancel after the actual ordered weak contractions.
It does not assume a constant scalar or set the connection to zero.

Using `B+D=j partial rho` and the full product rule yields the off-shell
identity

```
P_s= -i(3g^2+g'^2)A/24-i g'^2 B/12
     +d rho([(K bar L)phi N_L]_F-[bar L phi A_N N_L]_F)
     +i d partial_mu(rho j^mu).                         (S1)
```

It would lose an actual derivative term to move the norm derivative onto
only `N_L`, or to drop the `Dphi` term inside the current divergence.
Both mistakes are rejected by the independent implementation.

## Nonzero sourced equations and an off-shell certificate

Keep the actual interacting residuals

```
e_barL=K bar L-bar S_L,
e_R=A_N N_L-M N_R-S_R.
```

Substitution in (S1), without deleting the residuals, gives

```
P_s,red= -i(3g^2+g'^2)A/24-i g'^2 B/12
       +d rho([bar S_L phi N_L]_F
                        -[bar L phi(M N_R+S_R)]_F),
P_s-P_s,red=i d partial_mu(rho j^mu)
       +d rho([e_barL phi N_L]_F-[bar L phi e_R]_F).     (S2)
```

In the independent executable the spatial derivative jets are arbitrary.
The temporal Weyl jets are solved explicitly from these two equations,
with **independent nonzero source and residual generators**, rather than
declaring a residual symbol equal to the desired answer. The literal
scalar-current expression then satisfies the complete polynomial identity
(S2), including its total derivative, for each separate squared gauge
coupling and a third mixed rational choice. No scalar EOM, vacuum condition,
mass inverse or division by `rho` is used.

The optional extra local map follows directly from the two residual rows:

```
ell_alpha=d rho sum_i F_(alpha i)phi N_Li,
x_i=-d rho sum_alpha F*_(alpha i)phi^dagger L_alpha,
L_old=L-kappa ell,       N_old=N-kappa(x+x^c).
```

Its row contraction is exactly the residual in (S2). The executable checks
the adjoint of `x`, its Majorana-real completion, and rejection of an
additional sterile half. No scalar or gauge field is shifted by this
particular map. As elsewhere, this is an explicitly specified coordinate
change, and must be composed on the complete retained action before it is
combined with earlier maps. Other EFT coefficients require their own
inverse-scale counting.

## Yukawa orientation, scalar symmetry and conjugation order

For the truncated neutrino Yukawa sources,

```
bar S_L,alpha=bar N_Rj Y*_(alpha j)phi^dagger,
S_Ri=sum_beta Y_(beta i)X_beta^c,
X_alpha=phi^dagger L_alpha,       bar X_alpha=bar L_alpha phi,
```

the nonderivative sourced part is

```
+d rho^2 (Y^dagger F)_(ji)bar N_Rj N_Li
-d rho (F M)_(alpha i)bar X_alpha N_Ri
-d rho (F Y^T)_(alpha beta)bar X_alpha X_beta^c.         (S3)
```

All matrix products and signs agree with direct expansion of the source
slots. The test uses noncommuting complex `F,Y` and a nonzero real symmetric
retained `M`; no sterile field is integrated out.

To check charge and Grassmann conventions explicitly, the independent
Weyl representation gives

```
epsilon=[[0,1],[-1,0]],
C=i gamma^2 gamma^0=diag(epsilon,-epsilon),
N_R=-epsilon n*,       bar N_R=n^T epsilon,
X^c=-epsilon bar X^T.
```

The script computes `C` from the gamma matrices. Scalar spin contractions
therefore have an antisymmetric spin matrix and a symmetric flavor
polynomial:

```
bar N_Rj N_Li=bar N_Ri N_Lj,
bar X_alpha X_beta^c=bar X_beta X_alpha^c.
```

It verifies these exchanges on actual exterior generators, verifies the
sterile bilinear adjoint by reversing the entire odd-field order, and
checks that antisymmetric active and sterile scalar coefficient parts
vanish. Flavor-symmetric coefficients must be retained. The symmetric
part of `Y^dagger F` or `F Y^T` multiplies an ordered scalar sum; conversion
to independent off-diagonal flavor pairs doubles those off-diagonal
coefficients. There is no inferred scalar-operator or mixed-mass half.

If the charged-lepton source is retained, its row term is proportional to
`(Y_e^dagger F)_(ri)bar e_Rr(H^dagger phi)N_Li`. The single commuting
doublet obeys `H^dagger phi=H^dagger epsilon H*=0`; this removes that
specific contraction. It does not remove derivatives such as
`(Dphi)^dagger H`, nor charged-lepton spectator effects elsewhere.
The independent scalar executable concentrates on the explicit neutrino
source and arbitrary other sources; the charged-lepton orthogonality
statement is a separate algebraic observation, not an extra counted test.

## Reproduction and scope

Run

```
python -I -S -B -O independent_scalar_audit.py
```

The standard-library implementation uses exact rational real/imaginary
coefficients, commuting scalar/connection jets, and anticommuting Weyl
components. It imports no producer code or saved coefficients, and failures
use explicit exceptions. The saved
[`INDEPENDENT_SCALAR_AUDIT_RESULTS.json`](INDEPENDENT_SCALAR_AUDIT_RESULTS.json)
records **76 positive checks and 20 rejected mutations**, 96 total.
These are exact algebraic controls, not distinct physical calculations.

This completes an independent audit of the **selected scalar-current
contribution**. Other portions of the original packet, spectator transport,
or subsequent reductions may contribute to the same scalar coordinates.
The formulas do not give their total physical coefficient, full EFT basis
closure, finite matching, an evanescent prescription, canonical normalization
of an interacting theory, running or an observable.
