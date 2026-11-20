# Algebraic basis of the active and retained-mass strength sectors

15 October 2026. Input candidate
`34acf849172231282a0a6e7a07a81157b2d48799`. This note constructs a
complete algebraically independent basis **inside the declared field-content
sectors**, and maps the previously reduced selected source packet into it.
The input coefficient `-1/12` is not recalculated as a loop coefficient here.
No equation of motion is used to identify operators in these sectors.

## Conventions and precise sector

Use the spin, scalar, Grassmann, flavor and gauge conventions of
[`ELECTROWEAK_REDUCTION.md`](../interacting_reduction/ELECTROWEAK_REDUCTION.md):

```text
phi=tilde H=epsilon H*, epsilon12=+1, epsilon^2=-1,
eta=(+---), epsilon^{0123}=+1,
sigma^{mu nu}=i[gamma^mu,gamma^nu]/2,
L^c=C bar L^T, C=i gamma^2 gamma^0,
D=partial+i Acal, Acal=g' y_L B 1_2+g T^I W^I,
y_L=y_phi=-1/2, T^I=tau^I/2,
Fcal_mu_nu=g'y_L B_mu_nu 1_2+g T^I W_mu_nu^I,
W_mu_nu^I=partial_mu W_nu^I-partial_nu W_mu^I
           -g epsilon^{JKI}W_mu^J W_nu^K.
```

There are `n` active flavors, weak indices `a,b=1,2`, and one displayed
lepton-number-violating orientation. The two active fermion slots are
`bar L_alpha^a` and `L_beta^{c,b}`; both are built from the same conjugate
active Grassmann fields. The non-derivative active sector contains exactly
these two fermions, two identical commuting `phi` fields and **one covariant
electroweak field strength**. All counts are complex coefficient coordinates
per orientation; the Hermitian conjugate supplies the conjugate sector.
No extra overall `1/2` is included in the following operator definitions.

The lower-field-count retained-mass sector contains one `bar L`, one
retained `N_R`, one `phi` and one strength. This is a dimension-six field
monomial reached with coefficient `F M=O(Lambda^-3 M)` from the same
dimension-seven selected packet. A retained mass is counted as order
`Lambda^0`, rather than integrated out.

## Lorentz and Grassmann reduction

Without derivatives or additional Lorentz tensors, the two same-chirality
fermion slots couple to a strength only through the chiral antisymmetric
spin tensor. Scalar bilinears cannot absorb its indices, and the dual
tensor does not give an additional complex coordinate. Specifically, with
`tilde V_mu_nu=eta_mu_mu eta_nu_nu epsilon^{mu nu rho sigma}V_rho_sigma/2`,

```text
bar L sigma^{mu nu}L^c tilde V_mu_nu
   =-i bar L sigma^{mu nu}L^c V_mu_nu.                     (1)
```

The sign is for the displayed right-chiral `L^c` and these gamma/epsilon
conventions; the conjugate orientation has the conjugate relation.
The exact checker evaluates all six independent strength components.
Their chiral spin-matrix image has rank three, and adjoining all six dual
images leaves that rank unchanged. This is four-dimensional algebra;
it fixes no dimensional-regularization evanescent or finite prescription.

For the primary spin matrix, `P_R sigma^{mu nu}C P_R^T` is symmetric under
ordinary spin-index transposition. Exchanging the two Grassmann slots
therefore gives

```text
bar L_alpha^a sigma^{mu nu}L_beta^{c,b}
   =-bar L_beta^b sigma^{mu nu}L_alpha^{c,a}.              (2)
```

It is the combined `(flavor,weak)` slot that is exchanged, not flavor
alone. All later symmetry statements follow from (2).

## Complete weak invariant space

Define the quadratic scalar tensors

```text
P_B^{ab}=phi^a phi^b,
P^{I,ab}=(T^I phi)^a phi^b,
P_+^{I,ab}=[(T^I phi)^a phi^b+phi^a(T^I phi)^b]/2,
P_-^{I,ab}=[(T^I phi)^a phi^b-phi^a(T^I phi)^b]/2.       (3)
```

`P_B` and `P_+` are symmetric in the two weak slots; `P_-` is
antisymmetric. An equivalent singlet expression is
`P_-^{I,ab}=-epsilon^{ab}(phi^T epsilon T^I phi)/2`.
The scalar fields commute, so their own antisymmetric contraction
`phi^T epsilon phi` is zero; this does not set `P_-` to zero.

Representation theory gives completeness for arbitrary `n`. Two identical
scalars furnish `Sym^2(2)=3`. The two open fundamental indices furnish
`2 tensor 2=1+3`. For `B`, the overlap contains one copy of `3`, giving
one invariant. For `W`, `3 tensor 3=1+3+5` overlaps the open `1+3`
in two irreducible components, giving exactly two invariants. The
`P_+` and `P_-` tensors realize those components. The hypercharges cancel;
the color-singlet fields cannot absorb a single gluon adjoint.

The checker independently solves the **entire** weak covariance problem,
without substituting numerical scalar backgrounds. For `B`, it starts
with all 12 coefficients in a general two-index matrix quadratic in
`phi`. For `W`, it starts with all 36 coefficients in three such matrices.
For `U=exp(-i theta T^A)`, the coefficient equations are

```text
partial_phi P_B[T^A phi]-T^A P_B-P_B(T^A)^T=0,
partial_phi P^J[T^A phi]-T^A P^J-P^J(T^A)^T
       +i epsilon^{AJI}P^I=0.                             (4)
```

Their exact ranks are respectively **11 and 34**, leaving nullities
**one and two**. The displayed tensors solve every polynomial coefficient
equation and span those full nullspaces. Thus their completeness does
not rest on sampled scalar values or a selected neutral background.

## Operators and independent flavor coordinates

Keep the strengths outside the spin bilinear with the contracted Lorentz
indices understood:

```text
O_B^{alpha beta}
 =sum_ab P_B^{ab} bar L_alpha^a sigma^{mu nu}L_beta^{c,b} B_mu_nu,
O_W^{alpha beta}
 =sum_Iab P^{I,ab} bar L_alpha^a sigma^{mu nu}L_beta^{c,b} W_mu_nu^I,
O_W,A^{alpha beta}
 =sum_Iab P_+^{I,ab} bar L_alpha^a sigma^{mu nu}L_beta^{c,b} W_mu_nu^I,
O_W,S^{alpha beta}
 =sum_Iab P_-^{I,ab} bar L_alpha^a sigma^{mu nu}L_beta^{c,b} W_mu_nu^I.
```

The labels `A,S` refer to **flavor** symmetry; the `+,-` labels on the
scalar tensors refer to **weak-index** symmetry. Equation (2) gives

```text
O_B^{beta alpha}=-O_B^{alpha beta},
O_W,A^{alpha beta}=(O_W^{alpha beta}-O_W^{beta alpha})/2,
O_W,S^{alpha beta}=(O_W^{alpha beta}+O_W^{beta alpha})/2,
O_W=O_W,A+O_W,S.                                          (5)
```

A complete independent basis and its complex-coordinate counts are

| Strength | Independent operators | Flavor coordinates |
| --- | --- | ---: |
| `B` | `O_B^{alpha beta}`, `alpha<beta` | `n(n-1)/2` |
| `W`, weak symmetric | `O_W,A^{alpha beta}`, `alpha<beta` | `n(n-1)/2` |
| `W`, weak antisymmetric | `O_W,S^{alpha beta}`, `alpha<=beta` | `n(n+1)/2` |
| All `W` | equivalently unrestricted `O_W^{alpha beta}` | `n^2` |

All ordered `O_W^{alpha beta}` are independent. In contrast, the full
kernel of the ordered `B` coefficient map is the symmetric flavor
subspace. The basis therefore has total active dimension
`n(n-1)/2+n^2`; for `n=3` it has **3 B plus 9 W coordinates**.

These independence statements are also verified by expanding every
operator in the complete exterior algebra of the chiral fermion
components and every quadratic scalar/strength monomial for `n=1,2,3`.
The ordered ranks are `B: 0,1,3` and `W: 1,4,9`; the `B` kernel
dimensions are `1,3,6`. No additional relation is found or expected from
the general representation-theory proof above. In particular the
single-flavor weak-symmetric tensor class is absent, but the
weak-antisymmetric class survives through charged weak components.

## Exact coefficient transport of the selected source

The input from the source-retaining reduction is

```text
P_active=-(1/12)K_(alpha beta)
  bar L_alpha^a(Fcal_mu_nu)_ab phi^b sigma^{mu nu}X_beta^c,
X_beta^c=sum_d phi^d(L_beta^d)^c,
K=F Y_nu^T,  F=C_seed Y_nu*,
K_A=(K-K^T)/2, K_S=(K+K^T)/2.                             (6)
```

In an unrestricted ordered-index presentation this is exactly

```text
P_active=-(g'y_L/12)sum_(alpha beta)(K_A)_(alpha beta) O_B^{alpha beta}
         -(g/12)sum_(alpha beta)[
            (K_A)_(alpha beta) O_W,A^{alpha beta}
           +(K_S)_(alpha beta) O_W,S^{alpha beta}].        (7)
```

Equivalently the unrestricted raw `W` coefficient is `-g K/12`.
For the independent flavor-coordinate basis above, write
`P_active=sum_(alpha<beta)C_B O_B+C_W,A O_W,A`
`+sum_(alpha<=beta)C_W,S O_W,S`. The coefficients are

| Independent coordinate | Coefficient in `P_active` |
| --- | --- |
| `C_B^{alpha beta}`, `alpha<beta` | `-g'y_L(K_A)_(alpha beta)/6` |
| `C_W,A^{alpha beta}`, `alpha<beta` | `-g(K_A)_(alpha beta)/6` |
| `C_W,S^{alpha beta}`, `alpha<beta` | `-g(K_S)_(alpha beta)/6` |
| `C_W,S^{alpha alpha}` | `-g K_(alpha alpha)/12` |

The off-diagonal pair factor two and its absence on the symmetric diagonal
are essential. The checker proves (7) and this table on **every matrix
unit of the coefficient space**, not on a random choice of `K`.
The actual counterterm is `kappa(P_active+h.c.)`, with
`kappa=1/(16 pi^2 epsilon)`. There is no duplicate Majorana or active
flavor `1/2` to append to these explicit basis definitions.

An arbitrary complex `K` parametrizes these algebraic operator coordinates.
The selected symmetric seed instead gives
`K=C_seed H_Y^T`, `H_Y=Y_nu Y_nu^dagger`, `C_seed^T=C_seed`.
It generally occupies a restricted source image. The full sector count
does not imply that this particular seed produces every independent
coefficient freely. Similarly, the displayed correlation of the `B`
and weak-antisymmetric-flavor coefficient matrices is a property of this
selected source map, not a universal constraint on the full EFT.
The [seed-image derivation](SEED_FLAVOR_IMAGE.md) gives those additional
constraints and its rank for an arbitrary Yukawa rank.

## Retained-mass tensor sector

For `s` retained sterile flavors, define

```text
D_B^{alpha i}=bar L_alpha^a phi^a sigma^{mu nu}N_Ri B_mu_nu,
D_W^{alpha i}=bar L_alpha^a(T^I)_ab phi^b sigma^{mu nu}N_Ri W_mu_nu^I.
```

There is one invariant for each strength and no flavor-exchange constraint
between the distinct active and sterile slots. Completeness follows from
`2 tensor 3=2+4` for the weak insertion, and is also checked by the full
linear weak-covariance systems: the `B` system has four unknown
coefficients and rank three; the `W` system has 12 unknown coefficients
and rank 11. Their nullities are each one. Complete Grassmann component
ranks for `n=1,2,3` and `s=1,2` are `n s` per strength.

The retained-mass part of the input representative is therefore

```text
P_mass=-(g'y_L/12)(F M)_(alpha i)D_B^{alpha i}
       -(g/12)(F M)_(alpha i)D_W^{alpha i}.                (8)
```

Mass and flavor order is `F M`; `M F` is not an allowed reordering.
The Majorana relation supplies the conjugate sector and does not impose
an antisymmetry between unlike active/sterile flavor slots. This is a
retained-mass dimension-six monomial at the fixed inverse-scale order,
not a separately matched observable dipole coefficient.
The checker separately verifies the same `-i` dual relation for the
right-chiral `N_R` column in this mixed sector. Its count is **`n s` per
strength for arbitrary `n,s`**, with no mixed-operator `1/2`.

## Relation to named classes and scope

The counts agree with the known SMEFT `O_lHB` and unrestricted `O_lHW`
classes, but our definitions and internal map (6)-(8) are primary. To
compare with [Zhang, arXiv:2306.03008v2, Table 1 and Eq. (5)](https://arxiv.org/html/2306.03008v2#S2.T1),
one must first conjugate its `LL HH` orientation into the displayed
`bar L L^c phi phi` orientation, use `phi=epsilon H*` with the declared
epsilon sign, translate its `D=partial-i g A` convention to our all-Dplus
gauge fields and strengths, and convert its Pauli insertion to
`tau^I=2T^I`. Its tabulated strength operators have no explicit gauge
coupling; ours also have none, while the **source coefficients** in
(7)-(8) do. These facts justify identifying the field-content and flavor
classes; no coefficient or overall sign is imported from that table.
The separate literature record preserves the primary-source convention
evidence rather than treating agreement in counting as new mathematics.
See the [convention and literature map](LITERATURE_BASIS_MAP.md), which
also compares the older weak flavor-symmetric/antisymmetric definitions.

Run `python3 -I -S -B -O verify_active_tensor_basis.py`. The completed
standard-library exact run passes **166 positive checks** and rejects
**eight deliberate mutations**. The adjacent JSON includes the complete
weak nullspace bases, exhaustive flavor ranks and coefficient-map controls.

The result is algebraic sector independence, not closure of a complete
physical operator basis under EOM or RG evolution. Derivative classes,
other source descendants and field redefinitions can change these
coordinates when the full retained-sterile basis is reduced. No additional
loop is integrated, no scheme is fixed, and no physical Wilson prediction,
beta function, finite amplitude, novelty or observable is established by
the `-1/12` transport or these rank computations.
