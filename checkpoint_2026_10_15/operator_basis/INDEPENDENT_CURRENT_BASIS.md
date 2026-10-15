# Independent component audit of the nonderivative current channel

15 October 2026. Independent exact component and Grassmann-algebra derivation
from the selected gauge-current descendant in
[`CURRENT_DESCENDANTS.md`](../interacting_reduction/CURRENT_DESCENDANTS.md).
The adjacent checker imports no repository calculation and no producer
coefficient table. This is AI-assisted internal verification, not independent
human peer review or a new loop integration.

## Fixed sector and source normalization

Keep `phi=tilde H`, `F=C_seed Y_nu*`, `T^I=tau^I/2`, `y_L=-1/2`,
and the `D=partial+i A` reference action. Suppress only the original common
counterterm factor `kappa`; Hermitian conjugates are required. The input is

```
P_current=-(g'^2 y_L j_mu J_Y^mu+g^2 j_I,mu J_I^mu)/6,
j_mu=F_(alpha i)(bar L_alpha,a gamma_mu N_Li)phi_a,
j_I,mu=F_(alpha i)(bar L_alpha,a gamma_mu N_Li)(T^I phi)_a.
```

The spectator currents are `bar f_beta gamma_mu f_gamma`, initially with
`delta_(beta gamma)`, with a color sum for quarks. The five species are
`f=Q_L,L_L,u_R,d_R,e_R`. Their hypercharges are respectively
`1/6,-1/2,2/3,-1/3,-1`. Right-handed fields are used in that convention,
rather than silently replacing them by left charge conjugates.

The sector here has one scalar, four fermions, no derivatives, the displayed
chirality assignment, and one sterile left component. Scalar-current
descendants, curvature terms, other external-field sectors, and EOM relations
to derivative operators lie outside this component-basis calculation.

## Two weak contractions and their exact coefficients

For a fundamental spectator define the ordered monomials

```
D^f_(alpha beta gamma i)
  =(bar L_alpha,a gamma_mu N_Li)phi_a
       (bar f_beta,c gamma^mu f_gamma,c),
E^f_(alpha beta gamma i)
  =(bar L_alpha,a gamma_mu N_Li)phi_b
       (bar f_beta,b gamma^mu f_gamma,a).
```

The color index is summed inside each quark bilinear. There is no additional
factor three multiplying an operator already defined with that sum.
Explicit half-Pauli components give

```
sum_I (T^I)_ab (T^I)_cd
    =delta_ad delta_bc/2-delta_ab delta_cd/4.
```

Consequently the doublet contribution before the identical-lepton reduction is

```
sum_(alpha beta i) F_(alpha i)
   [(g^2/24+g'^2 y_f/12) D^f_(alpha beta beta i)
      -g^2 E^f_(alpha beta beta i)/12].                 (A1)
```

For a weak singlet only the direct contraction occurs, with coefficient
`g'^2 y_f F_(alpha i)delta_(beta gamma)/12`.

| Spectator | Direct coefficient multiplying `F_(alpha i)delta_(beta gamma)` | Exchange coefficient |
|---|---|---|
| `Q_L` | `(3g^2+g'^2)/72` | `-g^2/12` |
| `L_L`, before Fierz reduction | `(g^2-g'^2)/24` | `-g^2/12` |
| `u_R` | `g'^2/18` | absent |
| `d_R` | `-g'^2/36` | absent |
| `e_R` | `-g'^2/12` | absent |

For `Q_L` the exchanged contraction includes the off-diagonal charged weak
currents. Omitting it changes the component polynomial. This table describes
the specified representative of the sourced packet, with its common `kappa`
and conjugate channel; its entries are not physical beta functions.

## Identical barred leptons: Fierz swaps flavor labels with a plus sign

Use `(+---)` and two-component matrices
`sigma^mu=(1,tau^I)`, `bar sigma^mu=(1,-tau^I)`.
The checker computes all components of the contractions directly and finds

```
bar_sigma^mu_(a b) bar_sigma_mu,(c d)=2 epsilon_ac epsilon_bd,
bar_sigma^mu_(a b) sigma_mu,(c d)=2 delta_ad delta_bc,
epsilon_01=1, epsilon_10=-1.
```

These are the left-left and left-right spectator spinor channels. The
four-fermion word is expanded in its actual source order
`bar L_alpha, N_i, bar f_beta, f_gamma`, then sorted as an exterior product.
Repeated identical components vanish. This explicitly includes the
anticommuting-field signs, rather than applying a commuting spinor identity.

For `f=L` the result is

```
E^L_(alpha beta gamma i)=D^L_(beta alpha gamma i).       (A2)
```

Exchanging the two dotted barred-spinor indices changes the spinor tensor
sign. Exchanging the corresponding Grassmann fields supplies the second
minus sign. Their product is the plus sign in (A2). No flavor symmetry or
antisymmetry may therefore be imposed on `D^L` by deleting either ordering.

An unrestricted ordered-flavor direct basis has coefficient

```
C^L_(alpha beta gamma i)
   =(g^2-g'^2)F_(alpha i)delta_(beta gamma)/24
       -g^2 F_(beta i)delta_(alpha gamma)/12.           (A3)
```

The second term changes the flavor index of `F`. In the one-active-flavor
limit the coefficient is `-(g^2+g'^2)F/24`, and its operator is nonzero.
There is no extra identical-field half in this convention: the action is
defined by the displayed ordered polynomial with no such prefactor.

Equation (A2) is strictly four-dimensional. Its difference defines a possible
evanescent operator away from four dimensions; the present component proof
does not select its subtraction or its finite effects.

## Component witnesses prove the ordered families are independent

The ranks below are ranks of operator polynomials with unrestricted flavor
coefficients. They are not counts of independently generated seed parameters:
the selected coefficients in (A1)-(A3) occupy a restricted image controlled
by `F`, the couplings and the flavor deltas.

For the lepton family, select this component in each `D^L_(alpha beta gamma i)`:

```
phi_0,
bar L_alpha,(weak 0,dotted spin 0),
bar L_beta,(weak 1,dotted spin 1),
N_i,(spin 0),
L_gamma,(weak 1,spin 1).                              (A4)
```

After canonical Grassmann sorting it has coefficient `+2` or `-2`.
The weak labels identify which barred flavor is `alpha` and which is `beta`.
This component occurs in no other ordered member of the direct family.
It remains nonzero when `alpha=beta`, because the weak/spin components are
distinct. Different `gamma` and `i` are distinguished by the unbarred
generators. Thus these selected rows form a diagonal minor, proving
`dim(D^L)=n_L^3 n_N` for any positive flavor counts.

For the quark doublet choose color zero and the same spin components.
The direct witness has primary barred-lepton weak index zero and barred
spectator weak index one; the exchanged witness reverses those two weak
indices. Both have `phi_0` and unbarred-spectator weak index one. Each is
nonzero only in its matching direct/exchanged family member. Hence
`dim(D^Q,E^Q)=2 n_L n_N n_Q^2`.

For a right-handed singlet use `phi_0`, primary barred-lepton weak/spin
indices `(0,0)`, sterile spin zero, and spectator barred/unbarred spin zero,
again selecting one quark color if relevant. The explicit left-right tensor
has coefficient two. Flavor labels separate all columns, proving
`dim(D^f)=n_L n_N n_f^2` separately for `u_R,d_R,e_R`.

Within these fixed chiral field contents there are two dotted and two
undotted spinor indices, with one Lorentz singlet contraction in each pair.
The two weak invariants for distinct doublets, or the one invariant for a
singlet spectator, are therefore exhausted. The identical-lepton relation
(A2) replaces the second weak contraction by a distinct ordered-flavor
member, whose independence is certified by (A4). This is algebraic
completeness within the specified nonderivative component sector. It is
not closure under the full EFT equations of motion.

## Reproduction and actual coverage

Run

```
python -I -S -B -O independent_current_basis.py
```

The standard-library script uses exact rational real/imaginary arithmetic,
explicit Pauli and Weyl entries, and canonical exterior monomials. It keeps
each `F_(alpha i)` and each squared gauge coupling as a separate formal
symbol while comparing the literal source with (A1)-(A3). It stores explicit
witness rows and diagonal entries in
[`INDEPENDENT_CURRENT_BASIS_RESULTS.json`](INDEPENDENT_CURRENT_BASIS_RESULTS.json).

For three active/spectator flavors and two sterile flavors the component
minors have ranks `54,108,54,54,54` for `L,Q,u_R,d_R,e_R`: 324 members in
the union, separated also by external species. All diagonal entries are
`+2` or `-2`, and the exact nonzero determinants are recorded. The check
passes **23,544 positive component comparisons** and rejects **66 mutations**.
Most positive comparisons, **23,328**, are entries of these explicit square
rank matrices; this count must not be presented as a count of physical
calculations, observables, or loop diagrams.

The other checks cover all Pauli and Weyl contractions, 54 full Grassmann
Fierz exchanges, 108 weak-current reconstructions, all five complete symbolic
source maps, and the one-flavor limit. Wrong Fierz sign or flavor routing,
an extra identical-field half, an extra color factor, use of Pauli rather
than half-Pauli generators, a missing quark weak exchange, reversed
hypercharge sign, and omitted currents fail. These checks do not determine
other Wilson sectors, evanescent finite terms, finite matching, physical
running, or an observable prediction.
