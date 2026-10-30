# Algebraic basis for the selected nonderivative current descendants

15 October 2026. Input: the selected sourced representative and the current
definitions of
[`ELECTROWEAK_REDUCTION.md`](../interacting_reduction/ELECTROWEAK_REDUCTION.md)
and [`CURRENT_DESCENDANTS.md`](../interacting_reduction/CURRENT_DESCENDANTS.md).
This note supplies the **four-dimensional algebraic quotient** for its
four-fermion-plus-`phi` current sectors. It adds spinor Fierz and identical
field constraints beyond the previous weak-index contraction.

The result is complete for the specified nonderivative field content and
Lorentz/gauge contractions. It is not a complete dimension-seven EFT basis
modulo derivatives, integration by parts, equations of motion, other field
sectors, or higher-operator field changes. No new loop is integrated here.

## 1. Fix the fields, currents, and ordered convention

Use `phi=tilde H`, `y_L=y_phi=-1/2`, `T^I=tau^I/2`, and the original
ordered flavor coefficient `F=C_seed Y_nu*`. Define

    j_mu = sum_(alpha i a) F_(alpha i)
                 (bar L_alpha^a gamma_mu N_Li) phi^a,
    j_I,mu = sum_(alpha i a b) F_(alpha i)
                 (bar L_alpha^a gamma_mu N_Li) T^I_ab phi^b,
    P_current=−[g'^2 y_L j_mu J_Y^mu+g^2 j_I,mu J_I^mu]/6.       (B1)

The action is `kappa(P_current+P_current^dagger)`, with the same
`kappa=1/(16 pi^2 epsilon)` as the input. A retained spectator contributes

    J_Y^mu(f)=y_f sum_beta bar f_beta gamma^mu f_beta,
    J_I^mu(f_L)=sum_beta bar f_Lbeta gamma^mu T^I f_Lbeta.

The standard spectator charges are `y_Q=1/6`, `y_u=2/3`, `y_d=−1/3`,
and `y_e=−1`. A right-handed singlet has no weak current. The scalar
current and its derivatives remain in the input representative but are
outside this four-fermion sector.

For a doublet spectator define the ordered operators

    O0^f_(alpha beta gamma i)
       = (bar L_alpha^a gamma_mu N_Li) phi^a
                      (bar f_beta^c gamma^mu f_gamma^c),
    Oex^f_(alpha beta gamma i)
       = (bar L_alpha^a gamma_mu N_Li) phi^c
                      (bar f_beta^c gamma^mu f_gamma^a).       (B2)

All weak indices in (B2) are summed. The two currents are Grassmann even,
but each is expanded in its displayed fermion order when applying Fierz.
For a right-handed singlet there is just

    O0^f_(alpha beta gamma i)
       = (bar L_alpha^a gamma_mu N_Li) phi^a
                      (bar f_Rbeta gamma^mu f_Rgamma).       (B3)

Quark operators include one sum over the common color index of the two
quark fields. This is the unique color singlet because every other field
here is colorless. The color sum is included in the operator, so no extra
factor of three multiplies its coefficient. Flavor dimensions below count
generations, not color components.

No factor `1/2` is inserted into these ordered definitions. Coefficients in
a basis with additional normalization factors must be converted explicitly.

## 2. Independent weak and Lorentz counts

For a distinct doublet, the two invariant weak tensors are

    W0_abcd=delta_ab delta_cd,
    Wex_abcd=delta_ad delta_cb.

They span the two singlets in `(2bar tensor 2) tensor (2bar tensor 2)`.
An epsilon-pair contraction is their difference, not a third invariant.
Explicitly,

    sum_I T^I_ab T^I_cd=Wex_abcd/2−W0_abcd/4.                (B4)

The checker independently builds all infinitesimal weak-generator equations
on the full sixteen-component tensor space. Their rank is fourteen;
`W0,Wex` are independent null vectors. This determines the weak invariant
space before using (B4).

For the Lorentz count use two-component spinors, metric `(+---)`,
`sigma^mu=(1,tau)`, and `bar sigma^mu=(1,−tau)`. For left spectators the
component tensor is

    T_LL_(dot a b dot c d)
      =sum_mu eta_mu (bar sigma^mu)_(dot a b)
                         (bar sigma^mu)_(dot c d).

There are two dotted and two undotted Weyl slots. Each pair has one
antisymmetric invariant, hence their scalar tensor space is one-dimensional.
The right-current tensor is instead

    T_LR_(dot a b c dot d)
      =sum_mu eta_mu (bar sigma^mu)_(dot a b)
                         (sigma^mu)_(c dot d).

Its dotted fundamental/dual pair and undotted fundamental/dual pair each
have one delta invariant, again giving one Lorentz scalar. The checker
constructs separate traceless-generator equations for the dotted and
undotted groups on all sixteen components. Both systems have rank fifteen;
the displayed nonzero tensors span their kernels. This uses the complex
span of the two `SL(2)` Lie algebras, rather than imposing rotations alone.

Thus a distinct doublet has two weak structures with one Lorentz structure
each. A right singlet has one weak structure and one Lorentz structure.
These counts do not yet identify repeated lepton fields.

## 3. Identical barred leptons identify the two weak structures

For anticommuting four-dimensional left fields the chiral vector identity is

    (bar A gamma_mu P_L B)(bar C gamma^mu P_L D)
      =(bar A gamma_mu P_L D)(bar C gamma^mu P_L B).       (B5)

The sign is positive. In components the antisymmetric undotted spin
contraction changes sign when its slots are exchanged; the odd Grassmann
permutation supplies a second minus sign. Treating the spinors as commuting
variables would give the wrong identity for this action.

Apply (B5) to `Oex^L`. The two resulting bilinears commute because they
are even. Relabeling their weak sums then gives

    Oex^L_(alpha beta gamma i)=O0^L_(beta alpha gamma i).   (B6)

Equation (B6) removes the exchange structure as an additional basis
direction, while exchanging the two barred flavor labels. It does **not**
impose `O0_(alpha beta gamma i)=O0_(beta alpha gamma i)`. Distinct ordered
barred flavors still carry distinct weak tensors. The unrestricted ordered
basis is therefore

    {O0^L_(alpha beta gamma i)},
    dim = n_L^3 n_N.                                    (B7)

A constructive nonredundancy witness works at every flavor multiplicity.
For each ordered tuple, select `phi^1` and the component containing

    bar L_alpha^(weak 1,dot spin 1),
    bar L_beta^(weak 2,dot spin 2),
    N_i^(spin 1),  L_gamma^(weak 2,spin 2).

Its nonzero coefficient is `+2` or `−2` depending on canonical exterior
ordering. The scalar weak index fixes the first barred slot, and the
unbarred lepton's weak index fixes the second. No other ordered basis
tuple contributes that same Grassmann monomial. These components form a
diagonal nonsingular minor. The construction also works when `alpha=beta`,
because the two weak/spin components are different generators.

The exact exterior calculation additionally expands every weak and spinor
component and finds ranks `1`, `16`, and `27` at `(n_L,n_N)=(1,1),(2,2),
(3,1)`. Appending all `Oex` columns leaves these ranks unchanged. Thus
one active flavor does not make this sector vanish.

For distinct `Q`, separate component witnesses distinguish direct and
exchange contractions, so

    dim(Q sector)=2 n_L n_Q^2 n_N.                       (B8)

For example, with `phi^1`, the weak choices
`(bar L,bar Q,Q)=(1,2,2)` and `(2,1,2)` select the two structures separately.
Each right-singlet sector has dimension

    dim(f_R sector)=n_L n_f^2 n_N.                       (B9)

These statements are for complex coefficient coordinates before adding
the Hermitian conjugate. They use no mass relation or equation of motion.

## 4. Representative coefficients, including identical-flavor transport

For a distinct doublet, direct substitution of the currents in (B1) gives

    P_f=sum_(alpha beta gamma i) F_(alpha i)delta_beta,gamma
             [A_f O0^f_(alpha beta gamma i)
                    +B_f Oex^f_(alpha beta gamma i)],
    A_f=g^2/24+g'^2 y_f/12,       B_f=−g^2/12.           (B10)

The hypercharge sign follows from `y_L=−1/2` in (B1). For `Q`,

    A_Q=(3g^2+g'^2)/72,       B_Q=−g^2/12.              (B11)

For the identical lepton current one must use (B6), including its flavor
permutation. The independent ordered coefficient is

    P_L=sum_(alpha beta gamma i) C^L_(alpha beta gamma i)
                                     O0^L_(alpha beta gamma i),
    C^L_(alpha beta gamma i)
      =(g^2−g'^2)F_(alpha i)delta_beta,gamma/24
          −g^2 F_(beta i)delta_alpha,gamma/12.            (B12)

The second term is not an arbitrary coefficient symmetrization. It is
transported from the exchanged weak structure. With one active flavor,

    C^L=−(g^2+g'^2)F/24.                                (B13)

This nonzero check is a useful discriminator against both a wrong Fierz
sign and an erroneous flavor-antisymmetry rule.

For singlets the coefficients multiplying `O0^f` are as follows; each
entry includes `F_(alpha i)delta_beta,gamma`.

| Spectator | Vector-current coefficient | Scalar-exchange coefficient |
|---|---:|---:|
| `u_R` | `+g'^2/18` | `−g'^2/9` |
| `d_R` | `−g'^2/36` | `+g'^2/18` |
| `e_R` | `−g'^2/12` | `+g'^2/6` |

The scalar convention in the last column is explicitly

    S^f_(alpha beta gamma i)
       =(bar L_alpha^a f_Rgamma)phi^a(bar f_Rbeta N_Li),
    O0^f=−2 S^f.                                       (B14)

The negative sign in (B14) is the Grassmann sign. The component Lorentz
identity has factor `+2 delta_dot a,dot d delta_c,b`; reordering its
fermions into the displayed scalar products introduces a minus sign.
The checker verifies (B14) as an exterior polynomial, including all
sixteen spinor-component choices, for every tested flavor tuple.
Literature scalar representatives require their own charge-conjugation,
Higgs-epsilon, flavor-order, and normalization conversion before importing
these numerical coefficients.

## 5. What has and has not been reduced

Equations (B7)-(B14) provide an independent algebraic basis and coefficient
map for every retained `L,Q,u_R,d_R,e_R` current contribution in (B1).
The representative populates only a restricted subspace of these bases:
the source currents are diagonal in spectator flavor, and the map is
linear in `F`. Basis dimension is not the number of independent UV
parameters generated by this selected insertion. Additional sectors or
source terms can contribute to the same final operators and would need
to be combined before physical matching.

The Lorentz/Fierz quotient here is strictly four-dimensional. It is
appropriate as the declared representative conversion for the simple UV
pole, with dimension-dependent differences retained as a scheme issue.
It does not specify finite evanescent counterterms or authorize using
four-dimensional Fierz identities inside an unqualified `d`-dimensional
calculation. Derivative/EOM relations, complete retained-sterile basis
matching, parameter and field counterterms, and lower-order operator
products remain separate requirements for physical running.

This supplies neither a beta function nor a finite amplitude, observable,
new loop result, or claim of novel fundamental mathematics.

## 6. Exact reproduction and discriminating controls

Run `python -I -S -B -O verify_current_operator_basis.py` here. It imports
only the standard library and uses exact Gaussian rational arithmetic.
The weak and Lorentz invariant-space ranks are computed separately from
the representative map. Fermion monomials are sorted with their actual
exterior permutation signs, and repeated component generators vanish.
Identical barred leptons use the same generators from the outset.

All basis operators are expanded in weak and spinor components. The direct
source side uses the original Pauli-generator current sum before applying
weak completeness or Fierz. It is compared against (B11)-(B14) separately
in the `g^2` and `g'^2` directions and at mixed rational couplings, with
generic complex `F`. Treating `F` as arbitrary verifies the linear operator
identity; interpretation of the selected insertion keeps
`F=C_seed Y_nu*` in its original order.

The optimized isolated run passes **427 positive controls** and rejects
**112 mutations**, for **539 checks total**. Negative controls reverse
Fierz signs or hypercharge, substitute `tau` for `T=tau/2`, change the
transported exchange coefficient, impose a spurious barred-flavor symmetry, or
discard the one-flavor lepton sector. The adjacent JSON records actual
component ranks and counts. These are algebraic verification controls;
they are not additional physical observations or loop calculations.
