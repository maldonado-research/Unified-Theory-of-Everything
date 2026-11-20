# Independent audit of the local hypercharge lift

15 October 2026. This is separately authored AI-assisted algebra verification,
not external human peer review. The starting loop result is the reviewed
`R_H-R_L=G_J+D` at repository commit `4ac1f18`. No additional loop integral,
anomalous dimension, complete EFT basis, or observable is established here.
This note and its checker were authored independently in a local audit area.

## Direct local lift and conventions

Use (+---), `D=partial+i Gamma`, `D barL=partial barL-i barL Gamma`,
`sigma^{ab}=i[gamma^a,gamma^b]/2`, incoming `N` wave `exp(-ip.x)` and
outgoing `barL,phi=tildeH,B` waves `exp(+ir.x),exp(+ih.x),exp(+ik.x)`.
Thus `p=r+h+k`. Suppress the same right `P_L` projector as in the loop.

Define the following operators with no gauge coupling in their definitions:

    Qsigma = i(D_mu barL) sigma^{ab} gamma^mu phi N_L B_ab,
    Qdiv   = barL gamma^mu phi N_L partial^nu B_{nu mu},
    QT     = Qsigma/4 + Qdiv/3.

The row derivative in Qsigma acts only on barL. Direct Fourier substitution,
including both factors of i from the row derivative and the outgoing field
strength, gives the action kernels

    Qsigma -> [slash k,slash epsilon] slash r,
    Qdiv   -> (k.epsilon)slash k-k^2 slash epsilon.

Consequently QT has exactly the required D kernel. These are local,
dimension-seven, electroweak-invariant monomials. For the operator coefficient
one must still restore the coupling normalization: the four-point correction
to OJ is `g' y_L QT` in the action convention explained below.

The two displayed kernels are linearly independent as off-shell polynomials.
This statement does not assert that they exhaust the transverse kernel or
remain independent after EOM and integration by parts.

## Ordered row EOM and integration by parts

Set `K barL=-i(D_mu barL)gamma^mu=barS_L`. The identity

    [sigma^{ab},gamma^mu]=2i(g^{bmu}gamma^a-g^{amu}gamma^b)

gives, without using a free-field approximation,

    Qsigma = -barS_L sigma^{ab} phi N_L B_ab - 4 V_L,
    V_L    = (D^nu barL)gamma^mu phi N_L B_{mu nu}.

The derivative and gamma ordering are essential. In momentum space, moving
slash r from the right of the commutator to the left produces the vector
field-strength kernel. It cannot simply be annihilated by `ubar(r)slash r=0`.

Introduce

    QH = barL gamma^mu(D^nu phi)N_L B_{mu nu},
    QN = barL gamma^mu phi(partial^nu N_L)B_{mu nu}.

The total derivative of `barL gamma^mu phi N_L B_{mu nu}` gives
`V_L+QH+QN-Qdiv ~= 0`; the minus sign is from `B_mu nu=-B_nu mu`.
Therefore

    QT ~= QH+QN-(2/3)Qdiv-(1/4)barS_L sigma^{ab}phi N_L B_ab.   (A)

Here `~=` means equality modulo integration by parts and the declared sourced
reference-action EOM. It is not an identity of unrestricted off-shell vertices
after sources have been set to zero.

For the renormalizable sources in this specific hypercharge singlet channel,
the lepton-source term has two special cancellations:

* `bar e_R Y_e^dagger H^dagger sigma B phi N_L` vanishes because
  `H^dagger tildeH=0` for the same commuting scalar doublet.
* With flavor coefficient `F=C7 Ynu*` and `C7^T=C7`, the sterile source has
  coefficient `A=Ynu^dagger C7 Ynu*`, so `A^T=A`. The Majorana tensor
  `barN_i sigma^{ab}P_L N_j` is antisymmetric in i,j. Its contraction with A
  vanishes. An arbitrary mixed coefficient F need not have this property.

Thus these particular renormalizable source contributions can vanish after
their field and flavor identities have been checked. General EFT sources
cannot be discarded by this argument. A weak-generator insertion also
invalidates the simple `H^dagger tildeH=0` contraction.

## A second representation and the sterile source

Take `epsilon^{0123}=+1`, `gamma5=i gamma0 gamma1 gamma2 gamma3`,
`dualB^{mu nu}=epsilon^{mu nu rho sigma}B_{rho sigma}/2`, and define

    QHdual = barL gamma^mu(D^nu phi)N_L dualB_{mu nu},
    QAN    = barL phi sigma^{ab} B_ab (A N_L),    A=i slash partial.

The exact four-dimensional Clifford/Fourier identity is

    QT ~= (QH+i QHdual)/2-Qdiv/6-QAN/4,                        (B)

on this fixed chiral channel, modulo a total derivative, with covariant
completions understood. This can
be seen directly by using r=p-h-k in D and
`[slash k,slash epsilon]slash k=2[(k.epsilon)slash k-k^2 slash epsilon]`.
The rightmost pslash is the QAN kernel. The h term is the Higgs-derivative
and dual-Higgs-derivative pair. The independent checker fixes the dual sign;
using the opposite epsilon convention requires changing it.

The sterile EOM is `A N_L=M N_R+S_R^N`. Hence (B) includes

    -(1/4)barL phi sigma^{ab} B_ab M N_R
    -(1/4)barL phi sigma^{ab} B_ab S_R^N,

with the inherited flavor order. Replacing A N_L by only M N_R is not the
sourced reduction. The source is defined by that EOM and must be expanded
with the correctly charge-conjugated Yukawa channel if desired.

For the stated Majorana Yukawa action,
`S_R^N=Ynu^T phi^T L^c`. Let `X=phi^dagger L`. Then the sterile-source
term contains `(F Ynu^T)_{alpha gamma} barX_alpha sigma X_gamma^c`.
The spinor matrix sigma*C on this chirality is symmetric, so Grassmann
exchange makes this bilinear antisymmetric in active flavor. With
`Hflavor=Ynu Ynu^dagger`, its coefficient is therefore

    K=(C7 Hflavor^T-Hflavor C7)/2,    K^T=-K.

This active tensor generally survives, despite the sterile tensor cancellation
in the lepton-source route. For `Ynu=diag(1,2)` and off-diagonal unit C7,
`K_12=3/2`. The checker independently verifies the composite-spinor
antisymmetry, a complex unitary congruence, and the necessity of the displayed
transpose and complex conjugation.

If `U^dagger Hflavor U=diag(lambda_i)` and `Chat=U^dagger C7 U*`, then
`Khat_ij=(lambda_j-lambda_i)Chat_ij/2`. Thus its kernel consists precisely
of symmetric blocks within equal-eigenvalue subspaces, including null blocks.
For multiplicities d_a, the complex kernel dimension is
`sum_a d_a(d_a+1)/2`, and the rank is `(n^2-sum_a d_a^2)/2`. The Frobenius
norm is `sum_ij (lambda_j-lambda_i)^2 |Chat_ij|^2/4`. Exact controls cover
distinct, degenerate, entirely equal, and null spectra. These are conditional
linear-algebra facts about a source descendant, not a physical transition
moment or a new mathematical theorem.

The different coefficients of Qdiv in (A) and (B) are consistent because
their other operators differ. A contact coefficient cannot be moved between
these representations while retaining the other terms from the old one.

Liao and Ma, arXiv:1612.04527v1, Eq. (17), already give the independent
H-derivative B and dual-B operator family. Their Eq. (34) explicitly retains
a sourced Majorana tensor contact when reducing a derivative-on-L operator.
Our orientation is conjugate to their displayed one; equation-level matrix
translation, not an operator-name comparison, fixes phases and dual signs.

## Gauge EOM and coupling order

For `-B_ab B^{ab}/4`, the scalar kinetic term `|D H|^2`, and fermion kinetic
terms `bar f i slash D f`, variation gives

    partial_nu B^{nu mu} = +g' jY^mu,
    jY^mu = (i/2)H^dagger <->D^mu H
        + sum_flavors [ (1/6)barQ gamma^mu Q
                       +(2/3)baru_R gamma^mu u_R
                       -(1/3)bard_R gamma^mu d_R
                       -(1/2)barL gamma^mu L
                       -bare_R gamma^mu e_R ].

Color and weak indices are summed within the currents. Neutral sterile
fields do not add a hypercharge current. The scalar current contains its
own gauge completion: for general yH,
`i yH H^dagger <->D_mu H=i yH H^dagger <->partial_mu H
-2g' yH^2 B_mu H^dagger H` in the abelian part. The exact checker varies
the scalar kinetic term and confirms `delta Lmatter/delta B_mu=-g'jY^mu`.
Neither quark currents nor the Higgs current may be dropped.

The gauge EOM is the gauge-invariant classical-action equation; unrelated
gauge-fixing terms are not part of this operator reduction. Additional
reference-action interactions can change the sourced equation at the
power-counting order being reduced.

Since the action contains `g'yL QT`, reducing Qdiv adds another factor g'.
For example, (B) yields `-(g'^2 yL/6)(barL gamma^mu phi N_L)jY_mu`.
Those five-field and higher contacts occur in different external sectors and
at a different explicit power of g' from the one-B residue. Their presence
is an action identity, not evidence that all their loop matching is done.

## Charge-conjugated internal line and the action sign

The physical open-chain fermion is `chi=L^c`, of charge `-yL`, and the
scalar is phi=tildeH, of charge `yL`. With the declared D-plus convention,
their vertices are `-i g' q_chi gamma^nu` and
`-i g' q_phi(s+s')^nu`. The independent exact checks establish

    C gamma^{mu T} C^{-1}=-gamma^mu,
    C S^T(-q) C^{-1}=S(q),

so the charge flip must be accompanied by the proper momentum reversal.
Using charge yL for that internal conjugate fermion is the source of the
apparent sign mismatch between the old reduced integral labels and OJ.

Calibrated to the already declared positive-C7 seed vertex and zero-gauge
pole, the bubble's product of seed, fermion, Yukawa and scalar factors is
`i*i*(-i)*i=-1`. Adding a gauge vertex and a propagator multiplies by
`(-i)*i=+1` times g' times its physical charge. Thus the weighted loop
tensor is `yL(RH-RL)`, and

    iM_loop^UV = -i g'yL F (GJ+D)/(16 pi^2 epsilon),
    Lct        = +F[OJ+g'yL QT]/(16 pi^2 epsilon)+h.c.

This phase check uses the specified seed/Yukawa channel and its audited
zero-gauge normalization. It does not fix an arbitrary action whose seed
normalization, charge-conjugation rules, or flavor order differ. In particular
the literal `+J^T W1 Box J/2` has `C7=-W1` in this Fourier convention.

## Scope of the five-dimensional tensor classification

The independent review of `verify_transverse_space.py` finds its Ward map
correct under its premises. With generic independent p,r,k, constant scalar
coefficients, locality, degree two in momenta, and a fixed odd gamma chain,
there are 15 grade-one and three grade-three momentum structures. Adding
the separate m2 gamma^nu spurion gives 19 columns. Contraction with k has
rank 14, and the transverse kernel has dimension five:

    T_p, T_r, T_k, Gamma(k,nu,p), Gamma(k,nu,r),
    T_a^nu=a^nu slash k-(a.k)gamma^nu.

The independent checker constructs these columns as actual Weyl gamma
matrices at eight exact integer kinematic points. Stacking the entries
gives rank 19 for the columns and rank 14 for their contractions. Thus the
abstract column-independence premise has a concrete exact witness, including
the right chiral projection. D has coordinates `(0,1/2,1/3,0,1/2)`.

This is a complete kernel only within that declared kinematic ansatz. It
excludes extra spurions, other external fields, inverse momentum invariants,
and other chiral/mass sectors. No EOM or IBP quotient is taken. The compact
two-term expression for D is not the full five-dimensional kernel.

## Gauge completions and limitations

The operator definitions determine multi-gauge-field contacts. For example,
`iD_mu barL=i partial_mu barL+barL Gamma_mu` supplies Qsigma contacts
containing B_mu B_ab and W_mu B_ab. They cannot be deleted while continuing
to call the operator covariant. Whenever squared covariant derivatives are
reordered in the OJ reduction, `[D_mu,D_nu]=i F_mu nu` supplies the
corresponding curvature contacts, including the weak contribution. Replacing
covariant derivatives by commuting partial derivatives loses these terms.

The one-B calculation does not independently match all such additional
external-field sectors or preclude other gauge-invariant operators invisible
to it. Four-dimensional Clifford/dual identities here concern the retained
simple-pole tensor, not a choice of finite evanescent subtraction.

## Errors and failed shortcuts explicitly rejected

1. Treating the rightmost slash r as if it were next to the external row
   discards a nonzero field-strength term. An exact mutation rejects this.
2. My initial informal suggestion that the renormalizable charged-lepton
   source survives in this B singlet channel was wrong. Direct contraction
   gives Hdagger*tildeH=0; the correction is included above and in the checker.
3. The Majorana source cancels for the symmetric seed flavor tensor, but the
   arbitrary-F negative control is nonzero. General source removal fails.
4. Swapping B_mu nu and B_nu mu without its minus sign changes the IBP/gauge
   contact coefficient. A wrong-divergence-sign mutation is rejected.
5. The auxiliary producer symbol QE_h has its own minus normalization. With
   QE_h=-i g'yL QHdual, (B) contains `-QE_h/(2g'yL)`. A relayed plus sign
   for this symbol was inconsistent; the raw QH/QHdual formula is unambiguous.
6. A longitudinal Ward check alone neither fixes the ordered fermion chain
   nor validates an action-level overall sign. Charge conjugation, Fourier
   factors, and the zero-gauge normalization must be specified separately.

## Reproduction

Run `python3 -I -B -O verify_operator_critic.py` from the checker directory.
It imports SymPy, uses no producer or repository code, and checks exact
symbolic polynomials plus rational matrix-rank witnesses. Explicit exceptions
remain active with Python optimization. Its result is internal algebra
verification within the scope above.
