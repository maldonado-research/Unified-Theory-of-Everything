# Local normalization of the selected packet with weak interactions retained

15 October 2026. Continuation of the neutral calculation in
[`CANONICAL_NORMALIZATION.md`](../field_redefinition/CANONICAL_NORMALIZATION.md)
and the fully covariant action in
[`COVARIANT_ACTION_AND_WARD.md`](../nonabelian_gauge/COVARIANT_ACTION_AND_WARD.md).
The input is the actual six-placement `O_J` plus the displayed curvature packet,
not an arbitrary kinetic matrix. This supplies a concrete local change of
coordinates and its remaining interactions using standard field-redefinition
algebra. It does not supply a complete
canonical EFT basis or physical running.

## 1. A local covariant lift of the actual neutral normalization

Keep `D_mu=partial_mu+i A_mu`, `A_mu=g'y_L B_mu 1+g T^I W_mu^I`,
`T^I=tau^I/2`, and `F=C_seed Y_nu*`, in that flavor order. Both `L` and
`phi=tilde H` transform in the same weak representation; `N_L` is neutral.
Let `Phi=phi tensor 1_active`, with weak and flavor indices kept distinct,
and define

    P0=-partial^2,       P=P0+m_H^2,
    H=Phi F P,          H^ddagger=P F^dagger Phi^dagger,
    q=(L,N_L)^T,
    Z=[0 H; H^ddagger 0],
    K=diag(i slash D, i slash partial).

The symbol `ddagger` is the integrated Dirac adjoint of a differential
operator. Thus `H^ddagger L=P(F^dagger Phi^dagger L)`: derivatives act on
the **whole neutral composite**, including the scalar. It is not
`F^dagger Phi^dagger P L`. The connection does not appear explicitly in
`Z`, since `Phi^dagger L` and `N_L` are neutral. Nevertheless `H N_L`
transforms as `L`, so this is a covariant map between the two representation
spaces. It is not a constant mixing of inequivalent representations.

`Z^ddagger=Z`. Since `[F]=-3`, `[Phi]=1`, and `[P]=2`, `Z` has engineering
dimension zero. With real first-order bookkeeping parameter `kappa`, the
local perturbative substitution is

    q_old=(1-kappa Z/2) q_c,
    L_old=L_c-kappa Phi F P N_Lc/2,
    N_Lold=N_Lc-kappa P(F^dagger Phi^dagger L_c)/2.        (I1)

All right sides use the new fields. `N_Rold=(N_Lold)^c` completes the
Majorana field; no additional factor `1/2` belongs to that completion.
The inverse through this order changes the shift sign. No mass inverse
is used. For `phi=(u,0)` constant and zero gauge background, `H` reduces
exactly to the previous `u(p^2+m_H^2)F` neutral mixing, and the charged
lepton is unchanged on that background.

This fermion-only map is applied directly to the original covariant packet.
It includes no scalar shift, so it introduces no scalar-tadpole contraction.
It is distinct from the simultaneous source-reduction substitutions in
`FIELD_REDEFINITION.md`; combining those maps requires composing all their
shifts and retaining that earlier off-shell background caveat.

This is a specified extension of the neutral normalization to interacting
fields. Other covariant operators that vanish on that background could
be added to `Z`; their effect would be a different split of the remaining
interactions. We fix the displayed `Z` throughout.

The [source-reducing field map](ELECTROWEAK_REDUCTION.md) and (I1) are
separate maps of the same original packet. To relate their coordinates,
write `Xi_old=Xi_S+kappa v_S` and `Xi_old=Xi_C+kappa v_C` for the full
field vector and its conjugates. At first order their composition is
`Xi_S=Xi_C+kappa(v_C-v_S)`. It must be applied to the full action and
sources with consistent order counting. The calculation below starts from
the original packet; it does not show that the source-reduced representative
is already canonical or justify adding the two induced coefficient tables.

## 2. Direct reduction of all six placements

Write `A_N=i slash partial`. Integrate the row derivatives in the literal
six terms of `O_J` onto their complete right composites. With `F` suppressed
only for readability, their upper mixed differential kernel is

    Q_J N = (5/6) phi A_N P0 N
          +(1/3)(-D^2)(phi A_N N)
          +(1/6)(-D^2 phi) A_N N
          -(1/6)i slash D(phi P0 N)
          +(m_H^2/2) phi A_N N
          +(m_H^2/2)i slash D(phi N).                    (I2)

Each `phi ... N` means `Phi F ... N`. Applying the covariant product rule
to (I2) gives

    Q_J N = phi A_N P N
          -(D^2 phi) A_N N/2
          -(2/3)(D_mu phi) partial^mu(A_N N)
          -(i/6) gamma^mu(D_mu phi) P0 N
          +(i m_H^2/2) gamma^mu(D_mu phi) N.             (I3)

The upper block of the normalization action is

    Q_norm=({K,Z}/2)_LN,
    Q_norm N=phi A_N P N+(i/2)gamma^mu(D_mu phi) P N.    (I4)

Consequently the exact local difference is

    R_J N = -(D^2 phi) A_N N/2
            -(2/3)(D_mu phi) partial^mu(A_N N)
            -(2i/3)gamma^mu(D_mu phi) P0 N,
    Q_J=Q_norm+R_J.                                    (I5)

All `m_H^2` terms cancel in `R_J`. This follows from the literal packet,
without using any scalar or fermion equation of motion or imposing a
stationary background. The derivatives in (I5) remain part of the action.
In particular, on a varying scalar background `R_J` still contains
derivative fermion bilinears; (I1) does not canonically normalize the
entire interacting differential operator.

For the curvature terms use `Fcal_mu nu=partial_mu A_nu-partial_nu A_mu
+i[A_mu,A_nu]` and its adjoint covariant derivative `Dcal`. Row integration
by parts gives the full additional upper kernel

    R_D N= -(i/4) sigma^{mu nu} gamma^rho
                    D_rho(Fcal_mu nu phi N)
           +(1/3)gamma_mu (Dcal_nu Fcal^{nu mu}) phi N.  (I6)

The first derivative in (I6) acts on the curvature, scalar, and sterile
field. Equivalently its bracket expands as

    (Dcal_rho Fcal_mu nu) phi N
    +Fcal_mu nu(D_rho phi)N+Fcal_mu nu phi partial_rho N.

The order `Fcal phi` and the order `sigma gamma` are retained. Deleting
the middle term changes the interacting packet.

For `Q=Q_J+R_D`, the selected kinetic counterterm is therefore

    kappa bar q {K,Z} q/2
       +kappa [bar L(R_J+R_D)N_L+h.c.].

Substitution (I1) in the reference kinetic action supplies
`-kappa bar q {K,Z} q/2`. Thus that anticommutator cancels, while the
displayed interaction remainder, induced mass terms, and transformed
sources must all be kept.

## 3. Retained Majorana mass: transpose is a differential transpose

Use the same holomorphic left-Weyl convention as the previous checkpoint:

    L_mass=-1/2[q^T epsilon Mcal q+h.c.],
    Mcal=[0 B; B^t M*],       B=Phi* Y_nu*.

Here `t` denotes formal transpose with integration by parts, without
complex conjugation. It acts on differential operators as well as weak
and flavor indices; `P^t=P`. In particular

    Z^t=[0 Phi* F* P; P F^T Phi^T 0],
    Mcal_c=Mcal-kappa(Z^t Mcal+Mcal Z)/2.               (I7)

For clarity, the three induced coefficient blocks are the operators

    delta A =-[Phi* F* P Y_nu^dagger Phi^dagger
               +Phi* Y_nu* P F^dagger Phi^dagger]/2,
    delta B =-Phi* F* P M*/2,
    delta D =-[P F^T (Phi^T Phi*)Y_nu*
               +Y_nu^dagger(Phi^dagger Phi)F P]/2.     (I8)

For example the first `P` in `delta A` acts on
`Y_nu^dagger Phi^dagger L`; it acts on the scalar. In `delta D` the first
`P` acts on the scalar norm and the field to its right. On a constant
neutral background these reduce to the previous block shifts. They
retain the full nonzero symmetric `M*`, and are formally symmetric as a
mass operator. A numerical commuting test spinor is used only to check
this coefficient transpose property; it is not a replacement for
Grassmann statistics. Replacing `Z^t` by its Hermitian adjoint fails.

These are coordinate-induced coefficients, not independently generated
physical masses. The previously stated linear-seed counting still
applies: varying the original dimension-seven seed with (I1) is a product
of seed coefficients and lies beyond the selected `Lambda^-3` truncation.
Equations (I7)-(I8) display the declared neutrino/Majorana sector. Additional
retained renormalizable interactions, such as a charged-lepton Yukawa,
must also undergo (I1); their induced terms are not silently set to zero.
Other retained operators require their own consistent inverse-scale
counting; the `C6 C5` issue is not removed by this convention.

## 4. Sources and weak currents transform with the fields

For independent row sources `bar s_L L+bar s_N N_L+h.c.`, (I1) produces

    delta L_source= -kappa/2[
       bar s_L Phi F P N_L+bar s_N P(F^dagger Phi^dagger L)]
       +h.c.                                             (I9)

Integration by parts moves the second `P` onto the entire source row.
Dropping those derivatives changes source response. More general
derivative source vertices require the same operator composition, in
the displayed order.

At fixed scalar fields `Z` contains no explicit gauge connection. For an
arbitrary variation `delta A_mu=V_mu`,

    delta_A Q_norm N= -gamma^mu V_mu Phi F P N/2.        (I10)

The shift of the original weak current from substituting (I1) into
`bar L i slash D L` gives the opposite cross contact. Both are required
for cancellation of the normalization part of the interacting action.
For example, at `phi=(u,0)` the original charged interaction contains
`-g W^-_mu bar e_L gamma^mu nu_L/sqrt(2)+h.c.`. Its induced term is
`+kappa g u W^-_mu bar e_L gamma^mu F P N_L/(2sqrt(2))+h.c.`;
the normalization packet contributes its negative. The remainder (I5)
and curvature terms still supply their own current contacts. An apparent
sterile charged-current mixing read from the field substitution alone
would omit this cancellation.

An off-shell check stronger than the constant-background projection is
local covariance with `delta phi=-i alpha phi` and
`delta A_mu=partial_mu alpha+i[A_mu,alpha]`. Each of `Q_J N`,
`Q_norm N`, `R_J N`, and `R_D N` transforms as `-i alpha` times itself;
`P(F^dagger Phi^dagger L)` is invariant. The checker below differentiates
the actual operator expressions under nonconstant, noncommuting weak
transformations. This is an action-level current/Ward consistency test,
not a separate loop calculation.

## 5. What happens to the already checked higher-gauge contacts

Because `Z` has no explicit `A`, `Q_norm` has at most one gauge leg.
The retained Yukawa/Majorana terms and independent source rows in (I7)
and (I9) acquire no explicit gauge connection under this substitution.
At first order in this specified `kappa F` redefinition, the two- and
three-gauge local action kernels of the selected packet therefore remain
entirely in `R_J+R_D` and have exactly their previous coefficients.
The previously obtained photon-W+-W- local UV contact consequently is
not removed by this particular normalization.

This is a bounded local-action statement. It does not promote that
contact into an invariant Wilson coefficient, a scattering amplitude,
or an observable: further basis changes, all other coefficient sectors,
source and mass transport, external-leg conventions, finite terms, and
the physical matching calculation still matter. In particular it does
not determine a finite photon amplitude or a physical beta function.

## 6. Reproduction

Run `python -I -S -B -O verify_interacting_canonical.py` here. The checker
uses only the standard library, exact Gaussian rational arithmetic, and
finite Fourier sums with generic off-shell momenta. It implements local
products, derivatives on whole composites, weak matrices, Dirac matrices,
and a formal variation parameter. It imports no earlier calculation.

Three nonconstant-background trials use two active and two sterile
flavors, both weak components, noncommuting gauge matrices, complex
`F,Y_nu`, and a nonzero symmetric complex retained Majorana mass. The
checker varies `F` and `Y_nu` independently because these operator
identities hold for arbitrary fixed matrices; interpretation of the
selected loop coefficient restores `F=C_seed Y_nu*` from the input packet.
No additional loop coefficient is inferred from that algebraic freedom.
The checks include all six placements against (I4)+(I5), the curvature
product rule, local weak variations, differential adjoints, source
transport, actual mass substitutions in both slots, and degree-two and
degree-three gauge contacts. Negative controls omit individual scalar
derivatives, reverse their sign, omit transformed currents or sources,
replace the mass transpose, and discard higher-gauge terms. Counts in
`INTERACTING_CANONICAL_RESULTS.json` describe algebraic controls within
this derivation, not independent physical measurements or loop diagrams.
