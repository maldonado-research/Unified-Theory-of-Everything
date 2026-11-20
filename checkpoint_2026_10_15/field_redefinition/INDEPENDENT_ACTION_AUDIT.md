# Independent action-level audit

15 October 2026. Independent audit independently derived from Eq. (5) of
`SELECTED_PACKET_REDUCTION.md` and the literal `ACTION_SIGN.md` action.
No producer script or checker was imported. The covectors were sent to the
producer before inspecting its implementation. Their displayed expressions
agree term for term.

This audits an explicit first-order field redefinition for the already
calculated selected packet. It is not another loop calculation and does not
establish an observable, a full EFT basis, or a physical running coefficient.

## Variation conventions and literal sources

Keep the gauge-invariant part of the reference action and its positive real
diagonal retained mass matrix. Define `E_N=e_R+e_L`. The exact first variation,
after integrating derivatives off the variations, is

    δS_inv = ∫ [ e_barL δL + δbarL E_L
                 + δphi† e_phi + e_phi† δphi
                 + e_B^mu δB_mu + bar(E_N) δN ].

The fermionic ordering is part of this statement. Equivalently, a left odd
Euler derivative appears to the right of `δL`, whereas a right derivative
appears to its left. They cannot be interchanged without its Grassmann sign.

For the displayed reference Yukawa interaction, the sources in the residual
definitions are

    (bar S_L)_(alpha a) = bar N_Ri (Y†)_(i alpha) phi_a*,
    (S_phi)_a = -bar N_Ri (Y†)_(i alpha) L_(alpha a),
    S_L^N = Y† phi† L,       S_R = Y^T X^c.

In particular the scalar source is negative. Differentiating the literal
Hermitian-conjugate Yukawa term with respect to `phi†` gives
`-bar N_R Y† L`. With the source convention of Eq. (5),
`e_phi=-D²phi-(m_H²+2lambda rho)phi+S_phi` is exactly that Euler derivative.
The result does not permit inserting an independently chosen scalar source
sign into the same certificate.

The factor `1/2` in the Majorana reference kinetic/mass action is removed by
the two constrained variations. `δN` must itself be Majorana. There is no
additional `1/2` in `bar(E_N)δN`.

## Independent covectors

Let `κ` be the real pole prefactor and let all covectors below contain `F`.
Use `K0 Q=-i(partial_mu Q)gamma^mu` on a gauge-neutral row. Define

    Q_i = Σ_(alpha,a) F_(alpha i) bar L_(alpha a) phi_a,
    Q_Si = Σ_(alpha,a) F_(alpha i) (bar S_L)_(alpha a) phi_a,
    Q_KSi = Σ_(alpha,a) F_(alpha i) (K bar S_L)_(alpha a) phi_a,
    Q_phiSi = Σ_(alpha,a) F_(alpha i) bar L_(alpha a) (S_phi)_a,
    Q_sigmaBi = Σ_(alpha,a) F_(alpha i) bar L_(alpha a) phi_a
                  sigma^{mu nu} B_mu nu.

The sourced packet difference and its conjugate factorize, modulo a total
derivative and the already-declared abelian Bianchi identity, as

    P_F-P_red+h.c.
      = e_barL ℓ + t e_phi + c_Bmu e_B^mu
        + xbar_i e_Ri + ybar_i e_Li + h.c.,

where

    ℓ_(alpha a) = Σ_i F_(alpha i) [
       i gamma^mu(D_mu phi_a) A N_Li / 3
       + phi_a A² N_Li / 6 + m_H² phi_a N_Li / 2 ],
    t_a = Σ_(alpha,i) F_(alpha i) bar L_(alpha a) A N_Li / 6,
    c_Bmu = -g_y Σ_(alpha,i,a) F_(alpha i)
                    bar L_(alpha a) phi_a gamma_mu N_Li / 6,

    xbar = (5/6)(Q M²+K0²Q)
           +(2m_H²/3+lambda rho/3)Q
           +Q_KS/3-Q_phiS/6-K0 Q_S/6-g_y Q_sigmaB/12,
    ybar = (5/6)(K0 Q)M-Q_S M/6.

The lepton expression follows from

    ∫(K e_barL) phi A N_L
      = ∫e_barL i gamma^mu D_mu(phi A N_L).

The derivative on `phi` is essential, and combines the coefficient of
`phi A²N_L` into `1/3-1/6=1/6`. Likewise, `K0²Q` and `K0 Q_S` differentiate
their entire gauge-neutral composites. `Q_KS` differentiates only the
charged row source before multiplication by `phi`; it cannot silently be
identified with `K0 Q_S`. Every mass acts to the right on the sterile index.

The chiralities are consistent: `ℓ=P_Lℓ`, `xbar=xbar P_R`, and
`ybar=ybar P_L`. If `x=(xbar)^D` and `y=(ybar)^D` are Dirac adjoint columns,
then `x=P_Lx` and `y=P_Ry`.

## Explicit shifts and reality

For `Phi_old=Phi_new+δPhi`, the shifts cancelling the positive redundant
packet are

    δL = -κℓ,
    δphi† = -κt,               δphi = -κt†,
    δB_mu = -κ(c_Bmu+c_Bmu†),
    δN_L = -κ(x+y^c),          δN_R = -κ(y+x^c).

Thus `δN=-κ(x+y+x^c+y^c)` is manifestly Majorana. A second `1/2` would cancel
only half of the redundant action. Dropping the charge-conjugate part does
not preserve the Majorana constraint.

With `phi=epsilon H*`, `epsilon=[[0,1],[-1,0]]`, the scalar shifts are
constrained, not independent choices:

    δH = +κ epsilon t^T,       δH* = +κ epsilon t†.

These obey `δphi†=-(δH)^T epsilon=-κt`. The real gauge field requires the
sum `c_B+c_B†`; shifting by the complex current alone is invalid.

The resulting pullback is

    S_inv[Phi_old] + κ∫(P_F+h.c.)
      = S_inv[Phi_new] + κ∫(P_red+h.c.) + O(κ²)

at the declared inverse-heavy-scale order. Varying the original
dimension-seven seed produces products of EFT coefficients that require a
separate retained-order decision; it is not part of this first-order
renormalizable-action identity.

The literal action also contains gauge fixing. A claim about its full
gauge-fixed pullback must include

    δS_gf = -∫(partial·B)partial·δB / xi

and any corresponding BRST/ghost bookkeeping. The stated gauge-invariant
classical operator equivalence does not itself include those terms.

## Canonical-normalization sanity checks

At the unbroken zero-field origin, every shift has field degree at least
two. The selected packet therefore changes no quadratic Hessian at that
origin at order `F`. This does not establish canonical normalization about
a nonzero scalar background.

Independently projecting the six original `O_J` monomials onto the neutral
constant background `phi=(u,0)` with real `u`, and setting gauge fields to
zero, gives

    P_F|quadratic = u F bar nu_L A(-Box+m_H²)N_L+h.c.

The cubic-derivative coefficients add as `5/6+1/3-1/6=1`; the mass-squared
coefficients add as `1/2+1/2=1`. The scalar-derivative term vanishes for a
constant background. This off-diagonal mixed bilinear has no additional
Majorana half. A nonstationary background can still define this algebraic
projection; interpreting physical masses requires the vacuum conditions.

There is a further nonstationary-background subtlety if comparing this
projection to the field-redefined representative. At `phi=(u,0)`, the scalar
tadpole is `e_phi^bg=-(m_H²+2lambda u²)phi`. The bilinear scalar shift then
contributes to the quadratic fermionic action:

    δphi† e_phi^bg+h.c.
      = +κ u(m_H²+2lambda u²)F bar nu_L A N_L/6+h.c.

An action Hessian under a nonlinear coordinate change contains both the
Jacobian congruence and the equation-times-second-derivative term,
`H'_ij=J_ai H_ab J_bj+E_a f^a_,ij`. Keeping only linear fermion shifts misses
the displayed term on a nonstationary scalar background. It vanishes at a
stationary broken vacuum. The direct original-packet projection remains a
valid algebraic off-shell statement; comparing it with `P_red` requires this
tadpole term or the corresponding vacuum condition.

For left Weyl coordinates `q=(nu_L,N_L)`, distinguish the two equivalent
mass conventions. The four-component orientation

    -1/2[bar q_L mathcalM q_L^c+h.c.]

has `mathcalM=[[0,uY],[uY^T,M]]`. The left holomorphic convention

    -1/2[q^T epsilon M_holo q+h.c.]

used by the canonical-normalization calculation has
`M_holo=[[0,uY*],[uY†,M*]]`, with the same declared epsilon convention.
The original retained diagonal mass is real, so its conjugation is invisible
in this restricted mass basis. A generic complex Yukawa matrix exposes the
otherwise easy-to-miss conjugation.

## Executable controls and scope

`independent_action_audit.py` implements its own finite Grassmann exterior
algebra, Weyl gamma matrices, conjugation, left/right odd derivatives, and
exact rational Fourier tests. The optimized isolated run passed 33 checks
and rejected 21 deliberately wrong alternatives. These include a spurious
Majorana half, missing conjugate shift, reversed conjugate sign, confusing
left and right odd derivatives, the scalar source sign, the constrained
tilde-map sign, omitted derivatives on `phi`, wrong first/second integration
by parts, commuted flavor masses, and a complex gauge-field shift.

These are independent action-map lemmas and generic off-shell controls.
They do not independently recompute the full prior residual certificate,
the selected UV diagrams, a complete canonical operator basis, finite
matching, the Jacobian/quantum change of variables, or physical running.
The concrete advance supported is an explicit local, gauge-covariant,
chirality- and reality-preserving first-order classical field redefinition
for the declared selected packet, with its correlated source terms retained.

## Independent acceptance of the neutral canonical-coordinate note

The final canonical-coordinate derivation at `CANONICAL_NORMALIZATION.md`
was reviewed separately, including its equations (N6)-(N12). The following
claims agree with an independent matrix derivation under the note's stated
first-order, local-derivative, and invertibility assumptions:

* The left-holomorphic reference mass has `uY*` in its active-sterile block;
  its conjugate has `uY` in the antiholomorphic four-component orientation.
* With doubled coordinates `T_N=diag(T,T*)`, the rest-frame kernel
  `D=[[EK,-Mcal†],[-Mcal,EK^T]]` transforms as `T_N† D T_N`.
* Its determinant acquires `|det T|^4=1-2κ tr Z+O(κ²)`. Eliminating the
  conjugate equation gives `Q(s)=sK-Mcal†K^(-T)Mcal`, which transforms by
  `T† Q T` and acquires `|det T|²=1-κ tr Z+O(κ²)` in its determinant.
  These powers differ for the stated reason; neither requires a mass inverse.
* The full off-shell source-dressed kernel is invariant only when the source
  row transforms too: `W_c=W T_N`. A mass or pole comparison alone can miss
  omitted source transport.
* Static Schur reduction requires invertible retained sterile block `D_s`.
  Both `G_eff=G_a-G_s D_s^(-1)B^T` and the local source contact
  `C_J=G_s D_s^(-1)G_s^T` are necessary. The inverse active Schur block is
  used only for the explicitly invertible comparison benchmark.
* The neutral original-packet projection, derivative character of `Z(p²)`,
  perturbative interpretation of poles/rank, and nonstationary-background
  tadpole caveat are stated correctly.

No mathematical discrepancy was found in those identities. This acceptance
is limited to their declared selected neutral quadratic problem and its
coordinate/source transport. It does not complete the gauge-covariant full
EFT canonical basis, independent diagram matching, finite matching, or a
physical anomalous dimension.
