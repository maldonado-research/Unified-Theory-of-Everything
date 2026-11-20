# Composing the sourced and locally normalized action coordinates

16 October 2026. Input candidate `5f9889b` and the checked maps in
[`ELECTROWEAK_REDUCTION.md`](../../checkpoint_2026_10_15/interacting_reduction/ELECTROWEAK_REDUCTION.md)
and [`INTERACTING_CANONICAL.md`](../../checkpoint_2026_10_15/interacting_reduction/INTERACTING_CANONICAL.md).
This note gives their explicit coordinate bridge, then checks the actual
neutral quadratic actions independently by substitution. The scalar
equation-times-second-derivative term is retained off shell.

The full interacting bridge below follows from the previously certified
residual and the ordered first variation. The new component calculation
tests the declared neutral neutrino/Yukawa/Majorana/scalar action; it does
not re-expand every interacting electroweak vertex. All statements are
classical first-order action identities at the selected inverse-scale order.

## 1. Two maps from one original action

Use `D=partial+i Acal`, `phi=tilde H`, `F=C_seed Y_nu*`, `P0=−partial^2`,
`P=P0+m_H^2`, `rho=phi^dagger phi`, `T^I=tau^I/2`, `y_L=−1/2`, and
`kappa=1/(16 pi^2 epsilon)`. The retained sterile mass `M` is real positive
diagonal; constant real symmetric `M` also permits the ordered expressions
below. Variational interaction sources always mean those of the declared
gauge-invariant renormalizable action `S0`.

Write both substitutions with their new coordinates on the right:

    Xi_old=Xi_S+kappa vS[Xi_S],
    Xi_old=Xi_C+kappa vC[Xi_C].                            (F1)

The first is the source-reducing map. Its components are

    vS_L=−ell,       vS_phi*=−t,       vS_phi=−t^dagger,
    vS_B,mu=+g'y_L(j_mu+j_mu^dagger)/6,
    vS_W,mu^I=+g(j_mu^I+j_mu^{I dagger})/6,
    vS_NL=−(x+y^c),       vS_NR=−(y+x^c).                (F2)

Here `ell,t,x,y,j,jI` are exactly the fully indexed sourced coefficients
of the input map, including all retained interaction sources. In
particular `x` is left chiral, `y` right chiral, and neither source nor
Majorana completion is divided by another two.

The second map removes the chosen normalization anticommutator from the
original packet:

    vC_L=−Phi F P N_L/2,
    vC_NL=−P(F^dagger Phi^dagger L)/2,
    vC_NR=(vC_NL)^c,
    vC_phi=vC_phi*=vC_B=vC_W=0.                           (F3)

`Phi=phi tensor 1_active`; the derivative in the second line acts on
the complete neutral composite. This is the particular normalization
choice already specified in the input. It leaves interaction remainders
and induced mass/source terms; it is not full EFT canonicalization.

## 2. Explicit all-field bridge

Equating (F1) gives, at first order,

    Xi_S=Xi_C+kappa b[Xi_C],       b=vC−vS.              (F4)

Every coefficient on the right of the following formulas is evaluated on
`Xi_C`. The lepton bridge simplifies substantially:

    bL_alpha^a=(1/3)sum_i F_(alpha i)[
        i gamma^mu(D_mu phi)^a A N_Li−phi^a P0 N_Li],
    A=i slash partial.                                  (F5)

The `m_H^2` terms cancel, as do the corresponding parts of the two
derivative normalizations. Its explicit Dirac-adjoint row is

    bbarL_alpha^a=−(1/3)sum_i F*_(alpha i)[
        i overline(A N_Li)gamma^mu(D_mu phi)^{a*}
          +overline(P0 N_Li)(phi^a)*].                  (F6)

No integration by parts is hidden in (F6); it is the conjugate field
substitution itself. Scalar and gauge components are

    bphi_a*=+t_a,       bphi^a=+t_a^dagger,
    t_a=(1/6)sum_(alpha i)F_(alpha i)bar L_alpha^a A N_Li,
    bB_mu=−g'y_L(j_mu+j_mu^dagger)/6,
    bW_mu^I=−g(j_mu^I+j_mu^{I dagger})/6.                (F7)

The gauge shifts are real; the weak shift is adjoint covariant. For
`phi=epsilon H*`, the scalar constraint gives `bH=−epsilon t^T`.

For the sterile component one equivalent explicit form is

    bN_Li=x_i+y_i^c−P(F^dagger Phi^dagger L)_i/2,
    bN_Ri=(bN_Li)^c,
    bbarN_Li=overline(bN_Li),
    bbarN_Ri=overline((bN_Li)^c).                         (F8)

To make every derivative and source in (F8) visible, define only the
neutral composites

    U_i=sum_(alpha a)F*_(alpha i)(phi^a)* L_alpha^a,
    U_M2,i=sum_(alpha a)(F M^2)*_(alpha i)(phi^a)*L_alpha^a.

Then the full left component is

    bN_Li = (5/6)U_M2,i+(1/3)P0 U_i
             +(m_H^2/6+lambda rho/3)U_i+y_i^c
       +(i/6)sum_(alpha a)F*_(alpha i)[
           (phi^a)*gamma^mu(D_mu S_L,alpha)^a
           −(D_mu phi)^{a*}gamma^mu S_L,alpha^a]
       −(1/6)sum_(alpha a)F*_(alpha i)(S_phi^a)*L_alpha^a
       −(1/12)sum_(alpha a b)F*_(alpha i)(phi^b)*
              (Fcal_mu nu)_ba sigma^{mu nu}L_alpha^a,    (F9)

    y_i=(5i/6)sum_(alpha a)(F M)*_(alpha i)gamma^mu[
           (phi^a)*(D_mu L_alpha)^a+(D_mu phi)^{a*}L_alpha^a]
          −(1/6)sum_(alpha a)(F M)*_(alpha i)
                                      (phi^a)*S_L,alpha^a.

`P0 U` differentiates both scalar and lepton. Equation (F9) follows by
combining the covariant product derivatives in `x` into the ordinary
derivative of a neutral composite, then subtracting (F3). It keeps the
ordered mass products and the `ba` curvature index after conjugation.
The right component and both adjoints are fixed by (F8), not independent
choices. Thus the entire bridge preserves `N=N^c` without an additional
half. All components are local polynomial changes with a perturbative
inverse obtained by reversing the first-order bridge sign.

## 3. Full interacting action identity and its limits

Let `P_orig` denote the original selected packet plus its Hermitian conjugate.
With the same convention for the reduced representative, the existing
residual certificate says

    P_S=P_orig+delta_vS S0,
    P_C=P_orig+delta_vC S0.                              (F10)

The ordered first variation is the one derived from the literal action:
it includes `e_barL delta L+delta barL E_L`, both scalar slots, real
gauge variations, and the constrained Majorana variation. For example,
`bar E_N delta N` already includes the two Majorana slots that cancel
the half in the reference quadratic action. Use that graded variation
without commuting odd coefficients past odd Euler derivatives.

Substitute the actual bridge into the sourced action:

    (S0+kappa P_S)[Xi_C+kappa b]
       =S0[Xi_C]+kappa(P_S+delta_b S0)[Xi_C]
       =S0[Xi_C]+kappa(P_orig+delta_vC S0)[Xi_C]
       =(S0+kappa P_C)[Xi_C]+O(kappa^2 F^2).             (F11)

Derivatives of a shift functional first enter this composition at second
order, so every first-order coefficient may be evaluated on `Xi_C`.
Equation (F11) is not obtained by applying the original normalizer
unchanged to `P_S`. That would add `vC` to the source-reducing direction;
the bridge requires subtracting `vS`.

The prior residual establishes (F10) modulo its declared boundary and
Bianchi identities. Hence this is a full selected-packet classical
interacting bridge based on that certificate. It does not constitute
another independent component verification of all interactions. Nor does
it remove separate `C6 C5` contributions or other retained-operator
variations. The dimension-seven seed itself is excluded from `S0`; its
variation is a product of seed coefficients beyond this selected
`Lambda^-3` truncation at fixed retained-mass counting.

## 4. Actual neutral quadratic action, projected two independent ways

Set `phi=(u,0)` with real constant `u`, `B=W=0`, and zero fermion
background. Keep the literal neutrino interaction sources only:

    S_L=u Y N_R,        S_L^N=u Y^dagger nu_L,
    S_R=(S_L^N)^c,      bar S_L=overline(S_L),
    c_u=2m_H^2/3+lambda u^2/3,
    a_u=m_H^2+2lambda u^2.                               (F12)

No stationarity is imposed unless expressly stated. The reference
quadratic action is

    S0^(2)=integral[bar nu_L A nu_L
          +(1/2)bar N(A−M)N
          −u bar nu_L Y N_R−u bar N_R Y^dagger nu_L].    (F13)

The original six placements give

    P_orig^(2)=u sum F_(alpha i)bar nu_Lalpha
                               A(P0+m_H^2)N_Li+h.c.     (F14)

There is no extra Majorana half on the mixed kinetic term. Independently
projecting the reduced representative gives the explicit quadratic kernel

    R=M N_R+S_R,
    T=M^2 N_L+M S_L^N+A S_R,
    V=M^3 N_R+M^2 S_R+M A S_L^N+P0 S_R,

    P_S^(2)=u sum_(alpha i)F_(alpha i)[
         (5/6)bar nu_Lalpha V_i+c_u bar nu_Lalpha R_i
         +(1/3)(K bar S_Lalpha)R_i
         −(1/6)bar S_Lalpha T_i
         +(m_H^2/2)bar S_Lalpha N_Li]+h.c.,
    K bar psi=−i(partial_mu bar psi)gamma^mu.             (F15)

These are the five surviving types of terms from the actual sourced
representative, not a definition in terms of a transformed Hessian.
Terms containing `S_phi` there have fermion degree four, and gauge-current
products have no degree-two part on this background.

The fermion-linear source coefficients used for a second, direct action
substitution are

    ell=u F(P0/6+m_H^2/2)N_L,
    x_i=u sum_alpha[(5/6)(F M^2)*_(alpha i)nu_Lalpha
                +F*_(alpha i)(5P0/6+c_u)nu_Lalpha]
             +(u^2/6)sum_(alpha j)F*_(alpha i)Y_(alpha j)A N_Rj,
    y_i=(5u/6)sum_alpha(F M)*_(alpha i)A nu_Lalpha
             −(u^2/6)sum_(alpha j)(F M)*_(alpha i)Y_(alpha j)N_Rj.
                                                               (F16)

Use `vS_nu=−ell`, `vS_NL=−(x+y^c)`, and their conjugates. The scalar
shift is bilinear and must still be substituted in the **full potential**,
although it has no fermion-linear part. Equation (F5) yields the especially
simple neutral bridge `bnu=−u F P0 N_L/3`.

## 5. The nonlinear Hessian term is nonzero off shell

For `V=m_H^2 rho+lambda rho^2`, the scalar gradient is

    e_phi^background=−a_u phi.

Substituting `vS_phi*=−t` and its conjugate in `−V` produces

    P_tad^(2)=+a_u u(t_1+t_1^dagger)
             =+(a_u u/6)sum F bar nu_L A N_L+h.c.       (F17)

The bridge has the opposite scalar shift and produces `−P_tad^(2)`.
The normalizing map has no scalar shift. Gauge shifts begin at fermion
degree two, but their background Euler derivatives vanish at `B=W=0`,
so they contribute no neutral quadratic action here.

To state the Hessians without an unmanageable printed matrix, write the
integrated quadratic action as `theta^T H theta/2` in real exterior
component generators. Its matrix is antisymmetric. Let `H0` be (F13),
`Hp` the coefficient matrix of (F14), `Hs` that of the independently
projected (F15), and `Hc` the first-order coefficient after directly
substituting (F3) in the original action. Let `JS,JC` be the actual
fermion-linear Jacobians extracted from (F16) and (F3), and `Htad` the
matrix of (F17). Then direct substitution gives

    Hs=Hp+JS^T H0+H0 JS+Htad,
    Hc=Hp+JC^T H0+H0 JC,
    Hc=Hs+(JC−JS)^T H0+H0(JC−JS)−Htad.                (F18)

The checker evaluates each entire component matrix from both sides, rather
than defining `Hs` by its congruence. It also obtains (F17) from direct
substitution of the bilinear scalar shift into `−m_H^2 rho−lambda rho^2`.
This realizes the graded version of the classical chain rule

    H'_ij=J_ai H_ab J_bj+E_a f^a_,ij.

Keeping only a linear fermion Jacobian misses the last term away from
stationarity. Its disappearance requires `a_u=0` or a separately specified
external source whose own pullback is also retained.

For completeness, the normalized neutral kinetic perturbation cancels
exactly. Its remaining quadratic coefficient is the induced mass action.
In the declared holomorphic left-Weyl convention,

    M0=[0 uY*; uY^dagger M],
    Z=[0 u(s+m_H^2)F; u(s+m_H^2)F^dagger 0],
    delta M_C=−(Z^T M0+M0 Z)/2,       s=p^2.             (F19)

Its explicit blocks are

    delta A=−u^2(s+m_H^2)(F*Y^dagger+Y*F^dagger)/2,
    delta B=−u(s+m_H^2)F* M/2,
    delta D=−u^2(s+m_H^2)(F^T Y*+Y^dagger F)/2.

The real nondegenerate retained `M` is never inverted. These coordinate
coefficients, including their momentum dependence, are not physical mass
shifts read in isolation.

## 6. External sources also require the bridge

External probe rows are distinct from the variational interaction sources
`S_L,S_phi,S_R` in (F12). Let an external source functional be

    S_ext=integral[bar eta_L L+bar eta_N N_L+h.c.
                   +J_phi* phi+phi* J_phi+J_B^mu B_mu+J_W^{I mu}W_mu^I].

Under the bridge it includes the composite source vertices

    delta_b S_ext=kappa integral[
        bar eta_L bL+bar eta_N bN_L+h.c.
        +J_phi* t^dagger+t J_phi
        −g'y_L J_B^mu(j_mu+j_mu^dagger)/6
        −g J_W^{I mu}(j_mu^I+j_mu^{I dagger})/6].         (F20)

The fermion source rows have the chirality needed to pair their left
columns. Derivative operators in `bL,bN` act on their whole composites;
integration by parts transfers them to source rows as well. Coupling only
to unchanged elementary fields loses (F20). Scalar and gauge probes are
especially relevant: they multiply quadratic fermion composites and can
therefore change an off-shell fermion Hessian even without a fermion-linear
boson shift.

The benchmark Hessians (F13)-(F18) contain **no external background source**.
As a separate tested illustration, adding the fixed real linear source
`J0(phi+phi*)`, with `J0=a_u u` at the chosen background, cancels its
scalar gradient. Pulling back that source contributes
`−a_u u(t+t^dagger)`, precisely canceling (F17). Holding that source term
in elementary form after a field change would miss this cancellation.
It is a different, explicitly declared externally supported action, not
an implicit stationarity assumption for the benchmark.

## 7. Reproduction, controls, and remaining scope

Run `python -I -S -B -O verify_field_map_composition.py` here. It uses
only the standard library, exact Gaussian rational coefficients, explicit
Weyl gamma matrices, derivative Fourier modes, and a finite Grassmann
exterior algebra. Physical Weyl fields have independent opposite-momentum
modes and real/imaginary component generators; the sterile four-component
field is completed by explicit charge conjugation.

Three generic spatial-momentum trials use two active and two sterile
flavors, complex noncommuting `F,Y`, real retained masses `2,5`, real
`u=3/2`, and `lambda=2/5`. The scalar-gradient factors are `62/15`, `7/15`,
and `0`. The last case is expressly stationary. Every comparison evaluates
the full 32-component real fermionic quadratic matrix or its equivalent
exterior polynomial. No on-shell or mass-inverse shortcut is used.

The script compares the literal original placements, independent sourced
projection, direct kinetic/Yukawa/Majorana/scalar substitutions, actual
Jacobians and nonlinear Hessian term, normalizing kinetic cancellation,
and independent odd and even source probes. Negative controls omit the
tadpole, scalar/gauge composite sources, induced mass or external-source
transport; reverse the scalar bridge; use the sum of maps; apply the old
normalizer directly to the reduced action; divide the Majorana map by
two; delete its conjugate channel; or misorder the retained mass.

The full interacting identity uses the prior gauge-invariant residual.
Gauge fixing, ghost/BRST terms, the quantum functional measure, and the
ordinary Legendre effective action require additional bookkeeping; this
note does not assert that a nonlinear classical map makes an ordinary
quantum 1PI effective action transform by naive substitution. Full basis
matching, lower-order insertions, evanescent/finite prescriptions, and
physical running remain open. No loop, observable, novelty, or breakthrough
claim follows from this classical coordinate-completion calculation.
