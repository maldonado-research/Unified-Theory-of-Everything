# Sourced reduction of the selected hypercharge packet

15 October 2026. This note combines two previously calculated pieces of one
counterterm. It does not calculate an additional loop or a physical anomalous
dimension. Its new result is the explicit source-retaining reduction of their
sum and an off-shell equation-residual certificate for that reduction.

## Input and declared representative

Use the literal action, field orientation, Fourier conventions and positive
diagonal retained Majorana mass `M` of `ACTION_SIGN.md`. Write
`phi=tilde H`, `g_y=g' y_L`, `y_L=-1/2`, `rho=H^dagger H`,
`A=i slash partial`, and `K bar L=-i(D_mu bar L)gamma^mu`.
Only the U(1) connection is included in this calculation:
`[D_mu,D_nu]=i g_y B_mu nu` on `L` and `phi`.

For a flavor matrix `F`, define the unambiguous contraction

    [bar U V Z]_F = sum_(alpha,i) F_(alpha i) bar U_alpha V Z_i.

The field strength and spin matrices may occur inside `V`. In particular a
mass on `Z` becomes `F M`, never `M F`. The selected packet is

    P_F = F_(alpha i) (O_J + O_D)_(alpha i),
    O_D = g_y (Q_sigma/4 + Q_partialB/3).

Here the `Q_sigma,Q_partialB` definitions are those in `ACTION_SIGN.md`;
equivalently `O_D=(Q_Vr+Q_Er)/2+Q_Vk/3` in `OPERATOR_MAP.md`.
The actual cancelling action is

    L_ct = P_F/(16 pi^2 epsilon) + h.c.,   F=C_seed Y_nu*, C_seed^T=C_seed.

No extra factor `1/2` multiplies this mixed operator. The pole-to-action sign
and the one- and two-B matching are inputs already established in the existing
checkpoint. This note must not be counted as another integration of those
diagrams.

Use the same reference-action sources as the earlier notes:

    K bar L = bar S_L,
    -D^2 phi = (m_H^2+2 lambda rho) phi - S_phi,
    partial_nu B^{nu mu} = g' J_Y^mu,
    A N_L = M N_R + S_R,
    A N_R = M N_L + S_L^N.

The quartic potential `-lambda rho^2` is an explicitly declared optional
spectator extension. Set `lambda=0` for the displayed truncated action in
`ACTION_SIGN.md`. Other Yukawa and matter spectators may be included only
with their contributions to these same sources. The hypercharge current
contains the scalar current and every retained charged fermion current.
Sources from further EFT operators must be included if the retained order
requires them. No SU(2) result is asserted here.

The gauge EOM is that of the gauge-invariant classical action before gauge
fixing. Gauge fixing specifies loop propagators; its auxiliary terms are not
part of the gauge-invariant operator equivalence claimed here.

## Combined result

Define the sterile source polynomials

    R = M N_R + S_R,
    T = M^2 N_L + M S_L^N + A S_R,
    U = M^3 N_R + M^2 S_R + M A S_L^N + A^2 S_R.

All matrices in these expressions act on the sterile column in the displayed
order. A sourced representative of the full selected packet is

    P_F ~= (5/6)[bar L phi U]_F
          +(2m_H^2/3+lambda rho/3)[bar L phi R]_F
          +(1/3)[(K bar S_L) phi R]_F
          -(1/6)[bar L S_phi R]_F
          -(1/6)[bar S_L phi T]_F
          +(m_H^2/2)[bar S_L phi N_L]_F
          -(1/2)F_(alpha i)(Q_Vh+Q_Eh)_(alpha i)
          -(g'^2 y_L/6)[bar L phi gamma_mu N_L]_F J_Y^mu
          -(g_y/12)[bar L phi sigma^{mu nu} R]_F B_mu nu.       (1)

The equivalence is modulo a total derivative, the abelian Bianchi identity,
and the declared sourced equations at first order in `F`. It is not an
identity of arbitrary off-shell amplitudes after simply deleting equation
residuals. Equation (1) leaves derivatives on sources, including
`K bar S_L`, explicitly in the representative.

The coefficient of the explicit retained mass dipole in this representative
is therefore

    delta L_ct |explicit mass dipole
       = -g_y/(12*16 pi^2 epsilon)
         (F M)_(alpha i) bar L_alpha phi sigma^{mu nu} N_Ri B_mu nu
         + h.c.                                                (2)

The same coefficient multiplies the paired source
`[bar L phi sigma S_R]_F B`. The isolated transverse packet would give
`-1/4`; the coefficient of the selected sum is `-1/12`. The difference is
the lepton-curvature term in `O_J`. This is a coefficient in the declared
source-retaining representative, not a unique physical dipole coefficient.

For clarity, the complete sterile expansion in (1) consists of the following
ordered terms; this list does not silently set any source to zero:

    (5/6)[bar L phi M^3 N_R]_F
    (5/6)[bar L phi M^2 S_R]_F
    (5/6)[bar L phi M A S_L^N]_F
    (5/6)[bar L phi A^2 S_R]_F
    +(2m_H^2/3+lambda rho/3)[bar L phi M N_R]_F
    +(2m_H^2/3+lambda rho/3)[bar L phi S_R]_F
    +(1/3)[(K bar S_L) phi M N_R]_F
    +(1/3)[(K bar S_L) phi S_R]_F
    -(1/6)[bar L S_phi M N_R]_F
    -(1/6)[bar L S_phi S_R]_F
    -(1/6)[bar S_L phi M^2 N_L]_F
    -(1/6)[bar S_L phi M S_L^N]_F
    -(1/6)[bar S_L phi A S_R]_F
    +(m_H^2/2)[bar S_L phi N_L]_F
    -(1/2)F_(alpha i)(Q_Vh+Q_Eh)_(alpha i)
    -(g'^2 y_L/6)[bar L phi gamma_mu N_L]_F J_Y^mu
    -(g_y/12)[bar L phi sigma^{mu nu} M N_R]_F B_mu nu
    -(g_y/12)[bar L phi sigma^{mu nu} S_R]_F B_mu nu.           (3)

## Derivation of the combined coefficient

The row Dirac-square identity in these conventions is

    K^2 bar L = -D^2 bar L - (g_y/2) bar L B_mu nu sigma^{mu nu}.

Consequently the `O_L/3` term in `O_J` gives

    +(g_y/6) bar L phi sigma^{mu nu} A N_L B_mu nu.

The exact integration-by-parts identity for the already calculated
transverse lift is

    O_D = -(Q_Vh+Q_Eh)/2 - Q_Vk/6
          -(g_y/4) bar L phi sigma^{mu nu} A N_L B_mu nu
          modulo a total derivative.

It follows from the explicit Clifford identity and the abelian Bianchi
identity; it does not use the light-lepton or sterile EOM. With

    Q_Vk = g_y bar L phi gamma_mu N_L partial_nu B^{nu mu},

the current coefficient in (1) is negative. Adding the two sigma terms
before applying the sterile equation gives `1/6-1/4=-1/12`, including the
paired source after `A N_L` is replaced by `M N_R+S_R`.

Substitution of the scalar equation into `O_H/6` gives the independent
coefficient `m_H^2/2+m_H^2/6=2m_H^2/3`. The remaining row and scalar source
terms give (1) without further rearranging sources. The coefficient of the
explicit Yukawa-shaped monomial is the already known
`F[(5/6)M^3+(2/3)m_H^2 M]`; the new combined reduction does not present it
as a newly calculated result.

## Exact residual certificate

Define arbitrary off-shell equation residuals, without setting any to zero:

    e_barL = K bar L - bar S_L,
    e_phi = -D^2 phi - (m_H^2+2lambda rho)phi + S_phi,
    e_B^mu = partial_nu B^{nu mu} - g'J_Y^mu,
    e_R = A N_L - M N_R - S_R,
    e_L = A N_R - M N_L - S_L^N.

Direct differentiation, with constant `M`, gives the exact identities

    A N_L - R = e_R,
    A^2 N_L - T = M e_L + A e_R,
    A^3 N_L - U = M^2 e_R + M A e_L + A^2 e_R.              (4)

Let `P_red` denote the right side of (1). The complete difference, modulo
only a total derivative and the abelian Bianchi identity, is

    P_F-P_red = (1/3)[(K e_barL) phi A N_L]_F
               +(1/6)[bar L e_phi A N_L]_F
               -(1/6)[e_barL phi A^2 N_L]_F
               +(m_H^2/2)[e_barL phi N_L]_F
               -(g_y/6)[bar L phi gamma_mu N_L]_F e_B^mu
               +(5/6)[bar L phi (M^2 e_R+M A e_L+A^2 e_R)]_F
               +(2m_H^2/3+lambda rho/3)[bar L phi e_R]_F
               +(1/3)[(K bar S_L) phi e_R]_F
               -(1/6)[bar L S_phi e_R]_F
               -(1/6)[bar S_L phi (M e_L+A e_R)]_F
               -(g_y/12)[bar L phi sigma^{mu nu} e_R]_F B_mu nu.
                                                                  (5)

This certificate identifies every residual discarded in (1). In particular
it exposes the lost terms if the free sterile equation, free lepton equation,
or free Maxwell equation is substituted while their sources are nonzero.

## Flavor and retained-mass counting

For the renormalizable sterile source `S_R=Y_nu^T X^c`,
`X=phi^dagger L`, the explicit last source of (3) has the active matrix

    -(g_y/12) F Y_nu^T.

Because `bar X_alpha sigma X_beta^c` is antisymmetric under exchanging
the two Grassmann composite slots, its relevant ordered coefficient is

    -(g_y/12) K_F,   K_F=(F Y_nu^T-Y_nu F^T)/2.

For the selected symmetric seed this is
`K_F=(C_seed H_Y^T-H_Y C_seed)/2`, `H_Y=Y_nu Y_nu^dagger`.
That flavor kernel was already derived for the isolated transverse packet;
the new statement is its coefficient `-1/12` in the combined representative.
If summing only independent pairs `alpha<beta`, the coefficient becomes
`-(g_y/6)(K_F)_(alpha beta)`. The ordered-index and independent-pair
normalizations must not be mixed. This is still a source descendant, not a
physical transition moment.

The fixed mass basis is sufficient for (1)-(5). A real orthogonal change of
sterile basis keeps the same equations with real symmetric `M`; the ordered
products transform covariantly. For a general complex symmetric mass,
specify the chiral mass action first. With
`L_mass=-1/2(bar N_L M N_R+h.c.)`, the second equation instead has `M^dagger`:

    A N_L=M N_R+S_R,      A N_R=M^dagger N_L+S_L^N,
    T=M M^dagger N_L+M S_L^N+A S_R,
    U=M M^dagger M N_R+M M^dagger S_R+M A S_L^N+A^2 S_R.

The mass-dipole coefficient remains `F M`. Replacing these alternating
products by `M^2,M^3` before fixing a real mass basis is generally incorrect.
No inverse of any sterile mass is needed in the reduction.

Take `[F]=-3`, and count the retained masses and renormalizable couplings
as order `Lambda^0`. The explicit Yukawa-shaped coefficient scales as
`F M^3` or `F m_H^2 M`; the dimension-six dipole coefficient scales as
`g_y F M` and has mass dimension `-2`, while its inverse-heavy-scale order
remains `Lambda^-3`. The active field-strength source has coefficient
`g_y F Y_nu^T`, and the gauge-current packet has coefficient `g'^2 y_L F`.
They are correlated descendants of the selected packet despite their
different numbers of external fields and coupling powers.

A source from a dimension-`d` interaction with coefficient
`Lambda^{-(d-4)}` multiplies `F` at order `Lambda^{-(d-1)}`. For this isolated
packet and ordinary counting, renormalizable sources suffice through
`Lambda^-3`. This does not establish completeness of the full EFT at that
order: a redundant dimension-six coefficient multiplying a dimension-five
source is also `Lambda^-3` and is not supplied by this selected packet.
If a mass is assigned compensating heavy-scale growth, this counting must
be reconsidered before truncating any source.

## Representative dependence and limits

The same `O_D` can be written using the lepton-source route as

    O_D ~= Q_Vr -(g_y/4)bar S_L phi sigma N_L B
                 +(g'^2 y_L/3)(bar L phi gamma_mu N_L)J_Y^mu.

In that representative the explicit mass dipole visible from `O_J` is
`+g_y F M/6`; the derivative, source and current terms differ in a
correlated manner. This explicit alternative proves why (2) cannot be
advertised as a basis-invariant dipole prediction. Further reduction of
`K bar S_L`, scalar sources, and derivative composite sources can move
coefficients among the displayed structures.

There is a concrete lower-dimensional methodological warning in Di Zhang,
[2504.00792v2, Section 3](https://arxiv.org/html/2504.00792v2#S3): its
[Eq. (13)](https://arxiv.org/html/2504.00792v2#S3.EGx13) specifies field
redefinitions and [Eq. (14)](https://arxiv.org/html/2504.00792v2#S3.EGx14)
tracks correlated parameter shifts; Eqs. (16)-(17) show a flavor-coordinate
ambiguity from applying naive EOM replacements to different forms of a
redundant sterile derivative operator. That work treats renormalizable and
dimension-five running. Its coefficients are not dimension-seven inputs
here. The residual certificate (5) establishes an EOM equivalence class; a
fully indexed field redefinition, canonical normalization and correlated
parameter shifts are still needed before a physical coefficient vector can
be extracted. Neither the `-1/12` arithmetic nor the certificate closes that
remaining problem.

The result supplies neither a complete independent retained-sterile operator
basis nor all diagrams, parameter and field counterterms, evanescent
prescriptions, finite matching terms, a physical renormalization-group
coefficient, or a measurable dipole. The four-dimensional identities are
applied to the declared simple-pole representative. Finite terms from
dimension-dependent operator choices require an additional prescription.

`verify_selected_packet.py` is a standalone exact checker using SymPy and
no repository imports. It checks Clifford/IBP kernels, the row-curvature
sign from derivative jets, the full off-shell residual certificate on
independent nonzero source jets, flavor order, and deliberately wrong
coefficients/source omissions. Its arithmetic checks supplement the
displayed derivation; they do not replace an independent full EFT matching.

From a Python environment with SymPy installed, run

```sh
python -I -B -O verify_selected_packet.py
```

The optimized run is saved in `SELECTED_PACKET_RESULTS.json`. It includes
eight full residual trials with two active and three sterile flavors,
sixteen off-shell Clifford/IBP trials, the direct chiral Grassmann tensor
symmetry, and twenty rejected sign/order/source omissions. These counts
describe exact implementation controls, not independent loop calculations.
