# Independent electroweak source and field-shift audit

15 October 2026. Input commit `d822e9f4d70ad5b455436f40abc3c2b45909bbef`
(candidate.6). This derivation used the existing non-Abelian local action,
the hypercharge source certificate, and the prior Abelian field-redefinition
note. It was derived before reading the producer's electroweak reduction.
No new loop residue is calculated here.

The matrix extension of the source reduction is valid. There is an important
qualification: the non-Abelian sterile-source tensor cannot in general be
replaced by the flavor-antisymmetric part of `F Y_nu^T`. Weak-slot exchange
must be included. A charged component survives even for one active flavor.

## Conventions and independently derived identities

Use `D_mu=partial_mu+i A_mu`,
`A_mu=g' y_L B_mu 1_2+g T^I W_mu^I`, `y_L=-1/2`, `T^I=tau^I/2`,
and `[D_mu,D_nu]=i Fcal_mu nu`, with

    Fcal_mu nu = partial_mu A_nu-partial_nu A_mu+i[A_mu,A_nu],
    Dcal_rho Fcal_mu nu = partial_rho Fcal_mu nu+i[A_rho,Fcal_mu nu].

Throughout this note `F` without the suffix `cal` is the active-by-sterile
flavor coefficient. The gamma convention is `sigma^{mu nu}=i[gamma^mu,gamma^nu]/2`
and the metric is `(+---)`. `phi=tilde H`, `A=i slash partial` on neutral
sterile fields, and `K bar L=-i(D_mu bar L)gamma^mu`.

The commutator Jacobi identity gives
`Dcal_[rho Fcal_mu nu]=0`; replacing `Dcal` by `partial` is false. In components,

    W^K_mu nu = partial_mu W^K_nu-partial_nu W^K_mu
                -g epsilon^{IJK} W^I_mu W^J_nu,
    (Dadj_nu W^{nu mu})^K = partial_nu W^{K nu mu}
                -g epsilon^{IJK} W^I_nu W^{J nu mu}.

The two minus signs follow from the plus connection in `D`. Varying the
canonical kinetic action `-W^I_mu nu W^{I mu nu}/4` gives precisely
`(Dadj_nu W^{nu mu})^I`, with a positive coefficient in the first variation.
The matter action gives `-g J_I^mu`, so

    e_W^{I mu} = (Dadj_nu W^{nu mu})^I-g J_I^mu,
    e_B^mu = partial_nu B^{nu mu}-g' J_Y^mu,
    Dcal_nu Fcal^{nu mu}
       = g'y_L 1_2(e_B^mu+g'J_Y^mu)
         +g T^I(e_W^{I mu}+g J_I^mu).                         (A)

For retained scalar and fermion doublets, the relevant current contributions
are `i phi^dagger T^I <->D^mu phi` and `bar f gamma^mu T^I f`;
hypercharge weights replace `T^I` in the Abelian current. The scalar current
contains the connection. All retained charged fields contribute. These are
gauge-invariant classical variational equations before gauge fixing.

Define matrix monomials with the weak placement explicitly fixed:

    Q_hV = -bar L gamma^mu Fcal_mu nu D^nu phi N_L,
    Q_hE = -1/2 bar L G^{mu nu rho} Fcal_mu nu D_rho phi N_L,
    Q_div = bar L gamma_mu (Dcal_nu Fcal^{nu mu}) phi N_L,
    Q_AN = bar L sigma^{mu nu} Fcal_mu nu phi A N_L,
    O_D = i/4 (D_rho bar L) sigma^{mu nu}gamma^rho
                         Fcal_mu nu phi N_L + Q_div/3,

where `G^{mu nu rho}=gamma^[mu gamma^nu gamma^rho]` has unit-weight
antisymmetrization. Weak matrices commute with gamma matrices, but their
positions relative to weak fields and the connection cannot be exchanged.
Direct differentiation of the neutral boundary current gives

    O_D + (Q_hV+Q_hE)/2 + Q_div/6 + Q_AN/4
      = partial_rho [i/4 bar L sigma^{mu nu}gamma^rho
                             Fcal_mu nu phi N_L].              (B)

The `Dcal Fcal` term in this differentiation becomes `Q_div/2` by the
Clifford identity and covariant Bianchi identity. The connection terms cancel
between the row, curvature, and scalar factors of the neutral product.
This is an exact local identity, with no matter equations imposed.

The row Dirac-square identity is

    K^2 bar L = -D^2 bar L-1/2 bar L Fcal_mu nu sigma^{mu nu}.   (C)

Consequently the row-Laplacian term of `O_J`, whose coefficient is `1/3`,
supplies `+Q_AN/6`.
Combining it with (B) gives `1/6-1/4=-1/12`, for both the retained mass
and sterile-source tensor. No commutator or additional separately chosen
two-gauge coefficient remains to be added to this matrix identity.

## Reduced packet and residual implementation

Use the same complete variational sources as the input certificate and set

    R=M N_R+S_R,
    T=M^2 N_L+M S_L^N+A S_R,
    U=M^3 N_R+M^2 S_R+M A S_L^N+A^2 S_R,
    c=2m_H^2/3+lambda rho/3,   rho=phi^dagger phi.

Here `M` is constant real symmetric, including the positive diagonal input
convention; its flavor order cannot be reversed. Brackets denote contraction
with `F_alpha i` as in the input note. The source-retaining representative is

    P_red = 5/6 [bar L phi U]_F + c[bar L phi R]_F
        +1/3[(K bar S_L)phi R]_F-1/6[bar L S_phi R]_F
        -1/6[bar S_L phi T]_F+m_H^2/2[bar S_L phi N_L]_F
        -1/2 F_alpha i(Q_hV+Q_hE)_alpha i
        -1/6(g'^2 y_L j_mu J_Y^mu+g^2 j_mu^I J_I^mu)
        -1/12[bar L sigma^{mu nu}Fcal_mu nu phi R]_F,           (D)

where `j_mu=[bar L phi gamma_mu N_L]_F` and
`j_mu^I=[bar L T^I phi gamma_mu N_L]_F`. In particular, the SU(2) current
coefficient contains `g^2 T^I`, not `g T^I` or `g^2 tau^I`.

The full residual certificate is the prior Abelian certificate with full
covariant derivatives, with its gauge residual replaced by

    -1/6(g'y_L j_mu e_B^mu+g j_mu^I e_W^{I mu}),

and its last tensor residual replaced by

    -1/12[bar L sigma^{mu nu}Fcal_mu nu phi e_R]_F.

Every other residual and sterile polynomial is unchanged in form. This
follows from (B), (C), and the unchanged neutral sterile derivative algebra;
the sources themselves are those of the full declared electroweak action.

For completeness, it can be implemented without an inverse mass. Put
`kappa=1/(16 pi^2 epsilon)`. In an old-field-as-function-of-new-fields
substitution, use

    L_old=L-kappa ell,                 phi*_old=phi*-kappa t,
    B_mu,old=B_mu+kappa g'y_L(j_mu+j_mu^dagger)/6,
    W_mu,old^I=W_mu^I+kappa g(j_mu^I+(j_mu^I)^dagger)/6,
    N_old=N-kappa(w+w^c),              w=x+y.                 (E)

The conjugate scalar and lepton substitutions are required. The light
coefficients are

    ell_alpha^a=sum_i F_alpha i [i/3 gamma^mu(D_mu phi)^a A N_Li
                     +phi^a A^2 N_Li/6+m_H^2 phi^a N_Li/2],
    t_a=sum_alpha i F_alpha i bar L_alpha^a A N_Li/6.

The sterile coefficients are most unambiguously specified by their
integrated row covectors. Define neutral rows
`Omega_i=sum_alpha a F_alpha i bar L_alpha^a phi^a`,
`Xi_i=sum_alpha a F_alpha i bar S_L,alpha^a phi^a`, and
`K0=-i partial_mu(.) gamma^mu`. Then

    bar x_i=5/6[(Omega M^2)_i+K0^2 Omega_i]+c Omega_i
        +1/3 sum_alpha a F_alpha i(K bar S_L,alpha)^a phi^a
        -1/6 sum_alpha a F_alpha i bar L_alpha^a S_phi^a
        -K0 Xi_i/6
        -1/12 sum_alpha ab F_alpha i bar L_alpha^a
                              sigma^{mu nu}Fcal_mu nu,ab phi^b,
    bar y_i=5/6(K0 Omega)_j M_ji-Xi_j M_ji/6.                 (F)

Set `x=gamma^0 bar x^dagger`, `y=gamma^0 bar y^dagger`, and use
`w^c=C bar w^T`, `C=i gamma^2 gamma^0`. This explicitly fixes all
conjugations. In particular the curvature term of the column `x_i` is

    -1/12 sum_alpha ab F_alpha i^* (phi^a)^*
                         Fcal_mu nu,ab sigma^{mu nu}L_alpha^b. (G)

Transposing `Fcal` in (G) is wrong. The connection and curvature are
Hermitian for real gauge potentials, so their matrix adjoints reverse the
weak-slot order exactly as displayed.

The gauge shifts in (E) are real, and the W shift transforms in the real
adjoint representation. The Majorana shift preserves `N=N^c`. The
variation of the Majorana action is `overline(delta N)(e_R+e_L)` with no
extra `1/2`; both slots of the quadratic action contribute. For Grassmann
spinors and a Majorana residual `E`, `bar E w=overline(w^c)E`, so the
sterile residual plus its Hermitian conjugate is precisely
`overline(w+w^c)E`. A half-sized shift leaves half the packet behind.

Equations (A)-(G) establish the first-order identity for the declared
gauge-invariant `S_0+S_ct`. Gauge fixing must also be transformed if working
with a gauge-fixed generating functional. Varying other EFT terms produces
the separately counted `kappa F C_d` terms; these are not automatically
quadratic in `kappa`. No canonical normalization or complete EFT coefficient
vector is supplied by this identity.

## Grassmann exchange includes the weak slots

For `S_R=Y_nu^T X^c`, `X_alpha=phi^dagger L_alpha`, define
`Kfl=F Y_nu^T`. The source term associated with a fixed generator `T^I` is

    -g/12 (Kfl)_alpha beta
       (bar L_alpha T^I phi) sigma^{mu nu}X_beta^c W^I_mu nu. (H)

The weak tensor before Grassmann antisymmetrization is
`P_ab=(T^I phi)_a phi_b`. The spin coefficient `sigma^{mu nu}C` is
symmetric, including its same-chirality projected block. Anticommutation
therefore projects the *combined* active/weak coefficient:

    [Kfl tensor P]_antisym
      = (Kfl)_A tensor P_S+(Kfl)_S tensor P_A.                 (I)

For hypercharge `P=y_L phi phi^T` is symmetric and the old flavor-only
projection is valid. For a general SU(2) generator, `P_A` need not vanish.
Thus the symmetric flavor part cannot be deleted from (H).

Take `phi=(u,0)` with real `u`. At that scalar configuration the complete
sterile-source tensor is

    -u^2/12 (Kfl)_alpha beta [
       bar nu_Lalpha sigma^{mu nu}(nu_Lbeta)^c
                         (g'y_L B_mu nu+g W^3_mu nu/2)
       +g/sqrt(2) bar e_Lalpha sigma^{mu nu}(nu_Lbeta)^c
                                               W^-_mu nu],   (J)

where `W^-_mu nu=(W^1_mu nu+i W^2_mu nu)/sqrt(2)` denotes the full
nonlinear field-strength combination. The neutral term selects `(Kfl)_A`.
The charged term retains the full `Kfl`. For one active flavor a nonzero
`Kfl_11` gives a nonzero charged Grassmann polynomial, even though the
neutral same-flavor tensor vanishes. The checker verifies this with explicit
chiral gamma matrices and `T^1,T^2`; it rejects replacing (H) by flavor-only
antisymmetrization. Equation (J) is an algebraic component projection, not a
measurable decay rate, physical dipole prediction, or complete broken-phase
matching result.

## Independent reproduction and limitations

Run from a Python environment with SymPy 1.14.0:

```sh
python -I -B -O independent_ew_audit.py > INDEPENDENT_EW_AUDIT_RESULTS.json
```

The checker imports no repository module and uses no producer formulas or
results as computational inputs. It builds exact rational complex Pauli and
Dirac matrices and three deterministic generic connection/matter derivative
jets. It verifies covariant Bianchi, the IBP identity against the *ordinary*
derivative of its boundary current, the row square from covariant second
derivatives, and the nonlinear Yang-Mills Euler term by direct variation of
the curvature square. Separate controls cover the scalar current, real
adjoint shifts, curvature Dirac adjoints, current coupling/normalization,
the Majorana factor, and combined weak/flavor Grassmann exchange.

The saved optimized run passes **122 exact checks** and rejects **81
deliberate mutations or asserted-zero negative controls**. These numbers
count implementation controls, not diagrams or independent physical
calculations. The full arbitrary-source scalar algebra of the previous
certificate is reused analytically as stated above; this independent checker
targets its new matrix, variational, and Grassmann content. It does not claim
exhaustive polynomial coefficient coverage, full EFT completeness, a finite
matching prescription, an anomalous dimension, or an observable.
