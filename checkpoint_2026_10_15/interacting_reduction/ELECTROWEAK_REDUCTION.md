# Interacting electroweak reduction and indexed field substitutions

15 October 2026. This is action algebra applied to the selected matched
packet of candidate `d822e9f4d70ad5b455436f40abc3c2b45909bbef`. It adds no
loop integration. It extends the earlier hypercharge reduction to the
matrix-valued action whose two- and three-gauge vertices were explicitly
matched in the non-Abelian checkpoint.

## Reference action and ordering

Use the reference action and Grassmann conventions of
[`ACTION_SIGN.md`](../operator_map/ACTION_SIGN.md), adding the canonical
weak kinetic term `-W^I_mu_nu W^{I mu nu}/4`. The optional scalar quartic
is `-lambda (phi^dagger phi)^2`; set `lambda=0` for the original truncated
potential. Keep the real positive diagonal retained Majorana mass `M`;
constant real symmetric `M` also works with every displayed product in its
given order. The coefficient is `F=C_seed Y_nu*`, and
`kappa=1/(16 pi^2 epsilon)`. Weak, spin and flavor indices are distinct.

```text
phi=tilde H, rho=phi^dagger phi, y_L=y_phi=-1/2, T^I=tau^I/2,
D_mu=partial_mu+i Acal_mu,
Acal_mu=g' y_L B_mu 1_2+g T^I W_mu^I,
Fcal_mu_nu=partial_mu Acal_nu-partial_nu Acal_mu+i[Acal_mu,Acal_nu],
Dcal_rho Fcal_mu_nu=partial_rho Fcal_mu_nu+i[Acal_rho,Fcal_mu_nu],
A=i slash partial, K bar L=-i(D_mu bar L)gamma^mu,
sigma^{mu nu}=i[gamma^mu,gamma^nu]/2, eta=(+---).
```

In components,

```text
Fcal_mu_nu=g' y_L B_mu_nu 1_2+g T^I W_mu_nu^I,
W_mu_nu^I=partial_mu W_nu^I-partial_nu W_mu^I
           -g epsilon^{JKI}W_mu^J W_nu^K,
(D_ad,rho X)^I=partial_rho X^I-g epsilon^{JKI}W_rho^J X^K.
```

The minus nonlinear weak sign follows from **D=partial+i Acal**. Use the
row derivative `D_mu bar L=partial_mu bar L-i bar L Acal_mu`. In particular
`bar L Fcal phi` means `bar L_alpha^a (Fcal)_ab phi^b`; it never means
`bar L phi Fcal` or a trace of `Fcal`.

Let `[bar U V Z]_F=sum_(alpha i)F_(alpha i) bar U_alpha V Z_i`, retaining
the indicated weak contractions and spin ordering. The input counterterm is

```text
L_ct=kappa(P_F+P_F^dagger),
P_F=F_(alpha i)(O_J+O_D)_(alpha i),
O_D=i(D_rho bar L) sigma^{mu nu}gamma^rho Fcal_mu_nu phi N_L/4
    +bar L gamma_mu (Dcal_nu Fcal^{nu mu}) phi N_L/3.
```

`O_J` is exactly the six placements of
[`COVARIANT_ACTION_AND_WARD.md`](../nonabelian_gauge/COVARIANT_ACTION_AND_WARD.md),
with these full derivatives. No mixed-operator `1/2` is inserted.

## Gauge equations derived from the action

Differentiate the gauge-invariant renormalizable action `S_0`, including
every declared renormalizable spectator. Gauge fixing is not included in
these variational equations. Variation of the literal weak curvature gives

```text
delta W_mu_nu^I=(D_ad,mu delta W_nu)^I-(D_ad,nu delta W_mu)^I,
delta S_gauge=integral (D_ad,nu W^{nu mu})^I delta W_mu^I.
```

The fermion kinetic action contributes
`-g bar f gamma^mu T_f^I f delta W_mu^I`. The scalar kinetic action
contributes `-g i phi^dagger T^I <->D^mu phi delta W_mu^I`. Thus

```text
e_B^mu=partial_nu B^{nu mu}-g' J_Y^mu,
e_W^{I mu}=(D_ad,nu W^{nu mu})^I-g J_W^{I mu},
J_Y^mu=i y_phi phi^dagger <->D^mu phi
        +sum_retained_f y_f bar f gamma^mu f,
J_W^{I mu}=i phi^dagger T^I <->D^mu phi
           +sum_retained_f bar f gamma^mu T_f^I f.
```

`phi^dagger T^I <->D phi` denotes
`phi^dagger T^I D phi-(D phi)^dagger T^I phi`. A weak singlet contributes
zero to `J_W`; color and generations are summed. The scalar formula in
terms of `H` is identical for the weak current; the hypercharge formula
uses `y_H=+1/2`. The scalar connection terms are part of these currents.
An additional retained interaction contributes its actual variational
source, rather than an imported free current.

The matrix divergence therefore splits exactly as

```text
Dcal_nu Fcal^{nu mu}
 =g' y_L(e_B^mu+g'J_Y^mu)1_2
  +g T^I(e_W^{I mu}+g J_W^{I mu}).                         (1)
```

This derivation fixes both the weak adjoint sign and the current
normalization `T=tau/2`; neither is inferred from a replacement of `g_y B`.

## Covariant integration by parts before any sourced equation

Define the genuinely ordered matrix monomials

```text
Q_Vh^F=-bar L Fcal_mu_nu(D^nu phi)gamma^mu N_L,
Q_Eh^F=-bar L Fcal_mu_nu(D_rho phi)G^{mu nu rho}N_L/2,
Q_div^F=bar L gamma_mu(Dcal_nu Fcal^{nu mu})phi N_L,
Q_AN^F=bar L Fcal_mu_nu phi sigma^{mu nu}A N_L,
G^{mu nu rho}=gamma^{[mu}gamma^nu gamma^{rho]}.
```

Apply the ordinary derivative to the gauge-singlet product
`bar L sigma^{mu nu}gamma^rho Fcal_mu_nu phi N_L`. Its four terms are
the covariant derivative of the row, of `Fcal`, of `phi`, and the ordinary
derivative of `N`. The noncommuting connection terms cancel in that
order. Clifford algebra and `Dcal_[rho Fcal_mu_nu]=0` then give

```text
O_D=-(Q_Vh^F+Q_Eh^F)/2-Q_div^F/6-Q_AN^F/4
    modulo a total derivative.                              (2)
```

This uses the covariant Bianchi identity, including the adjoint connection.
Using `partial Fcal` in its place loses weak contact terms. Independently,
the row-square identity derived from two row covariant derivatives is

```text
K^2 bar L=-D^2 bar L-bar L Fcal_mu_nu sigma^{mu nu}/2.       (3)
```

Consequently `O_L/3` contributes `+Q_AN^F/6` to the selected packet, and
the combined coefficient is `1/6-1/4=-1/12` with the same ordered weak
matrix. Equations (1)-(3), rather than an Abelian substitution rule,
provide the non-Abelian extension.

## Source-retaining representative and exact residual

Keep all variational matter sources of `S_0`:

```text
e_barL=K bar L-bar S_L,
e_phi=-D^2 phi-(m_H^2+2lambda rho)phi+S_phi,
e_R=A N_L-M N_R-S_R,
e_L=A N_R-M N_L-S_L^N,
R=M N_R+S_R,
T=M^2 N_L+M S_L^N+A S_R,
U=M^3 N_R+M^2 S_R+M A S_L^N+A^2 S_R,
c=2m_H^2/3+lambda rho/3,
j_mu=sum_(alpha i a)F_(alpha i)bar L_alpha^a phi^a gamma_mu N_Li,
j_mu^I=sum_(alpha i a b)F_(alpha i)bar L_alpha^a T^I_ab phi^b gamma_mu N_Li.
```

The complete sourced representative of this selected packet is

```text
P_red=(5/6)[bar L phi U]_F+c[bar L phi R]_F
      +(1/3)[(K bar S_L)phi R]_F-(1/6)[bar L S_phi R]_F
      -(1/6)[bar S_L phi T]_F+(m_H^2/2)[bar S_L phi N_L]_F
      -(1/2)F_(alpha i)(Q_Vh^F+Q_Eh^F)_(alpha i)
      -(g'^2 y_L/6)j_mu J_Y^mu-(g^2/6)j_mu^I J_W^{I mu}
      -(1/12)[bar L Fcal_mu_nu phi sigma^{mu nu}R]_F.        (4)
```

Neither the weak current nor the source `S_R` in the last line is zero in
the interacting theory. Derivatives on sources are retained. The exact
difference, modulo only the total derivative and covariant Bianchi identity
of (2), is

```text
P_F-P_red=(1/3)[(K e_barL)phi A N_L]_F
         +(1/6)[bar L e_phi A N_L]_F
         -(1/6)[e_barL phi A^2 N_L]_F
         +(m_H^2/2)[e_barL phi N_L]_F
         -(g' y_L/6)j_mu e_B^mu-(g/6)j_mu^I e_W^{I mu}
         +(5/6)[bar L phi(M^2e_R+M A e_L+A^2e_R)]_F
         +c[bar L phi e_R]_F+(1/3)[(K bar S_L)phi e_R]_F
         -(1/6)[bar L S_phi e_R]_F
         -(1/6)[bar S_L phi(M e_L+A e_R)]_F
         -(1/12)[bar L Fcal_mu_nu phi sigma^{mu nu}e_R]_F.   (5)
```

This is an off-shell certificate with independent nonzero residuals, not
an assertion that arbitrary off-shell amplitudes are equal after deleting
those residuals. With `g=0`, (4)-(5) reduce to the earlier hypercharge
representative and residual certificate.

For the truncated Yukawa action the explicit sources are

```text
S_L,alpha^a=Y_(alpha j)phi^a N_Rj,
bar S_L,alpha^a=bar N_Rj Y*_(alpha j)(phi^a)*,
S_phi^a=-bar N_Rj Y*_(alpha j)L_alpha^a,
S_Li^N=sum_(alpha a)Y*_(alpha i)(phi^a)*L_alpha^a,
S_Ri=(S_Li^N)^c.
```

Additional Yukawa/spectator sources belong in the total sources of (4)-(5).
For example `H^dagger phi=0` does **not** imply
`H^dagger T^I phi=0`; weak insertions can prevent an Abelian source
simplification. The current terms contain all retained scalar and fermion
contacts, and their scalar connections provide additional gauge vertices.

## Fully indexed field map

All fields on the right are new fields. A first-order local substitution is

```text
L_old,alpha^a=L_alpha^a-kappa ell_alpha^a,
(phi^a)*_old=(phi^a)*-kappa t_a,   phi_old^a=phi^a-kappa t_a*,
B_old,mu=B_mu+kappa g'y_L(j_mu+j_mu^dagger)/6,
W_old,mu^I=W_mu^I+kappa g(j_mu^I+j_mu^{I dag})/6,
N_old,Li=N_Li-kappa(x_i+y_i^c),
N_old,Ri=N_Ri-kappa(y_i+x_i^c).                              (6)
```

The lepton adjoint is the conjugate of its displayed substitution. The
coefficients are

```text
ell_alpha^a=sum_i F_(alpha i)[
  (i/3)gamma^mu(D_mu phi)^a A N_Li
  +(1/6)phi^a A^2 N_Li+(m_H^2/2)phi^a N_Li],
t_a=(1/6)sum_(alpha i)F_(alpha i)bar L_alpha^a A N_Li.       (7)
```

In the scalar conjugate, the entire Grassmann-even bilinear is conjugated
and the fermion order reversed. Thus
`t_a*=(1/6)sum F* overline(A N_Li)L_alpha^a`.
For `phi=epsilon H*`, the corresponding map is
`H_old=H+kappa epsilon t^T`, with `epsilon^12=+1`.

The neutral chiral columns, with every weak and flavor slot fixed, are

```text
x_i=(5/6)sum_(alpha a)(F M^2)*_(alpha i)(phi^a)*L_alpha^a
    +c sum_(alpha a)F*_(alpha i)(phi^a)*L_alpha^a
    -(5/6)sum_(alpha a)F*_(alpha i)[
       (phi^a)*(D^2 L_alpha)^a
       +2(D_mu phi)^{a*}(D^mu L_alpha)^a
       +(D^2 phi)^{a*}L_alpha^a]
    +(i/6)sum_(alpha a)F*_(alpha i)(phi^a)*gamma^mu(D_mu S_L,alpha)^a
    -(i/6)sum_(alpha a)F*_(alpha i)(D_mu phi)^{a*}gamma^mu S_L,alpha^a
    -(1/6)sum_(alpha a)F*_(alpha i)(S_phi^a)*L_alpha^a
    -(1/12)sum_(alpha a b)F*_(alpha i)(phi^b)*
          (Fcal_mu_nu)_ba sigma^{mu nu}L_alpha^a,

y_i=(5i/6)sum_(alpha a)(F M)*_(alpha i)gamma^mu[
       (phi^a)*(D_mu L_alpha)^a+(D_mu phi)^{a*}L_alpha^a]
    -(1/6)sum_(alpha a)(F M)*_(alpha i)(phi^a)*S_L,alpha^a.   (8)
```

The repeated Lorentz indices in the final term of `x` are summed.
`x=P_L x`, `y=P_R y`. In particular the curvature index order is `ba`
after taking the adjoint of `bar L^a Fcal_ab phi^b`; it is not an ordinary
transpose without conjugation. `Fcal` is Hermitian in weak space for real
gauge fields. Mass products remain `(F M)*` and `(F M^2)*`.

For a compact derivation, let
`W_i=sum F_(alpha i)bar L_alpha phi`,
`V_i=sum F_(alpha i)bar S_L,alpha phi`, and
`K_0 W=-i(partial_mu W)gamma^mu`. Integration by parts in (5) yields

```text
bar x_i=(5/6)[(W M^2)_i+K_0^2W_i]+cW_i
        +(1/3)sum F_(alpha i)(K bar S_L,alpha)phi
        -(1/6)sum F_(alpha i)bar L_alpha S_phi
        -(1/6)K_0V_i
        -(1/12)sum F_(alpha i)bar L_alpha Fcal_mu_nu phi sigma^{mu nu},
bar y_i=(5/6)(K_0W)_j M_ji-(1/6)V_j M_ji,
P_F-P_red=e_barL ell+t e_phi
         -(g'y_L/6)j_mu e_B^mu-(g/6)j_mu^I e_W^{I mu}
         +bar x_i e_Ri+bar y_i e_Li.                        (9)
```

The covariant product rules for `K_0W` and `K_0^2W` hold because the row
and column connections cancel with their correct matrix order. Expanding
them and then taking the Dirac adjoint gives exactly (7)-(8).

The real gauge shifts in (6) cancel (9) and its conjugate. Under a weak
transformation `U`, the matrix
`K_mu^a_b=sum F_(alpha i)phi^a bar L_alpha^b gamma_mu N_Li`
transforms as `U K_mu U^dagger`; hence
`j_mu^I=tr(T^I K_mu)` transforms in the adjoint representation. Its real
part has the same transformation. The difference between old and new weak
connections is therefore adjoint covariant, as required for a connection
substitution. A real shift of `B` is invariant.

With `w=x+y`, the Majorana substitution is `N_old=N-kappa(w+w^c)`.
For an odd Majorana residual `E_N=e_R+e_L`, the Grassmann identity
`bar E_N w=overline(w^c)E_N` gives

```text
bar x e_R+bar y e_L+h.c.=overline(w+w^c)(e_R+e_L).
```

There is no additional `1/2`: variation of the two Majorana slots cancels
the `1/2` in the quadratic action. Thus, to first order,

```text
S_0[old]+kappa integral(P_F[old]+h.c.)
 =S_0[new]+kappa integral(P_red[new]+h.c.)+O((kappa F)^2).  (10)
```

The local inverse follows by reversing the first-order signs. Varying
another EFT operator of coefficient `C_d` adds `kappa F C_d`, which has
its own inverse-scale order; it is not automatically quadratic in kappa.

## Weak source tensor and its essential flavor limitation

Substituting `S_R=Y_nu^T X^c`, `X_alpha=phi^dagger L_alpha`, into the last
line of (4) gives

```text
-(FY_nu^T)_(alpha beta)
  bar L_alpha^a(Fcal_mu_nu)_ab phi^b sigma^{mu nu}X_beta^c/12. (11)
```

For the hypercharge identity matrix, this is proportional to
`bar X_alpha sigma X_beta^c`, so flavor antisymmetry alone applies. For
a weak insertion, the two slots contain **different weak composites**.
The tensor to antisymmetrize is instead

```text
K_(alpha a;beta c)=(FY_nu^T)_(alpha beta)(Fcal_mu_nu phi)^a phi^c,
K_relevant=(K_(alpha a;beta c)-K_(beta c;alpha a))/2.        (12)
```

The Lorentz indices in this statement stay fixed. The symmetric spin
contraction matrix associated with `sigma C` makes the combined weak/flavor
coefficient antisymmetric after exchanging the Grassmann slots. It does
not set (11) to zero for one active flavor. For example, with `phi=(u,0)`
and a `T^1` insertion, `(T^1 phi)^a phi^c` has a nonzero `21` entry and
a zero `12` entry. Its combined-slot antisymmetric part survives even
though every flavor-only antisymmetric matrix is zero in one flavor.
Using the earlier hypercharge flavor kernel for this full weak source
would therefore discard a real descendant of the selected action.

## Reproduction and scope

Run `python3 -I -S -B -O verify_electroweak_reduction.py`. It uses only the
standard library and exact Gaussian rational arithmetic. It differentiates
the canonical gauge, scalar and fermion actions, evaluates generic noncommuting
connection jets and covariant Bianchi/IBP identities, checks the entire
source residual with two active and three sterile flavors, and compares
the integrated rows with the indexed adjoint columns. It also tests real
adjoint weak shifts, Majorana completion and combined weak/flavor source
antisymmetry. Deliberate wrong signs, source omissions, matrix order,
generator normalization and Majorana half factors must fail. The adjacent
JSON records the actual categories, counts and negative controls.
The completed optimized run passes **1,410 positive exact checks** and
rejects **20 deliberate mutations**. These are implementation checks of
the displayed action algebra, not separate physical calculations.

This gives the complete reduction and field map **of the selected matched
packet**, not a complete retained-sterile Green basis or a full interacting
canonical normalization. Renormalizable sources suffice for this isolated
`F=O(Lambda^-3)` packet through that order under fixed retained-mass
counting. Separate dimension-six field changes acting on dimension-five
interactions can contribute at `Lambda^-3` and remain required. A general
complex symmetric mass requires the correctly alternating chiral products
`M M^dagger` and `M M^dagger M` rather than the real-basis powers used here.
Gauge-fixed quantum generating functionals also require transforming
their gauge fixing and ghosts. Four-dimensional simple-pole Clifford
identities do not fix finite evanescent prescriptions. No physical dipole,
finite amplitude, beta function, observable or novelty claim follows from
this action-coordinate calculation alone.
