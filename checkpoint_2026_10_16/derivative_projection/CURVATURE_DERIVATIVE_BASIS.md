# Chiral curvature derivatives and complete sourced placement transport

16 October 2026. Input: candidate
`4ed8560e019dfe74ed6c5178b0f990fe3ccb3008`. This note resolves the algebraic
Lorentz/weak basis and the selected dual channel of the remaining
`-(Q_Vh^F+Q_Eh^F)/2` term in the actual source-retaining representative.
It then gives an exact boundary/EOM certificate for moving its derivative
to the lepton row. The alternative placement has correlated tensor and
current sources; it does not establish that a physical operator is redundant.

## Conventions and the algebraic derivative-on-scalar basis

Keep the conventions of the
[`matrix source reduction`](../../checkpoint_2026_10_15/interacting_reduction/ELECTROWEAK_REDUCTION.md):

```text
eta=(+---), epsilon^{0123}=+1, gamma5=i gamma0 gamma1 gamma2 gamma3,
sigma^{mu nu}=i[gamma^mu,gamma^nu]/2,
G^{mu nu rho}=gamma^{[mu}gamma^nu gamma^{rho]},
phi=tilde H=epsilon H*, epsilon12=+1,
D=partial+i Acal, D_mu bar L=partial_mu bar L-i bar L Acal_mu,
Fcal_mu_nu=g'y_L B_mu_nu 1_2+g T^I W_mu_nu^I,
y_L=-1/2, T^I=tau^I/2,
W_mu_nu^I=partial_mu W_nu^I-partial_nu W_mu^I
           -g epsilon^{JKI}W_mu^J W_nu^K,
Dcal_rho Fcal_mu_nu=partial_rho Fcal_mu_nu+i[Acal_rho,Fcal_mu_nu].
```

Define `tilde V^{mu nu}=epsilon^{mu nu rho sigma}V_rho_sigma/2` and lower
its indices with `eta`; on two-forms `tilde tilde V=-V`. Write
`[bar U V Z]_F=sum_(alpha i)F_(alpha i)bar U_alpha V Z_i`, with
`F=C_seed Y_nu*`. Flavor, weak and spin matrix slots stay separate.

For a fixed active/sterile flavor pair, a complete **algebraic** basis
with the derivative on `phi` is

```text
H_B=bar L_alpha^a gamma^mu N_Li(D^nu phi)^a B_mu_nu,
H_Btilde=bar L_alpha^a gamma^mu N_Li(D^nu phi)^a tilde B_mu_nu,
H_W=bar L_alpha^a gamma^mu N_Li T^I_ab(D^nu phi)^b W_mu_nu^I,
H_Wtilde=bar L_alpha^a gamma^mu N_Li T^I_ab(D^nu phi)^b tilde W_mu_nu^I.
                                                               (1)
```

These definitions include no gauge coupling or flavor coefficient. They
have canonical dimension seven. There is no exchange symmetry between
the unlike `bar L` and retained `N_L` slots, and no mixed-operator `1/2`.
For `n` active and `s` sterile flavors the algebraic space in (1) has
`4 n s` complex coordinates per displayed orientation, before IBP or EOM.
The Hermitian conjugate supplies the conjugate orientation.

Completeness follows by separating Lorentz and weak representations. The
two vector slots are the chiral bilinear and `D phi`; their antisymmetric
product can contract with a strength or its dual. Those are the two
Lorentz scalars. Chirality does not identify them: the four matrices
`P_R gamma^mu P_L` span a four-dimensional vector space. The weak
intertwiner is unique for `B` and for `W`, since `2 tensor 3=2+4` contains
one output doublet. The exact checker also constructs the entire
coefficient systems without starting from these tensors: a general
`V^a w^b F^{rs}` has 96 coefficients, and all 576 infinitesimal Lorentz
equations have rank **94**, hence nullity **two**. The ordinary and dual
tensors span that nullspace. Weak covariance has rank three on four
coefficients for `B`, and rank 11 on 12 coefficients for `W`; each
nullity is one.

These statements concern the explicit derivative-on-`phi` placement.
They do **not** certify independence after all derivatives, source
operators and field-coordinate directions of a full EFT are admitted.

## Exact chiral triple-gamma sign and selected coefficient image

The literal four-dimensional Clifford identity is

```text
G^{mu nu rho}=i epsilon^{mu nu rho sigma}gamma_sigma gamma5.
```

Since `gamma5 N_L=-N_L`, the actual input operators obey

```text
Q_Vh^F=-H(Fcal),
Q_Eh^F=-i H(tilde Fcal),
H(V)=bar L V_mu_nu(D^nu phi)gamma^mu N_L,
H_plus=-(Q_Vh^F+Q_Eh^F)/2
      =[H(Fcal)+i H(tilde Fcal)]/2=H(Fcal_plus),
Fcal_plus=(Fcal+i tilde Fcal)/2,
tilde Fcal_plus=-i Fcal_plus.                             (2)
```

Thus one linear combination is selected inside the two-dimensional
Lorentz space; the space itself still has two directions. In the basis
(1), the ordered coefficient matrices in this selected subpacket are

| Algebraic coordinate | Coefficient |
| --- | --- |
| `H_B` | `g'y_L F/2` |
| `H_Btilde` | `i g'y_L F/2` |
| `H_W` | `g F/2` |
| `H_Wtilde` | `i g F/2` |

The common action factor is `kappa=1/(16 pi^2 epsilon)`, with the
Hermitian conjugate included once. In normalized plus/minus bases
`H_{V,plus/minus}=[H(V) +/- i H(tilde V)]/2`, the coefficients are
`g'y_L F` and `gF` in the plus channels and zero in the minus channels.
The Hermitian-conjugate sector has the conjugate chirality relation.
If `B` and `W` coefficients were separately free, projection onto the
plus channels would have rank two per flavor component. The actual
packet further correlates those two coefficients through the same `F`.

For a freely specified raw `n x s` matrix `F` and at least one nonzero
gauge coupling, this algebraic coefficient map has rank `n s` rather
than `4 n s`: dual coefficients are `i` times ordinary coefficients,
and `g C_B=g'y_L C_W`. This rank describes transport of raw `F`.
The particular symmetric-seed relation `F=C_seed Y_nu*` can restrict
its input image further; no arbitrary-seed reachability claim follows.

A concrete Fourier control fixes the circular-mode sign without importing
a helicity label. For the outgoing mode `exp(+i k.x)`, take
`k^mu=(1,0,0,1)` and `e^mu=(0,1,+i,0)` (normalization immaterial).
Its literal strength is `V_mu_nu=i(k_mu e_nu-k_nu e_mu)` and has
`tilde V=-i V`, so (2) retains it. The mode with `e_y=-i` has the
opposite dual eigenvalue and is killed by this plus projection. A
scattering external-state label additionally specifies whether its
polarization is conjugated. This is a local curvature-channel statement,
not a finite helicity amplitude or observable selection rule.

## Full placement relation before any equations of motion

Introduce the derivative-on-row monomials and the exact source-shaped
objects

```text
Q_Vr=-(D^nu bar L)Fcal_mu_nu phi gamma^mu N_L,
Q_Er=-(D_rho bar L)Fcal_mu_nu phi G^{mu nu rho}N_L/2,
Q_div=bar L gamma_mu(Dcal_nu Fcal^{nu mu})phi N_L,
Q_L=(K bar L)Fcal_mu_nu phi sigma^{mu nu}N_L,
Q_N=bar L Fcal_mu_nu phi sigma^{mu nu}A_N N_L,
K bar L=-i(D_mu bar L)gamma^mu,  A_N=i slash partial,
T^rho=bar L sigma^{mu nu}gamma^rho Fcal_mu_nu phi N_L.
```

Every expression in this section carries the common flavor contraction
`F_(alpha i)` when used in the packet. Direct Clifford algebra, before
any EOM or boundary manipulation, gives

```text
Q_Er=Q_Vr-Q_L/2,
O_D=(Q_Vr+Q_Er)/2+Q_div/3
   =Q_Vr-Q_L/4+Q_div/3.
```

The companion derivative-on-sterile placements make the IBP content
explicit. Define

```text
Q_Vp=bar L Fcal_mu_nu phi gamma^mu partial^nu N_L,
Q_Ep=bar L Fcal_mu_nu phi G^{mu nu rho}partial_rho N_L/2,
B_V=partial^nu[bar L Fcal_mu_nu phi gamma^mu N_L],
B_E=partial_rho[bar L Fcal_mu_nu phi G^{mu nu rho}N_L]/2.
```

Before any sourced equation is used, the two full product derivatives and
the sterile Clifford identity give

```text
Q_Vp=Q_Vr+Q_Vh^F+Q_div+B_V,
Q_Ep=Q_Er+Q_Eh^F+B_E,
Q_Vp+Q_Ep=-Q_N/2.
```

The second product rule uses the covariant Bianchi identity. These
relations retain both matter-derivative placements and the strength
divergence; none authorizes a free fermion or gauge equation. Their
boundaries satisfy `(B_V+B_E)/2=-i partial_rho T^rho/4`.

Covariant differentiation of all four factors in `T^rho`, using
`Dcal_[rho Fcal_mu_nu]=0`, gives the exact placement identity

```text
H_plus=Q_Vr-Q_L/4+Q_div/2+Q_N/4-i partial_rho T^rho/4.    (3)
```

Its adjoint-connection and row/curvature ordering are essential in the
non-Abelian sector. An ordinary `partial Fcal` Bianchi identity or a
reversed matrix insertion changes the boundary/weak contact packet.
Equation (3) demonstrates why the basis (1) cannot be called a complete
physical basis merely from its algebraic rank: derivative placements and
source-shaped operators are already connected by IBP and Dirac identities.

## Source-retaining alternative and its exact residual

Use the total variational sources of the same gauge-invariant
renormalizable reference action:

```text
e_barL=K bar L-bar S_L,
e_R=A_N N_L-M N_R-S_R,
e_B^mu=partial_nu B^{nu mu}-g' J_Y^mu,
e_W^{I mu}=(D_ad,nu W^{nu mu})^I-g J_I^mu,
j^mu=[bar L phi gamma^mu N_L]_F,
j_I^mu=[bar L T^I phi gamma^mu N_L]_F,
G_current=g'^2 y_L j_mu J_Y^mu+g^2 j_I,mu J_I^mu,
G_residual=g'y_L j_mu e_B^mu+g j_I,mu e_W^{I mu}.
```

Then `[Q_div]_F=G_current+G_residual`. In particular scalar currents and
all charged fermion currents survive this split. Let

```text
S_Ltensor=[bar S_L Fcal_mu_nu phi sigma^{mu nu}N_L]_F,
S_Ntensor=[bar L Fcal_mu_nu phi sigma^{mu nu}(M N_R+S_R)]_F,
P_alt=[Q_Vr]_F-S_Ltensor/4+G_current/2+S_Ntensor/4.
```

Equation (3) becomes the full off-shell certificate

```text
[H_plus]_F-P_alt
 =-[e_barL Fcal_mu_nu phi sigma^{mu nu}N_L]_F/4
  +G_residual/2
  +[bar L Fcal_mu_nu phi sigma^{mu nu}e_R]_F/4
  -i partial_rho[T^rho]_F/4.                             (4)
```

No source is set to zero. The retained mass is real positive diagonal as
in the input; no inverse mass is used. The remainder of the actual
`P_red` already contains `-G_current/6-S_Ntensor/12`. Thus its curvature
and gauge-current subpacket can equivalently be written, modulo (4), as

```text
[Q_Vr]_F-S_Ltensor/4+G_current/3+S_Ntensor/6.             (5)
```

The explicit mass tensor coefficient has moved from `-1/12` to `+1/6`,
and the current coefficient from `-1/6` to `+1/3`, accompanied by the
new row derivative and lepton tensor source. They are equivalent
representatives, not competing physical predictions. This is a transfer
of a selected derivative channel; it does **not** prove that an independent
operator direction of the full EFT has disappeared.

## Explicit induced neutrino, charged-lepton and scalar-current sources

For the retained neutrino and charged-lepton Yukawa interactions,

```text
bar S_L,alpha^a=bar N_Rj Y_nu*_(alpha j)(phi^a)*
               +bar e_Rr Y_e*_(alpha r)(H^a)*+other sources,
S_Ri=(Y_nu)_(beta i)X_beta^c+other sources,
X_beta^c=sum_a phi^a(L_beta^a)^c.
```

The new lepton tensor source in either (4) or (5) contains

```text
-(Y_nu^dagger F)_(ji)(phi^dagger Fcal_mu_nu phi)
        bar N_Rj sigma^{mu nu}N_Li/4
-(Y_e^dagger F)_(ri)(H^dagger Fcal_mu_nu phi)
        bar e_Rr sigma^{mu nu}N_Li/4.                    (6)
```

The sterile tensor bilinear is antisymmetric in `j,i`. Its first
coefficient therefore selects `(Y_nu^dagger F)_A`. For the selected
`F=C_seed Y_nu*`, `Y_nu^dagger C_seed Y_nu*` is symmetric, so that
particular neutrino-source tensor cancels by the Grassmann identity.
It need not cancel for an arbitrary `F` or another source.

In the charged-lepton term `H^dagger phi=0` removes its hypercharge
piece, whereas `H^dagger T^I phi` is generally nonzero. The surviving
charged weak tensor must be retained. The fact that an Abelian
charged-lepton source vanishes does not justify deleting it before the
matrix contraction. Other spectator contributions remain their actual
variational source tensors.

The sterile-source tensor in (4) supplies `+F M/4` times the retained-mass
tensor and `+F Y_nu^T/4` times the full active weak tensor. In the entire
subpacket (5), these ordered coefficients are instead `+F M/6` and
`+F Y_nu^T/6`. The prior
[`active tensor basis`](../../checkpoint_2026_10_15/operator_basis/ACTIVE_TENSOR_BASIS.md)
specifies their independent-flavor-pair normalization. In a weak active
source, antisymmetry acts on the combined weak/flavor slots; it cannot be
replaced by flavor antisymmetry alone.

No scalar equation or explicit scalar source `S_phi` is used in (3)-(6).
Scalar sources still occur through the gauge current. The current part
of (5) is `-2` times the previous selected `-G_current/6`. If its scalar
subpacket is subsequently reduced using the specified
[`scalar-current map`](../action_completion/SCALAR_CURRENT_REDUCTION.md),
that subpacket's two derivative remainders and its correlated
`rho^2(Y_nu^dagger F)bar N_R N_L`, `rho(F M)bar X N_R`, and
`rho(F Y_nu^T)bar X X^c` source terms all scale by the same `-2`.
This conditional transport does not include contributions from other
parts of the selected action, and it must not be counted again as an
independent loop contribution.

## An explicit optional placement field map

The residual in (4), together with its conjugate, is removed against the
same gauge-invariant `S_0` by the local first-order shifts below. All
right-hand fields are new fields; `kappa=1/(16 pi^2 epsilon)`.

```text
ell_alpha^a=-(1/4)sum_(i b)F_(alpha i)(Fcal_mu_nu)_ab phi^b
                                    sigma^{mu nu}N_Li,
x_i=+(1/4)sum_(alpha a b)F*_(alpha i)(phi^b)*(Fcal_mu_nu)_ba
                                    sigma^{mu nu}L_alpha^a,
L_old=L-kappa ell,
N_old=N-kappa(x+x^c),
B_old,mu=B_mu-kappa g'y_L(j_mu+j_mu^dagger)/2,
W_old,mu^I=W_mu^I-kappa g(j_I,mu+j_I,mu^dagger)/2.        (7)
```

The row lepton shift is the actual conjugate; the scalar is unshifted.
The adjoint weak shift is real, as in the original source map. The
curvature indices in `x` follow from its **Dirac adjoint**, using Hermitian
matrix curvature for real gauge fields. Its Majorana-real completion is
`x+x^c` without an extra half. Variation of `S_0` cancels the three
residuals in (4), yielding the alternative representative at first order.
The boundary term remains a boundary term. This is an optional placement
coordinate map, not an unannounced composition with the prior
normalization/source transformations. A gauge-fixed functional also
requires transforming its gauge fixing and ghost terms consistently.
Other retained EFT operators have their separate `kappa F C_d` variations
at the inverse-scale order appropriate to those coefficients.

## Reproduction and remaining gap

Run `python3 -I -S -B -O verify_curvature_derivative.py`. The standalone
standard-library checker covers all 64 triple-gamma components, complete
Lorentz/weak covariance coefficient systems, explicit circular Fourier
modes, raw-`F` coefficient-image constraints, eight generic noncommuting
connection/source IBP trials, full off-shell certificates, exact
Grassmann source expansions, and the indexed Majorana map. The adjacent
JSON records the completed check counts and deliberately rejected errors.
The completed optimized run passes **320 positive exact checks** and
rejects **12 deliberate mutations**.
These controls are action/tensor implementation checks, not independent
loop integrations or evidence of novelty.

The derivation establishes the complete algebraic basis in (1) and the
selected coefficient/dual image in it. It explicitly assesses the
IBP/fermion/gauge-EOM transport rather than assuming that algebraic
independence is physical independence. A complete relation analysis
including the rest of the derivative Green basis, all retained source
classes and lower-order field-coordinate contributions is still needed
before claiming a physical basis quotient or RG coefficient. No additional
redundant physical direction is claimed removed here. The four-dimensional
dual/Clifford identity does not prescribe finite evanescent terms; a
subtraction convention, complete loops/counterterms and finite matching
remain open. Nothing here determines a physical helicity amplitude,
dipole prediction, running coefficient or observable.
