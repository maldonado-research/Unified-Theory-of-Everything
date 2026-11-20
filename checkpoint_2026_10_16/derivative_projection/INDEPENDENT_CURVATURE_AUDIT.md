# Independent audit of the curvature–scalar derivative sector

16 October 2026. Independent derivation from the actual selected packet in
[`ELECTROWEAK_REDUCTION.md`](../../checkpoint_2026_10_15/interacting_reduction/ELECTROWEAK_REDUCTION.md),
with input candidate `4ed8560e`. The core signs, invariant kernels and
covariant-product equations were derived and implemented without reading
the producer executable. Formula comparisons were exchanged with its author
and the final prose was cross-checked. The exact executable imports no project
code or saved coefficients. This is AI-assisted internal verification, not external
peer review, another loop integration, or a full EFT basis calculation.

## Triple gamma, dual strength and the selected projection

Fix all signs explicitly:

```
eta=(+---), epsilon^0123=+1, gamma5=i gamma0 gamma1 gamma2 gamma3,
sigma^(mu nu)=i[gamma^mu,gamma^nu]/2,
G^(mu nu rho)=gamma^[mu gamma^nu gamma^rho],
(*F)_(mu nu)=eta_mu eta_nu epsilon^(mu nu rho sigma)F_(rho sigma)/2,
**F=-F.
```

The last display defines lower components of the dual; the Levi-Civita
symbol on its right has upper indices. Matrix multiplication gives

```
G^(mu nu rho)=i epsilon^(mu nu rho sigma)gamma_sigma gamma5.
```

The explicit chiral projector obeys `gamma5 N_L=-N_L`. With all weak
matrix order retained, define

```
H(F)=bar L gamma_mu F^(mu nu)(D_nu phi)N_L,
Q_Vh=−H(F),
Q_Eh=−bar L F_mu_nu(D_rho phi)G^(mu nu rho)N_L/2.
```

It follows, including the minus from interchanging the last two indices of
the dual, that

```
Q_Eh=−i H(*F),
H_plus=−(Q_Vh+Q_Eh)/2=[H(F)+i H(*F)]/2.                (D1)
```

Thus the selected combination couples to `(F+i*F)/2`, whose dual
eigenvalue is `−i`. This eigenvalue statement fixes the convention without
depending on competing names for self-duality. Its Hermitian conjugate must
also be retained. The executable verifies the full Clifford identity,
chirality, dual-square sign and the selected dual eigenvalue.

Here the physical weak matrix is
`Fcal=g'y_L B 1_2+g T^I W^I`, with `T^I=tau^I/2` and `y_L=-1/2`.
The flavor coefficient remains `F_flavor=C_seed Y_nu*`; it must not be
confused with a field-strength symbol in (D1). Equations below use `Fcal`
for that weak matrix, and `[ ... ]_F` for the original ordered flavor sum.

## Algebraic invariant rank in the fixed scalar-derivative placement

The four complex component products `bar L_(dotted s) N_(undotted t)`
are an invertible basis for the vector bilinear `bar L gamma^mu N_L`.
The independent Weyl component matrix has nonzero determinant. It is
therefore legitimate to classify a vector `b^mu`, a scalar-derivative
covector `u_nu`, and an antisymmetric covariant strength `F_mu_nu` first,
then return to actual Grassmann spinors.

The most general Lorentz tensor has `4×4×6=96` coefficients. The six
rotation/boost generators give 576 coefficient constraints. Exact rational
elimination gives rank 94 and nullity two. The ordinary and dual
contractions in (D1) are separately in this complete kernel.

The weak singlet-strength channel has four possible `bar L_a Dphi_b`
coefficients; its covariance constraints have rank three. The adjoint
strength channel has twelve coefficients; its constraints have rank eleven.
The unique weak contractions are respectively `bar L Dphi B` and
`bar L T^I Dphi W^I`. The checker separately verifies those invariant
vectors, using a traceless matrix basis for the adjoint during elimination
and the actual half-Pauli matrices for physical component witnesses.

Hence the fixed-placement algebraic rank is **four**:

```
H_B, H_*B, H_W, H_*W.                                (D2)
```

An actual chiral component minor supplements the invariant-kernel proof.
Choose weak index zero for `bar L` and `Dphi`, dotted/undotted spin zero,
and derivative component one. Use strength components `01` and `23` in
the B channel and again in the `W^3` channel. The four-by-four minor is
diagonal with entries `−1,1,−1/2,1/2`, with nonzero determinant `1/4`.
Thus the conclusion is also a component/exterior-polynomial statement,
not just a formal vector-tensor identity.

Three distinct ranks must be kept separate:

1. The unrestricted ordinary/dual B/W placement has rank four per
   unconstrained flavor pair.
2. Applying the chiral projector to independently variable B/W channels
   has rank two, one combination per channel.
3. The **actual selected packet at fixed couplings** has coefficient column
   `(g'y_L/2, i g'y_L/2, g/2, i g/2) F_flavor` in (D2). Its formal raw-F
   coefficient image has rank one per flavor pair when that column is
   nonzero. The seed/Yukawa constraints on `F_flavor` can further restrict
   its flavor image.

These are complex ranks for the displayed non-Hermitian channel, with its
conjugate understood. None of them proves independence after all EFT
integration-by-parts or interacting EOM relations. The next calculation
explicitly moves the scalar derivative into other placements.

## An independent route through the row derivative

Keep `K bar L=-i(D_mu bar L)gamma^mu` and `A_N=i slash partial`.
Define

```
Q_Vr=−(D^nu bar L)Fcal_mu_nu phi gamma^mu N_L,
Q_div=bar L gamma_mu(Dcal_nu Fcal^(nu mu))phi N_L,
Q_AN=bar L Fcal_mu_nu phi sigma^(mu nu)A_N N_L,
S^rho=bar L sigma^(mu nu)gamma^rho Fcal_mu_nu phi N_L.
```

The commutator of `sigma` and `gamma` fixes a pointwise row identity,
without integration by parts:

```
O_D=Q_Vr−(K bar L)Fcal_mu_nu phi sigma^(mu nu)N_L/4
                                      +Q_div/3.       (D3)
```

Independently differentiate the complete neutral product `S^rho`. Its
four terms differentiate the row, curvature, scalar and sterile field.
The adjoint Bianchi identity removes the antisymmetric curvature derivative.
Retaining the total derivative, the earlier scalar-placement identity is

```
O_D=H_plus−Q_div/6−Q_AN/4+i partial_rho S^rho/4.
```

Combining this with (D3) gives the exact alternative placement

```
H_plus=Q_Vr−(K bar L)Fcal_mu_nu phi sigma^(mu nu)N_L/4
             +Q_div/2+Q_AN/4−i partial_rho S^rho/4.    (D4)
```

The separate derivative-on-sterile placements also check directly. Define
`Q_Vp=bar L Fcal_mu_nu phi gamma^mu partial^nu N_L` and
`Q_Ep=bar L Fcal_mu_nu phi G^(mu nu rho)partial_rho N_L/2`, with
`Q_Er` the corresponding negative derivative-on-row triple-gamma term.
Let `B_V` and `B_E` be the ordinary derivatives of their complete vector
and triple-gamma products. Independently expanding every slot gives

```
Q_Er=Q_Vr−Q_L/2,              Q_L=(K bar L)Fcal phi sigma N_L,
Q_Vp+Q_Ep=−Q_AN/2,
Q_Vp=Q_Vr+Q_Vh+Q_div+B_V,
Q_Ep=Q_Er+Q_Eh+B_E,
(B_V+B_E)/2=−i partial_rho S^rho/4.
```

The independent checker evaluates all five identities on the same full
connection/fermion jets. In particular the absence of a strength term in
the triple-gamma product rule uses the full adjoint Bianchi identity;
the vector product rule retains its strength divergence.

In the executable, noncommuting Hermitian connection values, first
derivatives and symmetric second derivatives construct the curvature and
its full adjoint derivative. The Bianchi identity is checked on those
jets. Independent exterior components represent the row, sterile field
and every fermion derivative; scalar covariant jets are arbitrary. Both
routes (D3)-(D4) then agree coefficient by coefficient at two different
connection/scalar choices. No field is put on shell.

## Interacting sources and the exact residual

Use the actual variational equations of the retained reference action:

```
e_barL=K bar L−bar S_L,
e_R=A_N N_L−M N_R−S_R,
Dcal_nu Fcal^(nu mu)
 =g'y_L(e_B^mu+g'J_Y^mu)1_2+gT^I(e_W^(I mu)+gJ_I^mu).
```

With the ordered flavor sum restored, the source-retaining alternative is

```
H_alt=[Q_Vr]_F−[bar S_L Fcal_mu_nu phi sigma^(mu nu)N_L]_F/4
      +(g'^2 y_L j_mu J_Y^mu+g^2 j_I,mu J_I^mu)/2
      +[bar L Fcal_mu_nu phi sigma^(mu nu)(M N_R+S_R)]_F/4.
```

Its off-shell difference is

```
[H_plus]_F−H_alt=−i partial_rho[S^rho]_F/4
 −[e_barL Fcal_mu_nu phi sigma^(mu nu)N_L]_F/4
 +(g'y_L j_mu e_B^mu+g j_I,mu e_W^(I mu))/2
 +[bar L Fcal_mu_nu phi sigma^(mu nu)e_R]_F/4.          (D5)
```

The signs of all three residual types are fixed. The independent checker
tests matter-source replacement in independent nonzero source/residual
slots. It separately tests the gauge replacement with actual scalar and
quark-current polynomials and nonzero gauge residuals. This source-slot
test is separate from the generic connection-jet identity; it does not
claim arbitrary specified jets already solve the interacting gauge EOM.
Total currents and all additional retained sources belong in (D5).

The physical neutrino source is `S_R=Y_nu^T X^c`, with
`X=phi^dagger L`. The lepton row source includes both
`bar N_R Y_nu^dagger phi^dagger` and
`bar e_R Y_e^dagger H^dagger`. The latter contracts as
`bar e_R Y_e^dagger(H^dagger Fcal_mu_nu phi)sigma^(mu nu)N_L`.
It cannot be dropped using `H^dagger phi=0`: a weak-strength insertion
generically gives `H^dagger T^I phi!=0`. The pure B part does vanish by
that orthogonality. Exact one-flavor exterior examples retain the charged
weak source and the sterile-Yukawa source while checking the expected
one-flavor sterile tensor zero. Charge-conjugate components use
`C=i gamma^2 gamma^0=diag(epsilon,-epsilon)` and the whole Grassmann order.

The residual also identifies an optional local placement map. Its rows
are `ell=−F_flavor(Fcal phi)sigma N_L/4` and
`bar x=+F_flavor bar L Fcal phi sigma/4`; thus
`x=+F_flavor* phi^dagger Fcal sigma L/4`, with the actual weak adjoint
order. `L_old=L−kappa ell`, `N_old=N−kappa(x+x^c)` and real gauge shifts
`B_old=B−kappa g'y_L(j+j^dagger)/2`,
`W_old^I=W^I−kappa g(j_I+j_I^dagger)/2`
cancel the three residuals. No scalar shift or extra Majorana half is
required. This map is an explicitly optional coordinate choice; the
independent executable certifies (D5) and its source signs rather than
claiming a new fully reduced basis.

In particular, replacing this sector inside the earlier representative
transfers `+Q_AN/4` and `+Q_div/2` into its explicit tensor/current entries.
Together with their old `−1/12` and `−1/6` entries, those entries become
`+1/6` and `+1/3` **in this alternate placement**, with `Q_Vr`, row
sources, field transport and the boundary certificate retained. This is
an explicit demonstration that the isolated representative coefficients
are coordinate dependent. It is not an additional physical correction.

## Reproduction and limits

Run

```
python -I -S -B -O independent_curvature_audit.py
```

The standard-library implementation uses exact rational real/imaginary
arithmetic and explicit failure exceptions. It imports no producer code.
[`INDEPENDENT_CURVATURE_RESULTS.json`](INDEPENDENT_CURVATURE_RESULTS.json)
records **1,535 positive comparisons and 24 rejected adverse controls**,
1,559 total. Most positive counts are invariant-system rows and Bianchi
components, not separate physical calculations. Controls reject reversed
dual sign, a missing boundary/divergence/sterile derivative term, omitted
sources, a free gauge equation and a false weak-inserted orthogonality.

The four-structure rank is only the fixed scalar-derivative algebraic
placement. Equation (D5) exhibits the other placements and sources that a
full IBP/EOM quotient must retain; it does not close that quotient or
establish a smaller physical operator count. Evanescent finite terms,
other loop sectors and counterterms, full matching, physical running and
observables remain outside this calculation.
