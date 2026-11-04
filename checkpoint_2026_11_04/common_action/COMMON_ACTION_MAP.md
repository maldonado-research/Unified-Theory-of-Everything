# One action convention for the selected source and matter counterterms

4 November 2026 (Pacific). Input commit `41f9b709c0fafdefdc9e1b990caa5421ee99f7dc`.
This is classical first-order action algebra using earlier loop poles. It
does not integrate a new loop or complete a physical operator basis or beta
function. Its concrete new projection combines the source normalization,
matter normalization, scalar-background transport and known renormalizable
inputs without adding tables evaluated in different field coordinates.

## Fix the sign and the reference action first

Let `S0` contain the canonical retained-sterile kinetic terms, the scalar
potential `-V(rho)` with `V=m_H^2 rho+lambda rho^2`, and all four declared
Yukawa interactions. Keep independent probes for the lepton, sterile and
scalar fields. Use `phi=tilde H`, `rho=phi^dagger phi`, `D=partial+i Acal`,
`P0=-Box`, `P=P0+m_H^2` and `kappa=1/(16 pi^2 epsilon)`.

Distinguish the three inputs by their action sign:

```text
+kappa P : the earlier cancelling source counterterm, including its h.c.;
+kappa Q : the primitive matter LOOP pole, including its h.c.;
+kappa R : renormalizable cancelling action from field/parameter inputs.
```

Thus the selected cancelling action contains **`+P-Q+R`**, not `P+Q+R`.
`P` is the literal six-placement/curvature packet of the
[matched source action](../../checkpoint_2026_10_15/interacting_reduction/ELECTROWEAK_REDUCTION.md).
`Q` is the three previously integrated primitive matter poles in the
[matter ledger](../../checkpoint_2026_10_30/real_matter/REAL_MATTER_TRANSPORT.md).
The independently calculated
[renormalizable inputs](../../checkpoint_2026_10_16/yukawa_counterterms/YUKAWA_TWO_POINT.md)
fix a subset of `R`; unspecified additional proper counterterms remain
separate coefficients, not zeros.

Write

```text
F_s=C7 Y_nu*, K_m=F_s M Y_nu^dagger,
G=(K_m+K_m^dagger)/2, A=(K_m-K_m^dagger)/(2i),
g=Tr G, t=Tr A, W=G/2+iA, b=g+it,
H_L=Y_nu Y_nu^dagger, H_R=Y_nu^dagger Y_nu, T_Y=Tr H_R.
```

The source and matter flavor matrices have different sizes and dimensions:
`F_s` is active-by-sterile with dimension `-3`, whereas `W` is
active-by-active with dimension `-2`. They cannot be equated.

## A single chosen field map

Use the source packet's specified local normalizing map `vC`, rather than
its different sourced-EOM map `vS`:

```text
vC_L=-Phi F_s P N_L/2,
vC_NL=-P(F_s^dagger Phi^dagger L)/2, vC_NR=(vC_NL)^c,
vC_phi=vC_B=vC_W=0,  Phi=phi tensor 1_active.
```

Every derivative acts on its complete right composite. Let the matter map
from the preceding checkpoint be

```text
w_L=rho W L, w_barL=rho bar L W^dagger,
w_phi=rho b phi, w_phi^dagger=rho b* phi^dagger,
w_N=w_B=w_W=0.
```

With **old fields equal new fields plus the displayed change**, choose

```text
Xi_old=Xi+kappa V[Xi], V=vC-w.                         (1)
```

The tilde constraint fixes `V_H=-rho b* H`. All adjoint and Majorana
partners are fixed by conjugation. In particular, the scalar part of (1)
is nonzero, although `vC` itself leaves the scalar fixed.

For the graded first variation `delta_v S0`, define

```text
P_C=P+delta_vC S0, Q_M=Q+delta_w S0.
```

Then the common action identity is

```text
(S0+kappa[P-Q+R])[Xi+kappa V]
 =S0[Xi]+kappa[P_C-Q_M+R][Xi]+O(kappa^2).              (2)
```

The `R` action is already one loop, so substituting (1) into it only
changes second-loop-order terms. The source packet and primitive matter
poles likewise are not transformed a second time at this order. Equation
(2) follows from the literal first variation; it is not a sum of two
independently canonicalized coefficient tables.

If starting with the earlier sourced representative `P_S`, the necessary
bridge is `vC-vS-w`, since `P_S=P+delta_vS S0`. Applying `vC-w` directly
to `P_S` would count `vS` again. The
[all-field bridge](../../checkpoint_2026_10_16/action_completion/FIELD_MAP_COMPOSITION.md)
includes the scalar, gauge, Majorana and probe changes; a fermion-only
bridge is insufficient off shell.

For explicit interpretation of (2), `P_C` retains the previously derived
`R_J+R_D` interaction remainder, differential Majorana mass terms and
probe terms. `Q_M` retains the real scalar/lepton currents, primitive
four-lepton contact, transported potential, all four Yukawa changes and
probes of the preceding matter ledger. Their signs here are reversed.
For example the matter part of the cancelling scalar potential is

```text
-2[Re K3-g m_H^2]rho^2+4g lambda rho^3,
K3=Tr(C7 Y_nu* M M^dagger M Y_nu^dagger).              (3)
```

The source-normalizing `vC` gives zero charged-lepton and quark spectator
variation by `phi^dagger H=0` and its unchanged scalar/quark fields. Their
joint variations therefore equal minus their previously calculated `w`
variations. The neutrino interaction and Majorana action also contain
the nonzero `vC` terms; they are kept in `P_C`, not silently discarded.

## Renormalizable bare inputs and the tree seed

Let `z` denote the coefficients of the already calculated constant bare
field changes:

```text
z_L=-H_L L/4, z_phi=-T_Y phi/2,
z_NR=-H_R N_R/2, z_NL=-H_R^T N_L/2.
```

The known parameter pole coefficients are

```text
rY=H_L Y_nu/4+Y_nu H_R/2+T_Y Y_nu/2,
rM=(H_R^T M+M H_R)/2,
r_m2=m_H^2 T_Y-2Tr(H_R M^dagger M).                  (4)
```

Expanding `S0[Xi+kappa z,p+kappa r]` gives the corresponding subset of
`R`. Directly in the action, the pure-neutrino-Yukawa cubic coefficient
and sterile chirality-flipping mass coefficient cancel to zero, whereas
the scalar mass action contains `+2Tr(H_R M^dagger M)rho`. This does not
set the parameter poles (4) to zero. The field changes also act on the
quartic and spectator Yukawas; their proper parameter poles have not all
been calculated here and must remain separate. Equivalently, the chosen
common bare map is `Xi_0=Xi+kappa(z+V)` with the same parameter inputs.

The tree seed is **not** already a one-loop term. Its field variation
under `z` must be included at first-loop order:

```text
S7=(1/2)bar X C7 P0 X^c+h.c., X=phi^dagger L,
Delta C7_fields=-kappa r7,
r7=(H_L C7+C7 H_L^T)/4+T_Y C7.                       (5)
```

The supplied external-field-compensation part of the parameter pole is
`+kappa r7`. It cancels (5) **once** in that declared bare convention.
This identifies one subset of the parameter counterterm; other proper
and source/basis contributions remain undetermined. Substituting `V`,
which is itself linear in `C7`, into `S7` is quadratic in that seed and
is excluded from the selected linear-seed calculation. In contrast,
an additional inverse-`z` canonicalization would act on `S7` at order
`kappa C7 Y_nu^2` and cannot be discarded as `kappa^2`.

## A concrete common neutral quadratic projection

Set the **new** scalar coordinate to `phi=(u,0)` with real constant `u`,
and gauge/fermion backgrounds to zero. No scalar stationarity equation
is imposed. For a common real momentum invariant write

```text
q=(nu_L,N_L)^T, H=u(p^2+m_H^2)F_s,
Z=[0 H; H^dagger 0], D_s=Z+diag(u^2 G,0),
R_z=diag(-H_L/2,-H_R^T).
```

The literal selected cancelling kinetic action has matrix
`I+kappa(D_s+R_z)`: the matter contribution is **positive** `u^2 G`
because it comes from `-Q`. On this background (1) has fermion Jacobian

```text
q_old=(I+kappa J)q,
J=[-u^2 W  -H/2; -H^dagger/2  0],
J+J^dagger=-D_s.                                    (6)
```

Hence the transformed kinetic matrix is
`I+kappa R_z`, not `I`. This chosen map normalizes the seed-dependent
part and retains the standard renormalizable field counterterms.

Fix the holomorphic left-Weyl mass convention

```text
L_mass=-(q^T epsilon_spin Mcal q+h.c.)/2,
Mcal=[0 B; B^T M*], B=uY_nu*, M=M^T.
```

The scalar part of (1) gives `delta phi=-kappa u^3 b(1,0)` and thus an
additional holomorphic Dirac block `delta B_scalar=-u^3 b*Y_nu*`.
Including both fermion slots and this scalar-background change gives

```text
delta Mcal=J^T Mcal+Mcal J+[0 delta B_scalar;
                              delta B_scalar^T 0],
delta A=-(H*B^T+B H^dagger)/2,
delta B=-H*M*/2-u^3[W^T+b*I]Y_nu*,
delta D=-(H^T B+B^T H)/2.                            (7)
```

The transpose in the holomorphic mass and the conjugate scalar phase
are essential. Complex symmetric retained masses are allowed because
this chosen normalizing map does not use the real-mass EOM powers of
the alternative sourced representation. The known renormalizable bare
inputs cancel in the literal neutral mass action before (1); their
nonzero parameter poles remain in the bare-parameter relations.

For an independent source row `bar s q+h.c.`, the new term is
`kappa bar s J q+h.c.`. Independent scalar probes also acquire
`-kappa rho[b j^dagger phi+b*phi^dagger j]`. Omitting either (7)'s
background term or the source transport changes the declared action.
The one-active-flavor limit retains a radial neutrino Yukawa coefficient
`+3kappa g u^3Y_nu/2` in the cancelling action even when the earlier
joint phase alone left that Yukawa unchanged.

## An independent lower-order coefficient is still an input

If `S5=(1/2)bar X C5 X^c+h.c.` is independently retained, add its actual
variation `kappa delta_(z+V) S5` as appropriate to the chosen bare action.
It is not contained in the renormalizable `S0` used to define (2).
With `E=W^dagger+bI=G/2+gI-i(A-tI)`, the coordinate part is

```text
delta_V L5=-rho/2 bar X(E C5+C5 E^T)X^c
 -rho/4[overline(P N_L)F_s^dagger C5 X^c
                 +bar X C5 F_s*P N_R]+h.c.            (8)
```

The first term is the negative of the previous conditional active
dimension-seven insertion. The second is a retained-sterile contact
from `vC`; dropping it would not be the same field map. Its derivative
piece has canonical operator dimension eight, and its `m_H^2` piece
has dimension six with an explicit mass-squared coefficient. The
`z` part separately has the external-field coefficient form (5) with
`C7` replaced by `C5`; no total `C5` parameter pole is inferred.

Canonical dimension and inverse-heavy-scale order are distinct when
retained masses appear. In a strict truncation retaining only
`O(Lambda^-3)` with fixed masses and `C5=O(Lambda^-1)`, products
`C7 C5` are beyond that truncation. If canonical-dimension or mixed-input
bookkeeping instead retains them, (8) supplies the conditional terms
and their correlated retained-field partners. No independent `C5` is
generated or identified with a seesaw coefficient by this calculation.

## A ledger that prevents false coefficient addition

The entries below describe literal input classes before a complete EOM
quotient. `Y` degree counts explicit neutrino Yukawa factors and their
conjugates in the coefficient; retained masses are displayed separately.
The lepton-charge convention is `L,N_R: +1`, `N_L: -1`, `phi:0`, with
spurions `C7:+2`, `M:-2` and the opposite charges for their conjugates.

| Class | Canonical operator dimension | Coefficient | `Y` degree | Mass order | Spurion charge of coefficient |
| --- | --- | --- | --- | --- | --- |
| Original derivative seed | 7 | `C7` | 0 | 0 | `+2` |
| Derivative part of source packet | 7 | `F_s=C7Y*` | 1 | 0 | `+2` |
| Scalar-mass part of source packet | 5 | `m_H^2 F_s` | 1 | scalar mass squared | `+2` |
| Matter four-lepton and derivative poles | 6 | `C7Y*MY^dagger` | 2 | one `M` | 0 |
| Primitive scalar potential pole | 4 | `Tr(C7Y*MM^dagger MY^dagger)` | 2 | three masses | 0 |
| Field-compensation seed parameter input | 7 | `r7` | 2 | 0 | `+2` |
| Radially induced Yukawa contacts | 6 | `(W^dagger+bI)Y_nu` and spectators | 3 for neutrino | one retained `M` | 0 |
| Conditional active `C5` transport | 7 | `E C5+C5 E^T` | 2 before any `C5` model | one retained `M` | `+2` |
| Conditional derivative mixed-sterile `C5` contact | 8 | `F_s^dagger C5` or `C5 F_s*` | 1 before any `C5` model | 0 | 0 |

The source packet, its mass descendants and the physical operator
coefficients cannot be compared merely because a neutral projection has
the same number of fermions. In particular an LNC dimension-six matrix
containing `M` is not the LNV dimension-seven seed coefficient. The
common action (2), not a naive sum of projected matrices, is the starting
point for the remaining independent-basis reduction.

## Reproduction and remaining work

Run `python -I -B -O verify_common_action.py` in the root pinned symbolic
environment. Its adjacent receipt records literal exterior-field kinetic,
mass, probe, seed and conditional `C5` variations, complex noncommuting
flavor inputs and deliberate wrong-sign, omitted-background, wrong-adjoint
and double-map controls. It imports no project calculations and uses
explicit failures active under optimization. The general identities
follow from the displayed first variation and matrix algebra, rather
than numerical examples alone.

This assembles selected existing inputs in one action convention and
checks a definite common neutral projection. It does not complete every
source-sector coefficient, remaining coupling-sector counterterms,
derivative/EOM or evanescent basis, quantum source map, gauge-loop
completion of the new matter poles, physical running, finite matching or
an observable. No sterile field is integrated out, and no new law,
unification, empirical confirmation or external novelty is asserted.
