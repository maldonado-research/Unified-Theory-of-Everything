# Primitive pole transport, flavor image and retained scalar sources

23 October 2026, Pacific. Input commit: `f1f9a892e15545c70e4fa05827ea7fefc5a438eb`.

This note transports the independently calculated [four-lepton pole](FOUR_LEPTON_POLE.md) and [four-scalar pole](FOUR_SCALAR_POLE.md). It performs no additional loop integration. It fixes their normalized flavor image, counterterm-to-RG weights and an optional local Higgs phase change with interaction and probe sources retained. These are specified projections in unreduced Green coordinates, not a complete physical operator basis or physical running vector.

## Symmetric-pair image of the four-lepton pole

Let

```text
C=C^T, R=Y* M Y^dagger=R^T,
B_(alpha u)(beta v)=l_alpha^{u T} epsilon_spin l_beta^v,
epsilon_spin=[[0,1],[-1,0]].
```

The odd-spinor pair is symmetric under exchange of its two complete flavor/weak labels. With barred components conjugate to the unbarred ones, `B^dagger=-bar B`. The primary ordered loop result is

```text
L_loop,4L^UV=(kappa/2) C_alpha beta R_gamma delta
                    sum_(u,v) bar B_(alpha u)(beta v)
                              B_(gamma u)(delta v) + h.c.,
kappa=1/(16 pi^2 epsilon).
```

Choose real orthonormal symmetric flavor matrices `E_I` and weak matrices `e_A`, orthonormal in the Frobenius metric. Their off-diagonal entries are `1/sqrt(2)`. Define actual pair fields

```text
q_IA = sum_(alpha beta u v) (E_I)_alpha beta (e_A)_uv
                                      B_(alpha u)(beta v),
c_I=Tr(E_I^T C), r_I=Tr(E_I^T R), v_I=r_I*,
H=c v^dagger+v c^dagger.
```

The weak index `A=1,2,3` is the symmetric doublet-pair triplet. Contracting a symmetric flavor coefficient into `B` makes its remaining weak pair symmetric as well; the antisymmetric-flavor/antisymmetric-weak channel is absent from this primitive contribution. The fully normalized ordered action is

```text
L_loop,4L^UV=-(kappa/2) sum_(I J A) H_IJ q_IA^dagger q_JA,
L_ct,4L=+(kappa/2) sum_(I J A) H_IJ q_IA^dagger q_JA.      (1)
```

The executable checks (1) by expanding every actual Grassmann component for two active flavors and both weak components. This fixes the simultaneous flavor/weak off-diagonal normalization; it is not inferred from separate informal square-root factors. Its general derivation follows from completeness of the two symmetric-pair bases and `bar B=-B^dagger`.

The image of `H` lies in `span(c,v)`, because `H x=c(v^dagger x)+v(c^dagger x)`. Consequently its rank is at most two per weak triplet block, or six across the three weak components. Let `a=v^dagger c`, `N_c=c^dagger c`, and `N_v=v^dagger v`. The two possibly nonzero eigenvalues are

```text
h_plus/minus=Re(a) +/- sqrt(N_c N_v-[Im(a)]^2).
```

Their product is `|a|^2-N_c N_v`. The Gram determinant

```text
N_c N_v-|a|^2=sum_(I<J)|c_I v_J-c_J v_I|^2
```

is positive for independent sources, so that case has one positive and one negative eigenvalue. Proportional, vanishing and one-dimensional pair spaces can have lower rank. These are algebraic properties of the single primitive pole. Other action contributions can populate its kernel or change its signature. It is neither a positivity bound nor a rank prediction for a total physical four-lepton coefficient.

Under the preceding passive flavor convention,

```text
Y'=V^dagger Y U*, M'=U^dagger M U*, C'=V^dagger C V*,
R'=V^T R V.
```

Thus `c` and `v` transform in the same unitary symmetric-pair representation and `H` by unitary similarity. The actual exterior expansion, complex-mass ordering and that covariance are checked independently of any mass inverse.

## Bare-mu powers and projected RG entries

The [topology certificate](SEED_TOPOLOGY.md) proves that the original `bar L bar L phi phi` proper four-point vertex is absent at one seed and total pure-Yukawa degree two. It also distinguishes external-field attachments from the surviving reducible earlier three-point source bridge. In the fixed original unreduced Green coordinates, the target seed coefficient therefore has the already calculated parameter pole

```text
r7=(H_L C+C H_L^T)/4+T C,
delta C|target Green=kappa r7,
H_L=Y Y^dagger, T=Tr(Y^dagger Y),
16 pi^2 beta_C|target Green,pureY=2r7.                    (2)
```

Equation (2) closes this specified projection, with all other higher-operator inputs held as declared. It does not remove the earlier `C Y*` three-point counterterms or their sourced field/EOM transport into a physical basis.

For the generated four-field operators, explicitly choose the parameter conventions

```text
L_4L=sum_(I J A) D_IJ q_IA^dagger q_JA,
L_scalar=-a4 rho^2+d O+d* O^dagger,
rho=phi^dagger phi,
O=(phi^dagger partial_mu phi)^2
       -rho partial_mu phi^dagger partial^mu phi.
```

Here `a4` names the coefficient of `rho^2` in this displayed convention; it must be converted before comparison with any earlier quartic normalization. The scalar-loop inputs are `K1=Tr(C Y* M Y^dagger)` and `K3=Tr(C Y* M M^dagger M Y^dagger)`. The primitive parameter residues in these precise coordinates are

```text
delta D=kappa H/2,
delta d=-kappa K1,
delta a4=kappa(K3+K3*).                              (3)
```

The last sign follows from `L_scalar=-a4 rho^2`: the cancelling action is negative the computed positive scalar potential pole. There is no additional identical-pair half in the declared `D` convention.

In `d=4-2epsilon`, each of these four-field recipient parameters has a bare prefactor `mu^(2epsilon)`. A four-fermion operator and the displayed four-scalar derivative operator have dimension `6-4epsilon`, while `rho^2` has dimension `4-4epsilon`. The four-dimensional parameter dimensions are respectively `-2,-2,0`. The original seed has the same explicit bare prefactor, whereas `Y0=mu^epsilon(...)` and the dimensionful `M0` has no explicit mu power.

Treat conjugate parameters independently during differentiation. Every residue in (2)-(3) is linear in a weight-two seed coefficient and quadratic in weight-one Yukawas; mass powers have weight zero. Its canonical weight is four. Holding each bare recipient fixed gives `(4-2)r=2r`, hence the following **primitive projected MS inputs**:

```text
16 pi^2 beta_D|this pole=H,
16 pi^2 beta_d|this pole=-2 K1,
16 pi^2 beta_a4|this pole=2(K3+K3*).                   (4)
```

These are statements about the specified unreduced action coefficients and the calculated monomial, not total physical dimension-six or Higgs beta functions. The mass factors allow a dimension-seven insertion to feed lower-dimensional action terms while the sterile field remains present. No threshold expansion, sterile removal, fit or prediction of the input masses/couplings has occurred. Finite terms and evanescent definitions remain separate; the four-lepton flavor image refers to the primary physical-spinor ordered projection, not a completed `d`-dimensional Fierz basis.

## A complex derivative coefficient and an explicit phase change

The scalar pole must keep its conjugate derivative placement. Write

```text
a_mu=phi^dagger partial_mu phi,
b_mu=(partial_mu phi^dagger)phi=a_mu*,
J_mu=i(a_mu-b_mu), partial_mu rho=a_mu+b_mu,
A=rho partial_mu phi^dagger partial^mu phi.
```

Then its Hermitian derivative packet is exactly

```text
K1(a_mu a^mu-A)+K1*(b_mu b^mu-A)
 = Re(K1){[(partial rho)^2-J_mu J^mu]/2-2A}
       +Im(K1)(partial_mu rho)J^mu.                     (5)
```

Replacing the packet by twice its real coefficient times one placement would lose the second term. It cannot simply be erased by a free scalar-current equation in the interacting theory.

For real constant `alpha` of dimension `-2`, perform the first-order field change

```text
delta phi=i alpha rho phi, delta phi^dagger=-i alpha rho phi^dagger.
```

Direct differentiation gives

```text
delta rho=0,
delta L_kin=-alpha(partial_mu rho)J^mu.                 (6)
```

The scalar mass and potential are unchanged by (6). If this change is used to cancel the displayed imaginary loop packet, its value is `alpha=kappa Im(K1)`. At first loop order, linear in `C`, its action on the loop itself is higher order; the original action's sources must still be varied.

For the declared neutrino interaction `L_Y=-(V_Y+V_Y^dagger)`, with `V_Y=bar L Y phi N_R`, and an optional independent scalar probe `L_probe=j^dagger phi+phi^dagger j`, the literal induced terms are

```text
delta L_Y=-i alpha rho V_Y+i alpha rho V_Y^dagger,
delta L_probe=i alpha rho(j^dagger phi-phi^dagger j).     (7)
```

Equations (6)-(7) are one source-retaining action-coordinate example, not a completed physical reduction. Additional interactions in a larger reference action must transform as well; spectator Yukawas are not silently omitted from a claimed full Standard Model map. The zero-gauge calculation does not independently evaluate the external-gauge vertices of a covariant lift.

The checker substitutes `phi=x+i y` with independent first derivatives, verifies (5)-(7) component by component, and rejects the lost imaginary, Yukawa and probe terms. The real local phase vector has zero pointwise divergence; this check is not a general quantum effective-action/source theorem.

## Reproduction

Run `python -I -B -O verify_green_mixing.py` in the pinned SymPy environment. The explicit-guard standalone program imports no repository modules. [Its receipt](GREEN_MIXING_RESULTS.json) records 71 positive checks and 11 rejected mutations. These are algebraic checks using independently integrated poles, not another loop integration or independent human peer review. Complete mixed matter poles, the original sourced reductions, lower-order insertions, dimensional/evanescent conventions and physical coefficient collection are still required before predicting an observable.
