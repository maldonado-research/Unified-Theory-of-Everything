# Correlated phase transport of the matter-pole packet

26 October 2026 (Pacific). This is a first-order action-coordinate calculation using the preceding independently calculated mixed and four-Higgs poles. It does not evaluate a new loop or establish a complete physical operator map, finite matching, or observable.

## Inputs and the correlated cancellation

Use the literal conventions and two primitive pole inputs in [the mixed-loop note](../mixed_yukawa/MIXED_SEED_POLE.md) and [the four-Higgs note](../../checkpoint_2026_10_23/seed_yukawa/FOUR_SCALAR_POLE.md):

```text
phi=tilde H, rho=phi^dagger phi,
R=Y_nu* M Y_nu^dagger, K=C7 R,
G=(K+K^dagger)/2, A=(K-K^dagger)/(2i), t=Tr A=Im Tr K,
kappa=1/(16 pi^2 epsilon),
J_phi^mu=i[phi^dagger D^mu phi-(D^mu phi)^dagger phi].
```

`G,A` are Hermitian. `K` is an arbitrary generally nonnormal complex flavor matrix determined by the declared seed, Yukawa and symmetric retained mass. At zero external gauge fields, the summed displayed one-loop action includes

```text
L_mixed^UV = (kappa/2)[-rho bar L G i<->slash partial L
                        +J_phi,mu bar L G gamma^mu L]
             +kappa(partial_mu rho) bar L A gamma^mu L,
L_scalar^UV includes +kappa t (partial_mu rho) J_phi^mu.
```

The remainder of the scalar pole contains the real derivative packet and the mass-cubic potential pole. The separately calculated four-lepton pole is also retained. The displayed current relation `t=Im Tr(C7 Y_nu* M Y_nu^dagger)` is algebraic; it does not imply physical cancellation.

Consider the joint, conjugate-safe first-order field substitution

```text
delta L=+i kappa rho A L,     delta bar L=-i kappa rho bar L A,
delta phi=+i kappa rho t phi, delta phi^dagger=-i kappa rho t phi^dagger.
```

With `epsilon=i sigma^2` and `phi=epsilon H*`, the Standard Model doublet obeys `H=-epsilon phi*`, so `delta H=-i kappa rho t H`. Taking `A,t` constant and flavor generators to commute with gauge matrices, direct product rules give

```text
delta L_kin,L   =-kappa(partial_mu rho) bar L A gamma^mu L,
delta L_kin,phi =-kappa t(partial_mu rho)J_phi^mu,
delta rho=0.
```

Thus, at first order in the pole and linear in `C7`, these kinetic variations cancel both displayed imaginary-current terms. The real scalar and `G`-current terms, the scalar potential pole, the four-lepton pole and all source variations below remain. This is a chosen classical action-coordinate map. The covariant product rules are an algebraic extension of the zero-gauge action; they are not an independent integration of the loop's external-gauge vertices. The transformation of the full quantum 1PI functional and its composite sources is a separate problem.

## Interactions and independent probes

With the stated Standard Model reference interactions

```text
L_Y= -bar L Y_nu phi N_R -bar L Y_e H e_R
     -bar Q Y_u phi u_R -bar Q Y_d H d_R + h.c.,
```

and no shifts of the singlet or quark fields, the action change is

```text
delta L_Y = +i kappa rho bar L (A-t I)Y_nu phi N_R
            +i kappa rho bar L (A+t I)Y_e H e_R
            -i kappa rho t bar Q Y_u phi u_R
            +i kappa rho t bar Q Y_d H d_R + h.c.                 (1)
```

These interactions are not set to zero. The scalar potential and mass terms depend on `rho` and have zero first-order variation under this phase. For independent sources `L_probe=bar eta L+bar L eta+j^dagger phi+phi^dagger j`, the same transformation yields

```text
delta L_probe = i kappa rho (bar eta A L-bar L A eta)
               +i kappa rho t(j^dagger phi-phi^dagger j).         (2)
```

Other interactions and composite sources in a fuller action must also be transformed. Equations (1)-(2) exhibit the correlated terms from this particular map, not a complete physical operator coefficient.

## An explicit lower-order insertion

If the reference action independently contains a symmetric dimension-five input

```text
L5=(1/2) bar X C5 X^c+h.c., X=phi^dagger L, C5=C5^T,
```

then `delta X=i kappa rho phi^dagger (A-t I)L`. Writing `B=A-t I`, the full conjugate field variation gives

```text
delta L5=-(i kappa rho/2) bar X (B C5+C5 B^T) X^c+h.c.         (3)
```

The transpose is on the charge-conjugated flavor slot. Since `C5` has mass dimension `-1` and `A` has mass dimension `-2`, this is a dimension-seven contribution. It is conditional on the separately supplied `C5`; the calculation does not claim that the present loop generated it. This term and any other lower-order insertions must be included consistently before deriving a physical dimension-seven coefficient.

The exact condition for the neutrino Yukawa term in (1) to vanish is

```text
(A-t I) Y_nu=0.
```

Because `A` is Hermitian, the image of `Y_nu` must lie in its eigenspace with eigenvalue `t`. For one active flavor this holds identically, so `X` and the declared neutrino Yukawa are invariant under this joint phase. Charged-lepton, quark and probe terms can still change. For more than one active flavor with full-row-rank `Y_nu`, the condition forces `A=0`: a right inverse gives `A=t I`, and tracing gives `(n_L-1)t=0`. For rank-deficient `Y_nu`, nonzero `A` can satisfy the condition; for example `Y_nu=(1,0)^T`, `M=1`, `C7=diag(i,2)` gives `A=diag(1,0)` and `t=1`. These are exact conditional flavor statements, not physical selection rules for the total EFT.

## Checks and remaining work

Run `python -I -B -O verify_combined_phase.py` in the pinned SymPy environment. Its [receipt](COMBINED_PHASE_RESULTS.json) records 16 positive checks and nine deliberate errors rejected. The standalone checker varies actual scalar first jets, flavor-matrix Yukawa interactions, independent probe sources and the literal two-flavor/weak/spinor Grassmann form of (3); it tests transpose versus adjoint and concrete one-flavor, full-rank and singular-rank examples. The general identities follow from the displayed product rule and matrix argument, not sampling alone.

This map transports a specified part of the primitive loop action. The real current and scalar packets, four-lepton contact, earlier three-point source counterterms, other mixed-dimension terms and all other interactions must still be collected in a fixed, complete retained-sterile basis. External-gauge loop vertices, evanescent and subtraction definitions, quantum source transport, finite matching and an observable are not resolved here. No sterile field is integrated out, and no new fundamental law, empirical confirmation or complete physical running is asserted.
