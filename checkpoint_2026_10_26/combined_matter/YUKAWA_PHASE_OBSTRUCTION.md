# When the selected joint phase leaves the neutrino source unchanged

26 October 2026 (Pacific). Input: the independently evaluated [mixed pole](../mixed_yukawa/MIXED_SEED_POLE.md) and [four-scalar pole](../../checkpoint_2026_10_23/seed_yukawa/FOUR_SCALAR_POLE.md), from input commit `9b1e18450c0c1e18fb6f917105368b05ca123df2`. This is an exact flavor-algebra statement for one chosen first-order action-coordinate change. It is not a no-go theorem for other field redefinitions, a complete physical basis map, or a loop or observable prediction.

## The two correlated phase terms

Set

```text
R=Y* M Y^dagger, K=C R,
A=(K-K^dagger)/(2i), a=Tr A=Im Tr K,
rho=phi^dagger phi, phi=tilde H.
```

Here `C=C^T` is the original seed matrix, `M=M^T` is the constant retained Majorana mass, and `Y` has active rows and sterile columns. `A` is Hermitian and `a` is real. Both may be nonzero. The mixed lepton/Higgs pole contains `kappa (partial_mu rho) bar L A gamma^mu L`; the scalar pole contains `kappa a (partial_mu rho) J_phi^mu`. With the preceding literal conventions, the selected first-order changes

```text
delta L=i kappa rho A L,
delta phi=i kappa rho a phi
```

cancel those two displayed terms through the kinetic actions. The other terms, gauge and probe sources, interaction variations and their Hermitian conjugates remain. In particular the neutrino interaction generates

```text
delta L_Ynu=+i kappa rho bar L (A-a I)Y phi N_R+h.c.     (1)
```

The question addressed here is exactly when the coefficient `(A-a I)Y` in this *particular coordinate change* vanishes. It does not decide whether a different action-coordinate choice can be made or whether a physical observable vanishes.

## A rank-sensitive obstruction

Let `P` be the orthogonal projector onto the active-flavor image of `Y`, and `r=rank(Y)`. The proof works for rectangular, complex and singular `Y`, and uses no inverse of `M`.

Since `Y^dagger(1-P)=0`,

```text
R(1-P)=0, K(1-P)=0,
(1-P)A(1-P)=0.                                        (2)
```

Equation (1) vanishes exactly when `(A-a I)Y=0`, or equivalently

```text
A P=a P.                                               (3)
```

Because `A=A^dagger` and `a` is real, (3) also gives `P A=a P`; its off-diagonal blocks across `P` vanish. Together with (2), this forces

```text
A=a P,        a=Tr A=r a.                              (4)
```

The reverse implication is immediate: if `A=aP`, then `(A-aI)Y=a(P-I)Y=0`. Thus:

- If `rank(Y)>=2`, the selected joint phase leaves the neutrino interaction unchanged **if and only if `A=0`**. In this case the two displayed imaginary-current coefficients already vanish.
- If `rank(Y)=1`, it is unchanged **if and only if `A=aP`**. This allows a nonzero phase acting on the one coupled active direction.
- If `rank(Y)=0`, then `Y=0`, hence `R=K=A=0`.

This sharpens the earlier one-flavor example: it is the rank of the coupled active subspace, rather than only the number of listed flavors, that controls this particular neutrino-source cancellation.

The criterion is basis covariant. For passive unitary rotations

```text
Y'=V^dagger Y U*, M'=U^dagger M U*, C'=V^dagger C V*,
K'=V^dagger K V, A'=V^dagger A V, P'=V^dagger P V,
(A'-a I)Y'=V^dagger (A-a I)Y U*.
```

The trace and rank are unchanged. `P` is a conceptual orthogonal support projector; neither the statement nor the proof requires a sterile mass inverse, positive masses, or invertible Yukawas. A construction from an independent column basis may use a nonsingular *column Gram matrix*, which is not a sterile-mass inverse.

## The charged-lepton spectator does not disappear automatically

The phase of the ordinary Higgs doublet is opposite to that of `phi=tilde H`: `delta H=-i kappa rho a H`. For the reference interaction `L_Ye=-bar L Y_e H e_R+h.c.`, direct variation gives

```text
delta L_Ye=+i kappa rho bar L (A+a I)Y_e H e_R+h.c.     (5)
```

Thus the criterion for (1) says nothing by itself about (5). In particular, with one active flavor `A=a`, the neutrino source in (1) vanishes while the charged-lepton source in (5) is proportional to `2a Y_e` and can be nonzero. External probes and other interactions must also be varied in any larger reference action. A source-free result for the displayed neutrino term is not a complete Standard Model or quantum effective-action map.

## Reproduction

Run `python -I -B -O verify_yukawa_phase_obstruction.py` in the pinned SymPy environment. The standalone explicit-guard checker imports no repository modules. Its [receipt](YUKAWA_PHASE_OBSTRUCTION_RESULTS.json) records exact projector, trace, source and flavor-covariance checks, singular/zero-mass and rank 0–3 witnesses, and deliberate rejected omissions. Finite examples check implementation; equations (2)–(4) supply the general proof. The result does not classify the full derivative/EOM quotient, calculate new loops, determine a total physical Wilson coefficient, or establish an observable or unification.
