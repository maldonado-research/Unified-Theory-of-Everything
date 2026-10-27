# Correlated matter-pole field changes and flavor sources

26 October 2026 (Pacific). Provisional candidate.14, continuing research commit `9b1e18450c0c1e18fb6f917105368b05ca123df2`.

This continuation combines specified terms from the previously evaluated mixed lepton–Higgs and four-Higgs UV poles. It constructs a correlated first-order field change, tracks the induced interaction and probe sources, and proves a flavor-rank condition for that particular change. These are classical action-coordinate results using existing pole inputs; no new loop integral or physical observable is calculated here.

## The selected action map

Use the preceding literal conventions and define

```text
phi=tilde H, rho=phi^dagger phi,
R=Y_nu* M Y_nu^dagger, K=C7 R,
A=(K-K^dagger)/(2i), t=Tr A=Im Tr K,
kappa=1/(16 pi^2 epsilon).
```

The earlier mixed and four-Higgs loop packets contain, respectively,

```text
+kappa (partial_mu rho) bar L A gamma^mu L,
+kappa t (partial_mu rho) J_phi^mu,
J_phi^mu=i[phi^dagger D^mu phi-(D^mu phi)^dagger phi].
```

At first order in the displayed pole and linear in `C7`, the correlated changes

```text
delta L=i kappa rho A L,     delta phi=i kappa rho t phi,
```

with conjugate changes of the adjoints, cancel these two terms through the lepton and scalar kinetic actions. Since `phi=epsilon H*` and `epsilon^2=-1`, the inverse is `H=-epsilon phi*` and `delta H=-i kappa rho t H`. The source terms induced by these changes remain. In particular, the neutrino Yukawa variation is

```text
delta L_Ynu=+i kappa rho bar L (A-t I) Y_nu phi N_R+h.c.
```

The charged-lepton variation instead contains `(A+t I)Y_e`; up- and down-quark terms and independent lepton/scalar probes are also explicitly transported. The [producer derivation](COMBINED_PHASE_TRANSPORT.md) and [independent action audit](INDEPENDENT_COMBINED_AUDIT.md) give the complete declared expressions. The real-current and scalar packets, four-lepton contact and earlier three-point source terms remain. Covariant product rules do not themselves evaluate additional external-gauge loop vertices.

If a symmetric dimension-five coefficient `C5` is independently present, the field change also induces a dimension-seven contact proportional to `(A-t I)C5+C5(A-t I)^T`; its literal half and charge-conjugate transpose are checked. This is a conditional lower-order insertion, not a new seed-loop coefficient.

## A precise flavor condition

Let `P` be the projector onto the image of `Y_nu`, and `r=rank(Y_nu)`. Because `R(1-P)=0`, one has `(1-P)A(1-P)=0`. The neutrino Yukawa variation vanishes exactly when

```text
(A-t I)Y_nu=0, equivalently AP=tP.
```

Hermiticity then forces `A=tP`, and its trace gives `t=r t`. Thus:

- At rank zero, `A=0`.
- At rank one, `A=tP` permits a nonzero special case.
- At rank at least two, the variation vanishes if and only if `A=0`.

The [flavor proof and checker](YUKAWA_PHASE_OBSTRUCTION.md) retain singular complex Yukawa and Majorana matrices and active/sterile unitary rotations. Even in the one-flavor neutrino-source-free case, charged-lepton, quark and probe sources can change. The statement is about this selected coordinate change, not a total physical selection rule or an observable prediction.

## Reproduce and assess

- [Combined action map](COMBINED_PHASE_TRANSPORT.md), [checker](verify_combined_phase.py), [receipt](COMBINED_PHASE_RESULTS.json).
- [Independent action audit](INDEPENDENT_COMBINED_AUDIT.md), [checker](independent_combined_audit.py), [receipt](INDEPENDENT_COMBINED_RESULTS.json).
- [Flavor obstruction](YUKAWA_PHASE_OBSTRUCTION.md), [checker](verify_yukawa_phase_obstruction.py), [receipt](YUKAWA_PHASE_OBSTRUCTION_RESULTS.json).
- [Bounded literature assessment](LITERATURE_COMBINED.md) and [provenance](LITERATURE_PROVENANCE.json).

The three standalone checkers use explicit failure guards and no repository imports. The [full replay](../../verification/COMBINED_MATTER_REPLAY_LINUX.json) includes them alongside the preserved earlier checks. Internal AI-assisted derivations and implementations do not constitute independent human peer review.

The complete retained-sterile operator/source basis, the original three-point source counterterms, other lower-order insertions, gauge and evanescent definitions, the quantum effective-action/source map, finite matching and an observable remain open. The flavor condition does not establish a complete physical running or unification. Earlier science and receipts are preserved; research PR #2 precedes website PR #58, and both remain open for review. No new Zenodo version or DOI is made.
