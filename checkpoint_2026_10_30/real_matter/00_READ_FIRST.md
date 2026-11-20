# Real matter normalization and the remaining scalar geometry

30 October 2026 (UTC). Provisional GitHub candidate.15, using research input commit `5094584ebf402736b0bc5b597f6769f640178fd1`.

The preceding checkpoint transported selected imaginary-current poles through correlated lepton and Higgs phases. This continuation also normalizes selected **real** kinetic terms and records the interactions that the field changes induce. It uses the previously integrated matter poles; it does not evaluate another loop.

The result makes the remaining calculation more concrete. Some derivative terms move into Yukawa interactions, the scalar potential and independent probe sources. A scalar derivative interaction survives, and its nonzero field-space curvature explains why local scalar coordinates alone cannot remove the entire selected interaction. This is an application of established field-space geometry, not new fundamental mathematics or a spacetime-gravity result.

## The selected action map

Use the conventions of the [mixed-loop calculation](../../checkpoint_2026_10_26/mixed_yukawa/MIXED_SEED_POLE.md), [four-scalar pole](../../checkpoint_2026_10_23/seed_yukawa/FOUR_SCALAR_POLE.md) and [correlated phase map](../../checkpoint_2026_10_26/combined_matter/COMBINED_PHASE_TRANSPORT.md). Write

```text
phi=tilde H, rho=phi^dagger phi, kappa=1/(16 pi^2 epsilon),
K=C7 Y_nu* M Y_nu^dagger,
G=(K+K^dagger)/2, A=(K-K^dagger)/(2i),
g=Tr G, t=Tr A,
J_phi,mu=i[phi^dagger D_mu phi-(D_mu phi)^dagger phi].
```

In the convention `old field = new field + delta field`, combine the earlier phase changes with

```text
delta L=kappa rho (G/2+iA)L,
delta phi=kappa rho(g+it)phi.
```

The adjoint fields are transformed by Hermitian conjugation; `H=-epsilon phi*` implies `delta H=kappa rho(g-it)H`. Direct substitution in the reference kinetic action cancels the selected lepton normalization term, the scalar `rho |D phi|^2` term and both previously identified phase-current terms. Among the residual derivative terms are

```text
+(kappa/2) J_phi,mu bar L G gamma^mu L
+(kappa g/2)[3(partial_mu rho)^2-J_phi,mu J_phi^mu].
```

The [action derivation and source ledger](REAL_MATTER_TRANSPORT.md) retains all four declared Yukawa interactions, independent lepton/scalar probes and the scalar potential. The primitive four-lepton pole also remains. An independently supplied dimension-five coefficient has its own conditional induced dimension-seven term; that input is not inferred from the seed loop. A field change is not permission to erase its induced interactions.

## A precise limit on scalar normalization

Express the complex doublet using four real components, `phi=(q1+i q2,q3+i q4)/sqrt(2)`. The residual scalar kinetic metric is

```text
metric_ij=delta_ij+kappa g[3 q_i q_j-(Jq)_i(Jq)_j],
```

where `J` is the constant real complex structure. The [geometry calculation](SCALAR_GEOMETRY.md) computes the same linearized Riemann tensor before and after the radial field change. With its explicit sign convention, the sectional curvature for orthonormal real vectors `X,Y` is

```text
3 kappa g [1+(X^T J Y)^2].
```

Consequently the selected metric is not flat when `g` is nonzero. A local, invertible, derivative-free reparametrization of the scalar fields cannot flatten it throughout a neighborhood. Normal coordinates can still normalize the metric at one point. Derivative-dependent field changes and equations of motion can redistribute operators with their sources retained.

This is the coefficient of a formal UV-pole contribution in a specified truncated scalar sector. It is not a finite renormalized curvature, a curvature of spacetime, a positivity theorem or an obstruction to completing the full EFT. Other independent operator inputs can change the total metric.

## Reproduce and assess

| Material | Scope |
| --- | --- |
| [Action transport](REAL_MATTER_TRANSPORT.md) · [checker](verify_real_matter.py) · [receipt](REAL_MATTER_RESULTS.json) | Selected real and phase normalization, residual terms and interaction/source variations |
| [Independent action audit](INDEPENDENT_REAL_AUDIT.md) · [checker](independent_real_audit.py) · [receipt](INDEPENDENT_REAL_RESULTS.json) | Separately authored scalar jets, Grassmann kinetic/Yukawa/probe and conditional lower-order insertion checks |
| [Scalar geometry](SCALAR_GEOMETRY.md) · [checker](verify_scalar_geometry.py) · [receipt](SCALAR_GEOMETRY_RESULTS.json) | Linearized curvature before/after the field change and its restricted coordinate obstruction |
| [Literature comparison](LITERATURE_REAL_MATTER.md) · [provenance](LITERATURE_PROVENANCE.json) | Bounded current searches, primary-source comparisons and coverage limits |

Run the three adjacent programs in the root pinned Python environment. They use explicit failure checks and are intended to run with `python -I -B -O`; the root publication runner itself must run without optimization. The [full replay receipt](../../verification/REAL_MATTER_REPLAY_LINUX.json) records the extended checks on a disposable copy. Computational checks and AI-assisted independent derivations do not constitute independent human peer review or empirical validation.

## The next physical step

The selected primitive matter-pole representative is now more explicit. The earlier three-point source counterterms, renormalizable field and parameter inputs, other interaction sectors and independent lower-order coefficients must still be collected in one retained-sterile operator/source convention. Complete the derivative/EOM quotient, external-gauge matter amplitudes, subtraction and evanescent definitions, and redundant-coefficient differentiation before claiming a physical running vector. Quantum composite-source transport, finite matching and one observable remain open.

Research PR #2 precedes the source-pinned website update in PR #58; both remain open for review. No merge, live deployment, release or new Zenodo version is made. DOI `10.5281/zenodo.22319172` continues to identify only the unchanged N00AK-r1 baseline.
