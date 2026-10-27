# Mixed lepton–Higgs pole and interaction sources

26 October 2026 (Pacific). Provisional candidate.13; input commit `8999d3de1a0e3dc2576fe1586295426da687fc67`.

This continuation computes the remaining mixed two-lepton/two-Higgs proper triangle at one original derivative-seed insertion and two neutrino-Yukawa vertices. The retained Majorana mass produces a momentum-dependent local UV pole. The calculation is for the declared zero-external-gauge-field sector; the full physical operator map, gauge completion, finite matching and observables remain open.

## Result

Use `d=4-2 epsilon`, `kappa=1/(16 pi^2 epsilon)`, `phi=tilde H`, the original symmetric seed coefficient `C7`, and the symmetric retained Majorana mass `M`. Define

```text
R = Y* M Y^dagger,       K = C7 R,       rho = phi^dagger phi.
```

The [Feynman-parameter calculation](MIXED_SEED_POLE.md) and the [independent angular/Wick audit](INDEPENDENT_MIXED_AUDIT.md) give the zero-gauge primitive loop action

```text
L_loop^UV = (i kappa/2) K_ab phi^v
  [(partial_mu bar L_a^u) phi_dagger^v
   - bar L_a^u (partial_mu phi_dagger^v)] gamma^mu L_b^u + h.c.
```

The displayed weak indices pair each lepton current and scalar bilinear; no free weak trace is appended. The cancelling counterterm is the negative of this action. The same-composite derivative placement is UV finite. The opposite placement has a local momentum-linear pole: in the producer's all-incoming Fourier convention, the stripped vector residue is `(p_barL-p_phi_dagger)/2`. The two routes initially used opposite loop-momentum directions for the ordered fermion propagator; the independent audit records the corresponding sign correction. Retaining the derivative on the full scalar/lepton composite is essential.

The prior four-lepton and four-Higgs poles and this mixed pole cover the three proper external-matter sectors identified by the earlier [topology certificate](../../checkpoint_2026_10_23/seed_yukawa/SEED_TOPOLOGY.md), at the stated seed insertion and Yukawa order. This statement concerns those primitive zero-gauge graphs. Reducible sterile-source bridges and their counterterms remain separate, and other interactions and operator insertions are outside this calculation.

## Interaction sources and field changes

The independently checked [source identity](MIXED_SOURCE_IDENTITY.md) uses the coefficient-free representative

```text
O(K)=i K_ab[(D_mu phi)^dagger phi (bar L_a gamma^mu L_b)
            -rho (D_mu bar L_a) gamma^mu L_b].
```

At zero gauge fields the computed mixed pole is `-kappa O(K)/2+h.c.`. The identity keeps the boundary term, the scalar current, the full declared lepton interaction sources, the Hermitian conjugate and a conjugate-safe first-order field change. Its covariant product-rule identity does not itself evaluate the loop's external-gauge vertices. A quantum effective-action/source map, dimensional and evanescent conventions, lower-order insertions and collection of all contributions into a physical basis remain necessary.

## Reproduction and literature

- [Loop derivation](MIXED_SEED_POLE.md), [checker](verify_mixed_seed_pole.py), and [receipt](MIXED_SEED_RESULTS.json).
- [Independent audit](INDEPENDENT_MIXED_AUDIT.md), [checker](independent_mixed_audit.py), and [receipt](INDEPENDENT_MIXED_RESULTS.json).
- [Source identity](MIXED_SOURCE_IDENTITY.md), [checker](verify_mixed_source_identity.py), and [receipt](MIXED_SOURCE_IDENTITY_RESULTS.json).
- [Bounded literature comparison](LITERATURE_MIXED_YUKAWA.md) and [provenance](LITERATURE_PROVENANCE.json). Two paper searches screened eight distinct records. A broader web search was credit-blocked; the sources do not provide this retained-sterile coefficient or establish novelty.

The three standalone programs use explicit failure guards and no project imports. The [full replay receipt](../../verification/MIXED_YUKAWA_REPLAY_LINUX.json) records the extended check. These are internal AI-assisted calculations, not independent human peer review.

The next priority is to collect the three matter poles and earlier three-point source counterterms into a consistently normalized retained-sterile operator basis, retaining lower-order insertions and external sources. Only after that closure can one assess physical running, finite threshold matching and a declared observable. No new fundamental law or unification is claimed.

Research PR #2 is updated before the source-pinned website PR #58; both remain open for review. No merge, website deployment, GitHub release, new Zenodo version or DOI is made. The existing DOI `10.5281/zenodo.22319172` continues to identify the unchanged N00AK-r1 baseline.
