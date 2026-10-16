# Reproducing the focused research package

## Environment

The publication checks use Python 3.12.14, NumPy 2.3.5, SciPy 1.18.1, SymPy 1.14.0 and mpmath 1.3.0. Install the root `requirements.txt` in a fresh virtual environment. The checkpoint's original dependency pins are retained as historical source metadata; they are not the environment used for this GitHub replay. Matplotlib is only needed to regenerate plots and is not required by the focused verification commands.

## Commands

```sh
python verify.py
python verify.py --full
```

The runner verifies the repository payload manifest, copies the scientific package into a temporary directory, verifies the baseline's original manifest, and executes:

1. N00AK-r1 tree and conditional checks in normal Python: 700 tree cases, 500 conditional cases and frozen controls.
2. N00AK-r1 optimized replay and deterministic semantic comparison. Its checks use explicit exceptions and stay active under optimization.
3. Routed ultraviolet seed: 815 exact rational controls.
4. Separate tensor-routing comparison.
5. Exact rational nonzero weak-pulse bound.
6. Scalar continuum and pole-bound controls.
7. Local Green representative, sourced sterile tower and transverse gauge-ambiguity controls: 605 exact checks.
8. Sourced sterile contact and complex-flavor Schur-invariance controls: 640 exact checks.
9. One-hypercharge-leg seed residues, Ward contraction and nonzero transverse remainder.
10. Separately authored large-momentum angular reduction of both gauge attachments, including altered-formula controls.
11. Literal-action phases, charge conjugation, seed normalization and counterterm signs: 891 controls.
12. Local transverse operators and sourced contacts: 158 checks and 13 deliberate-error controls.
13. Exact polynomial Ward-kernel classification: 19 coefficients, rank 14 and transverse dimension five.
14. Seven two-hypercharge-leg diagrams from exact Feynman-simplex moments: 2232 checks, including every local polynomial coefficient and 60 mutation rejections.
15. Independently authored angular extraction of those seven diagrams: 2483 checks, both Ward descents, Bose symmetry and six mutation classes.
16. Exact unintegrated two-B vertices, charge-conjugated chains, phases and Ward shift identities: 338 positive checks and eight mutation classes at 24 rational points.
17. One-W weak-index contractions and actual crossed-triangle integrations: 695 checks at 24 generic rational points.
18. Independent one-W angular integration and exhaustive local-polynomial coefficient controls: 1461 checks, including neutral-photon, Z and charged projections.
19. Combined hypercharge source packet: 142 exact identities, including eight full off-shell residual trials with independent sources, and 20 rejected mutations. This stage requires SymPy.
20. Indexed field substitutions realizing that selected residual certificate: 192 exact checks and 14 rejected mistakes, including integrated adjoints, source signs, Majorana reality, six direct sterile-residual Fourier reconstructions and constant-background checks.
21. Separately implemented action audit using exterior Grassmann algebra: 33 exact checks and 21 rejected mistakes. This audits action-map lemmas rather than recomputing the complete prior certificate.
22. Actual neutral quadratic projection and canonical mass/source transport: 82 exact checks and 24 rejected mistakes, with doubled/squared pole kernels, external response and rank-safe Schur domains.
23. Explicit dimension-six normalization × dimension-five sterile-scalar benchmark: 124 controls including 21 rejected omissions/adjoint errors. Its optional input coefficients are not generated seed-loop results.
24. Separate component Grassmann closure review: seven controls, including the Majorana factor, a nonzero contact with singular mass, an incorrect-dagger rejection and local spectral normalization.
25. Ordered two-W/mixed-BW residues from seven matter diagrams: 7,168 controls, separating 2,304 inherited reproductions, 4,768 other positive controls and 96 mutation rejections; 256 exhaustive coefficients and 12 generic points.
26. Separate two-gauge angular integration: 4,285 controls covering literal weak indices, all generator pairs, exhaustive coefficients and both non-Abelian Ward contractions.
27. Independently differentiated matrix action and two-/three-gauge Ward recursion: 642 controls including 56 rejected nonlinear-sign/order/omission errors.
28. Actual 33-diagram three-gauge UV integration: 5,805 total controls, including seven mutation rejections, all 64 polarization coefficients and a direct neutral charged-weak projection.
29. Independent exact 24-cell cubature integration: 10,407 low-level and coefficient controls, literal color contractions, six free words and deliberate sign/order omissions.
30. Reproducible cross-route comparison of all 2,112 individual graph/basis values and 64 complete color/Clifford tensors. The comparison imports both adjacent core routes; each core integration remains standalone.
31. Matrix electroweak source reduction and indexed field changes: 1,410 positive exact controls and 20 rejected mistakes, including covariant Bianchi/IBP, gauge-action variations, complete selected residuals, adjoint weak shifts and source ordering.
32. Independent non-Abelian action/Grassmann audit: 122 exact controls and 81 negative controls, including a surviving one-flavor charged source and combined weak/flavor exchange.
33. Ordered fermion/scalar weak-current contractions and exact neutral/charged gauge-coordinate projections: 61 positive controls and 20 rejected mistakes.
34. Local differential normalization of the actual six-placement packet: 85 positive controls and 33 rejected mistakes, including the explicit interacting remainder, differentiated composite adjoints and induced mass/source/current transport on three nonconstant backgrounds.

35. Active/retained-mass tensor invariant spaces, chiral dual relation, Grassmann ranks and exact source maps: 166 positive checks and eight rejected mistakes.
36. Independent tensor invariant/nullspace and full polynomial quotient: 80 exact controls and 25 negative controls, including unrestricted source-image and neutral/charged projections.
37. Current-sector full weak/Lorentz invariant spaces, identical-lepton flavor Fierz quotient and spectator maps: 427 positive checks and 112 rejected mistakes.
38. Independent current exterior expansion and constructive nonsingular minors: 23,544 positive comparisons and 66 rejected mistakes. Of those positive comparisons, 23,328 are entries of explicit rank matrices, not separate physical calculations.
39. Symmetric-seed image, supported symmetric-channel reconstruction and sharp relative source bound: 908 positive controls and ten rejected mistakes, with 45 exact singular/degenerate spectra and complex-unitary rotations.

40. Explicit all-field source-to-normalizing bridge and actual off-shell neutral action substitution: 55 positive checks and 31 rejected mistakes, including the scalar-potential Hessian and probe-source transport.
41. Independent neutral Grassmann action comparison: 59 exact checks and 18 rejected mistakes, including non-diagonal real masses and the stationary scalar control.
42. Scalar gauge-current invariant spaces, IBP/complete fermion divergence, interaction-source reduction and additional local field map: 75 positive exact controls and 12 rejected mistakes.
43. Independent scalar invariant/component, sourced-divergence and Yukawa/Majorana audit: 76 positive comparisons and 20 rejected mistakes.
44. Charged-lepton/quark spectator Yukawa transport, actual scalar/lepton substitutions, source contraction and bridge cancellation: 54 positive checks and 33 rejected mistakes.

The full mode adds:

45. Six pulse cases with midpoint refinement and DOP853 comparison.
46. A separate Hadamard-rotated Radau review, sign and zero-momentum controls, and weak-pulse comparison.

The baseline input-manifest check is an additional stage. The wrapper reports each stage and exits nonzero on any failure. Temporary outputs are not written over shipped receipts.

Counting the root payload check gives forty-six quick stages and forty-eight full stages. `verification/PUBLICATION_REPLAY.json` preserves the earlier September replay; `verification/OPERATOR_REPLAY_LINUX.json` records the October 1 edition. The root and numerical checkpoint commands must run without `-O`. The scientific derivations and limits are documented in `checkpoint_2026_10_01/`, `checkpoint_2026_10_15/` and `checkpoint_2026_10_16/`; case counts are implementation controls, not proof by sampling or proof of a complete operator basis. The two-B implementations also exhaust the 256 coefficients of their declared degree-one local momentum polynomial, bilinear in polarizations; this completeness statement relies on the proved UV degree bound.

An additional independent symbolic review is reproducible in the root pinned environment with `python -B checkpoint_2026_10_01/independent_clifford_checks.py`. Its four exact control groups include a noncommuting gauge connection and a physical-polarization trace. This optional review uses pinned SymPy/mpmath and is separate from the root runner's forty-eight full stages. The same SymPy/mpmath pins are now root requirements because the selected source-packet stage requires them; the historical review requirements file remains unchanged.

The same pinned symbolic environment supports the optional general tensor-residue audit: `python -B checkpoint_2026_10_15/tensor_engine.py --self-test` (73 exact checks, including 18 independent angular-expansion comparisons and three wrong-formula controls). This reusable UV-pole helper does not calculate finite threshold terms or distinguish the full integral's infrared poles. It is separate from the root runner's forty-eight full stages.

## Publication correction to replay behavior

The original checkpoint wrapper invoked every script with `-O`. In `fermion_portal/independent_review.py`, that would remove its scientific assertions. The publication wrapper runs the checkpoint scripts normally, and both it and the independent review reject optimized execution. Formula code, test values and tolerances are otherwise unchanged.

This is an implementation correction to the checking path, not a new physical result. The unchanged baseline's optimized comparison remains valid because that baseline uses explicit guards.

## What passing means

A passing check supports the named algebra, numerical tolerances and toy dynamics under their declared premises. Floating-point agreement is not a certified enclosure of every error. The scalar parent theorem, full loop matching, source energy budget, unseen experimental data and complete unification are not verified by these scripts.

The time-tail bound, numerical refinement and finite momentum interval address different omissions. They must not be merged into one universal accuracy certificate. Results can vary slightly across compatible platforms; explicit scientific tolerances control acceptance.

For a third derivation of the one-B pole, run `python -B checkpoint_2026_10_15/gauge_critic.py` in the pinned symbolic environment. It compares full shifted-numerator simplex integrals to the reduced residues, then checks source anchors, transverse terms, on-shell rows and deliberate omissions. Run without `-O`; the script rejects optimization because its scientific assertions must remain enabled. This symbolic check is optional and separate from the forty-eight-stage full root replay.

The action/source continuation adds an optional independent audit using the same pinned SymPy environment: `python -B checkpoint_2026_10_15/operator_map/verify_operator_critic.py`. Its 45 exact controls independently check the action charges, Fourier lift, dual sign, 19-column independence and Ward rank, total currents, and flavor spectral/null-block statements. The original hypercharge receipt remains unchanged; `verification/OPERATOR_MAP_REPLAY_LINUX.json` records the extended seventeen-stage replay.

The two-B continuation is documented in `checkpoint_2026_10_15/two_gauge_leg/`. Its three standard-library scripts are included in both root modes and also run under `python3 -I -S -B -O`; they use explicit failures. The seven-diagram integration checks a previously predicted counterterm contact. The vertex audit instead checks the unintegrated Ward hierarchy including a loop-momentum shift; it is not a third integration. The new `verification/TWO_GAUGE_REPLAY_LINUX.json` records the twenty-stage replay. Earlier receipts retain their historical stage counts. The conditional three-B Ward consequence is not an integrated three-B result.

The electroweak/source packet in `checkpoint_2026_10_15/electroweak_source/` adds two standard-library stages and one SymPy stage to both root modes. Its weak scripts also run with `python3 -I -S -B -O`; its source checker uses `python3 -I -B -O` in the pinned environment. The latter has explicit failures, not removable assertions. `verification/ELECTROWEAK_SOURCE_REPLAY_LINUX.json` records the extended full replay. The one-W pole is independently integrated; the combined source reduction instead certifies an EOM equivalence of previously matched operators. Neutral-photon cancellation concerns the selected simple pole after the constant neutral-Higgs projection, not a finite amplitude or a physical dipole beta function.

The field-change continuation in `checkpoint_2026_10_15/field_redefinition/` adds five SymPy stages. Their explicit failures remain active under `python -I -B -O`; use that isolated command only for these standalone symbolic checkers, while retaining normal execution for the root runner and numerical historical checks. `verification/FIELD_REDEFINITION_REPLAY_LINUX.json` records the extended twenty-eight-stage full replay. The new calculations implement a selected first-order action change and its neutral quadratic normalization, plus a concrete optional lower-order insertion witness. They do not establish a complete physical coefficient vector, new loop diagrams, finite matching, a quantum Jacobian prescription or an observable. The scalar-tadpole background condition and Schur invertibility domains are part of the derivations.

The non-Abelian continuation in `checkpoint_2026_10_15/nonabelian_gauge/` adds six standard-library stages. All also run with `python3 -I -S -B -O`; the comparison harness imports only its two adjacent independently written core scripts. `verification/NONABELIAN_GAUGE_REPLAY_LINUX.json` records the thirty-four-stage full replay. These stages compare actual two-/three-gauge integrations with a separate Fourier action and Ward hierarchy. The degree bound closes external gauge multiplicity only for the declared single-seed/source primitive family at lowest gauge order. A nonzero neutral mixed photon–charged-weak local pole is not a finite amplitude or an observable prediction. No complete external-matter operator basis, physical RGE or evanescent/finite matching prescription is certified.

The interacting reduction in `checkpoint_2026_10_15/interacting_reduction/` adds four stages. The electroweak and differential-normalization scripts use the standard library and support `python3 -I -S -B -O`; the independent audit and current projection require pinned SymPy and support `python3 -I -B -O`. They implement the sourced matrix field changes, combined weak/flavor exchange, weak-current contractions and a specified local normalization of the actual packet, not new loop integrals. `verification/INTERACTING_REDUCTION_REPLAY_LINUX.json` records the thirty-eight-stage full replay. The scalar-derivative interaction remainder, composite-source derivatives, differential Majorana transpose and higher-gauge contacts are retained. Full basis closure, other external-matter diagrams and coupling sectors, mixed-order insertions, quantum/subtraction prescriptions, physical running and finite observables are outside these checks.

The operator-sector continuation in `checkpoint_2026_10_15/operator_basis/` adds five standalone stages. Three use the standard library with `python -I -S -B -O`; the separate tensor audit and seed-image checker require pinned SymPy with `python -I -B -O`. All use explicit failures and no project imports. `verification/OPERATOR_BASIS_REPLAY_LINUX.json` records the forty-three-stage full replay. Invariant nullspaces, full exterior expansions and constructive minors establish algebraic independence only for the stated nonderivative chiral field sectors. The coefficient tables transport explicit representative contributions, not totals after all derivative/source reductions. The flavor-image bound concerns raw selected projections, not physical moments or complete coefficient vectors. Full EFT EOM/subtraction closure, other loop sectors and physical running remain open.


The action completion in `checkpoint_2026_10_16/action_completion/` adds five stages. All support isolated optimized execution with explicit failures; the scalar-current and composition producer scripts and independent scalar audit use only the standard library; the independent composition and spectator-Yukawa scripts use pinned SymPy. The direct new coefficient calculations cover the specified scalar-current/spectator sectors and neutral quadratic action. The full first-order bridge is proved using the sourced action residual; every interacting vertex is not separately recomputed. `verification/ACTION_COMPLETION_REPLAY_LINUX.json` records the extended forty-eight-stage full replay. Earlier science and receipts retain their original dates and stage counts. The off-shell scalar-tadpole contribution is essential, and no stationary background is silently imposed. These classical action checks do not compute the transformed quantum 1PI effective action, a finite amplitude, total physical Wilson coefficients or running.
