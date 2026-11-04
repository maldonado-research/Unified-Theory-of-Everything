# A common counterterm action and sector grading

4 November 2026 (Pacific). Provisional GitHub candidate.16. Input research commit `41f9b709c0fafdefdc9e1b990caa5421ee99f7dc`.

This continuation brings the earlier selected three-point counterterm, the primitive matter poles and the known pure-neutrino-Yukawa renormalization inputs into a specified classical action convention. A neutral quadratic calculation checks the combined kinetic, mass and probe terms directly. A separate grading ledger tracks canonical dimension, inverse-scale order, lepton spurions and dimensional-regularization weights.

These are action algebra and consistency calculations using previously integrated poles. No new loop, complete physical coefficient vector, finite matching or observable is calculated.

## The common sign and coordinate convention

The earlier [three-point packet](../../checkpoint_2026_10_15/interacting_reduction/INTERACTING_CANONICAL.md) was written as a **counterterm**, `+kappa P`. The later [primitive matter results](../../checkpoint_2026_10_30/real_matter/REAL_MATTER_TRANSPORT.md) were written as **loop poles**, `+kappa Q`, whose cancelling action is `-kappa Q`. They therefore enter the selected cancelling action as

```text
+kappa P-kappa Q+the declared renormalizable counterterms.
```

Let `vC` be the earlier local normalization direction for `P`, and `w` the October 30 matter normalization direction for `Q`. In the convention `old field = new field + kappa v`, the common seed-dependent substitution is

```text
v=vC-w,
P_C=P+delta_vC S0,
Q_M=Q+delta_w S0,
P-Q+delta_(vC-w) S0=P_C-Q_M.
```

Here `S0` is the declared renormalizable reference action. The ordinary renormalizable counterterms remain explicit in their bare-field convention. The [common-action calculation](COMMON_ACTION_MAP.md) gives the field, interaction, source and neutral mass terms; its [independent audit](INDEPENDENT_COMMON_AUDIT.md) checks the selected action by a separately authored route. This construction does not add the source-reduced and locally normalized tables, which belong to different field coordinates.

The matter map includes a nonlinear scalar change. On a constant neutral scalar background it changes the background entering the Dirac mass block as well as the fermion fields. Keeping only the fermion Jacobian misses that contribution. The independently held probe sources change with the fields too.

## Ordinary field factors still act on the tree seed

An ordinary Yukawa field counterterm is first order in the loop expansion and independent of `C7`. Acting on the tree derivative seed, it consequently gives a retained term linear in `C7`. With the earlier conventions,

```text
H_L=Y_nu Y_nu^dagger, T_Y=Tr(Y_nu^dagger Y_nu),
r7=(H_L C7+C7 H_L^T)/4+T_Y C7,
Delta C7_fields=-kappa r7.
```

The separately specified external-field part of the parameter counterterm is `+kappa r7`, counted once. Any further inverse ordinary field normalization must also transport the tree seed. This is distinct from applying a seed-dependent one-loop shift to an already one-loop pole, or to the tree seed itself at quadratic order in `C7`. The known field-compensation term does not determine the total physical seed counterterm.

If an independent dimension-five interaction is included, the full joint substitution also produces a retained-sterile term from `vC`, in addition to the matter-map contribution. These terms have different canonical dimensions. Their inclusion requires an explicit choice of higher-operator truncation; canonical dimension seven alone does not establish order `Lambda^-3` at fixed retained masses.

## A graded ledger before RG projection

For a zero-external-gauge monomial with `n_s` scalar fields, `n_f` fermion fields and `n_D` derivatives, its four-dimensional canonical dimension and bare coefficient exponent are

```text
D4=n_s+3 n_f/2+n_D,
w_out=n_s+n_f-2,       coefficient_0=mu^(w_out epsilon)(...).
```

The [sector-grading calculation](SECTOR_GRADING.md) derives the graph identity `w_in-w_out=2 L` from vertex and internal-line counting. The selected one-loop residues all have weight difference two, although their recipient exponents, canonical dimensions and lepton spurion charges differ. That common factor does not permit adding coefficients of different operators or prove graph existence.

For complex symmetric retained masses, `M M^dagger M` has mass degree three but the same net lepton-spurion charge as `M`. Replacing it by literal `M^3` changes the charge and generally the ordered matrix product. The ledger retains that distinction and separates canonical dimension from fixed-mass inverse-scale counting.

## Reproduce and continue

| Material | Scope |
| --- | --- |
| [Common action](COMMON_ACTION_MAP.md) · [checker](verify_common_action.py) · [receipt](COMMON_ACTION_RESULTS.json) | Common sign/field convention and selected neutral kinetic, mass and source transport |
| [Independent audit](INDEPENDENT_COMMON_AUDIT.md) · [checker](independent_common_audit.py) · [receipt](INDEPENDENT_COMMON_RESULTS.json) | Separately derived action and tree-seed normalization controls |
| [Sector grading](SECTOR_GRADING.md) · [checker](verify_sector_grading.py) · [receipt](SECTOR_GRADING_RESULTS.json) | Dimension, spurion, scale-weight and truncation ledger |
| [Literature comparison](LITERATURE_COMMON_ACTION.md) · [provenance](LITERATURE_PROVENANCE.json) | Bounded searches, primary-source comparisons and retrieval limits |

Run the three standalone checkers with isolated optimized Python in the root pinned environment. The root publication runner itself runs without optimization; the [full replay receipt](../../verification/COMMON_ACTION_REPLAY_LINUX.json) records its disposable-copy execution. Checks and independent AI-assisted authorship do not constitute human peer review or empirical validation.

The common convention supplies a starting point for collecting total coefficients. The full retained-sterile derivative/EOM and evanescent basis, all required interaction sectors and external-gauge matter amplitudes, redundant-coefficient differentiation, quantum composite-source transport, finite thresholds and one observable remain open. The neutral projection does not complete the interacting physical map.

Research PR #2 precedes the source-pinned website update in PR #58. Both remain open for review; no merge, live deployment, release or new Zenodo version is made. DOI `10.5281/zenodo.22319172` continues to identify only the unchanged N00AK-r1 baseline.
