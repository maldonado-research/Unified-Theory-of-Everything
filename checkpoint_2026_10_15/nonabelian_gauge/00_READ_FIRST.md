# Non-Abelian two- and three-gauge matching

15 October 2026. Provisional continuation prepared with AI assistance.
Input: candidate `5c77bc5243b32da7b97a0ab86ea8571cb1cff5c0` and the
[literal seed action](../operator_map/ACTION_SIGN.md).

This checkpoint fills the missing two-W/mixed-BW comparison and explicitly
integrates the next external-gauge multiplicity in the selected
one-loop seed/conjugate-Yukawa 1PI matter sector.

The [two-gauge calculation](TWO_WEAK_DERIVATION.md) retains the ordered
weak generators on every one of the seven matter-loop topologies. Its
result is

    R(X1,X2) = {X1,X2} R2/2 + [X1,X2] C2/2.

`R2` is the existing Abelian polynomial; `C2` is the evaluated ordered
combination invisible in the Abelian sum. Mixed BW has commuting
generators and reduces to `y_L T^I R2`, with no extra factor two.
A [separate angular integration](INDEPENDENT_NONABELIAN_AUDIT.md) agrees
on all coefficients of the declared local polynomial.

The [three-gauge calculation](THREE_GAUGE_DERIVATION.md) evaluates all
33 required matter-loop graphs: 24 ordinary-insertion diagrams and nine
with a scalar seagull. Its noncommuting cubic residue matches the
covariant action. A separately written
[24-cell cubature integration](INDEPENDENT_CUBATURE_AUDIT.md) agrees on
2,112 individual diagram evaluations and all 64 complete polarization
tensors, retaining six formal noncommuting color words. The adjacent
comparison harness reproduces that agreement.

An [independent action/Fourier audit](COVARIANT_ACTION_AND_WARD.md)
differentiates the full matrix curvature and its adjoint derivative.
For `D=partial+i A`, the required curvature is

    Fcal_mu nu=partial_mu A_nu-partial_nu A_mu+i[A_mu,A_nu].

With Hermitian `T^I=tau^I/2`, its weak component has the nonlinear
`-g epsilon^IJK W_mu^J W_nu^K` sign. The two-/three-gauge Ward identities
include lower-point commutator terms. A triple-gauge tree attached through
one gauge propagator is reducible and is excluded from these 1PI kernels.

The [hierarchy and photon projection](GAUGE_HIERARCHY_AND_PHOTON_PROJECTION.md)
make two consequences precise. Primitive diagrams with four or more
external gauge fields are UV finite by degree `3-n`, at the stated
lowest gauge-coupling order. Thus the selected external-gauge hierarchy
is now checked through its final potentially divergent multiplicity.
Also, a neutral `A W^+ W^-` local UV contact can remain nonzero even
though the selected one-photon pole cancels. Its full charged-basis
coefficient is calculated and independently checked; it is not a finite
photon amplitude or an observable prediction.

Run the six verification stages with Python's standard library:

```sh
python3 -I -S -B -O verify_two_weak.py --json
python3 -I -S -B -O verify_two_gauge_angular.py
python3 -I -S -B -O verify_covariant_action.py
python3 -I -S -B -O verify_three_gauge.py
python3 -I -S -B -O verify_three_gauge_cubature.py
python3 -I -S -B -O compare_three_gauge_routes.py
```

The two-gauge producer reports 7,168 controls, explicitly separating
2,304 inherited Lorentz/Abelian reproductions, 4,768 other positive
controls and 96 mutation rejections. The separate angular route reports
4,285 controls. The action audit reports 642 controls including 56
rejected mistakes. The three-gauge producer reports 5,805 checks;
the cubature route reports 10,407 low-level and coefficient controls.
The cross-route comparison separately records 2,112 individual and
64 full-tensor agreements. These counts are implementation controls,
not counts of independent physical predictions or measurements.

The root runner executes all six stages on disposable copies, preserving
the earlier scientific packages and receipts. Normal and optimized
standalone failures use explicit exceptions. The root/numerical runner
still requires normal Python execution.

The [primary-method note](METHOD_NOTE.md) checks the background Ward
functional, 1PI-versus-connected distinction and convention translation.
Two bounded 2025–2026 searches screened six distinct records without
finding another directly usable benchmark. Full downloaded papers and
private references are excluded; safe retrieval and search metadata
are retained in the adjacent JSON files.

A complete retained-sterile operator basis, other external-matter
amplitudes, full interacting source reduction/canonical normalization,
lower-order insertions, correlated counterterms, a subtraction/evanescent
prescription, physical running and finite matching remain open. This
checkpoint does not establish new fundamental mathematics, external
novelty, unification or an experimentally tested prediction. No new
Zenodo version is assigned; its earlier DOI retains the baseline scope.
