# Hypercharge-leg seed checkpoint — 15 October 2026

**Ricardo Maldonado · provisional calculation prepared with AI assistance**

The October 1 checkpoint showed that the zero-gauge three-point result cannot determine every transverse gauge contact. This checkpoint computes the two internal hypercharge attachments of that declared derivative seed, at exactly one seed insertion, one conjugate sterile Yukawa and first order in the hypercharge coupling. It resolves a transverse part of this restricted off-shell four-point UV pole.

The [two-gauge-leg continuation](two_gauge_leg/00_READ_FIRST.md) now computes all seven two-B loop diagrams at second hypercharge order and matches the action vertex predicted by the [literal action/sign audit](operator_map/ACTION_SIGN.md). The one-B result and its historical derivation remain below.

It does not supply a complete dimension-seven anomalous-dimension matrix, a physical beta coefficient, a finite threshold, an observable or a unified theory. The previously withdrawn coefficient remains withdrawn.

## Computed one-B object

Let p enter on the sterile leg, with r, h and k leaving on the lepton, scalar and hypercharge legs; p=r+h+k. Retain the public positive-momentum seed convention, Minkowski metric (+---), signed scalar mass m_H² and a massless internal lepton. The separate charged-line residues are combined in the Ward normalization R=R_H-R_L, with the charge factor stripped. The formula below defines a reduced integral tensor. The [action/source continuation](operator_map/00_READ_FIRST.md) now fixes its full pole/counterterm sign in a literal action, correcting the earlier unconjugated internal charge labels.

Let T(p,h) denote the previously computed bubble polynomial with lepton momentum p-h. Then k_nu R^nu=T(p,h+k)-T(p,h). Gauging the existing local representative gives, in this same normalization,

    G_J^nu = [-(2r+k)^nu/3 + (2h+k)^nu/6] slash(p)
             + (p²/6-m_H²/2) gamma^nu.

The new four-point residue is R^nu=G_J^nu+D^nu, where

    D^nu = [slash(k),gamma^nu] slash(r)/4
           + [k^nu slash(k)-k² gamma^nu]/3.

Both terms of D are exactly transverse. The scalar mass cancels from D. This is a coefficient statement in a declared local tensor representation of two loop attachments, not a physical running coefficient. The commutator expression is retained without a four-dimensional dual-epsilon/gamma5 replacement.

For a free massless outgoing lepton and a real transverse gauge wave, its matrix element reduces to

    ubar(r) D(epsilon) = ubar(r)[(r.epsilon) slash(k)-(r.k) slash(epsilon)].

This can be nonzero. Agreement with the Ward identity therefore does not imply that G_J alone matches the loop, or that the extra information disappears on shell.

## Computed two-B object

For `p=r+h+k1+k2`, the seven diagrams comprise two fermion orderings, two scalar orderings, two mixed assignments and one scalar seagull. The ordered propagators and literal vertices fix their reduced weights as `FF12+FF21+SS12+SS21-FS12-FS21-2SG`. Their evaluated sum is

    R_BB = (epsilon1.epsilon2) slash(p)
           - ([slash(k1),slash(epsilon1)] slash(epsilon2)
              + [slash(k2),slash(epsilon2)] slash(epsilon1))/4.

With `F=C_(7,seed)Y_nu*` and `y_L=-1/2`, the full loop pole is `-i F (g' y_L)^2 R_BB P_L/(16 pi^2 epsilon)`. It is cancelled by the positive counterterm action already fixed in the action/sign audit. The two-B residue is independent of the scalar mass; its degree-one UV polynomial cannot receive the mass terms from the propagator expansion. No external EOM or on-shell restriction is imposed in this comparison.

Read [TWO_B_DERIVATION.md](two_gauge_leg/TWO_B_DERIVATION.md), [INDEPENDENT_TWO_B_AUDIT.md](two_gauge_leg/INDEPENDENT_TWO_B_AUDIT.md) and [VERTEX_WARD_AUDIT.md](two_gauge_leg/VERTEX_WARD_AUDIT.md) for the integration routes, scalar-seagull sign and Ward controls. The [locality and source-method note](two_gauge_leg/SOURCE_METHOD.md) establishes uniqueness only for the specified local degree-one polynomial with fixed lower-point Ward kernels. Its three-B vanishing statement is conditional on the Ward hierarchy and degree bound; no three-B loop integration has been performed.

## Evidence and remaining work

Read [GAUGE_LEG_DERIVATION.md](GAUGE_LEG_DERIVATION.md) and the separately authored [INDEPENDENT_GAUGE_AUDIT.md](INDEPENDENT_GAUGE_AUDIT.md) for routing, pole reduction and exact controls. Validation compares parameter/denominator reduction with a separately derived large-momentum angular expansion. Deliberately altered coefficients, omitted attachments and omitted transverse terms must fail their stated checks. These are internal AI-assisted mathematical checks; no independent human peer review or empirical validation is claimed.

The source audit identifies the primary papers and records its bounded successful search coverage. No private raw research files, personal documents, historical transcripts or downloaded paper full texts are published here. Earlier source-calibration rounds are not relabeled as new results.

The [action/source continuation](operator_map/00_READ_FIRST.md) maps the one-B tensors into a declared five-direction off-shell space and derives their sourced descendants. The new two-B integration checks that chosen lift at the next gauge order. The next priority is to match the sourced descendants and correlated counterterms consistently in a complete retained-sterile basis, then derive physical running. Other coupling sectors and finite thresholds remain open. The source ledgers record bounded literature comparisons. Author review and the remaining basis/subtraction work precede any decision to archive a new Zenodo version. The baseline DOI is unchanged.

## Reproduce

From the repository root, use the pinned Python 3.12 environment:

```sh
python -B verify.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python -B verify.py --full --output /workspace/toe-validation/two-gauge-full-replay.json
```

The root runner checks the public payload and uses disposable copies. The original one-B standard-library stages are `verify_gauge_leg.py` (318 matrix equalities, three scalar equalities and six wrong-formula controls) and `verify_gauge_asymptotic.py` (24 generic kinematic choices, 1,536 matrix entries, omitted-attachment/transverse controls and four Ward-blind chain reversals). All are exact arithmetic. A Ward-only check accepts the deliberately reversed open chain, demonstrating why the full tensor check is necessary.

The two-B continuation adds three standard-library stages: `two_gauge_leg/verify_two_b.py`, `two_gauge_leg/verify_two_b_angular.py` and `two_gauge_leg/verify_two_b_vertices.py`. They check the actual loop residues, an independently authored angular extraction, and the literal vertices and Ward identities. The root runner contains eighteen quick stages and twenty full stages. Its completed full replay is recorded in `verification/TWO_GAUGE_REPLAY_LINUX.json`; the two-B directory separately retains the focused check outputs.

With the existing optional pinned SymPy/mpmath review dependencies, run `python -B checkpoint_2026_10_15/tensor_engine.py --self-test` for 73 general integral controls. See [TENSOR_ENGINE.md](TENSOR_ENGINE.md) and [LITERATURE_SCOPE.md](LITERATURE_SCOPE.md) for conventions and source coverage.

A further independent symbolic route is `python -B checkpoint_2026_10_15/gauge_critic.py`, documented in [GAUGE_CRITIC.md](GAUGE_CRITIC.md). It shifts the full numerators and integrates simplex moments before comparison. It requires the same pinned SymPy review environment and explicitly rejects optimized execution, preserving its scientific assertions.
