# One-hypercharge-leg seed checkpoint — 15 October 2026

**Ricardo Maldonado · provisional calculation prepared with AI assistance**

The October 1 checkpoint showed that the zero-gauge three-point result cannot determine every transverse gauge contact. This checkpoint computes the two internal hypercharge attachments of that declared derivative seed, at exactly one seed insertion, one conjugate sterile Yukawa and first order in the hypercharge coupling. It resolves a transverse part of this restricted off-shell four-point UV pole.

It does not supply a complete dimension-seven anomalous-dimension matrix, a physical beta coefficient, a finite threshold, an observable or a unified theory. The previously withdrawn coefficient remains withdrawn.

## Computed object

Let p enter on the sterile leg, with r, h and k leaving on the lepton, scalar and hypercharge legs; p=r+h+k. Retain the public positive-momentum seed convention, Minkowski metric (+---), signed scalar mass m_H² and a massless internal lepton. The separate charged-line residues are combined in the Ward normalization R=R_H-R_L, with the charge factor stripped. Full Feynman amplitude and counterterm signs require the separate coupling/vertex map; the formula below defines a reduced integral tensor.

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

## Evidence and remaining work

Read [GAUGE_LEG_DERIVATION.md](GAUGE_LEG_DERIVATION.md) and the separately authored [INDEPENDENT_GAUGE_AUDIT.md](INDEPENDENT_GAUGE_AUDIT.md) for routing, pole reduction and exact controls. Validation compares parameter/denominator reduction with a separately derived large-momentum angular expansion. Deliberately altered coefficients, omitted attachments and omitted transverse terms must fail their stated checks. These are internal AI-assisted mathematical checks; no independent human peer review or empirical validation is claimed.

The source audit identifies the primary papers and makes the fresh-access limitation explicit. No private raw research files, personal documents, historical transcripts or downloaded paper full texts are published here. Earlier source-calibration rounds are not relabeled as new results.

Next: map the transverse tensors into a declared independent retained-sterile Green basis, retain their sourced EOM descendants, determine all other counterterms at the chosen accuracy, and match the complete invariant combination with finite thresholds. Broader coupling sectors and a fresh literature comparison remain open. Author review and the remaining convention/basis work precede any decision to archive a new Zenodo version. The baseline DOI is unchanged.

## Reproduce

From the repository root, use the pinned Python 3.12 environment:

```sh
python -B verify.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python -B verify.py --full --output /workspace/toe-validation/hypercharge-full-replay.json
```

The root runner checks the public payload and uses disposable copies. The two new standard-library stages are `verify_gauge_leg.py` (318 matrix equalities, three scalar equalities and six wrong-formula controls) and `verify_gauge_asymptotic.py` (24 generic kinematic choices, 1,536 matrix entries, omitted-attachment/transverse controls and four Ward-blind chain reversals). All are exact arithmetic. A Ward-only check accepts the deliberately reversed open chain, demonstrating why the full tensor check is necessary.

With the existing optional pinned SymPy/mpmath review dependencies, run `python -B checkpoint_2026_10_15/tensor_engine.py --self-test` for 73 general integral controls. See [TENSOR_ENGINE.md](TENSOR_ENGINE.md) and [LITERATURE_SCOPE.md](LITERATURE_SCOPE.md) for conventions and source coverage.

A further independent symbolic route is `python -B checkpoint_2026_10_15/gauge_critic.py`, documented in [GAUGE_CRITIC.md](GAUGE_CRITIC.md). It shifts the full numerators and integrates simplex moments before comparison. It requires the same pinned SymPy review environment and explicitly rejects optimized execution, preserving its scientific assertions.
