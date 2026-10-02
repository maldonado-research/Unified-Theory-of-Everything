# TOE operator-reduction checkpoint — 1 October 2026

**Ricardo Maldonado · internal mathematical audit · prepared with AI assistance**

This checkpoint advances one step beyond the September 6 routed seed calculation: it supplies a convention-fixed local operator representative and derives the interaction contacts produced by a sourced sterile-field equation. It does not complete dimension-seven matching, establish a physical beta function, or connect the separate research branches into a unified action.

## What was established under the declared assumptions

1. The six structures of the existing off-shell seed polynomial have a local, electroweak-covariant dimension-seven representative. Derivative placement and Fourier signs are explicit. These six monomials are a representative, not a complete or independent operator basis.
2. Applying the sourced lepton and Higgs equations retains curvature, Yukawa-source and scalar-potential descendants. In that representative, the scalar mass contributes a coefficient `1/2 + 1/6 = 2/3` to its derivative-sterile term. This coefficient is not physical running.
3. Reducing the cubic sterile derivative produces a Yukawa-shaped term **and three source terms**. A two-component calculation makes the associated Weinberg and derivative contacts explicit.
4. The correlated Yukawa and Weinberg shifts cancel in the static combination `C5 + Y M^-1 Y^T` for the removable sterile-derivative monomial, with arbitrary complex flavor coefficients. This is a consequence of a consistent change of variables. It cannot determine the loop coefficient or establish cancellation of the complete seed packet.
5. An explicit gauge-field-strength operator is invisible to the zero-gauge three-point restriction, obeys the gauge Ward test, and has a nonzero on-shell four-point amplitude. The three-point polynomial therefore cannot determine its full gauge completion.

Read [GREEN_REPRESENTATIVE.md](GREEN_REPRESENTATIVE.md) for the operator map and light/sterile source tower; [SOURCED_STERILE_REDUCTION.md](SOURCED_STERILE_REDUCTION.md) fixes the two-component action and correlated contact signs. [AUDIT.md](AUDIT.md) records a separate skeptical review. The exact arithmetic checks support these algebraic statements; their scope and wrong-formula controls are recorded alongside the scripts and results.

## Corrections and limits

The N00AK-r1 withdrawal remains in force. Neither `5/9`, `5/6`, nor the scalar descendant `2/3` supplies a replacement physical beta coefficient. Free on-shell projection removes terms that the sourced equations retain. A basis-change cancellation is not evidence that every diagram or physical matching contribution cancels.

The unchanged September 4 baseline and September 6 scientific sources remain separately dated. The older DOI `10.5281/zenodo.22319172` identifies N00AK-r1 only. This checkpoint introduces no observational prediction, new particle detection, novel fundamental mathematics, experimental validation, or independent human peer review.

No new literature discovery is claimed. The cited earlier primary-source links were unavailable from the current research environment; the prior dated literature ledger has not been represented as a fresh web survey. No private raw research files, personal documents or historical transcripts are included in this checkpoint.

## Reproduce

From the repository root, with its pinned Python 3.12 virtual environment:

```sh
python -B verify.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python -B verify.py --full
```

The root runner checks payload integrity and executes the calculations in a temporary copy. Keep fresh replay receipts outside the repository with `--output /workspace/toe-validation/full-replay.json`. The new exact-arithmetic scripts also run individually without third-party packages; use the root runner to preserve the shipped numerical outputs. Do not optimize the root/checkpoint replay with `-O`.

## Next decisive calculation

Declare the complete sequential EFT action, power counting and operator basis. Compute the additional gauge-leg and multi-field Green functions, reduce every source descendant consistently, and assemble mass, field, Yukawa and contact counterterms in a declared subtraction scheme. Then compare the complete physical matching combination, including finite thresholds. Report a surviving correction or cancellation with the same standard of evidence.
