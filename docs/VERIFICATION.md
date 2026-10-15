# Reproducing the focused research package

## Environment

The publication checks use Python 3.12.14, NumPy 2.3.5 and SciPy 1.18.1. Install the root `requirements.txt` in a fresh virtual environment. The checkpoint's original dependency pins are retained as historical source metadata; they are not the environment used for this GitHub replay. Matplotlib is only needed to regenerate plots and is not required by the focused verification commands.

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

The full mode adds:

14. Six pulse cases with midpoint refinement and DOP853 comparison.
15. A separate Hadamard-rotated Radau review, sign and zero-momentum controls, and weak-pulse comparison.

The baseline input-manifest check is an additional stage. The wrapper reports each stage and exits nonzero on any failure. Temporary outputs are not written over shipped receipts.

Counting the root payload check gives fifteen quick stages and seventeen full stages. `verification/PUBLICATION_REPLAY.json` preserves the earlier September replay; `verification/OPERATOR_REPLAY_LINUX.json` records the October 1 edition. The root and numerical checkpoint commands must run without `-O`. The scientific derivations and limits are documented in `checkpoint_2026_10_01/` and `checkpoint_2026_10_15/`; case counts are implementation controls, not proof by sampling or proof of a complete operator basis.

An additional independent symbolic review is reproducible with `python -m pip install -r checkpoint_2026_10_01/requirements-review.txt` followed by `python -B checkpoint_2026_10_01/independent_clifford_checks.py`. Its four exact control groups include a noncommuting gauge connection and a physical-polarization trace. This optional review uses pinned SymPy/mpmath and is separate from the root runner's seventeen full stages; the cloud setup includes its dependencies.

The same pinned symbolic environment supports the optional general tensor-residue audit: `python -B checkpoint_2026_10_15/tensor_engine.py --self-test` (73 exact checks, including 18 independent angular-expansion comparisons and three wrong-formula controls). This reusable UV-pole helper does not calculate finite threshold terms or distinguish the full integral's infrared poles. It is separate from the root runner's seventeen full stages.

## Publication correction to replay behavior

The original checkpoint wrapper invoked every script with `-O`. In `fermion_portal/independent_review.py`, that would remove its scientific assertions. The publication wrapper runs the checkpoint scripts normally, and both it and the independent review reject optimized execution. Formula code, test values and tolerances are otherwise unchanged.

This is an implementation correction to the checking path, not a new physical result. The unchanged baseline's optimized comparison remains valid because that baseline uses explicit guards.

## What passing means

A passing check supports the named algebra, numerical tolerances and toy dynamics under their declared premises. Floating-point agreement is not a certified enclosure of every error. The scalar parent theorem, full loop matching, source energy budget, unseen experimental data and complete unification are not verified by these scripts.

The time-tail bound, numerical refinement and finite momentum interval address different omissions. They must not be merged into one universal accuracy certificate. Results can vary slightly across compatible platforms; explicit scientific tolerances control acceptance.

For a third gauge-leg derivation, run `python -B checkpoint_2026_10_15/gauge_critic.py` in the pinned symbolic environment. It compares full shifted-numerator simplex integrals to the reduced residues, then checks source anchors, transverse terms, on-shell rows and deliberate omissions. Run without `-O`; the script rejects optimization because its scientific assertions must remain enabled. This symbolic check is optional and separate from the seventeen-stage full root replay.

The action/source continuation adds an optional independent audit using the same pinned SymPy environment: `python -B checkpoint_2026_10_15/operator_map/verify_operator_critic.py`. Its 45 exact controls independently check the action charges, Fourier lift, dual sign, 19-column independence and Ward rank, total currents, and flavor spectral/null-block statements. The original hypercharge receipt remains unchanged; `verification/OPERATOR_MAP_REPLAY_LINUX.json` records the extended seventeen-stage replay.
