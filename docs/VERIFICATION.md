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

The full mode adds:

7. Six pulse cases with midpoint refinement and DOP853 comparison.
8. A separate Hadamard-rotated Radau review, sign and zero-momentum controls, and weak-pulse comparison.

The baseline input-manifest check is an additional stage. The wrapper reports each stage and exits nonzero on any failure. Temporary outputs are not written over shipped receipts.

## Publication correction to replay behavior

The original checkpoint wrapper invoked every script with `-O`. In `fermion_portal/independent_review.py`, that would remove its scientific assertions. The publication wrapper runs the checkpoint scripts normally, and both it and the independent review reject optimized execution. Formula code, test values and tolerances are otherwise unchanged.

This is an implementation correction to the checking path, not a new physical result. The unchanged baseline's optimized comparison remains valid because that baseline uses explicit guards.

## What passing means

A passing check supports the named algebra, numerical tolerances and toy dynamics under their declared premises. Floating-point agreement is not a certified enclosure of every error. The scalar parent theorem, full loop matching, source energy budget, unseen experimental data and complete unification are not verified by these scripts.

The time-tail bound, numerical refinement and finite momentum interval address different omissions. They must not be merged into one universal accuracy certificate. Results can vary slightly across compatible platforms; explicit scientific tolerances control acceptance.
