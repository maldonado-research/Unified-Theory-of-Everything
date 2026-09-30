# Verification scope

## Public revision

The included results in `code/results/` record a fresh replay under Python 3.12.14 and NumPy 2.3.5. Both normal and optimized runs pass 700 complex tree cases and 500 conditional matrix cases. The optimized run records exact agreement of deterministic results and input hashes with the normal run. Runtime fields, including completion time and optimization level, are intentionally separate.

The public runner checks four frozen tree-kernel errors, equal-mass handling, exact zero-residue deletion, excluded poles and the bound domain. Conditional tests cover local orders 0 through 5 and a 25-term comparison to the rational ansatz below the removed poles. Common coefficients 0, 1, 3, 7 and 1+2i all satisfy the assumed paired relation. Four hundred mismatched-contact controls reject changing only one partner. No physical loop coefficient is established by these tests.

Use `code/README.md` for execution. Across runtimes, acceptance is determined by the stated numerical tolerances; last digits need not agree with archived results. The checker raises exceptions and exits nonzero on a failed condition, including in optimized mode.

## Original local release audit

The original COMPLETE and PHONE archive hashes match the supplied release ledger. Both original archive validators pass. Rebuilding the unchanged original payload reproduces both original ZIPs byte-for-byte.

Separate fresh execution of the original core (700 tree plus 500 conditional cases), cleanroom (96 trials), and structural checker (240 matrix trials) passes in normal and optimized modes, with exact parity within each pair in the tested runtime. The exact-rational finite-pole audit for one through six nodes passes. Five deliberately wrong scientific constructions are rejected and one correct recurrence control passes. Fifteen malformed-container cases are rejected and the canonical control passes.

Regeneration of the original numerical outputs changes 8 of 11 generated artifact bytes on the newer runtime, including small numerical differences and propagated hashes; all fresh numerical thresholds pass. The nested-parent external replay was not executed. Neither complete historical generated-byte replay nor independent component-loop reproduction is claimed.

These original audit observations are provenance, not a claim that the focused public package reruns every historical container or parent workflow. Historical notes, raw receipts and nested ZIPs are omitted from this public package to keep the current scientific status explicit and exclude conversation exports.
