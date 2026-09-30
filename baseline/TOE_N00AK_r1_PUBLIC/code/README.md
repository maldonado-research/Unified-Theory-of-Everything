# N00AK-r1 numerical replay

This focused package checks exact tree-resolvent identities and clearly labeled conditional matrix constructions. It does not compute a loop integral, validate the proposed dimension-seven pole, or establish a physical beta function. Read the release scope correction before interpreting any conditional output.

Requires Python 3.10 or newer and NumPy. A fresh environment may be used:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python verify_public.py --output results/replay_normal.json
python -O verify_public.py --output results/replay_optimized.json --compare results/replay_normal.json
```

All checks raise explicit exceptions and remain active under `-O`; a failed check exits nonzero. The second command compares deterministic results and source hashes exactly while keeping timestamps and optimization mode separate. Across Python/NumPy versions tiny numerical changes can occur; the explicit tolerances are the portable acceptance criterion. Included results were produced with the Python and NumPy versions recorded in each result.

The deterministic replay uses seed 20260904 and includes:

- 700 complex tree cases: independent spectral-sum evaluation, geometric remainders, spectral-norm bounds, known-node reconstruction, moment recurrence, and complex symmetry. Additional controls cover equal masses, exactly canceled combined residues, and excluded poles.
- 500 conditional matrix cases: local orders 0–5, finite-ratio rational algebra, 25-term convergence under a strict mass hierarchy, and common coefficients 0, 1, 3, 7, and 1+2i. Four hundred mismatched-contact controls demonstrate that changing only one partner breaks the assumed relation.
- Four frozen tree error cases from the hash-locked Yukawa CSV, plus a same-input coefficient diagnostic. Cancellation occurs with many common coefficients; it cannot determine the value three.

The conditional functions take a coefficient as an input and compute the two partners independently from their formulas. They are **not** labeled or reported as physical beta functions. Their nonzero numerical norms quantify the supplied matrices only.

For the tree formulas, masses are positive Takagi masses, Yukawa arrays use one column per mass, and the spectral sum is evaluated away from its poles. Geometric-series and remainder-bound claims require |z| below the lowest mass squared. Known-node reconstruction assumes the pole locations are supplied; it does not infer unknown masses. Exact equal-mass grouping and exact zero-residue deletion use no arbitrary closeness tolerance. Random tests support reproducibility; the analytic proofs and their scope are given in the release note.

`PROVENANCE.json` identifies original research files and hashes, the unchanged frozen input, and the adaptation boundary. No parent archives, private chat exports, or superseded theorem certificates are required to replay this package.
