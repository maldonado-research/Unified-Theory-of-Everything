# TOE-GitHub-2026.09.30

This first GitHub edition makes Ricardo Maldonado's Unified Theory of Everything research program inspectable and reproducible. It preserves the published TOE-N00AK-r1 baseline and adds curated scientific material from the 6 September 2026 checkpoint.

## Included research

- Exact selected Majorana tree-kernel identities, approximation bounds and known-pole reconstruction.
- The published dimension-seven scope correction and a later evaluated off-shell seed-loop polynomial.
- Conditional scalar continuum and positive-pole comparisons.
- A prescribed mass-pulse example with a rigorously positive weak-pulse occupation and separate numerical solver checks.

## Reproducibility

The publication copy passed all ten runner stages with Python 3.12.14, NumPy 2.3.5 and SciPy 1.18.1. This includes baseline hashes, normal/optimized baseline comparison, 815 exact routed controls, weak-pulse and scalar checks, the six-case pulse benchmark and a separate rotated-basis Radau review. The checkpoint checking path now keeps scientific assertions enabled and rejects optimized execution that would remove them.

Run `python verify.py` for focused checks or `python verify.py --full` for the numerical pulse comparisons. The repository includes tested dependencies, a public manifest, checksums, source provenance, research-status notes and citation metadata.

## Scope and citation

This is a research and reproducibility release. It does not establish complete unification, empirical confirmation, an independently human-reviewed result or the completed physical dimension-seven beta function. The scalar premises remain inherited, and the pulse background remains prescribed. The five-ninths projection is not a replacement rule for the old physical coefficient.

The earlier DOI https://doi.org/10.5281/zenodo.22319172 identifies only the unchanged N00AK-r1 baseline. Cite this larger GitHub edition by its release tag or commit. Author-owned materials preserve CC BY 4.0; third-party dependencies retain their own licenses. AI assistance is disclosed.

The public adaptation removes private inventories, cloud identifiers, account paths and unrelated-project provenance. Original scientific files and earlier public records were preserved.
