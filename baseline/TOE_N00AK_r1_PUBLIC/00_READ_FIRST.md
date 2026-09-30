# TOE-N00AK-r1 public research release

Prepared 2026-09-04. Author: Ricardo Maldonado.

Start with `TOE_N00AK_r1_RELEASE_NOTE.pdf` for the overview and `TOE_N00AK_r1_TECHNICAL_NOTE.md` for the mathematical statements and derivation audit.

## Current claim status

- **Established mathematical identities:** finite-pole tree kernel, local expansion and exact remainder, known-node reconstruction, and tree approximation errors.
- **Conditional algebra only:** dimension-seven, higher-odd and finite-ratio Schur cancellation from the proposed correlated counterterms. Neither the proposed loop coefficient 3 nor a complete dimension-seven quantum-field-theory result is established.
- **Open:** a fully routed off-shell dimension-seven loop calculation, the complete operator/EOM packet, hard and finite matching, and validated physical observables.

This correction supersedes the interpretation of the local August 17 N00AK theorem, technical report, summary, and audit receipts. It does not retract the valid tree identities. The local N00AK predecessor remains a local source release. The preceding published record for this version series is 10.5281/zenodo.20241403, version `v2026.05.16-math-audit`; its title and version were verified on Zenodo on 4 September 2026. The all-versions DOI is 10.5281/zenodo.17872787.

## Reproduce

Use Python 3.10 or later and NumPy. See `code/README.md` for the exact command, tested versions, and result semantics. The code and benchmark data are self-contained; no historical ZIP or parent replay is required. A PASS certifies only the named checks and does not establish a loop coefficient or experimental discovery.

`SHA256SUMS.txt` and `PUBLIC_MANIFEST.json` describe the files in this package. The checksum file excludes itself and the manifest; the manifest lists each other payload file and its hash.

## Provenance and privacy

This package is a focused revision derived from the August 17 TOE-N00AK source release. The original input archive hashes and benchmark-data hash are recorded in `PROVENANCE.json`. The source archives remain local and unchanged. Nested historical parent ZIPs, chat exports, screenshots, and private task context are omitted.

Drafting, code adaptation, mathematical review and reproducibility checks used OpenAI Codex with GPT-6 Astra. These are AI-assisted checks within the project, not independent expert peer review. No unverified affiliation, ORCID, or new physical prediction is asserted.
