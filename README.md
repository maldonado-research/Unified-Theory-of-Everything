# Unified Theory of Everything — research program

[Readable project overview](https://maldonado-research.github.io/projects/unified-theory-of-everything/) · [All research projects](https://maldonado-research.github.io/)


**Ricardo Maldonado · hypothesis and reproducible research · GitHub edition 30 September 2026**

Could the universe's many rules come from something deeper? This project explores that question through specific calculations we can inspect and reproduce. The current work concerns selected neutrino kernels, a corrected quantum-loop argument, a conditional scalar continuum and a changing-mass toy model.

The program remains unfinished. This repository makes its limited results, corrections and open questions public; it does not establish a complete theory of everything or an experimentally confirmed new law.

## Start here

- [Plain-language overview](docs/OVERVIEW.md): what the research is asking and what has been calculated.
- [Research status](docs/RESEARCH_STATUS.md): supported results, assumptions and unresolved physical steps.
- [Published N00AK-r1 technical note](baseline/TOE_N00AK_r1_PUBLIC/TOE_N00AK_r1_TECHNICAL_NOTE.md): the unchanged 4 September baseline and scope correction.
- [6 September checkpoint, public adaptation](checkpoint_2026_09_06/00_READ_FIRST.md): the later seed-loop result and supporting diagnostics.
- [Verification guide](docs/VERIFICATION.md): dependencies, quick checks, full replay and their limits.
- [Next calculations](docs/ROADMAP.md): a concrete path toward a more complete physical test.

## What the package contains

| Branch | Result within its declared assumptions | Important limit |
| --- | --- | --- |
| Heavy-neutrino tree kernel | Exact finite-pole expression, local moments, remainder and known-node reconstruction | Selected tree object; masses and couplings are inputs, and unknown poles require more information |
| Quantum correction | Evaluated off-shell ultraviolet polynomial for one derivative seed | Full operator basis, gauge completion, sourced-equation reduction and physical running remain unfinished |
| Scalar continuum | Finite-moment criterion, nonanalytic low-energy behavior and a positive single-pole lower bound | Inherited spectral premises; no completed source-to-observable susceptibility |
| Changing-mass diagnostic | A rigorously positive weak-pulse occupation and solver comparisons | Externally prescribed background; no complete cosmic energy account or matter–antimatter asymmetry |

The earlier physical dimension-seven coefficient and its all-odd extension were withdrawn in the r1 correction. Their cancellation survives as conditional algebra. The later seed calculation provides a missing ingredient without reinstating the withdrawn physical claim. In particular, the five-ninths projected ratio is not a rule for replacing the old beta-function coefficient.

## Reproduce the focused checks

Use Python 3.12 and a disposable virtual environment:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python verify.py
python verify.py --full
```

The quick command verifies the public payload and runs tree/conditional, routed-seed, weak-pulse and scalar controls. The full command adds the six-case pulse benchmark and the separate rotated-basis Radau comparison. The runner works on a temporary copy and leaves the shipped files intact. Do not add `-O` to checkpoint replay commands: the numerical review uses scientific assertions, and the publication copy now rejects optimization explicitly.

Checks support reproducibility of the declared formulas and toy models. They are not experimental validation, a proof of the parent scalar assumptions or independent human peer review. [The publication replay receipt](verification/PUBLICATION_REPLAY.json) records the tested environment and results.

## Dates, citations and reuse

The [Zenodo N00AK-r1 record](https://zenodo.org/records/22319172), DOI **10.5281/zenodo.22319172**, identifies the unchanged baseline published on 4 September 2026. The larger GitHub edition is a distinct publication dated 30 September and contains public adaptations of the 6 September local checkpoint. The existing DOI must not be cited as if it identifies all the added GitHub material. Use the repository release tag or commit when citing this edition; [CITATION.cff](CITATION.cff) supplies citation metadata.

The author-owned archive is offered under [CC BY 4.0](LICENSE), preserving the baseline's license. Third-party dependencies and cited papers retain their own licenses. AI-assisted derivation, editorial preparation and computational reviews are disclosed; they do not constitute independent human peer review. The [provenance record](checkpoint_2026_09_06/PROVENANCE.json) identifies original source hashes and the limited publication edits.

For a specific correction or reproduction question, see [CONTRIBUTING.md](CONTRIBUTING.md). The author's different research programs retain separate physical assumptions; no merger is established by this repository.
