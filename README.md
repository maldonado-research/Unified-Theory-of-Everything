# Unified Theory of Everything — research program

[Readable project overview](https://maldonado-research.github.io/projects/unified-theory-of-everything/) · [All research projects](https://maldonado-research.github.io/)


**Ricardo Maldonado · hypothesis and reproducible research · two-gauge-leg candidate 15 October 2026**

Could the universe's many rules come from something deeper? This project explores that question through specific calculations we can inspect and reproduce. The current work concerns selected neutrino kernels, a corrected quantum-loop argument, a conditional scalar continuum and a changing-mass toy model.

The program remains unfinished. This repository makes its limited results, corrections and open questions public; it does not establish a complete theory of everything or an experimentally confirmed new law.

## Start here

- [Plain-language overview](docs/OVERVIEW.md): what the research is asking and what has been calculated.
- [Research status](docs/RESEARCH_STATUS.md): supported results, assumptions and unresolved physical steps.
- [Published N00AK-r1 technical note](baseline/TOE_N00AK_r1_PUBLIC/TOE_N00AK_r1_TECHNICAL_NOTE.md): the unchanged 4 September baseline and scope correction.
- [6 September checkpoint, public adaptation](checkpoint_2026_09_06/00_READ_FIRST.md): the later seed-loop result and supporting diagnostics.
- [1 October operator audit](checkpoint_2026_10_01/00_READ_FIRST.md): a local seed-operator representative, sourced interaction contacts and separate exact checks.
- [15 October gauge-leg candidate](checkpoint_2026_10_15/00_READ_FIRST.md): two internal hypercharge attachments and their explicitly computed transverse UV residue.
- [Action and source continuation](checkpoint_2026_10_15/operator_map/00_READ_FIRST.md): physical charge-conjugation signs, local operators, sourced currents and an exact flavor-contact criterion.
- [Two-hypercharge-leg calculation](checkpoint_2026_10_15/two_gauge_leg/00_READ_FIRST.md): seven evaluated loop diagrams match the predicted two-B action vertex, with independent integration and Ward controls.
- [Verification guide](docs/VERIFICATION.md): dependencies, quick checks, full replay and their limits.
- [Next calculations](docs/ROADMAP.md): a concrete path toward a more complete physical test.

## What the package contains

| Branch | Result within its declared assumptions | Important limit |
| --- | --- | --- |
| Heavy-neutrino tree kernel | Exact finite-pole expression, local moments, remainder and known-node reconstruction | Selected tree object; masses and couplings are inputs, and unknown poles require more information |
| Quantum correction | Evaluated off-shell ultraviolet polynomial for one derivative seed | Full operator basis, gauge completion, sourced-equation reduction and physical running remain unfinished |
| Operator-reduction audit | Local gauge-covariant representative, explicit sourced contacts and a basis-change cancellation | A complete basis, remaining sourced matching and physical loop running remain unfinished |
| One hypercharge leg | A transverse residue beyond the earlier representative for two declared seed-loop attachments, with a literal action/sign map | Restricted coupling sector; complete operator matching and physical running remain unfinished |
| Two hypercharge legs | Seven loop diagrams reproduce the mass-independent two-B kernel of the chosen action; locality and degree bounds fix its restricted Ward completion uniquely | No full EFT basis or physical running; the three-B statement is conditional and has no explicit loop integration |
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

The quick command verifies the public payload and runs tree/conditional, routed-seed, weak-pulse, scalar, operator/contact and independently authored gauge-leg/action/source controls, including three two-B stages (eighteen stages). The full command adds the six-case pulse benchmark and the separate rotated-basis Radau comparison (twenty stages). The runner works on a temporary copy and leaves the shipped files intact. Do not add `-O` to checkpoint replay commands: the numerical review uses scientific assertions, and the publication copy rejects optimization explicitly.

Checks support reproducibility of the declared formulas and toy models. They are not experimental validation, a proof of the parent scalar assumptions or independent human peer review. [The earlier publication receipt](verification/PUBLICATION_REPLAY.json) retains the September replay; [the operator-edition receipt](verification/OPERATOR_REPLAY_LINUX.json) retains the October 1 full replay. The gauge-leg checkpoint records its own validation separately. The [two-gauge-leg full replay](verification/TWO_GAUGE_REPLAY_LINUX.json) records all twenty stages passing; stage counts alone do not establish the scientific claims beyond their declared scope.

## Dates, citations and reuse

The [Zenodo N00AK-r1 record](https://zenodo.org/records/22319172), DOI **10.5281/zenodo.22319172**, identifies the unchanged baseline published on 4 September 2026. The September 30 GitHub edition added public adaptations of the September 6 checkpoint; the October 1 edition added the operator audit. The October 15 branch is provisional candidate.3, extending the one-B action/source map through explicit two-B loop matching. The existing DOI must not be cited as if it identifies later GitHub material. Use the reviewed repository commit when citing a candidate; [CITATION.cff](CITATION.cff) supplies citation metadata. A [Zenodo update description](docs/ZENODO_OPERATOR_UPDATE.md) remains specific to the earlier operator edition; no new deposit or DOI is claimed.

The author-owned archive is offered under [CC BY 4.0](LICENSE), preserving the baseline's license. Third-party dependencies and cited papers retain their own licenses. AI-assisted derivation, editorial preparation and computational reviews are disclosed; they do not constitute independent human peer review. The [provenance record](checkpoint_2026_09_06/PROVENANCE.json) identifies original source hashes and the limited publication edits.

For a specific correction or reproduction question, see [CONTRIBUTING.md](CONTRIBUTING.md). The author's different research programs retain separate physical assumptions; no merger is established by this repository.
