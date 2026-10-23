# Unified Theory of Everything — research program

[Readable project overview](https://maldonado-research.github.io/projects/unified-theory-of-everything/) · [All research projects](https://maldonado-research.github.io/)


**Ricardo Maldonado · hypothesis and reproducible research · seed/Yukawa matter-loop candidate 23 October 2026 (Pacific)**

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
- [Electroweak and selected-source continuation](checkpoint_2026_10_15/electroweak_source/00_READ_FIRST.md): a calculated one-W UV pole, its restricted neutral-photon cancellation, and an explicit reduction retaining the sources of the combined hypercharge packet.
- [Explicit field changes and neutral normalization](checkpoint_2026_10_15/field_redefinition/00_READ_FIRST.md): indexed lepton/Higgs/Majorana/gauge shifts, the actual neutral quadratic projection, correlated mass/source normalization, and a concrete lower-order insertion check.
- [Non-Abelian gauge completion](checkpoint_2026_10_15/nonabelian_gauge/00_READ_FIRST.md): evaluated two-W/mixed-BW poles, 33 three-gauge diagrams, independent cubature and Ward checks, and the scope of the neutral photon cancellation.
- [Interacting electroweak reduction](checkpoint_2026_10_15/interacting_reduction/00_READ_FIRST.md): source-retaining matrix reduction, indexed real weak-field changes, combined weak/flavor exchange, and a local differential normalization with its interaction remainder.
- [Algebraic operator sectors](checkpoint_2026_10_15/operator_basis/00_READ_FIRST.md): independent tensor/current coordinates, exact flavor/Fierz coefficient transport, and a symmetric-seed source-image constraint with a sharp algebraic norm bound.
- [Selected action-map and scalar-current completion](checkpoint_2026_10_16/action_completion/00_READ_FIRST.md): an explicit bridge between the field maps, an off-shell scalar-tadpole Hessian check, source-retaining scalar-current reduction, and charged-lepton/quark Yukawa transport.
- [Derivative placement and source projection](checkpoint_2026_10_16/derivative_projection/00_READ_FIRST.md): curvature-derivative invariants and an equivalent row-derivative representative, complete neutrino-source active field-content projection, and general complex-Majorana covariance.
- [One-loop Yukawa inputs](checkpoint_2026_10_16/yukawa_counterterms/00_READ_FIRST.md): independently integrated field/mass poles, proper renormalizable vertex UV analysis, known restricted parameter running, and the derivative seed's external-field counterterm input.
- [Derivative-seed matter loops](checkpoint_2026_10_23/seed_yukawa/00_READ_FIRST.md): exhaustive proper topology, independently integrated four-lepton/four-Higgs poles, primitive lower-dimensional mixing and a source-retaining scalar phase change.
- [Verification guide](docs/VERIFICATION.md): dependencies, quick checks, full replay and their limits.
- [Next calculations](docs/ROADMAP.md): a concrete path toward a more complete physical test.

## What the package contains

| Branch | Result within its declared assumptions | Important limit |
| --- | --- | --- |
| Heavy-neutrino tree kernel | Exact finite-pole expression, local moments, remainder and known-node reconstruction | Selected tree object; masses and couplings are inputs, and unknown poles require more information |
| Quantum correction | Evaluated off-shell ultraviolet polynomial for one derivative seed | Full physical operator basis, sourced matching and physical running remain unfinished |
| Operator-reduction audit | Local gauge-covariant representative, explicit sourced contacts and a basis-change cancellation | A complete basis, remaining sourced matching and physical loop running remain unfinished |
| One hypercharge leg | A transverse residue beyond the earlier representative for two declared seed-loop attachments, with a literal action/sign map | Restricted coupling sector; complete operator matching and physical running remain unfinished |
| Two hypercharge legs | Seven loop diagrams reproduce the mass-independent two-B kernel of the chosen action; locality and degree bounds fix its restricted Ward completion uniquely | No full EFT basis or physical running; its earlier conditional three-B result is now checked in the separate three-gauge calculation |
| Electroweak and selected-source continuation | One-W pole and selected neutral-photon UV cancellation; combined hypercharge packet reduced with an exact off-shell EOM-residual certificate | Cancellation concerns the selected local pole; the source representative does not determine a physical dipole |
| Field changes and neutral normalization | Explicit field substitutions realize that residual certificate; the selected neutral kinetic mixing is normalized with its induced mass/source terms; a dimension-six × dimension-five contact is checked | Selected classical action and quadratic background sector; full interacting normalization, complete EFT closure and physical running remain open |
| Non-Abelian gauge completion | Ordered two-W/mixed-BW residues and 33 integrated three-gauge diagrams match the full curvature; degree counting closes this restricted lowest-order external-gauge hierarchy | Other matter sectors, full interacting physical reduction and running remain open; a mixed photon–charged-weak local contact is not a finite observable |
| Interacting electroweak reduction | The selected matrix action has explicit source-retaining field changes and current descendants; a covariant local normalization retains scalar derivatives, mass and source terms | Selected action coordinates; other loop sectors, complete coefficient basis, subtraction prescription and physical running remain open |
| Algebraic operator sectors | Independent nonderivative tensor/current bases and explicit source maps; fixed-Yukawa symmetric-seed image, rank and relative norm bound | Only the declared algebraic sectors and explicit contributions; other terms can feed the same operators, and full derivative/EOM matching remains open |
| Selected action completion | Explicit source-to-normalizing coordinate bridge, directly checked off-shell neutral action including scalar tadpoles; scalar-current source reduction and spectator-Yukawa transport | Classical first-order packet and specified sectors; full derivative/EOM coefficients, quantum source transport and physical running remain open |
| Derivative placement and source projection | Exact curvature/source coefficient transport, complete active projection of the declared neutrino-source packet, and covariant complex-Majorana recurrences | Field-content projection is not sterile elimination; the placement basis and selected coefficient tables do not close the physical EFT |
| Pure-neutrino-Yukawa loop inputs | Independently integrated lepton/sterile/Higgs field and mass poles; proper vertex UV control and known restricted renormalizable running; seed external-field contribution | Gauge, other interactions and proper higher-operator mixing remain outside this calculation; the seed contribution is not its full physical running |
| Seed/Yukawa matter loops | No proper original seed-shaped four-point graph at this order; independently integrated Majorana-mass four-lepton, scalar derivative and potential poles, plus normalized primitive flavor image | Mixed lepton/Higgs pole, source/EOM transport and total physical coefficients remain open; retained fields are not integrated out |
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

The quick command verifies the public payload and runs tree/conditional, routed-seed, weak-pulse, scalar, operator/contact and independently authored gauge-leg/action/source/normalization controls in sixty-two stages. The full command adds the six-case pulse benchmark and the separate rotated-basis Radau comparison for sixty-four stages. The runner reports each stage, works on a temporary copy and leaves the shipped files intact. Do not add `-O` to checkpoint replay commands: the numerical review uses scientific assertions, and the publication copy rejects optimization explicitly.

Checks support reproducibility of the declared formulas and toy models. They are not experimental validation, a proof of the parent scalar assumptions or independent human peer review. [The earlier publication receipt](verification/PUBLICATION_REPLAY.json) retains the September replay; [the operator-edition receipt](verification/OPERATOR_REPLAY_LINUX.json) retains the October 1 full replay. The [two-gauge-leg full replay](verification/TWO_GAUGE_REPLAY_LINUX.json) preserves its twenty-stage result. The [electroweak/source full replay](verification/ELECTROWEAK_SOURCE_REPLAY_LINUX.json) preserves its twenty-three stages. The [field-redefinition replay](verification/FIELD_REDEFINITION_REPLAY_LINUX.json) preserves its twenty-eight stages. The [non-Abelian replay](verification/NONABELIAN_GAUGE_REPLAY_LINUX.json) preserves its thirty-four stages. The [interacting-reduction replay](verification/INTERACTING_REDUCTION_REPLAY_LINUX.json) preserves its thirty-eight stages. The [operator-sector replay](verification/OPERATOR_BASIS_REPLAY_LINUX.json) preserves its forty-three-stage result. The [action-completion replay](verification/ACTION_COMPLETION_REPLAY_LINUX.json) preserves its forty-eight-stage result. The [derivative/source replay](verification/DERIVATIVE_PROJECTION_REPLAY_LINUX.json) preserves its fifty-three-stage full check. The [Yukawa-counterterm replay](verification/YUKAWA_COUNTERTERM_REPLAY_LINUX.json) preserves its fifty-eight-stage full check. The [seed/Yukawa matter replay](verification/SEED_YUKAWA_REPLAY_LINUX.json) records the extended sixty-four-stage full check; stage counts alone do not establish scientific claims beyond their declared scope.

## Dates, citations and reuse

The [Zenodo N00AK-r1 record](https://zenodo.org/records/22319172), DOI **10.5281/zenodo.22319172**, identifies the unchanged baseline published on 4 September 2026. The September 30 GitHub edition added public adaptations of the September 6 checkpoint; the October 1 edition added the operator audit. The current branch is provisional candidate.12, dated October 23 Pacific, adding an exhaustive derivative-seed/Yukawa proper-topology certificate and independently integrated mass-dependent four-lepton/four-Higgs poles. Primitive lower-dimensional mixing, flavor image and a source-retaining scalar phase change are explicit; the mixed matter pole and complete physical running remain unfinished. The existing DOI must not be cited as if it identifies later GitHub material. Use the reviewed repository commit when citing a candidate; [CITATION.cff](CITATION.cff) supplies citation metadata. A [Zenodo update description](docs/ZENODO_OPERATOR_UPDATE.md) remains specific to the earlier operator edition; no new deposit or DOI is claimed.

The author-owned archive is offered under [CC BY 4.0](LICENSE), preserving the baseline's license. Third-party dependencies and cited papers retain their own licenses. AI-assisted derivation, editorial preparation and computational reviews are disclosed; they do not constitute independent human peer review. The [provenance record](checkpoint_2026_09_06/PROVENANCE.json) identifies original source hashes and the limited publication edits.

For a specific correction or reproduction question, see [CONTRIBUTING.md](CONTRIBUTING.md). The author's different research programs retain separate physical assumptions; no merger is established by this repository.
