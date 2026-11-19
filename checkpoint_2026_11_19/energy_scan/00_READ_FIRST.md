# What an energy scan can test when the lower-dimensional coupling is unknown

19 November 2026. Provisional candidate.19, following the [on-shell scattering checkpoint](../on_shell/00_READ_FIRST.md) at commit `e0f2556190ca906f1349a514be774de9e2951221`.

**An independent dimension-five neutrino coupling can erase the earlier scattering cut at a chosen energy. A common coupling cannot erase it at two distinct energies if the derivative slope is nonzero.** Across four settings a separate test works with only the muon-initiated process, even after fitting complex local four-lepton terms. These statements concern leading amplitudes in a specified effective theory with massless left-handed external leptons.

## How the cancellation changes

Use equal weak components in μ⁻μ⁻→e⁻e⁻, and its inverse where indicated. The scalar squared mass is a, the sterile masses are √a and 4√a, and the coherent sterile mass/Yukawa weights are +w₀ and −w₀. This cancels the derivative-seed ultraviolet pole. Both Yukawa couplings attach to the incoming muon flavor in the direct process. The independent dimension-five coefficient C₅ and two-derivative coefficient C₇ attach to the electron pair. The scalar decay into a coupled sterile and lepton has zero phase space at the lighter mass and is closed at the heavier mass.

The [lower-operator calculation](LOWER_OPERATOR.md) gives the complete scalar-triangle contribution of C₅, including its finite real part, in the reference action coordinates. In the two-scalar cut the electron-side factor is proportional to `(s/2−a) C₇−C₅`. A single collision energy therefore leaves a choice of complex C₅ that hides the cut. In particular, at `s=8a` the cancellation occurs at `C₅=3a C₇`. At that point a finite real bubble contribution remains; the absolute dispersive amplitude also depends on local four-lepton matching.

## Three scan tests

Let `ξ=s/a` and `g₇=a C₇ w₀`. Write `σ̄` for the arithmetic mean of the two reciprocal exclusive cross sections and `σ_dir` for μ⁻μ⁻→e⁻e⁻ alone. Exposure averages cross sections across runs. The benchmarks fix incident helicities and use the identical-final-particle phase space. Each lower bound below is on `a⟨σ⟩/|g₇|²` for nonzero g₇.

| Exposure and process | Unknown contributions allowed | Lower bound |
| --- | --- | ---: |
| Equal exposure at ξ=6 and 10, both reciprocal directions | One arbitrary complex C₅ and same-order Hermitian local four-lepton amplitudes | 3.78076 × 10⁻⁶ |
| Equal exposure in the two bins [5,7] and [9,11], uniform in ξ, both directions | Same | 3.87456 × 10⁻⁶ |
| Equal exposure at ξ=5,6,9,11, direct process only | One arbitrary complex C₅ plus independent complex constant and linear terms in the reduced local four-lepton amplitude | 1.33764 × 10⁻⁹ |

For the first two rows the unknown C₅ can adjust the constant term in the cut. Minimizing over it leaves the weighted variance of `ξ/2−1`. For the third row, any three of the four direct amplitudes can be tuned to vanish in the stated model. The four rows together contain a component that survives all three complex adjustments. An interval calculation also guarantees a coefficient above **8.61 × 10⁻¹⁰** if each of the four beam settings has any exposure law supported within ±0.001 in ξ. [Reciprocal scans](RATE_SCAN.md) · [Direct scan](SINGLE_DIRECTION.md) · [Rational certificate](CERTIFICATE.md)

The quantity these tests probe is the total linear slope of the on-shell lepton–scalar contact. With only the two specified seeds that slope fixes C₇ in the declared normalization. Other allowed two-derivative interactions can give the same projected slope and can cancel a nonzero literal C₇ at all energies. Higher-dimensional energy shapes, other cuts and additional light states can change the bounds. In strict EFT power counting, dimension-five by dimension-nine interference can enter at the order of a dimension-seven square. The benchmarks therefore assume a controlled truncation of the *amplitude* or the explicitly restricted model. Their numerical sizes require the masses and couplings, and an experimental test would require flux, resolution, backgrounds and event statistics. [Local contact classification and counterexamples](CONTACT_SCOPE.md)

## Reproduction and sources

Five programs accompany this checkpoint. In the pinned repository environment run each with `python -I -B -O` from this folder.

- [Lower-operator integral](verify_lower_operator.py) and [results](LOWER_OPERATOR_RESULTS.json)
- [Reciprocal energy scans](verify_rate_scan.py) and [results](RATE_SCAN_RESULTS.json)
- [On-shell contact projection](verify_contact_projection.py) and [results](CONTACT_PROJECTION_RESULTS.json)
- [Single-direction projection and quadratures](verify_single_direction.py) and [results](SINGLE_DIRECTION_RESULTS.json). This check reuses the lower-operator Feynman-parameter routine for its second integration route.
- [Independent rational enclosures](verify_scan_certificate.py) and [endpoints](SCAN_CERTIFICATE_RESULTS.json)

The [source review](SOURCES.md) compares contact bases, EFT truncation and weighted-design methods, and records changes in the author's public Zenodo and HDblast catalogs. The [root replay](../../verification/ENERGY_SCAN_REPLAY_LINUX.json) checks the full package. The N00AK-r1 baseline has its own [archival record](https://doi.org/10.5281/zenodo.22319172).
