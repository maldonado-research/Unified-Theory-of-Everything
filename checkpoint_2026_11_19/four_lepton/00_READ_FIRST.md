# A finite four-lepton response after cancellation of the ultraviolet pole

19 November 2026. Provisional GitHub candidate.17, based on research commit `4e4630ae9a4a617453e899c5f075156fc264e3ff`.

**Two unequal sterile masses can cancel the logarithmic pole exactly and still produce a finite, momentum-dependent four-lepton vertex contribution.** This result concerns a particular one-loop graph with a derivative neutrino interaction, two Yukawa couplings and retained Majorana masses. The October calculation fixed its ultraviolet residue. Here the entire Euclidean momentum difference of that graph is evaluated on a specified off-shell slice. [Original pole](../../checkpoint_2026_10_23/seed_yukawa/FOUR_LEPTON_POLE.md) · [Common action convention](../../checkpoint_2026_11_04/common_action/00_READ_FIRST.md)

## A test with opposite pole weights

Take one active flavor, a real unit seed coefficient in the chosen units, a scalar of mass one, two sterile masses 1 and 4 and a conjugated Yukawa row `(1,i/2)`. Their pole weights are exactly `1` and `−1`. To test the finite contribution, set the two barred momenta to zero and the two unbarred momenta to `+p,−p`. Subtract the vertex at zero momentum. A local dimension-six subtraction constant drops out. Both orientations of the seed must be evaluated at these same labelled external momenta.

| Reduced contribution | Pole weight sum | Subtracted value at Euclidean p²=1 |
| --- | ---: | ---: |
| Seed on the two zero-momentum legs | 0 | −0.06184365542 |
| Conjugate seed on the two legs carrying momentum | 0 | +0.36950480197 |
| Sum in the one-flavor pair tensor | 0 | **+0.30766114655** |

The entries omit the common spin tensor, the units carried by the seed and mass, and the loop factor `1/(32π²)`. Direct angular integration and Feynman parameters give the same functions. The difference between the two rows traces to the derivative on the whole composite at the seed. Reversing its placement at fixed momenta changes the numerator. [Integral and conventions](FINITE_TRIANGLE.md) · [Independent calculation](INDEPENDENT_FINITE.md)

For the ordered soft-seed graph, the pole depends only on the sum of its flavor tensors over sterile masses. To cancel the finite difference for an interval of momenta, **each group of equal positive masses must cancel separately**. Each distinct mass produces a distinct first branch point. Approaching the lightest one isolates that group's coefficient, and the argument repeats for heavier groups. For two opposite real pole weights in the one-flavor example, the paired slope is also strictly ordered by sterile mass. These statements hold for the specified integral and fixed field coordinates. [Proof and counterexample](MASS_RESOLUTION.md) · [Numerical values](MASS_RESOLUTION_RESULTS.json)

## Interpretation and reproduction

The scalar squared mass is positive, the momenta are Euclidean, and the external fermion legs are off shell. Independent local derivative operators can affect this same momentum dependence near zero. The calculation supplies a benchmark for matching and a check on any future claim of a finite cancellation. A prediction for a rate still requires the remaining contributions, physical external states and matched operator coefficients. The neutral-photon calculation is a separate channel.

The [one-channel evanescent convention](EVANESCENT_CONVENTION.md) gives the finite local conversion between two dimensional Fierz prescriptions and shows why it cancels in this reduced momentum difference. The [literature note](LITERATURE.md) compares published matching and scalar-integral conventions. The [archive survey](PUBLICATION_SURVEY.md) screens the author's Zenodo records and identifies a phase-sensitive test from HDblast. No numerical coupling or observable is transferred from that cosmological model.

Four standalone scripts accompany the calculations. Run each from this directory with the repository's pinned Python environment.

```sh
python -I -B -O verify_finite_triangle.py
python -I -B -O independent_finite.py
python -I -B -O verify_fierz_scheme.py
python -I -B -O verify_mass_resolution.py
```

The [root verification guide](../../docs/VERIFICATION.md) describes the full replay. The [public receipt](../../verification/FINITE_FOUR_LEPTON_REPLAY_LINUX.json) records the tested release. The original research remains under review. The archived N00AK-r1 DOI identifies the September baseline. [Zenodo baseline](https://doi.org/10.5281/zenodo.22319172)
