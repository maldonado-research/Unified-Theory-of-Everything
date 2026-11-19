# A full four-lepton rate in a renormalizable singlet model

19 November 2026 · Candidate.20 · Follow-up to the [energy-scan calculation](../energy_scan/00_READ_FIRST.md) at commit `109872358a2f85c50493953d759a02ccd6c67bda`.

**Three evaluations of the loop integrals give the complete leading rate for polarized $\mu^-_L\mu^-_L\to e^-_Le^-_L$ in a specified four-singlet model.** The preceding checkpoint left the analytic four-lepton interaction free. Here the masses and Yukawa couplings determine it. The calculation also quantifies how accurately the first two derivative terms reproduce the physical two-scalar cut.

## What is calculated

The scalar has mass $m>0$ in an unbroken, gaugeless theory with massless external leptons. A pair of Majorana singlets with masses $m,4m$ couples to the muon doublet. A pair with masses $M,2M$ couples to the electron doublet. For this benchmark $M=10m$. Denote the first column's Yukawa in each block by $y_\mu$ and $y_e$. The relative **squared** Yukawa couplings within the two blocks are $(1,-1/4)$ and $(1,-2)$ respectively. These phases can be realized with ordinary complex Yukawa couplings. The electron pair cancels its tree-level dimension-five interaction and leaves a dimension-seven coefficient $3y_e^2/(4M^3)$ in the original fields. The muon pair cancels a different mass-weighted sum. [Action and matching conventions](MATCHING.md)

| Squared energy $s/m^2$ | Full $m^2\sigma/\lvert y_e y_\mu\rvert^4$ | Lower bound from the $J=0$ scalar cut, same units |
|---:|---:|---:|
| 6 | $4.00925370\times10^{-11}$ | $7.01354913\times10^{-12}$ |
| 10 | $7.17400709\times10^{-11}$ | $3.48420371\times10^{-11}$ |

The full rate integrates over all scattering angles, includes the identical-electron factor and uses fixed incident helicities. The scalar cut accounts for about 17.5% and 48.6% of the full rate at these two settings. The remaining contribution comes from the real part after removing the common Yukawa phase. The transition converts two units of muon flavor into electron flavor with zero net change in total external lepton number.

A dispersion integral and a direct integration of the four-propagator Feynman integral agree on both components of the amplitude and the integrated rates. A third evaluation with QCDLoop 2.1.0 uses its published scalar box functions. It reproduces all ten tabulated amplitude and rate comparisons across both energies within a relative difference of $3\times10^{-15}$. Independently sewing the two physical scalar-production trees checks the imaginary part and the normalization of the bound. The bound applies to the direct process because the chosen flavor texture gives all of its diagrams the same overall coupling phase. [Dispersion calculation](DISPERSION.md) · [Direct box integration](BOX_CROSSCHECK.md) · [On-shell check](CUT.md) · [QCDLoop comparison](QCDLOOP.md)

For $4<s/m^2\le11$, keeping only the matched dimension-seven electron interaction incurs at most **11.46%** relative error in the absorptive amplitude. Including the matched dimension-nine term reduces this to **1.09%**. Both estimates use the exact result as the denominator and hold at every scattering angle. The full rate uses the unexpanded propagators. [Derivation of the error bounds](DISPERSION.md)

## What the benchmark establishes

The specified renormalizable action is an example in which zero-momentum cancellation coexists with a nonzero calculable flavor-changing rate. Prior work already describes the inverse-mass moment mechanism and the finite-momentum Majorana box topology. The contribution here is the matched two-pair benchmark, its normalized rates and the controlled range of its derivative approximation. The mass and coupling choices are assumptions. The broader unification hypothesis supplies no derivation of these inputs in this checkpoint. The tree cancellation can shift under radiative corrections. A comparison to experiment requires the broken electroweak spectrum, gauge interactions, physical masses and mixings, and experimental constraints. [Literature and public portfolio comparison](SOURCES.md)

## Reproduction

Each program accepts `--output PATH` or prints its result to standard output. The full repository replay regenerates and compares all four records in the pinned environment. The dispersion program also compares its results with the independent direct-box and QCDLoop records. The external library was run separately with the [pinned C++ reproduction commands](QCDLOOP.md). [Full repository receipt](../../verification/UV_COMPLETION_REPLAY_LINUX.json)

- [Tree matching and remainders](check_matching.py) · [Matching receipt](MATCHING.json)
- [Two-scalar cut and partial waves](check_cut.py) · [Cut receipt](CUT.json)
- [Upper-rim dispersion integral and rational bounds](check_bounds.py) · [Dispersion receipt](BOUNDS.json)
- [Direct simplex integration with grid refinement](check_box_direct.py) · [Direct receipt](BOX_DIRECT.json)

The next matching test is to compute the analytic four-lepton contact below threshold across several heavy masses, then use it to predict the amplitude at a second energy without refitting. The radiative stability of the assumed Yukawa ratios is a separate open calculation.
