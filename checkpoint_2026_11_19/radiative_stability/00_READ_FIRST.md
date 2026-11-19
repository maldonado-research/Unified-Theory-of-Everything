# What a scalar self-interaction does to the four-lepton rate

19 November 2026 · Candidate.21 · Follow-up to the [four-singlet calculation](../uv_completion/00_READ_FIRST.md) at commit `38205a427ae24bdaac7893079f2c3169b0201e38`.

**The scalar self-interaction lifts the electron-pair cancellation and changes the predicted rate of $\mu^-_L\mu^-_L\to e^-_Le^-_L$.** At the first order where this interaction contributes, the ratio of the cross sections at $s=10m^2$ and $s=6m^2$ is $1.78936222+0.72504369\,\lambda$. Here $m$ is the fixed scalar pole mass and $\lambda$ is defined by $V\supset+\lambda(H^\dagger H)^2$. The linear coefficient comes from a two-loop amplitude. The value of $\lambda$ remains an input.

## Does the zero survive running and matching?

The benchmark has two electron-coupled Majorana singlets of masses $M,2M$ with relative squared Yukawas $1,-2$. At tree level their inverse-mass sum vanishes. The muon block has masses $m,4m$ and relative squared Yukawas $1,-1/4$. The numerical example sets $M=10m$. Gauge couplings vanish and external leptons are massless in the unbroken theory.

- **Pure Yukawa interactions.** With this flavor texture the electron zero survives one-loop evolution above the heavy thresholds, sequential evolution at leading dimension five, and the corresponding pure-Yukawa one-loop matching onto the electron dimension-five interaction. An exact matrix calculation and an independent diagram count check these statements. Radiative generation of a scalar quartic from Yukawa couplings enters the quartic effect at a further order if its input value is zero. [Running](RUNNING.md) · [Matching audit](MATCHING_AUDIT.md)
- **Scalar self-interaction.** In the original field coordinates, matching at leading order in inverse heavy mass gives $\Delta C_{5,ee}=2\lambda y_e^2\ln4/(16\pi^2 M)$. Direct integration at finite $m$ gives the zero-external-momentum four-point coefficient $0.0168504162\,\lambda y_e^2/M$ at $M=10m$, about 4.03% below its leading heavy-mass limit. This latter coefficient includes light-scalar effects. A published one-loop matching formula gives the same sign and normalization for the leading heavy threshold. [Finite-mass calculation](QUARTIC.md) · [Published comparisons](SOURCES.md)

## What changes in scattering?

At fixed scalar pole mass the connected two-loop amplitude proportional to $\lambda y_e^2(y_\mu^*)^2$ reduces to two triangles joined at the quartic vertex. The remaining quartic insertions shift the scalar mass and cancel against its mass counterterm. This fixes the complete coefficient at that coupling order. Let $r=m^2\sigma/|y_e y_\mu|^4$. The table gives $r=r_0+\lambda r_1$ through the displayed orders.

| Squared energy $s/m^2$ | Leading $r_0$ | Quartic coefficient $r_1$ |
|---:|---:|---:|
| 6 | $4.00925370\times10^{-11}$ | $-4.82330547\times10^{-11}$ |
| 8 | $5.61821543\times10^{-11}$ | $-5.30510123\times10^{-11}$ |
| 10 | $7.17400709\times10^{-11}$ | $-5.72375648\times10^{-11}$ |

A small positive quartic reduces all three cross sections. It increases the endpoint ratio because the fractional reduction is largest at the lower energy. Checks include the physical two-scalar cut and an external QCDLoop calculation of all sixteen complex triangle reference values. The full next-order rate still requires the separate diagrams with extra Yukawa couplings. [Two-loop derivation, rate and scan](RESCATTERING.md) · [QCDLoop comparison](QCDLOOP.md)

## Can this effect be distinguished from a change in the Yukawas?

A real 1% change in the **squared Yukawa coupling of the second electron-coupled singlet** relative to its tuned value, at fixed masses and $\lambda=0$, changes the endpoint ratio to 1.48888 for a positive shift and 2.39097 for a negative shift. An energy scan therefore needs control of the coupling ratio. Using the reciprocal means of the forward and inverse rates at all three energies, a weighted ratio cancels a common normalization and the leading sensitivity to complex detuning. The resulting fractional response to $\lambda$ has magnitude only $0.00078244$ per unit quartic in the logarithmic statistic. This makes these three energies a poor choice for isolating the scalar interaction. [Detuning calculation](DETUNING.md) · [Definition and limitations of the statistic](RESCATTERING.md)

Three forward rates and one reverse rate at a shared energy give a stronger consistency check. Their signed difference fixes the imaginary part of the detuning up to normalization. Removing its predicted contribution from the forward rates gives three combinations proportional to the overall normalization, the real coupling change and its squared magnitude. These four numbers obey an exact quadratic identity for every finite complex change of the single squared Yukawa at leading order. The quartic violates this identity at first order with a small calculable coefficient. This supplies a model test even when the unknown Yukawa change exceeds the range of a linear approximation. [Exact rate relation and sensitivity](CLOSURE.md)

Neither the mass ratios nor the couplings have been derived from a unifying principle in this calculation. Testing the benchmark against data would require a physical electroweak spectrum, gauge interactions and constraints on the masses and Yukawas. The mechanism behind the dimension-five correction already appears in the seesaw literature. The contribution here is the matched finite-mass and rate calculation for this model and the quantitative test of how well an energy scan distinguishes its parameters. [Sources and project inventory](SOURCES.md)

## Reproduction

In the pinned repository environment run `python -B verify.py --full` from the repository root. Each of the following Python programs also accepts `--output PATH` and otherwise prints its JSON result. The optional external C++ calculation has separate pinned build instructions. [Full replay](../../verification/RADIATIVE_STABILITY_REPLAY_LINUX.json)

- [Matrix evolution](check_running.py) · [Results](RUNNING.json)
- [Quartic matching](check_quartic.py) · [Results](QUARTIC.json)
- [Independent diagram audit](check_audit.py) · [Results](AUDIT.json)
- [Complex Yukawa detuning](check_detuning.py) · [Results](DETUNING.json)
- [Scalar rescattering and three-energy response](check_rescattering.py) · [Results](RESCATTERING.json)
- [Exact finite-detuning relation](check_closure.py) · [Results](CLOSURE.json)
