# The first three-particle channel in the four-singlet benchmark

19 November 2026 · Candidate.22 · Follows the [scalar-quartic and rate consistency calculation](../radiative_stability/00_READ_FIRST.md) at commit `695b29ad39a9e39e630f3acc4c0a9223b650b9ea`.

**A charged scalar, a muon and the lightest singlet can propagate as a real intermediate state in $\mu^-_L\mu^-_L\to e^-_Le^-_L$.** We calculated its contribution to the imaginary part of the two-loop amplitude. We also obtained the finite real and imaginary contributions from a lepton–singlet loop inserted into either scalar line, and the renormalization-scale logarithm required of the complete two-loop Yukawa amplitude. These quantities answer different questions. The remaining finite diagrams and two-particle cuts are needed for the full next-order rate.

## The physical three-particle channel

The scalar has pole mass $m$. Majorana singlets of masses $m,4m$ couple to the muon doublet and singlets of masses $10m,20m$ to the electron doublet. Their relative squared Yukawas are $(1,-1/4)$ and $(1,-2)$ respectively. External charged leptons are massless, their helicities are fixed, and gauge interactions are switched off. These are the model inputs used for the earlier [complete leading box](../uv_completion/00_READ_FIRST.md).

When $s>4m^2$, enough energy is available to produce the intermediate charged scalar, a massless muon and the singlet of mass $m$. The next such channel opens at $s=25m^2$. The electron-initiated production tree has only the diagrams where a virtual scalar produces the muon–singlet pair. Identical incoming muons in the reciprocal tree allow additional exchange diagrams. The physical discontinuity sums the overlap of both full production trees over three-body phase space.

The additional exchange diagrams increase the magnitude of the angle-averaged ($J=0$) three-particle cut by 9.03%, 15.33% and 20.40% at the three benchmark energies. Divide the amplitude by $\gamma |y_\mu|^2$ with $\gamma=y_e^2(y_\mu^*)^2$ and use the external-spin phase of the leading-box calculation. The negative sign then specifies an interference amplitude between two different flavors.

| $s/m^2$ | Scalar-line part of the $J=0$ cut | Full three-particle $J=0$ cut |
|---:|---:|---:|
| 6 | $-3.99620179\times10^{-7}$ | $-4.35689215\times10^{-7}$ |
| 8 | $-1.26379608\times10^{-6}$ | $-1.45757749\times10^{-6}$ |
| 10 | $-2.45946981\times10^{-6}$ | $-2.96120030\times10^{-6}$ |

The three-particle phase-space calculation and a spin-summed integration check the signed exchange contribution. MadGraph5 separately generates the production helicity amplitudes and verifies the complete $J=0$ cut at $s=8m^2$ by phase-space integration. [Helicity-amplitude calculation](EXACT_MG.md) The result covers the complete three-particle channel for $4<s/m^2<25$ at this order in the Yukawas. The heavier thresholds quoted here are fixed-order propagator thresholds. At higher energies, unstable heavy singlets would require a treatment of their physical decay products.

Near the first threshold the three-particle contribution starts as $(s/m^2-4)^{3/2}$. The virtual scalar approaches its pole at the endpoint where the massless muon is soft. The decay probability vanishes there, and the combination remains integrable. Counting the three particles as featureless phase space would give the wrong threshold power. For $0<s/m^2-4\le0.01$ a separate spin-norm bound limits the relative change due to the additional exchange diagrams to $(2/5)(s/m^2-4)$. For $\Delta=s/m^2-4\downarrow0$, the signed exchange-to-scalar-line ratio is $7\Delta/120+O(\Delta^2)$. The two-particle channel also opens at this energy and has its own contributions. [Three-body phase space](THREE_BODY.md) · [Uniform threshold bound](THRESHOLD_BOUND.md) · [Complete three-particle cut](FULL_THREE_BODY.md)

## A finite part of the two-loop rate

A lepton and a singlet inside the scalar propagator give a continuous distribution of invariant masses. The spectral density is positive, and the charged scalar keeps unit pole residue in the stated subtraction prescription. Integrating one-loop boxes over that distribution gives a finite contribution at all scattering angles. A direct integration in Feynman parameters independently reproduces it. The heavier electron-coupled states contribute to the real amplitude even below their production thresholds. [Propagator and decay width](SPECTRUM.md) · [Direct integration](SIMPLEX.md) · [QCDLoop comparison](EXTERNAL.md)

For $r=m^2\sigma/|y_e y_\mu|^4$, the selected contribution has the form $\Delta r_{\rm lines}=c_\mu |y_\mu|^2+c_e |y_e|^2$. Here $y_e,y_\mu$ are the first Yukawas of their respective blocks; both are independent inputs.

| $s/m^2$ | Leading $r$ | $c_\mu$ | $c_e$ |
|---:|---:|---:|---:|
| 6 | $4.00925\times10^{-11}$ | $1.89542\times10^{-12}$ | $3.85911\times10^{-14}$ |
| 8 | $5.61822\times10^{-11}$ | $2.78740\times10^{-12}$ | $4.84449\times10^{-14}$ |
| 10 | $7.17401\times10^{-11}$ | $3.69883\times10^{-12}$ | $5.58344\times10^{-14}$ |

These terms contribute $+0.0076630\,|y_\mu|^2-0.00032971\,|y_e|^2$ to the earlier endpoint rate ratio $\sigma(10m^2)/\sigma(6m^2)$. The four-rate consistency test from the previous checkpoint also responds to these terms. A measured departure from that relation would therefore need an analysis of Yukawa effects as well as the scalar quartic. [Rate and consistency-test projections](PROJECTION.md)

## What is fixed about the remaining amplitude?

At a fixed scalar pole mass with all four singlets active, running the renormalized sterile mass matrix and Yukawas changes the entire leading rate. Scale independence of the physical result determines the compensating explicit logarithm in the complete two-loop amplitude. For example, the required coefficient of $\ln(\mu/\mu_*)$ in the **fractional** endpoint ratio is $0.0000356549\,|y_\mu|^2-0.000919096\,|y_e|^2$. This coefficient is a check on future finite calculations. The explicit and implicit scale dependence cancel through the stated order. [Derivation and independent matrix finite difference](RG_LOG.md)

A calculation of the remaining finite Yukawa terms must treat the lepton–scalar production vertices, sterile propagators, external residues and counterterms consistently. In particular, adding the complete three-particle cut on top of the scalar-line answer would double count the scalar-line share of that cut. The broader unification hypothesis still needs to determine or constrain the mass and coupling inputs. Prior results on subtracted dispersion, scalar decay widths and phase-space factorization are documented in the [literature and public archive review](SOURCES.md).

## Reproduction

The offline repository replay regenerates the seven numerical checks and compares their records. Reproduction with the two external programs uses [QCDLoop 2.1.0](EXTERNAL.md) and [MadGraph5 3.5.7](EXACT_MG.md). [Full repository replay](../../verification/YUKAWA_CONTINUUM_REPLAY_LINUX.json)

- [Scalar spectral density](check_spectrum.py) · [Results](SPECTRUM.json)
- [Three-body cut of the propagator insertions](check_threebody.py) · [Results](THREE_BODY.json)
- [Threshold bound independent of the exchange sign](check_threshold.py) · [Results](THRESHOLD_BOUND.json)
- [Full three-particle cut](check_full_threebody.py) · [Results](FULL_THREE_BODY.json)
- [Full-amplitude scale logarithm](check_rg.py) · [Results](RG_LOG.json)
- [Finite scalar-line integral](check_continuum.py) · [Results](CONTINUUM.json)
- [Rate and four-rate consistency test](check_projection.py) · [Results](PROJECTION.json)
