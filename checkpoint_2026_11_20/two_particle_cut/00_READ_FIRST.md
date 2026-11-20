# Two-particle unitarity and the threshold cusp

This checkpoint completes the finite pure-Yukawa absorptive part of polarized $\mu^-_L\mu^-_L\to e^-_Le^-_L$ at two loops in the four-singlet benchmark for $4<s/m^2<25$. The scalar pair contributes with opposite sign to the lightest three-particle state in the muon-Yukawa coefficient. At $s=8m^2$ the latter cancels 39.67% of the former. The fraction reaches one half at $s/m^2=9.84984$. All values use the stated Yukawa texture and $\overline{\rm MS}$ parameters at $\mu=m$, with $m$ the scalar pole mass. [Combined derivation, normalization and energy scan](ABSORPTIVE.md) · [Numerical results](ABSORPTIVE.json)

The fixed-order discontinuity predicts the leading square-root cusp and the next nonanalytic term in the **real** amplitude approaching the pair threshold from below. For $\eta=\sqrt{4m^2/s-1}$, the coefficient of the $\gamma|g|^2$ contribution to the real angular mean contains $-1.31099848\,10^{-6}\eta+2.29468998\,10^{-6}\eta^3$. The three-body state first contributes at the cubic power. The smooth real terms require additional information. [Threshold derivation and checks](CUSP.md) · [Cusp coefficients](CUSP.json)

Converting the lightest sterile mass to a pole mass reveals a limit on the common-threshold interpretation. At the chosen scale its pole lies below the scalar by $m|g|^2/(16\pi^2)$ at first order. This opens a narrow charged-scalar decay to a muon and the lightest sterile. The coincident cusp describes the fixed-order expansion. Resolving the two nearby thresholds requires including the mass shift and, at the width scale, the scalar lifetime. [Pole mass and partial width](POLE_THRESHOLD.md)

| Calculation | What can be checked |
|---|---|
| [Majorana self-energy](STERILE.md) | Internal mixing, both weak components, signs from inverting the propagator, and a change of sterile basis |
| [Proper Yukawa vertex](VERTEX.md) | A massive spacelike triangle, both vertex insertions, and its discontinuity in the external scalar invariant |
| [Crossed lepton–scalar channel](INDEPENDENT.md) | Reconstruction of the vertex by a second dispersion relation and an independent calculation of the two-particle angular mean |
| [External scalar and lepton residues](EXTERNAL.md) | Finite normalization factors and cancellation of the explicit scale logarithm against the running leading amplitude |
| [Two- and three-particle cut](ABSORPTIVE.md) | Angular sewing, the effect of their cancellation on energy dependence, and the absorptive contribution to the rate |
| [Threshold cusp](CUSP.md) | The first two nonanalytic terms and their continuation below threshold |
| [Pole-mass conversion](POLE_THRESHOLD.md) | Position of the lightest sterile pole, phase space for scalar decay and the range of the common-threshold expansion |
| [Literature and research inventory](SOURCES.md) | Primary-source comparisons, the author’s public Zenodo records and the HDblast context |

The calculation assumes the unbroken gaugeless theory with massless external charged leptons. Connecting the benchmark to physical masses and couplings and deriving its parameters from the unification hypothesis are separate research questions. The next amplitude calculation is the remaining finite dispersive contribution. Reconstructing it away from threshold requires accounting for higher-energy and crossed-channel cuts and fixing any subtraction terms, or independently integrating the remaining two-loop diagrams. [Earlier finite dispersive subset](../../checkpoint_2026_11_19/yukawa_continuum/SIMPLEX.md) · [Verification instructions](../../docs/VERIFICATION.md)
