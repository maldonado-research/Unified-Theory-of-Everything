# An on-shell four-lepton cut and what it reveals about sterile masses

19 November 2026. Provisional GitHub candidate.18, following the [off-shell four-lepton benchmark](../four_lepton/00_READ_FIRST.md) at commit `c5124b723043fa2c7fc73ebbd852035691721f95`.

**The finite effect survives putting the external leptons on shell.** Two sterile masses with equal and opposite ultraviolet pole weights give a nonzero discontinuity in a four-lepton scattering amplitude. The discontinuity is the difference between the two boundary values across the energy threshold where a pair of scalar lines can go on shell. A local four-lepton counterterm is polynomial in momenta and drops out of that difference.

The result uses the unbroken, gaugeless theory with zero scalar background, massless external left-handed doublets and a common positive scalar squared mass a. The selected order has one derivative interaction and two neutrino Yukawa vertices. Counting the field slots exhausts the connected contractions at this monomial. Only the two scalar edges can carry the physical s-channel cut. The sterile propagator exchanged across the tree subamplitude stays spacelike. We evaluate the cut both by integrating the loop and by sewing the on-shell tree amplitudes, including exchange of identical particles. [Loop and dispersion calculation](ON_SHELL_TRIANGLE.md) · [Independent two-body cut](INDEPENDENT_CUT.md) · [Scope of the channel](CHANNEL_SCOPE.md)

## A conditional bound on two reciprocal rates

For two identical-component left-handed lepton beams, take the flavor-changing channel μ⁻μ⁻→e⁻e⁻ and its inverse. The selected flavor texture has no lower-order amplitude for this transition. Let A be the anti-Hermitian part of the leading amplitude, fixed here by the two-scalar cut. For each pair of external states the average squared amplitude of the two reciprocal processes is the sum of two nonnegative terms, one of which is |A|². A Hermitian local four-lepton interaction added to this order can change the other term. The identical-final-state factor and the massless flux give a lower bound on the mean exclusive cross section, `σ_mean ≥ |A|²/(32πs)`, for this isotropic cut magnitude.

For scalar squared mass a, sterile masses √a and 4√a, opposite pole weights of magnitude |w₀|, and `s/a=8`, define the dimensionless coupling combination `g=a C_ee w₀`. Then **`a σ_mean ≥ 3.79085 × 10⁻⁵ |g|²`** at the square of the leading amplitude, in natural units. The two-particle scalar decay to a sterile and a lepton has zero phase space at the lighter mass and is closed at the heavier mass. The bound applies to this restricted weak-coupling calculation below its EFT cutoff. A numerical cross section in physical units requires the scalar mass and the couplings. Extra absorptive channels change A. The bound on the mean permits a cancellation in either direction individually. [Rate derivation and checker](INVERSE_RATE_BOUND.md) · [Independent normalization audit](INVERSE_RATE_AUDIT.md)

## Two cases with different physical implications

**Sterile masses 1 and 4, scalar mass 1.** Choose the derivative coupling on the outgoing e-flavor pair and Yukawa couplings on the incoming μ-flavor pair. The weights of the two mass groups can then be +1 and −1. Their ultraviolet sum is zero. For the common weak component, write K for the mean of the two direct/exchange kinematic kernels, with the overall spin, Wick multiplicity and loop normalization stripped. The discontinuity `ρ=Disc K/(2πi)` at s=8 for the weighted sum is **+0.54854587396**. It opens at s=4 with coefficient `lim ρ/β=15/34`, where `β=√(1−4/s)`. At s=8 the once-subtracted upper-rim value of the same sum is **1.18998510707 + 1.72330768778 i**. Its real part is relative to the chosen subtraction at s=0. The stated masses admit no open two-body scalar decay through the selected Yukawa interactions. [Conventions and numerical receipt](ON_SHELL_TRIANGLE_RESULTS.json)

**Sterile masses 1/2, 1 and 2, scalar mass 1.** With weights +1, −2 and +1, the equal-component discontinuity cancels **at every energy** above threshold. Both scalar lines crossing the cut have the same mass in every term. Their angular integration gives an identity between the kernels for squared masses b, a²/b and a. An external pair with distinct weak components retains an additional term odd in the scattering-angle cosine. At s=9 and `cosθ=1/2` the two weak routings give **−0.04010098542 and +0.04010098542**, including phase space and stripping the common factor `W/(16π)`, where W is the ordered external spinor contraction. These numbers describe signed contributions to the cut. [Exact identity, proof and numerical checks](MASS_IDENTIFIABILITY.md)

There are two ways to distinguish the mass groups in this calculation. If every coupled sterile mass is at least the scalar mass, the equal-component cut on an open energy interval already separates the positive-mass group weights. Without that restriction, the angle-even and angle-odd cut functions together separate them. States with exactly equal masses still enter only through the sum of their weights. This conclusion assumes access to the signed or complex amplitude functions. Inferring them from scattering data requires an observable and an uncertainty analysis.

The three-mass cancellation necessarily puts one coupled sterile below the scalar mass. That opens scalar decay into a sterile plus a massless lepton. The three-mass result uses free scalar propagators at the declared perturbative order. A physical line shape for that spectrum needs the decay products and a consistent treatment of the width. These examples establish a result in a specified neutrino effective theory. Extending it to a physical electroweak scattering rate also requires electroweak symmetry breaking, the remaining interactions and matched local coefficients.

## Sources and verification

The [literature review](SOURCES.md) compares the calculation with published cutting rules, a Majorana exchange amplitude, and a 17 November preprint on normal versus anomalous triangle thresholds. It rechecks the [same-day publication survey](../four_lepton/PUBLICATION_SURVEY.md) against the author's Zenodo catalog and the HDblast site and repository.

The four standalone checks can be run from this directory in the pinned environment.

```sh
python -I -B -O verify_on_shell.py
python -I -B -O independent_cut.py
python -I -B -O verify_mass_identifiability.py
python -I -B -O verify_inverse_rate.py
```

The [replay guide](../../docs/VERIFICATION.md) and [full receipt](../../verification/ON_SHELL_CUT_REPLAY_LINUX.json) cover the repository-wide checks. The September N00AK-r1 record remains the archived baseline. [Zenodo](https://doi.org/10.5281/zenodo.22319172)
