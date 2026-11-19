# A lower bound on the mean of the two reciprocal rates

Consider the flavor-isolated process `μ⁻_L μ⁻_L → e⁻_L e⁻_L` and its inverse in the unbroken theory. Take all four external leptons massless, with the same charged weak component. Only the electron seed and muon Yukawa row are present. The scalar-pair cut then exhausts the discontinuity of the leading flavor-changing amplitude at order `C Y²`. The two-mass example below has a stable scalar in this sector. [Channel and cut conventions](ON_SHELL_TRIANGLE.md) · [Completeness at this order](CHANNEL_SCOPE.md)

Let $M_{fi}$ and $M_{if}$ be the upper-rim amplitudes with initial and final two-particle states interchanged, including their momenta and helicities. Use the same state phases in both amplitudes. Define

$$
H_{fi}=\frac{M_{fi}+M_{if}^*}{2},\qquad
D_{fi}=\frac{M_{fi}-M_{if}^*}{2i}.
$$

Hermitian analyticity identifies $D_{fi}$ with the discontinuity divided by $2i$. At each scattering angle the parallelogram identity gives

$$
\frac{|M_{fi}|^2+|M_{if}|^2}{2}
=|H_{fi}|^2+|D_{fi}|^2 \geq |D_{fi}|^2 .
$$

This identity permits arbitrary complex flavor couplings. A finite Hermitian local four-lepton term of the same amplitude order adds to $H_{fi}$. Its reverse matrix element is the complex conjugate of its direct one, so it drops out of $D_{fi}$. The bound is on the mean of the two rates. Destructive interference can reduce one direction separately.

With fixed initial helicities, the massless two-body flux and final-state measure give $d\sigma/d\Omega=|M|^2/(128\pi^2s)$ for either process. This expression integrates over the full sphere and includes the $1/2!$ for the two identical outgoing leptons. The incoming beam labels have fixed momenta. A spatial rotation puts the inverse initial momenta along the beam axis while preserving the scattering angle and the squared matrix element.

Write $w_i=m_i(Y^*_{\mu i})^2$ in the Takagi basis and $\rho_i=\operatorname{Disc} K(s;a,m_i^2)/(2\pi i)$. The same-component angular sum removes the scattering-angle dependence of the cut magnitude. The ordered spin factor has magnitude $s$, giving $|D_{fi}|=s|C_{ee}|\,|\sum_i w_i\rho_i|/(8\pi)$. Integration gives

$$
\frac{\sigma_{\mu\mu\to ee}+\sigma_{ee\to\mu\mu}}2
\ \geq\ \frac{|D_{fi}|^2}{32\pi s}
=\frac{s|C_{ee}|^2}{2048\pi^3}
\left|\sum_i w_i\rho_i\right|^2 .
$$

For $a=m_0^2$, sterile masses $(m_0,4m_0)$ and weights $(w_0,-w_0)$, the last squared factor is $|w_0|^2(0.548545873956202579\ldots)^2$ at $s=8m_0^2$. Here $w_0$ has dimension one and $C_{ee}$ dimension minus three, yielding the cross-section dimension minus two. This is the contribution obtained by squaring the leading restricted amplitude, since neither direction has a lower-order flavor-changing tree. Higher-order amplitudes or other permitted on-shell intermediate channels in an enlarged model can change the total discontinuity coherently. A prediction for electroweak scattering would require those amplitudes and its physical particle spectrum. [Benchmark and subtraction](ON_SHELL_TRIANGLE.md)
