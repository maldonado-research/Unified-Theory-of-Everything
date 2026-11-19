# A conditional lower bound for reciprocal lepton-pair scattering

Consider $\mu^-_L\mu^-_L\to e^-_Le^-_L$ and its inverse, in the gaugeless unbroken theory with massless external leptons and common scalar squared mass $a>0$. Retain $C_{ee}$ as the only seed entry and the muon row as the only neutrino Yukawa coupling. At amplitude order $C Y^2$, the two-scalar cut exhausts the discontinuity in this projection. The seed acts on the electron pair in each process. The direct process has just the displayed seed orientation, while the inverse has just the conjugate one. For the two-mass benchmark below, the known discontinuity gives a strictly positive lower bound on the **arithmetic mean of the two cross sections computed from this leading amplitude**. Any Hermitian local four-lepton tree amplitude admitted at this order may be included. [Cut normalization](ON_SHELL_TRIANGLE.md) · [Contraction and channel count](CHANNEL_SCOPE.md)

## The reciprocal amplitude identity

Let $M_{fi}$ be the upper-rim amplitude for initial state $i$ and final state $f$. Form $M_{if}$ by swapping the complete normalized states, including momenta and helicities, using consistent phases. With Hermitian analyticity,

$$
 A_{fi}=\frac{M_{fi}-M_{if}^*}{2i}
       =\frac{\operatorname{Disc}_s M_{fi}}{2i},\qquad
 H_{fi}=\frac{M_{fi}+M_{if}^*}{2}.
$$

Both $A$ and $H$ transform as Hermitian matrices on the two-state space. The discontinuity can therefore have a complex entry between the two different flavors. At fixed kinematics $M_{fi}=H_{fi}+iA_{fi}$ and $M_{if}^*=H_{fi}-iA_{fi}$. The cross terms in their squared magnitudes have opposite signs, giving the exact excess above the bound,

$$
 \frac{|M_{fi}|^2+|M_{if}|^2}{2}-|A_{fi}|^2=|H_{fi}|^2. \tag{1}
$$

A Hermitian local tree contribution with direct value $L$ and inverse value $L^*$ shifts $H_{fi}$ by $L$. Its momentum polynomial has zero discontinuity at the scalar threshold. This statement permits complex couplings and arbitrary such polynomials. At a chosen energy the bound is attained when $H_{fi}=0$. In this same-component channel the reduced primitive is angle independent, so a local term of the same spin structure can enforce this equality at one energy. A different tuning $H_{fi}=-iA_{fi}$ sets the direct transition to zero and leaves $|M_{if}|=2|A_{fi}|$. The identity constrains the reciprocal mean even when one member vanishes.

## Flux, helicity and identical-particle counting

The massless two-body flux in the center of mass is $2s$. Write $k=\sqrt{s}/2$ for the magnitude of the final momentum. Spatial momentum conservation fixes the second final momentum. For the first, the radial measure gives $k^2/(2k)^2$ and $\delta(\sqrt{s}-2k)$ gives the Jacobian $1/2$. Restoring $(2\pi)^{4-6}$ yields $d\Phi_2=d\Omega/(32\pi^2)$. We integrate the full direction sphere for one outgoing momentum, so we include $1/2!$ for the two identical final leptons. Both initial helicities and beam momenta are fixed. This gives $d\sigma/d\Omega=|M|^2/(128\pi^2s)$ in each direction. Rotation invariance allows the inverse initial axis to serve as the beam axis and preserves the scattering angle between the two pair directions. Rotations may rephase the helicity states and leave each squared magnitude invariant.

For sterile weights $w_j=m_j(Y_{\mu j}^*)^2$ in a positive Takagi mass basis and $\rho_j=\operatorname{Disc} K(s;a,m_j^2)/(2\pi i)$, the direct amplitude is $2 C_{ee} W\sum_j w_j K_j/(16\pi^2)$. Here $W=-(u_3^T\epsilon u_4)^*(u_1^T\epsilon u_2)$, $|W|=s$, and reversing the state gives $W^*$. The inverse upper-rim amplitude conjugates $C_{ee}Ww_j$ and uses the same $K_j(s+i0)$. Consequently
$A_{fi}=C_{ee}W\sum_j w_j\rho_j/(8\pi)$. The identical-weak-component sum makes its magnitude independent of the scattering angle. Equation (1) integrates to

$$
 \overline{\sigma} \equiv \tfrac12(\sigma_{\mu\mu\to ee}+\sigma_{ee\to\mu\mu})
 =\frac1{128\pi^2s}\int d\Omega\,|H_{fi}|^2
 +\underbrace{\frac{s|C_{ee}|^2}{2048\pi^3}
                 \left|\sum_j w_j\rho_j\right|^2}_{\displaystyle \sigma_{\min}} .
 \tag{2}
$$

The scalar-intermediate-state count is already included in $\rho_j$. The identical-final-*lepton* factor enters when converting the external amplitude to a cross section. There is one occupied external orientation per flavor direction. In the separate one-active-flavor elastic example, both seed orientations contribute to a single amplitude. Doubling the present flavor-changing amplitude would introduce a factor of four in (2). [General weak and flavor cut](INDEPENDENT_CUT.md)

## Stable-scalar benchmark and dimensions

Take $m_0=\sqrt a$, sterile masses $(m_0,4m_0)$, and $(w_1,w_2)=(w_0,-w_0)$ with $w_0\ne0$. The primitive UV pole vanishes exactly. The scalar is stable against decay to a massless lepton and either of these sterile states in this restricted sector. Define $\xi=s/a$, $\beta=\sqrt{1-4/\xi}$, $g=a C_{ee}w_0$, and $\Delta\rho=\rho(s;a,a)-\rho(s;a,16a)$. With $[C_{ee}]=-3$, $[w_0]=1$ and $[a]=2$ in mass dimensions, $g$ and $A_{fi}$ are dimensionless and $\sigma_{\min}$ has dimension $-2$. The normalized results are

$$
 |A_{fi}|=\frac{\xi |g|\Delta\rho}{8\pi},\qquad
 a\sigma_{\min}=\frac{\xi |g|^2(\Delta\rho)^2}{2048\pi^3}.
$$

The transfer-angle representation shows that $\rho$ strictly decreases with sterile squared mass at each $s>4a$. Its threshold expansion gives $\Delta\rho=(15/34)\beta+O(\beta^3)$. The lower bound opens quadratically,

$$
 \frac{a\sigma_{\min}}{|g|^2}=
 \frac{225}{591872\pi^3}\,\beta^2 + O(\beta^4).
$$

This is the threshold law for the bound. A dispersive contribution to the actual mean can remain nonzero at threshold. At $\xi=8$ the gap is $0.548545873956202579$, giving $|A_{fi}|/|g|=0.174607574705586826$ and $a\sigma_{\min}/|g|^2=0.0000379084757248466264$. The checker integrates the difference of the two on-shell sterile transfer denominators, verifies the reciprocal helicity tensors, and integrates two distinct squared amplitudes with an angle-dependent complex contact. [Executable](verify_inverse_rate.py) · [Receipt](INVERSE_RATE_RESULTS.json)

Equation (2) applies to the stated amplitude and Hermitian analytic contact contributions at the same order. In this restricted flavor sector the transition first appears at $CY^2$, which allows its square to be interpreted as the leading cross section when the couplings are counted accordingly. Additional interactions may contribute new cuts or supply lower-order transitions that interfere with higher amplitudes. A lower bound for a more general measured cross section requires its full discontinuity. The signed interference term in the unitarity cut of the four-lepton amplitude need not be positive.
