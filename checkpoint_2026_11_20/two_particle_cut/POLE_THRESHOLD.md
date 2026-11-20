# The light sterile pole and the opening of charged-scalar decay

The benchmark identifies the lightest **running** sterile mass $M_1(m)$ with the scalar pole mass $m$. In the unbroken gaugeless four-singlet theory this gives a light sterile pole below the scalar. Write $G=|g|^2$ and $\alpha_N=G/(16\pi^2)$. To first order in the squared Yukawas,
$$
 M_{1,\mathrm{pole}}=m(1-\alpha_N)+O(mY^4).
$$
The charged scalar can therefore decay to $\mu^-N_1$. Its leading two-body partial width is $\Gamma_{\mu N_1}/m=G^3/(1024\pi^5)$ under this input prescription. These statements use a joint perturbative expansion in the Yukawas at fixed ratios. Once the decay opens, $m$ denotes the real scalar resonance mass to this accuracy. [Model inputs](../../checkpoint_2026_11_19/uv_completion/MATCHING.md) · [Numerical checks](POLE_THRESHOLD.json)

## Pole sign and scale dependence

Set $z=p^2$, $a=m^2$ and $t=z/a$. On the real sheet $0<t<1$, the inverse-kinetic bubble of the sterile fields is proportional to
$$
 f(z;\mu)=-2\int_0^1(1-x)\ln\frac{(1-x)(a-xz)}{\mu^2}\,dx
 =\ln\frac{\mu^2}{a}+2-\frac1t
   -\frac{(1-t)^2}{t^2}\ln(1-t).
$$
The factor two sums the internal weak doublet. The integrable endpoint gives $f(a;\mu)=1+\ln(\mu^2/a)$, with a finite first derivative in $z$. For either spin at rest, put $H=Y^\dagger Y$, $h=fH/(16\pi^2)$. The inverse propagator in the two chiral sectors is
$$
 \begin{pmatrix}
 -M^\dagger & \sqrt z\,(1+h)\\
 \sqrt z\,(1+h^T) &-M
 \end{pmatrix}.
$$
The plus sign of $h$ follows directly from the proper inverse-kernel loop and its kinetic subtraction. The full Majorana mass coefficient has zero proper one-loop correction in these coordinates. [Spacelike propagator derivation](STERILE.md) · [Mass and field counterterms](../../checkpoint_2026_10_16/yukawa_counterterms/YUKAWA_TWO_POINT.md)

The light root is nondegenerate in sterile flavor. Both off-diagonal insertions and changes in its eigenvector first affect its eigenvalue at second order in $h$ when the reference mass matrix is diagonal. The light pole condition to first order is $\sqrt z[1+G f(a;\mu)/(16\pi^2)]=M_1(\mu)$. Consequently $M_{1,\mathrm{pole}}=M_1(\mu)[1-G f(a;\mu)/(16\pi^2)]+O(mY^4)$ near the benchmark. The one-loop mass flow gives $d M_1/d\ln\mu=2GM_1/(16\pi^2)$ at the reference point, while $\partial f/\partial\ln\mu=2$. These derivatives cancel in the pole mass at this order. The equality $M_1(m)=m$ is imposed only at the reference scale. [Full mass-matrix running](../../checkpoint_2026_11_19/radiative_stability/RUNNING.md)

## Open channel and threshold scales

At invariant mass $m$, the fixed charged scalar selects one active weak component. The spin sum at tree level is $G(m^2-M_{1,\mathrm{pole}}^2)$ and the daughter momentum has magnitude $(m^2-M_{1,\mathrm{pole}}^2)/(2m)$. These give the first width induced by the splitting,
$$
 \frac{\Gamma_{\mu N_1}}m=
 \frac G{16\pi}\left(1-\frac{M_{1,\mathrm{pole}}^2}{m^2}\right)^2
 +\text{higher-order corrections}
 =\frac{G^3}{1024\pi^5}+\text{higher-order corrections}.
$$
The spin-sum and phase-space zeros each supply a power of the small mass gap. The width-to-gap ratio from this channel is $G^2/(64\pi^3)$ at leading order. [Scalar two-point optical theorem](../../checkpoint_2026_11_19/yukawa_continuum/SPECTRUM.md)

In this field content, $\phi\to-\phi$, $N_i\to-N_i$ with all active leptons even is an exact symmetry of the Yukawa and quadratic action and of a scalar quartic. With the perturbative spectrum above, $N_1$ is the lightest odd particle and is stable. This symmetry forbids scalar decay to active leptons alone. A simultaneous sign flip of the electron lepton row and its two sterile fields also forbids the open two-body assignment $e^-N_1$ in the block-diagonal benchmark. The other sterile columns are heavier than $m$, and the other components of the scalar multiplet are degenerate in the unbroken theory. A final state with one $N_1$ and additional active leptons requires an odd number of active fermions by fermion parity, beginning with one or three. The displayed coefficient specifies the two-body partial width. Multilepton partial widths have not been integrated here.

For $\sqrt{s}$ there are three distinct scales to this accuracy.

| Configuration | Energy divided by $m$ | Interpretation |
|---|---:|---|
| Two scalar resonances | $2$ | Narrow-resonance pair energy |
| One scalar resonance, $\mu^-$ and $N_1$ | $2-\alpha_N+O(Y^4)$ | Narrow-resonance production energy |
| $\mu^-\mu^-N_1N_1$ | $2-2\alpha_N+O(Y^4)$ | Threshold of this stable four-particle channel |

The fixed-order cuts through sixth degree in Yukawas use free masses in the cuts and a scalar pole mass prescription. Their two- and three-particle endpoints coincide at $2m$. The pole calculation resolves a separation of order $m\alpha_N$. Within that interval, expanding the moving three-particle endpoint at fixed $s$ is nonuniform. The much smaller partial width sets a resonance energy scale of order $mG^3/(1024\pi^5)$ from this mode. For nonzero width a stable-state unitarity sum contains the decay products of a scalar resonance. The two resonant configurations in the table then feed the common stable four-particle state. Higher Yukawa orders also allow massless two-lepton intermediate states through pairs of loop-induced four-lepton amplitudes. The table concerns the resonances and their daughter channel. The [fixed-order threshold cusp](CUSP.md) applies to the specified perturbative coefficients.

## Retuning the running mass to impose pole equality

An alternative benchmark holds the other three sterile masses and all Yukawas fixed and sets $M_1(m)=m(1+\alpha_N)$ to enforce $M_{1,\mathrm{pole}}=m+O(mY^4)$. This is a finite change of physical inputs. At fixed $s$ it adds the corresponding mass derivative of the leading box to the sixth-degree coefficient. For the reduced $J=0$ cut, use $m=1$ and $\beta=\sqrt{1-4/s}$. Differentiation of $M_1/(M_1^2+q)$ at $M_1=1$ gives the shift per $G$ in a single muon-side tree,
$$
 \delta F_\mu(q)=\frac{q-1}{16\pi^2(1+q)^2}.
$$
The positive transfer $q$ is uniform between $(1-\beta)/(1+\beta)$ and $(1+\beta)/(1-\beta)$. Its angular mean yields
$$
 \Delta D_{0,\mu}^{(6)}
 =\frac{s\beta}{8\pi}\overline F_e\,
   \frac{(1-\beta^2)(\operatorname{atanh}\beta-\beta)}
        {32\pi^2\beta}
 =\frac{\overline F_e(\operatorname{atanh}\beta-\beta)}{64\pi^3}.
$$
The full amplitude change is $\gamma G\,\Delta D_{0,\mu}^{(6)}$. The shifts at $s/m^2=6,8,10$ are $-5.94630821204\times10^{-8}$, $-1.88407110331\times10^{-7}$ and $-3.64692049048\times10^{-7}$. The tree mass derivative vanishes at $q=1$, giving zero shift in the linear muon cusp coefficient. The cubic coefficient shifts by $F_e(1)/(192\pi^3)=-1.24424143976\times10^{-7}$. This comparison expands at fixed energy before taking the threshold limit of the coefficients. The leading three-particle cut already carries the extra Yukawa pair, so its mass derivative first enters one order later at fixed energy. [Angular sewing and normalization](../../checkpoint_2026_11_19/uv_completion/CUT.md)

Run `python -I -B -O check_pole_threshold.py --output POLE_THRESHOLD.json`. The checker evaluates the original bubble at three scales, brackets the real zero of the full mixed chiral determinant below $m$, and checks the mass gap and spin-trace width as the coupling decreases. Advancing the full complex symmetric mass matrix to two other scales leaves the pole invariant through first order. Higher powers generated by the numerical determinant serve only to check convergence of its first-order coefficient. Independent symmetric finite differences of the leading angular integral verify all three finite-input shifts. [Checker](check_pole_threshold.py) · [Receipt](POLE_THRESHOLD.json)
