# How many charged daughters can the heavy sterile channel supply?

Prepare the heavy muon-coupled sterile species from the neutral wall’s instantaneous readout on the bounce slice. In a subsequent local two-body decay model with $|g|=1$, the full angular tower gives **0.00371–0.00442 charged scalars** and the same number of prompt charged leptons in expectation on the entire sphere by $t=10.5017$. A separate preparation late in expansion gives **0.00143–0.00213 of each** by turnaround. The massive daughters in the late trial retain mean energy of above $2.1068m$ at turnaround under free redshift. The counts use numerically resolved group populations and an analytic bound on the unresolved remainder. [Numerical receipt](CARRIERS.json) · [HDblast charged gas](https://github.com/maldonado-research/HDblast/blob/030a4a0b4c2e65fdbce15fbe76b458a7addffdc9/research/HDBLAST_CHARGED_MATTER_HEATING_20261120/coupled/COUPLED.md)

## Decay and charge bookkeeping

Use $\hbar=c=1$, curvature length $L$, $mL=1$ and $M_2=4m$. The interaction selects $h=Y_{\mu2}=-ig/2$. There are four tree channels, namely $\mu^-\phi^+$, $\mu^+\phi^-$ and the two neutral-component assignments with a neutrino or antineutrino. The charged-conjugate pair of channels carries one half of the width. Every charged event creates two oppositely charged daughters, of which exactly one is massive. The electron row couples to the separate $E_1,E_2$ columns. The exact odd parity stabilizes the lightest odd pole $N_1$ against single-particle decay. It allows reactions involving a pair of odd particles. [Matched model](../../checkpoint_2026_11_19/uv_completion/MATCHING.md) · [Pole ordering and exact parity](../two_particle_cut/POLE_THRESHOLD.md)

For one sign and weak component the spin-averaged squared amplitude is $|h|^2 p_N\cdot p_\ell=|h|^2(M_2^2-m^2)/2$. Multiplying by $p_*/(8\pi M_2^2)$, with $p_*=(M_2^2-m^2)/(2M_2)$, and summing the four channels gives

$$
\Gamma_2=\frac{|h|^2 M_2}{8\pi}\left(1-\frac{m^2}{M_2^2}\right)^2,
\qquad \Gamma_2 L=\frac{225|g|^2}{2048\pi}=0.03497057|g|^2 .
$$

The result uses the vacuum continuum phase space and the leading masses. In the parent rest frame $p_*=15m/8$, $E_\phi^*=17m/8$ and $E_\mu^*=15m/8$. Thus the scalar speed is $15/17$. An isotropic distribution with this one speed has $P_\phi/\rho_\phi=75/289$. The leading partial width for the later scalar decay $\phi^\pm\to\mu^\pm N_1$ is $\Gamma_{\mu N_1}L=|g|^6/(1024\pi^5)$. This transfers the scalar’s charge to a second active daughter. At $|g|=1$ the probability from that channel is bounded by $3.08\times10^{-5}$ during the longer trial below, even for a scalar born immediately and at rest. Additional partial widths were not computed in the pole calculation. [Scalar width](../two_particle_cut/POLE_THRESHOLD.md) · [Heavy Majorana width conventions, Pilaftsis](https://arxiv.org/abs/hep-ph/9707235)

## Two independent preparations and their clocks

At $t_s$ take an incoherent on-shell ensemble with the free spectator readout’s mean $N_{2n}=2(n+1)(n+2)|\beta_n(t_s)|^2$ for each retained shell. Take the increasing-cutoff limit after forming the decay observables below. Propagate this cohort as dilute, unpolarized particles with independent vacuum decays on the prescribed neutral radius. Omit further production, inverse processes, environmental forces, collisions and backreaction during this second stage. For $k=n+3/2$ and dimensionless wall time, $E_nL=\sqrt{(M_2L)^2+k^2/a^2}$, and the elapsed parent proper time in units of $L$ is $s_n(t)=\int_{t_s}^t M_2/E_n(u)\,du$. The expected number of prompt charged scalars equals the number of prompt charged leptons,

$$
N_{\phi^\pm}=N_{\mu^\pm}=\tfrac12\sum_n N_{2n}\,[1-e^{-(\Gamma_2L)s_n(t_f)}]
\le \tfrac12 N_2(t_s)\,[1-e^{-(\Gamma_2L)(t_f-t_s)}].
$$

Here each count sums both signs. The right side evaluated at $|g|=1$ bounds the entire range $0\le |g|\le1$ in this model. The underlying spherical quantum calculation supplies one-particle expectations on a slice. An interacting field calculation on the compact wall could have discrete final states, finite-time interference and inverse reactions that change the inferred conversion. Even within a kinetic description, continuous on-shell injection requires integrating every birth rate against its own subsequent decay probability. Free instantaneous occupation can decrease through coherent evolution, so its time derivative does not give that rate. [Mode occupations at each slice](MODES.json)

The following entries retain shells 0 through 32.

| Prepared neutral-wall interval (time in units of $L$) | Prepared $N_2$ | Maximum charged events per parent from the rest clock | Charged events using shell clocks | Scalar / lepton energy sums at birth, in $1/L$ |
|---|---:|---:|---:|---:|
| Bounce 0.867705 to 10.5017 | 0.0296677 | 0.143012 | 0.00333162 | 0.00951477 / 0.00839538 |
| 11.03911 to turnaround 12.42163 | 0.0492802 | 0.0235987 | 0.00103822 | 0.00247779 / 0.00218629 |

The first endpoint copies the clock time chosen for mobility onset in the separate charged-fluid history. Its local metric and density on the neutral trajectory are those of the spectator experiment. The second preparation occurs after that onset. The two rows count different prepared cohorts and cannot be added. The energy entries sum daughter energies locally at their individual creation times. In an unpolarized ensemble a boost gives mean daughter energy $x_i E_n$ at birth, where $x_\phi=17/32$, $x_\mu=15/32$. Thus the time integral for each entry is $x_i (\Gamma_2 L)(M_2 L)\sum_n N_{2n}\int e^{-(\Gamma_2L)s_n(u)}du/2$, in units of $1/L$. [Clock and energy integration](check_carriers.py)

For a freely propagating massive daughter born at $a_b$, set $r=a_b/a_f$. Its final energy is $E_f=\sqrt{r^2 E_b^2+(1-r^2)m^2}$, convex in $E_b$ during expansion. The angle average of $E_b$ is at least $E_\phi^*$. Monotone expansion gives $r\ge a_s/a_f$, so the full prompt scalar birth cohort, treated as persistent test particles, has $\langle E_f\rangle/m\ge\sqrt{1+(15a_s/8a_f)^2}$. The late interval expands only from 9.63691 to 9.74401 and yields 2.10684. The massless charged leptons have an isotropic pressure of one third their energy density and cannot carry a nonzero cold rest density. Neither a thermalization rate nor electromagnetic energy transfer follows from the vacuum decay width.

## Conversion across the full angular tower

For fixed endpoints the charged decay fraction $f_n=(1-e^{-(\Gamma_2 L)s_n(t_f)})/2$ decreases with shell number because larger $k$ shortens the parent clock. A group $l\le n\le h$ with positive total population $C$ therefore contributes between $Cf_h$ and $Cf_l$ charged scalars. Apply this to the numerical groups 33–128, 129–512 and 513–2048. Above them, $f_{2049}N_{n\ge2049}^{\rm upper}$ bounds the entire additional yield. This last contribution is at most $2.07\times10^{-7}$ for the bounce preparation and $1.66\times10^{-6}$ for the expanding preparation. [Computed group populations and mesh tests](BAND.md) · [Analytic remainder](TAIL.md)

| Independent preparation | All-shell charged scalars, or prompt charged leptons | All-shell scalar birth energy $LE_\phi$ | All-shell lepton birth energy $LE_\mu$ |
|---|---:|---:|---:|
| Bounce to reference onset time | 0.0037135–0.0044168 | 0.014767–0.015021 | 0.013029–0.013254 |
| Expanding slice to turnaround | 0.0014363–0.0021259 | 0.0067186–0.0069403 | 0.0059281–0.0061238 |

For energy, each parent’s expected contribution from both daughters of its charged branch is $U_n=(\Gamma_2 L)(M_2L)\int e^{-(\Gamma_2L)s_n(t)}dt/2$. Shorter clocks increase the survival factor, so $U_n$ **increases** with $n$ and is bounded by $(\Gamma_2L)(M_2L)\Delta t/2$. Applying the group endpoints in this order, bounding the ultraviolet contribution by this constant times the tail number ceiling, and multiplying by each daughter’s energy fraction gives the table. Time dilation suppresses the decay probability at high momentum, leaving a finite expected energy emitted during these finite intervals even though the bare instantaneous parent energy diverges on the expanding preparation slice. There the result is a limit of cutoff decay observables. Preparing an ensemble with a finite total energy would require a physical state or readout prescription.

The table uses $|g|=1$. Its number upper endpoints also bound the weaker couplings $0\le|g|\le1$. The group inequalities apply to exact group populations, and the final tail bound applies to the specified interpolated geometry. The resolved occupations remain numerical. The high-band mesh comparison gives less than $3\times10^{-9}$ in multiplicity-weighted absolute probability change, and the endpoint decay probabilities change by less than $10^{-9}$ under a second clock grid. These convergence checks quantify the resolved numerical terms independently of the group widths shown. [Per-group conversion intervals](CARRIERS.json) · [Mode solver checks](MODES.md)

## Necessary number and charge normalization

The HDblast starting seed has total cold rest density $D_i=10^{-4}$ in tension units at radius $a_i=2$. Its conserved individual-carrier count is $N_{\rm seed}=2\pi^2a_i^3 D_i S_0/(mL)=0.01579137 S_0$ for a common mass $mL=1$, where $S_0=\sigma L^4/(\hbar c)$ in ordinary units. This normalization already includes every stream and both signs. If the massive scalar daughters of the chosen one-species preparation alone must supply that count in expectation, even perfect survival and use of every scalar require $S_0\le0.21098$ for the first low-band row or $S_0\le0.065746$ for the second. Including the full angular tower relaxes these necessary upper bounds to **0.27970** and **0.13463** respectively. Filling all 26,180 available parent one-particle states in shells 0–32 at the initial slice of each trial would relax the two rest-clock ceilings only to $2.3710\times10^5$ and $3.9124\times10^4$, respectively. The fluid benchmark has its precursor from the beginning and assigns separate pairs to three independent Abelian fields. A one-species number budget is only one necessary condition for that assignment. [HDblast normalization](https://github.com/maldonado-research/HDblast/blob/030a4a0b4c2e65fdbce15fbe76b458a7addffdc9/research/HDBLAST_CHARGED_MATTER_HEATING_20261120/synthesis/BRIDGE.md)

Requiring the same scalar to reproduce the projected $Q=|q|/(\sqrt3\,mL)=10^8$ fixes its dimensionless single-field charge at $|q|=1.73205\times10^8$. Under the HDblast local wall kinetic approximation, $\alpha_{\rm wall}=q^2/(4\pi S_0 K)$, so $\alpha_{\rm wall}<1$ requires $S_0>2.3873\times10^{18}$ at $K_i=0.001$ or $S_0>1.8086\times10^{11}$ even at the late $K=13200.04$. Even at the late perturbativity threshold the required count exceeds the full-tower upper yield by factors greater than $6.46\times10^{11}$ and $1.34\times10^{12}$ for the two trials. The charge is an added matching assumption for the gaugeless Yukawa benchmark. Gauge dynamics, an external supply of carriers, annihilation, a changed mass spectrum or different couplings require separate calculations. [HDblast local gauge normalization](https://github.com/maldonado-research/HDblast/blob/030a4a0b4c2e65fdbce15fbe76b458a7addffdc9/research/HDBLAST_CHARGED_MATTER_HEATING_20261120/synthesis/BRIDGE.md)

After generating the mode, band and tail receipts, run `python -I -B -O check_carriers.py --output CARRIERS.json`. The check independently compares the spin-trace width with the closed form, brackets each shell’s clock and birth energy between their endpoint limits, and repeats the integration on a refined time grid. [Checker](check_carriers.py) · [Receipt](CARRIERS.json)
