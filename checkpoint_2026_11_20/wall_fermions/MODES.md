# Two Majorana spectra on the neutral closed wall

Starting in the instantaneous vacuum at wall proper time zero, one free Majorana species of mass $M=1/L$ has an expected **0.097367822 quanta** in shells $0\le n\le32$ at the first turnaround. For $M=4/L$ the same calculation gives **0.045241176 quanta**. Both numbers count the entire spatial three-sphere. The lowest shell alone accounts for 0.00216051 and 0.0000142900 quanta, respectively. The full shell distributions matter especially for the heavier species. [Mode probabilities and totals](MODES.json)

## Background and state

Assume the wall supports an intrinsic closed FLRW metric $ds^2=L^2[dt^2-a^2(t)d\Omega_3^2]$ with $c=\hbar=1$. The input supplies $t$, $a$, and the dimensionless proper-time Hubble rate $H$ from the **neutral parent** of HDblast's charged-matter study. The calculation uses the wall proper time column $t$. The archive also carries a separate displaced source clock. Trace fields outside these three columns enter only through their historical effect on $a$. Setting the TOE scalar mass parameter to $mL=1$ identifies our two free reference masses with that model's $N_1,N_2$ mass choices $(m,4m)$ at leading order. HDblast supplies no measurement of the dimensionless product $mL$. [Pinned neutral background trace](https://github.com/maldonado-research/HDblast/blob/030a4a0b4c2e65fdbce15fbe76b458a7addffdc9/research/HDBLAST_CHARGED_MATTER_HEATING_20261120/background/trace.csv.gz) · [TOE singlet action and mass ordering](../../checkpoint_2026_11_19/uv_completion/MATCHING.md)

The 29,558 archived rows start at $a=2,H=0$ and initially contract. Locate the first bounce and the subsequent turnaround from the sign changes of a shape-preserving cubic interpolant through $H$. Use a cubic Hermite radius through the archived radii with slopes $aH$. Insert both roots at the radius of that initial Hermite interpolant and impose zero slope at each root. Truncate the evolution at the turnaround. This gives the radii and slopes required for the endpoint projectors and the tail analysis. The slice at $t=11.03911$ samples the still-expanding neutral solution at a time chosen for comparison with the separate charged run. It has $a=9.63690865$ and $H=0.01617707$ on the neutral parent. [Geometry data](trace_geometry.csv.gz) · [Charged-run interval](https://github.com/maldonado-research/HDblast/blob/030a4a0b4c2e65fdbce15fbe76b458a7addffdc9/research/HDBLAST_CHARGED_MATTER_HEATING_20261120/coupled/COUPLED.md)

For an intrinsic, minimally coupled Majorana field of constant mass $\mu=ML$, rescale the field by $a^{3/2}$. Write $k_n=n+3/2$ and use one of the two equivalent spatial spinor sectors. The two-component evolution, momentum and positive energy in wall proper time are

$$
i\partial_t u_n=\left(\mu\sigma_z+p_n\sigma_x\right)u_n,\qquad
p_n=k_n/a,\qquad \omega_n=\sqrt{\mu^2+p_n^2}.
$$

Define $\theta_n=\tfrac12\arctan2(p_n,\mu)$, $u_-=(-\sin\theta_n,\cos\theta_n)^T$ and $u_+=(\cos\theta_n,\sin\theta_n)^T$. Each mode starts in $u_-(0)$. Its occupation on a readout slice is $P_n(t)=|u_+(t)^\dagger u_n(t)|^2$. Both the initialization and readout use this instantaneous diagonalization. Changing the initial state or the readout basis can change the result. In particular a stationary point of $a$ fixes its first derivative on one slice and leaves the evolution on neighboring slices dynamical. [TOE fermion mode convention](../../checkpoint_2026_09_06/fermion_portal/INDEPENDENT_NUMERICAL_REVIEW.md) · [Closed-space basis and contemporary analytic comparison](SOURCES.md)

For each sign of the two-component Dirac operator on the unit sphere the eigenvalue magnitude $k_n$ has multiplicity $(n+1)(n+2)$. The two sectors give equal occupation since reversing the sign of $p_n$ is a constant unitary conjugation. The positive-energy Majorana quasiparticles therefore carry $g_n=2(n+1)(n+2)$ independent states per shell, including both spin sectors. Charge conjugation fixes the negative-energy partners. A flat-space normalization check gives $g_n\,dn/(2\pi^2a^3)\to p^2 dp/\pi^2$, the phase-space measure for two spin states. The per-species total and excitation energy are

$$
N_{\le32}=\sum_{n=0}^{32} g_n P_n,\qquad
E_{\le32}L=\sum_{n=0}^{32}g_n\omega_nP_n.
$$

The number density in $L^{-3}$ and energy density in $L^{-4}$ follow by dividing these dimensionless sums by $2\pi^2a^3$. Finite-band excitation energy is defined by the stated particle basis. A covariantly renormalized stress tensor requires a separate prescription.

## Evolution and shell dependence

| Neutral-parent slice | $t$ | $a$ | $N_{\le32}$, $\mu=1$ | $N_{\le32}$, $\mu=4$ |
|---|---:|---:|---:|---:|
| Stationary initial slice | 0 | 2 | 0 | 0 |
| Bounce, minimum radius | 0.86770530 | 1.88094812 | 0.011851781 | 0.029667744 |
| Maximum archived expansion rate | 4.41550 | 4.71522376 | 1.1278049 | 1.6133861 |
| Expanding readout | 11.03911 | 9.63690865 | 0.10652915 | 0.049280197 |
| Turnaround, maximum radius | 12.42163330 | 9.74401259 | 0.097367822 | 0.045241176 |

The initial vacuum attains the minimum possible number, zero. Searching 4,097 equally spaced proper-time samples and refining all sampled local maxima gives peak totals of 1.3808335 at $t=4.65321$ for $\mu=1$ and 1.7135556 at $t=4.51837$ for $\mu=4$. Coherent evolution in a time-dependent particle basis allows these totals to decrease. The difference between the two late readouts includes interference with the amplitude already present. A fresh-vacuum production experiment on the short late segment would require initialization at its own starting slice.

| Turnaround quantity | $\mu=1$ | $\mu=4$ |
|---|---:|---:|
| $P_0$ | $5.4012715\times10^{-4}$ | $3.5724986\times10^{-6}$ |
| Shell with largest $P_n$ | 0 | 5 |
| Shell with largest $g_nP_n$ | 4 | 14 |
| $N_{\le8}$ | 0.07423777 | 0.00716379 |
| $N_{\le16}$ | 0.09199903 | 0.02341114 |
| $N_{\le32}$ | 0.09736782 | 0.04524118 |
| $E_{\le32}L$ | 0.13102701 | 0.20270956 |
| Number density $\times L^3$ | $5.3317798\times10^{-6}$ | $2.4773686\times10^{-6}$ |
| Excitation energy density $\times L^4$ | $7.1749284\times10^{-6}$ | $1.1100205\times10^{-5}$ |

At turnaround, shells 17–32 contribute 5.5% of the reported light-species number and 48.3% of the heavy-species number. The heavy species has fewer quanta over this band and carries more energy. Extrapolating either species' total to all shells requires control of the remaining ultraviolet contribution. [Spectral receipt](MODES.json) · [Ultraviolet analysis](TAIL.md)

## Numerical and geometry controls

The checker multiplies exact unitary exponentials for piecewise constant midpoint Hamiltonians on the archived mesh, split at all readout slices. It uses the mesh at spacings $h,h/2,h/4$. Richardson extrapolation of the $h/2,h/4$ probabilities cancels the leading quadratic error. Comparison with the $h,h/2$ extrapolation checks every probability within $2\times10^{-13}$ absolute and each weighted number within $5\times10^{-10}$. A separate DOP853 integration uses relative tolerance $10^{-11}$, absolute tolerance $2\times10^{-13}$ and maximum step $0.018$. It agrees within $8\times10^{-12}$ per occupation and $2\times10^{-8}$ per band number. These are numerical acceptance thresholds for the archived calculation. [Checker](check_modes.py)

A frozen radius at $a=2$ and a massless field on the evolving radius each have time-independent eigenspaces. Direct integration gives occupations below $10^{-22}$ in both controls. An exponential scale factor permits a Whittaker-function check of the same ODE and readout at two shell indices and three masses. The separate tanh mass-profile calibration tests nonzero mixing. [External analytic reference](SOURCES.md) · [Calibration checker](check_calibration.py)

Replacing the Hermite radius with direct PCHIP through $a$ changes either late number by less than $10^{-7}$. A different physical sensitivity check replaces $a(t)$ with $2+(1+\delta)[a(t)-2]$ at the same readout times, preserving the initial radius and the stationary times. For $\delta=-0.01,+0.01$, the light turnaround number ranges from 0.094857655 to 0.099918416 and the heavy number from 0.044274038 to 0.046215503. These deliberately deformed free-field backgrounds illustrate sensitivity to the radial history. The calculation contains no metric or matter backreaction, collisions, thermal population, or lepton-charge source. [Perturbed runs](MODES.json)

From this directory, run `python -I -B -O check_modes.py --output MODES.json` with NumPy and SciPy. The checker validates the little-endian binary64 values of the three input columns by SHA-256. Each receipt array lists $P_n$ in increasing shell order.
