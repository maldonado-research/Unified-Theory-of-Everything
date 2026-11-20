# How many cold charged fermions fit on the HDblast wall?

HDblast prescribes a cold population before switching on charge mobility. A microscopic interpretation must supply enough low-energy states for that population. Here we test three free, minimally coupled, massive Dirac species on the round wall, with one species per Abelian field. Pauli exclusion and the prescribed magnetic field give a necessary condition independent of how the particles were made. A separate calculation follows neutral TOE singlets and their charged decay products. [Gas simulation](https://github.com/maldonado-research/HDblast/blob/030a4a0b4c2e65fdbce15fbe76b458a7addffdc9/research/HDBLAST_CHARGED_MATTER_HEATING_20261120/coupled/COUPLED.md) · [Singlet decay calculation](CARRIERS.md)

## Counting the carriers required by the fluid

Write $\mu=m_{\rm carrier} L$, set $c=\hbar=1$ and let $a$ be the radius in units of the bulk curvature length $L$. The volume of the three-sphere is $2\pi^2(aL)^3$. The prescribed initial cold rest density is $D_i=10^{-4}$ in wall-tension units at $a_i=2$. With $S_0=\sigma L^4/(\hbar c)$ in ordinary units, the required number of individual carriers, summed over all streams and both signs, is

$$
N_{\rm req}=2\pi^2 a_i^3 D_i S_0/\mu.
$$

The primary charged-fluid simulation assumes a response $Q=|q|/(\sqrt3\mu)=10^8$. The local wall kinetic term gives $\alpha=q^2/(4\pi S_0K)$ when that term dominates the gauge response. At the first mobility step $t=\tau/L=10.5017$, the archived row has $a=9.55326542219306$, $K=13200.04$, $B=-9.634492309871888\times10^{-9}$ and zero initial drift and temperature. Eliminating $S_0$ gives

$$
N_{\rm req}=\mathcal A\mu/\alpha,\qquad
\mathcal A=\frac{3\pi a_i^3D_i Q^2}{2K}=2.855984666\ldots\times10^9 .
$$

The condition $\alpha<1$ thus gives a conservative lower threshold. A smaller coupling raises the required population. Both $Q$ and $K$ belong to this particular simulation. The curvature length, carrier mass and local coupling have no measured values in this comparison. [HDblast normalization](https://github.com/maldonado-research/HDblast/blob/030a4a0b4c2e65fdbce15fbe76b458a7addffdc9/research/HDBLAST_CHARGED_MATTER_HEATING_20261120/synthesis/BRIDGE.md) · [Extracted parameters](HD_INPUTS.json)

## Spectrum in the instantaneous magnetic field

On the unit sphere a homogeneous Hopf connection of amplitude $A$ has physical magnetic amplitude $2A/a^2$ in the HDblast convention. Each of the three fields carries the equal-amplitude projection $1/\sqrt3$. Write $s=qA/\sqrt3$ and $S=|s|=d\mu$, with $d=Qa^2|B|/2=43.9645393\ldots$ at the switch. The homogeneous electric field governs subsequent evolution. In temporal gauge the spatial connection fixes the instantaneous energy used below. [HDblast Hopf connection](https://github.com/maldonado-research/HDblast/blob/030a4a0b4c2e65fdbce15fbe76b458a7addffdc9/research/HDBLAST_CHARGED_MATTER_HEATING_20261120/quantum/PAIRS.md)

In a left-invariant orthonormal frame the spatial Dirac operator can be written $D_s=2\boldsymbol\sigma\cdot\mathbf J+3/2-s\sigma_3$. The spin-$j$ representation of the orbital generators occurs $\ell=2j+1$ times in functions on $SU(2)$. The two extremal spin states give edge eigenvalues $\ell+1/2\pm s$. For $M=-\ell/2+1,\ldots,\ell/2-1$, the paired block has diagonal elements $1/2\pm(2M-s)$ and off-diagonal magnitude $\sqrt{\ell^2-4M^2}$. Its two eigenvalues are

$$
\lambda_\pm=\tfrac12\pm\sqrt{\ell^2+s^2-4Ms}.
$$

Every listed eigenvalue has multiplicity $\ell$ before counting species and charge. Reversal of charge exchanges the two edges and reverses $M$. With three Dirac species and both charge signs the weight is $6\ell$. At $s=0$ each sign of $n+3/2$ occurs $(n+1)(n+2)$ times per two-component spatial operator. Direct diagonalization of the orbital and Pauli matrices and a free-sphere count check the normalization. The massive positive energy is $E/m=\sqrt{1+(\lambda/(a\mu))^2}$. [Spectrum and matrix check](check_bridge.py) · [Free spinor harmonics](https://arxiv.org/pdf/1603.00037v2)

**Exactly zero kinetic energy.** Write the squared interior radical as $R^2=(\ell-S)^2+4Sh$, with $h=1,\ldots,\ell-1$. For $\ell\ge2$ this is at least 4, so an interior eigenvalue cannot vanish. Only one edge can have $\lambda=0$ when $S=\ell+1/2$. Even at that tuned mass, its capacity is $6\ell<6d\mu$. Dividing by $N_{\rm req}$ at $\alpha<1$ gives a fraction below $9.24\times10^{-8}$ at every positive mass. This test interprets the ideal cold fluid literally as massive free fermions with vanishing microscopic kinetic energy.

## Allowing a mean kinetic energy of 10% of rest energy

Let $e_i=E_i/m-1$. To allow an ensemble with dispersion, impose $\sum e_i f_i\le\delta\sum f_i$ with $0\le f_i\le1$ and $\delta=0.1$. This 10% allowance is a sensitivity test chosen here. Filling states in order of increasing energy maximizes the count, with partial occupancy of the final degenerate level allowed. Fractions denote expectation values. Equal populations of opposite signs give a charge-neutral ensemble.

| Snapshot | $\mu$ | Maximum expected count | Required count for $\alpha<1$ |
|---|---:|---:|---:|
| Mobility switch | 1 | 4,014.37 | $>2.85598\times10^9$ |
| Mobility switch | 4 | 68,007.95 | $>1.14239\times10^{10}$ |
| Nearby magnetic zero | 1 | 677.24 | $>2.85598\times10^9$ |
| Nearby magnetic zero | 4 | 46,608.42 | $>1.14239\times10^{10}$ |

The zero-field comparison uses $t=10.50164818220685$, $a=9.553254051209148$ and the same late $K$. Even the small simulated field at the mobility step appreciably changes the number of low-energy states because $Q$ is large. Counts for an alternative preparation with charge already present at the initial radius appear in the numerical receipt. [Counts and marginal energies](BRIDGE.json)

### Exclusion over continuous mass intervals

Let $C(y)$ count individual states with $e_i\le y$. For every $z>\delta$ a feasible ensemble satisfies

$$
(z-\delta)N\le\sum_i(z-e_i)_+=\int_0^z C(y)\,dy .
$$

For small masses, at least half the occupied weight must have $e_i\le0.2$. The operator displacement $\|D_s-D_0\|=S$ bounds the displacement of singular values by $S$. With $p=(a\sqrt{0.44}+d)\mu$ and $J=\max(0,\lfloor p-3/2\rfloor+1)$ the free count gives $C(0.2)\le4J(J+1)(J+2)$. It vanishes for $p<3/2$ and is at most $4(p+1/2)^3$ thereafter. In this latter range $8(p+1/2)^3/\mu$ increases with $\mu$. The value at $\mu=52$ lies below $\mathcal A$, excluding all $0<\mu\le52$.

For the rest take $z=0.163$, $c=a\sqrt{y(2+y)}$, $w=c\mu$ and $t_*=w+1/2<S$ for $\mu\ge52$, $y\le z$. The edges give at most $6(S+w)(2w+1)$ states. Both interior branches together give at most

$$
\frac3S\sum_{\ell\ge1}\ell\,[t_*^2-(\ell-S)^2]_+ .
$$

Indeed, each branch requires $R\le t_*$, and counting positive integers $h$ in the expression for $R^2$ gives the bound. The nonnegative unimodal summand has sum at most its integral plus its maximum. Its integral is $4S t_*^3/3$, and its maximum is at most $(S+t_*)t_*^2$. Combining edge and interior counts, for $\mu\ge\mu_0$ and $u=c+1/(2\mu_0)$,

$$
 C(y)\le4u^3\mu^3+
 \left[6(d+c)(2c+1/\mu_0)+3(1+u/d)u^2\right]\mu^2.
$$

The coefficients increase with $y$. A right-endpoint sum with 4,096 bins therefore bounds the integral from above. The checker uses rational bounds $a<9.554$, $43.96<d<43.97$, $2.855\times10^9<\mathcal A<2.856\times10^9$ and rounds square roots outward with integer arithmetic. After division by $\mu$, the upper bound increases on each interval. Testing the upper endpoints of $[52,1900]$ and $[1900,1970]$ excludes **all masses through 1,970** under $\alpha<1$. For $\alpha\le0.1$ the required number grows by at least a factor of ten and the interval $[1970,6240]$ is excluded as well. At $\mu=2000$ the same inequality requires $1/\alpha<1.03$.

### Explicit populations at larger masses

Occupy a common fraction of all states below $z_*=0.1625$. Before scaling the count, the mean kinetic energy per rest energy is $z_*-\int_0^{z_*} C(y)\,dy/C(z_*)$. For $|\lambda|\le w$, the interior branches require respectively $0\le R\le w-1/2$ and $\max(0,1/2-w)\le R\le1/2+w$. The checker counts the allowed integer $h$ values in each interval and both edges. It uses $S$ rounded down to $0.001$, enlarges and shrinks the momentum cutoff to enclose the shift and square-root errors, and integrates the lower counts with 384 left-endpoint bins.

| $\mu$ | Trial $\alpha$ | Count in the eligible states, lower to upper | Upper mean kinetic/rest |
|---:|---:|---:|---:|
| 2,000 | $100/101$ | 5,811,869,612,934 to 5,811,873,274,998 | 0.0997860 |
| 6,500 | $10/101$ | 199,512,966,694,884 to 199,513,018,735,572 | 0.0997874 |

In each case the lower count already exceeds the number required at the chosen coupling. Uniformly thinning both charge signs to the target preserves the mean. The exclusion and the explicit populations bracket the infimum of admissible masses without assuming shell counts are monotonic. At 2,000 the required coupling is close to unity, where a free description deserves particular caution. The examples establish instantaneous capacity. Their production, cooling, conversion to mobile charge, and the interacting distribution through time are open problems. Other particle statistics or dynamics and a different gauge response change this test. [Executable bounds and witness](check_bridge.py) · [Numerical receipt](BRIDGE.json)
