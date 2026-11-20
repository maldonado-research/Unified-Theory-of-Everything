# Exact tests of the closed-wall fermion calculation

The same proper-time code that propagates the archived wall history reproduces the exact occupation of all 33 retained shells on a soluble test geometry, for both benchmark masses. A second test measures the factor of two in the Majorana count by taking a large-radius limit. These tests concern a free, minimally coupled spectator on a specified metric. [Checker](check_calibration.py) · [Numerical receipt](CALIBRATION.json)

## A smooth expansion with known scattering data

Write $d\eta=dt/a(t)$ and $M=ma$. The rescaled spinor on the round unit three-sphere has a two-component Hamiltonian $H_k=k\sigma_1+M(\eta)\sigma_3$ for each signed intrinsic Dirac eigenvalue $k=\pm(n+3/2)$. Conjugation by $\sigma_3$ interchanges the signs of $k$ and preserves the transition probability. For the calibration choose

$$
M(\eta)=\tfrac12(M_f+M_i)+\tfrac12(M_f-M_i)\tanh(\eta/T),
\qquad \omega_j=\sqrt{k^2+M_j^2}.
$$

The asymptotically static ends define the in and out states. The upper spinor component satisfies $u''+[k^2+M^2+iM']u=0$. For $z=(1+\tanh(\eta/T))/2$, take out the factor $z^{-iT\omega_i/2}(1-z)^{iT\omega_f/2}$. The regular incoming solution remaining at $z=0$ is a Gauss hypergeometric function. Define $\Delta M=M_f-M_i$, $\Delta\omega=\omega_f-\omega_i$, $x=T(\Delta\omega-\Delta M)/2$, $y=T(\Delta\omega+\Delta M)/2$. Its parameters are $ix,1+iy,1-iT\omega_i$. Continuation to $z=1$ gives the upper-component coefficient of the negative-frequency wave,

$$
B=\frac{\Gamma(1-iT\omega_i)\Gamma(iT\omega_f)}
        {\Gamma(ix)\Gamma(1+iy)} .
$$

Converting the upper-component coefficient into a probability requires the norm of each asymptotic eigenspinor. The incoming positive-energy spinor has upper-component norm squared $(\omega_i+M_i)/(2\omega_i)$, and the outgoing negative-energy spinor has $(\omega_f-M_f)/(2\omega_f)$. Taking the squared ratio times $|B|^2$, using $|\Gamma(iv)|^2=\pi/[v\sinh(\pi v)]$, gives

$$
P_k=
\frac{\cosh(\pi T\Delta M)-\cosh(\pi T\Delta\omega)}
     {2\sinh(\pi T\omega_i)\sinh(\pi T\omega_f)} .
$$

Unitarity gives the same modulus for a filled negative-energy input projected on the positive-energy output, the convention in the wall solver. Alves and Camilo give the tanh profile, scalar-component equation, spinor conversion and gamma-function connection coefficient. The spatial spectrum changes the allowed $k$ on the sphere. [Alves and Camilo, §4, Eqs. (48)–(54) and (62)–(66)](https://arxiv.org/abs/1712.01400)

The proper-time test sets $a_i=0.5,\ a_f=2,\ T=0.12$. To feed this exact profile to the archive solver, the checker inverts $t=\int_{-18T}^\eta a(s)\,ds$ and evaluates $a(\eta(t))$. At $|\eta|=18T$ each endpoint mass is within $e^{-36}|\Delta M|$ of its asymptote. For $m=1,4$, the exact occupations summed through $n=32$ are respectively $2.3227166997$ and $5.5841493196$. The numerical check compares each individual shell using an absolute allowance of $2\times10^{-14}$ plus a relative allowance of $3\times10^{-7}$. An independent conformal-time propagation starts in the positive-energy state and tests both spatial signs under rapid expansion, slow expansion, rapid contraction, static massive evolution, and massless evolution. [Archive propagator](check_modes.py) · [Calibration results](CALIBRATION.json)

| $M_i$ | $M_f$ | $T$ | $P_{n=0}$ | Expected Majorana count on the entire sphere | Out conformal excitation energy |
|---:|---:|---:|---:|---:|---:|
| 0.5 | 2 | 0.12 | 0.07530463 | 2.32271671 | 11.23209231 |
| 0.5 | 2 | 0.60 | 0.00509112 | 0.02370569 | 0.06180611 |
| 2 | 0.5 | 0.12 | 0.07530463 | 2.32271671 | 10.07654414 |

The counts and energies use $2(n+1)(n+2)P_k$, multiplied by $\omega_f$ for energy. Physical out energy is the conformal energy divided by $a_f$. The sum includes $n=0,\ldots,128$, with an analytic bound on the remainder. Indeed, for $k\ge k_0>0$, each $P_k$ is at most $Ce^{-2\pi T k}$, where $C=2[\cosh(\pi T\Delta M)-1]/[1-e^{-2\pi T k_0}]^2$. For first omitted shell $j=129$, take $k_0=j+3/2$. Consecutive upper-envelope terms have number ratio at most $q=e^{-2\pi T}(j+3)/(j+1)$ and energy ratio at most $q(1+1/k_0)$. Summing the two resulting geometric series bounds the omitted number by $4.0\times10^{-39}$ and the omitted conformal energy by $5.3\times10^{-37}$ across these three examples. This exponential bound uses the exact tanh spectrum and applies to these calibration profiles. [Bounds and assumptions for the archived trace](TAIL.md)

## What the rotation bound can establish

With the Hamiltonian angle $\theta=\operatorname{atan2}(|k|,M)$, the sudden limit $T\to0$ gives

$$
P_k=\frac12\left(1-\frac{k^2+M_iM_f}{\omega_i\omega_f}\right)
=\sin^2[(\theta_f-\theta_i)/2].
$$

The instantaneous eigenbasis sweeps an angle $|\theta_f-\theta_i|$ during a monotone step. A sudden step saturates the corresponding half-angle transition bound, while the slow example in the table has the same angular change and much smaller occupation. The bound by itself gives $P_k=O(k^{-2})$ at large $k$. Multiplication by the $O(k^2)$ degeneracy shows why this bound alone cannot control the total number.

For a return to the starting radius, consider a step $M_0\to M_1$, a conformal plateau of duration $L$, and the reverse step. Exponentiating the plateau Hamiltonian gives

$$
P_k=\frac{k^2(M_1-M_0)^2}{\omega_0^2\omega_1^2}
        \sin^2(\omega_1L).
$$

For $k=1.5,\ M_0=0.5,\ M_1=2$, the lowest shell has $P=0.324$ at $\omega_1L=\pi/2$ and exactly zero at $\omega_1L=\pi$. The endpoint radii and total angle swept are identical. The example isolates interference between the two changes in radius in a single shell. The archive calculation integrates its resolved, smooth history. [Matrix-exponential check](check_calibration.py)

## The Majorana multiplicity and its flat limit

The intrinsic two-component Dirac operator has eigenvalues $\pm(n+3/2)$ on the unit three-sphere, each with $d_n=(n+1)(n+2)$ spinor harmonics. For a four-component spacetime field, the two time-dependent components associated with each intrinsic harmonic give one positive- and one negative-energy eigenvector. Majorana reality identifies the negative-energy coefficients with conjugate creation operators and leaves one independent annihilation operator for each positive-energy mode. The two signs of the spatial operator therefore give $2d_n$ physical single-particle states. An expected *pair* count would be half the expected particle count in a vacuum production experiment. [Camporesi and Higuchi, §3.1, Eqs. (3.67)–(3.68)](https://arxiv.org/abs/gr-qc/9505009)

To test the normalization independently, take $a_f=R,\ a_i=0.4R,\ m=1,\ T=0.35/R$. Then $p=(n+3/2)/R$ and the common physical distribution $P(p)$ is the tanh formula with momentum $p$, effective endpoint masses $0.4,1$, and width $0.35$, all in final-radius coordinates. With $V=2\pi^2R^3$ and $\Delta p=1/R$,

$$
\frac{1}{V}\sum_n2d_nP(p_n)=\frac{1}{\pi^2}\sum_{n\ge0}\Delta p
\left(p_n^2-\frac{1}{4R^2}\right)P(p_n)
\longrightarrow \frac{1}{\pi^2}\int_0^\infty p^2P(p)\,dp .
$$

This is precisely the flat-space two-spin measure $2\,d^3p/(2\pi)^3$. The smooth even $P(p)$ in this example allows midpoint summation after including $p=1/(2R)$, whose weight vanishes. Odd derivatives at $p=0$ vanish, so the leading difference between the flat and spherical density is $\int_0^\infty P(p)\,dp/(4\pi^2R^2)$. The code independently integrates both continuum coefficients and sums on spheres of radii $4,8,16,32,64$. At $R=64$ it gives $0.004326741281$ compared with the flat density $0.0043269759724$, and the coefficient of $R^{-2}$ agrees with $0.00096129626534$. [Flat-limit table](CALIBRATION.json)

Run the calibration from this directory with

    python -I -B -O check_calibration.py
