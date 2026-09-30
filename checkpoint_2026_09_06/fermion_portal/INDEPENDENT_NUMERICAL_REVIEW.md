# Independent review of the prescribed fermion mass-pulse benchmark

Date: 6 September 2026. Reviewed benchmark.py and the completed BENCHMARK_RESULTS.json. This is a secondary numerical applicability example; the user's selected task is the TOE off-shell loop calculation. The pulse has not been derived from or calibrated to the TOE model.

**Outcome:** the Hamiltonian, vacuum projection, one-species Majorana counting and stated exact-arithmetic state-error estimates are consistent. Independent numerical controls pass. The results establish particle production for the specified external mass histories; they do not establish reheating, leptogenesis, an observational fit, or a fundamental discovery.

The independent implementation is independent_review.py; its output is INDEPENDENT_REVIEW_RESULTS.json. It does not import or mutate the benchmark. The reviewed benchmark-results SHA-256 is:

58eee78f9700ed885507ac9d4c4c4cdf533b22b9c21f376a6708d495721c6387

## 1. Physical interpretation and conventions

For a spatially homogeneous real fermion mass in flat spacetime, a fixed-momentum helicity mode can be represented by a two-component Hermitian Hamiltonian. A constant unitary basis rotation gives the implemented form

$$H_k(t)=m(t)\sigma_z+k\sigma_x,\qquad
m(t)=1-a\,\operatorname{sech}^2(t/\tau).$$

This is a fermionic two-level system. A scalar oscillator with the same mass profile would obey a different equation and different Bogoliubov normalization. The implemented unitary evolution preserves the fermionic occupation bound.

With the asymptotic mass equal to one, put $\theta=\arctan(k)$ and $E=\sqrt{1+k^2}$. The code uses

$$u_-=\begin{pmatrix}-\sin(\theta/2)\\\cos(\theta/2)\end{pmatrix},
\quad u_+=\begin{pmatrix}\cos(\theta/2)\\\sin(\theta/2)\end{pmatrix}.$$

Direct multiplication gives $H_0u_\pm=\pm Eu_\pm$ and $u_+^\dagger u_-=0$. Starting in $u_-$ and evaluating $|u_+^\dagger Uu_-|^2$ is the correct asymptotic transition probability for this representation. Finite start/end times are accounted for by the omitted-tail estimate below. Free asymptotic phases do not change the probability.

At $k=0$, all Hamiltonians commute and are diagonal in a time-independent basis. Because the initial and final signed masses are the same, the final production probability is exactly zero even when the mass crosses zero twice. Intermediate changes of instantaneous energy labels do not imply final particles.

For a real scalar mass history the two helicities have equal occupation. One Majorana species has two physical spin states, so

$$n_{[0,K]}=\frac{1}{\pi^2}\int_0^K dk\,k^2 n_k,\qquad
\rho_{[0,K]}=\frac{1}{\pi^2}\int_0^K dk\,k^2\sqrt{1+k^2}\,n_k.$$

These match the code's finite-band normalization. They are out-particle contributions, not an unrenormalized vacuum energy. Two degenerate Majorana species, used in the cancellation counterexample below, would double these total densities; the benchmark itself has one species.

## 2. Error bounds

The following are valid for the implemented cases $a\geq0$, $\tau>0$ and exact arithmetic.

For two Hermitian Hamiltonians, Duhamel's identity and unitarity imply

$$\|U-\widetilde U\|\leq\int dt\,\|H(t)-\widetilde H(t)\|.$$

On a cell of width $\Delta t$, replacing $m(t)$ with its midpoint value gives

$$\int_{\rm cell}|m(t)-m(t_{\rm mid})|dt
\leq\frac{\Delta t}{2}\int_{\rm cell}|\dot m(t)|dt.$$

The variation of the pulse over the full line is $2a$. Summing cells therefore yields

$$\epsilon_{\rm disc}\leq a\Delta t.$$

This deliberately loose first-order global bound is compatible with the much smaller observed midpoint-refinement differences; it does not claim the observed error saturates the bound.

The total omitted interaction on the two tails outside $|t|\leq L\tau$ is

$$\epsilon_{\rm tail}
\leq2a\tau(1-\tanh L)
=\frac{4a\tau}{e^{2L}+1}
\leq4a\tau e^{-2L}.$$

For a normalized exact target state and a state approximation at distance at most $\epsilon$, a projection probability differs by at most $2\epsilon+\epsilon^2$. Thus the reported probability bound follows from $\epsilon=\epsilon_{\rm disc}+\epsilon_{\rm tail}$. Exact normalization of both states would permit a slightly tighter bound, but the implemented one is valid.

These arguments exclude floating-point roundoff. A small norm error alone cannot rigorously bound phase error or the accumulated floating-point error. The code correctly labels the analytic estimate as excluding rounding; it must not be described as a certified machine enclosure. Likewise, Simpson refinement is a useful diagnostic, not an error certificate for the momentum integral or the uncomputed ultraviolet tail.

## 3. Independent numerical controls

The independent calculation uses the Hadamard-rotated Hamiltonian

$$H'_k=k\sigma_z+m(t)\sigma_x,$$

obtains its asymptotic eigenstates from a general eigensolver, and evolves four real components with SciPy's **implicit Radau** solver. This differs from the benchmark's explicit unitary midpoint products and DOP853 controls. The independent integration uses $L=14$, relative tolerance $5\times10^{-12}$ and absolute tolerance $5\times10^{-13}$.

| Case | $k/M$ | Benchmark probability | Independent Radau probability | Absolute difference |
|---|---:|---:|---:|---:|
| Constant mass | 0.600 | $5.75\times10^{-29}$ | $2.15\times10^{-30}$ | $5.53\times10^{-29}$ |
| Shallow, fast | 0.600 | 0.0407539957 | 0.0407539972 | $1.49\times10^{-9}$ |
| Crossing, fast | 0.550 | 0.9117087136 | 0.9117087059 | $7.74\times10^{-9}$ |
| Crossing, slow | 0.250 | 0.7248309238 | 0.7248308724 | $5.14\times10^{-8}$ |
| Shallow, slow | 0.275 | 0.00160304109 | 0.00160304137 | $2.80\times10^{-10}$ |

Selected-state norm errors were below $1.1\times10^{-13}$. The $k\leftrightarrow-k$ helicity/reflection control agreed to $6.7\times10^{-16}$. The zero-momentum crossing control returned $2.7\times10^{-34}$, consistent with its exact zero.

The benchmark's largest sampled DOP853/midpoint difference is $5.15\times10^{-8}$; the largest change when doubling midpoint steps is $1.58\times10^{-7}$. The combined cutoff/tolerance control changes probabilities by at most $1.19\times10^{-10}$. None of these empirical comparisons is a proof about untested momenta.

The finite-band number densities in the nonconstant cases are approximately $0.00228062$, $0.0343984$, $0.00175734$ and $0.00000642078$ in units of $M^3$, in the case order shown in the JSON. These are diagnostics over $0\leq k/M\leq8$. In particular, a value such as 0.9117 is a mode occupation, not the fraction of cosmological energy converted to particles.

## 4. Independent weak-pulse analytic control

In the interaction picture of $H_0$, the first-order transition amplitude for this real mass dip has magnitude

$$|A_1(k)|=\frac{2\pi a k\tau^2}{\sinh(\pi E\tau)},\qquad E=\sqrt{1+k^2}.$$

This follows by projecting the perturbation $-a\,\operatorname{sech}^2(t/\tau)\sigma_z$ between $u_-$ and $u_+$ and taking its Fourier transform at frequency $2E$.

One further exact Duhamel iteration, using unit norm of the evolved state, bounds the remainder after the first-order term by

$$|A-A_1|\leq\frac12\left(\int dt\,\|V(t)\|\right)^2
=2a^2\tau^2.$$

For $a=0.01$, $\tau=1$ and $k=0.5$, the first-order amplitude magnitude is approximately $0.00187563$, while this rigorous analytic remainder bound is $0.0002$. Consequently the exact infinite-time probability is strictly positive and lies between approximately

$$2.8077\times10^{-6}<n_{0.5}<4.3083\times10^{-6}.$$

The independent numerical probability is $3.6605750\times10^{-6}$, inside the interval. The $L=14$ tail state bound is $2.77\times10^{-14}$. The displayed decimal evaluations are ordinary floating arithmetic; the strict positivity has a large margin and follows from the analytic expression, not merely the positive numerical result.

## 5. Exact seesaw cancellation does not measure particle production

Take two degenerate singlets with

$$M_N(t)=m(t)I_2,\qquad Y=(y,\;iy),\qquad m(t)>0.$$

Then the tree Weinberg coefficient cancels identically:

$$C_5(t)=YM_N(t)^{-1}Y^T
=\frac{yy^T}{m(t)}(1+i^2)=0.$$

The same cancellation holds for the degenerate tree-level lepton-number-violating finite-pole kernel and its odd inverse-mass moments. It does **not** remove the mass-pulse dynamics of either singlet. The weak-pulse calculation above therefore supplies a controlled example with $C_5=0$ and nonzero particle production. No zero-mass crossing is required.

This cancellation pair is equivalent to a lepton-number-conserving Dirac system: defining $N_\pm=(N_1\pm iN_2)/\sqrt2$ turns the degenerate mass term into an off-diagonal Dirac mass and leaves only one linear combination coupled by $Y$. The symmetric model creates particles and antiparticles without a net lepton asymmetry. In particular, $\operatorname{Im}[(Y^\dagger Y)_{12}^2]=0$. Thus it is **not a leptogenesis example**.

This is a known symmetry structure combined with a controlled production calculation, not a claim of new fundamental mathematics. It disproves the proposed implication “a vanishing tree LNV coefficient forces zero particle production” within the stated model. It does not prove any cosmological TOE realization.

## 6. Scope and remaining limitations

No blocking defect was identified in the benchmark within its stated scope. Appropriate retained limitations are:

- The mass is externally prescribed; the scalar dynamics and its available energy have not been solved.
- Backreaction, expansion, particle collisions, decay, thermalization, CP violation and washout are absent.
- A free-particle occupation spectrum is not a cosmological abundance or a likelihood fit.
- The finite momentum range and floating arithmetic remain explicit limitations.
- Slow evolution can suppress production or rearrange interference; these five examples do not prove a universal monotonic relation with pulse duration.
- This benchmark does not resolve the TOE dimension-seven mixed-loop matching or its equation-of-motion reduction.
