# An interval certificate for the single-direction energy scan

For the mass ratio and four settings in the [single-direction calculation](SINGLE_DIRECTION.md), the profiled coefficient multiplying `|g₇|²` lies strictly between **1.33764501629 × 10⁻⁹ and 1.33764501631 × 10⁻⁹**. The following enclosure uses rational arithmetic and explicit bounds on infinite series. The same calculation gives a conservative lower bound when the actual energies vary over specified windows.

## Enclosing the complex triangle

Use `ξ=s/a`, `B=1−4/ξ`, `β=√B`, and `k=15/17`. At sterile squared masses `a,16a` the subtracted mass combination of triangle masters is `f=a(T_a−T_16a)`, using the *unsubtracted* master T for each mass. Its spectral density is `L(β)/ξ` with `L(u)=2 atanh(ku)`. Changing the dispersion variable to u, the intermediate two-scalar speed, gives the upper-rim formula

```text
ξ Re f = 2 PV ∫₀¹ u L(u)/(u²−B) du,
ξ Im f = π L(β).
```

Expanding L gives the real moments `J_r=PV ∫₀¹ u^(2r+2)/(u²−B) du`. The series for L and its derivative converge uniformly on [0,1] since k<1. To pass the series through the principal value, subtract the numerator value at the simple pole and use uniform convergence of its derivative on a neighborhood of β. The remainder is an ordinary integral and the subtracted constant has an elementary principal value. Polynomial division then yields

```text
J₀ = 1−β atanh β,             J_(r+1) = 1/(2r+3)+B J_r,
ξ Re f = 4 Σ k^(2r+1) J_r/(2r+1).
```

For the points and windows certified here `0<B<2/3`. The positive series for `β atanh β` is at most `B/(1−B)<2`, so `|J₀|≤1`. The recurrence gives `|J_r|≤1/(1−B)`. After retaining r=0 through N−1, the omitted real sum has absolute value bounded by

```text
4 k^(2N+1) / [(2N+1)(1−k²)(1−B)].
```

For the initial value and the imaginary part, the series `Σ z^j/(2j+1)` at `z=B` or `z=k²B` has a positive tail bounded by `z^N/[(2N+1)(1−z)]`. The checker uses N=180. It encloses √B with an integer square root and π with the alternating series for `16 atan(1/5)−4 atan(1/239)`. Rational lower and upper endpoints are rounded outward. For a window the calculation uses the entire input interval in all recurrences, with its upper B endpoint in the tail estimate. This gives a rectangular enclosure for both parts of each physical upper-rim value f.

## Eliminating the three complex nuisance coefficients

At four energies let `X_i=(f_i,1,ξ_i)`, `y_i=(ξ_i/2−1) f_i`. The rate statistic averages the four **squared amplitudes** with weights `v_i=ξ_i/4` before the common factor `1/(2048π⁵)`. Define the cofactor `l_i=(−1)^i det X_{−i}` with zero-based row indices. Laplace expansion gives `Σ l_i X_i=0`. Cauchy–Schwarz gives the minimum over all common complex nuisance coefficients c,

```text
min_c Σ v_i |y_i−X_i c|² = |Σ l_i y_i|² / Σ (|l_i|²/v_i).
```

For rank-three X the minimum attains equality because the weighted residual space has complex dimension one. Conjugating the cofactor to form the Hermitian left-null vector makes its action on the data the bilinear sum displayed here. Every three-row minor is bounded away from zero in the indicated cases. For tighter interval bounds the code removes a fixed complex constant from f and a fixed real constant from n before forming their product. The resulting change in `n f` lies in the nuisance span, and the cofactors annihilate it exactly.

The code bounds the numerator and denominator using the enclosed complex triangles. A strictly positive lower bound on the squared numerator certifies a residual after optimizing over all three complex coefficients. [Checker](verify_scan_certificate.py) · [Exact rational endpoints](SCAN_CERTIFICATE_RESULTS.json)

## What survives a spread in energy

Give each of the settings `5,6,9,11` one quarter of the exposure. Allow the actual ξ around each setting to have any probability distribution supported within **±0.001** of its stated value. Re-evaluating the interval recurrences and cofactors on those four windows gives a conservative coefficient **above 8.61 × 10⁻¹⁰** for every allowed quadruple of energies. For any fixed nuisance coefficients, apply the inequality to one draw from each window and integrate against the product of the four exposure laws. The resulting left side is the equal-weight mean cross section over those laws. This argument works without a specified within-window profile. The halfwidth is in squared-energy units ξ. Wider energy responses and uncertainties in the masses or couplings require separate bounds.

All results assume the finite-dimensional amplitude model and mass spectrum of the [single-direction derivation](SINGLE_DIRECTION.md). Positivity of this algebraic residual does not quantify sensitivity to event counts or backgrounds.
