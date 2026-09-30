# TOE scalar continuum and finite-pole matching: a conditional spectral result

**Date:** 2026-09-06. **Status:** derived conditional mathematics and independently checked toy model. No experimental discovery, new particle-production result, or novelty claim.

The main result is that a finite static response and its first derivative do not determine whether a model contains a gapless continuum or a discrete pole. A positive one-pole approximation can match both quantities exactly and still have the wrong nonanalytic response and spectral discontinuity. This provides a concrete test for extending TOE v294’s scalar branch into a low-energy effective description.

This note concerns the TOE scalar spectral branch. Other projects are reference material. It does not assert a connection to a higher-dimensional blast, or identify the scalar measure with the complex-symmetric Majorana kernel in N00AK.

## 1. Precise inherited assumptions

The packaged inherited source is [TOE_v294_THEOREM.md](../inputs/TOE_v294_THEOREM.md), especially equations (2.3), (2.4), (7.2), (7.3) and (8.1). The source supplies, conditionally on its locked scalar operator and endpoint trace:

- a positive self-adjoint scalar operator with spectrum `[0,infinity)`;
- no point spectrum, hence no discrete zero mode or other discrete spectral atom for this operator;
- a positive boundary spectral measure `d rho`, with finite total mass;
- near the lower endpoint, the density

\[
r(\lambda):=\frac{d\rho}{d\lambda}(\lambda)
=C\lambda^{\beta_{\rm spec}}[1+o(1)],\qquad C>0,
\]

\[
\boxed{\beta_{\rm spec}=\nu=1+\frac ak>1.}
\]

Thus the inherited spectrum is **gapless**. The source’s asserted absence of point spectrum is an input to this note, not a result re-proved here. A separate atom at zero would invalidate the finite inverse moments used below. A gap would change the expansion problem.

Notation matters: `r(lambda)` is the density, while `d rho(lambda)` is the measure. The v294 notation `rho'(lambda)` means this density; it is not a derivative to take again. Also, v294 already uses `beta=2/mu_b` for an endpoint coefficient. That endpoint parameter is unrelated to `beta_spec`. The Euclidean variable `s=Q^2>0` below is unrelated to the dimensionless endpoint coordinate `(1-tau)/2` also called `s` in some v294 formulas.

Define

\[
G_E(s)=\int_0^\infty\frac{d\rho(\lambda)}{\lambda+s}.
\]

This equals v294’s `G(-s)`. Spectral `lambda` and Euclidean `s` both have units of mass squared in a dimensioned application.

## 2. Exactly how far a local expansion exists

Let

\[
\mu_n=\int_0^\infty\frac{d\rho(\lambda)}{\lambda^{n+1}},
\qquad n=0,1,2,\ldots.
\]

Finite total measure controls the ultraviolet part of these inverse moments. At the lower endpoint the integrand scales as `lambda^(beta_spec-n-1)`. Therefore

\[
\boxed{\mu_n<\infty\quad\Longleftrightarrow\quad n<\beta_{\rm spec}.}
\]

This is v294’s inverse-moment condition with the index shifted: `mu_n=I_(n+1)`. When finite, `G_E^(n)(0+)=(-1)^n n! mu_n`.

For noninteger `beta_spec`, put `m=floor(beta_spec)`. The exact resolvent identity gives

\[
G_E(s)=\sum_{n=0}^{m}(-1)^n\mu_n s^n+R_m(s),
\]

\[
R_m(s)=(-1)^{m+1}s^{m+1}
\int_0^\infty\frac{d\rho(\lambda)}{
\lambda^{m+1}(\lambda+s)}.
\]

Substitute `lambda=s t` in the threshold part and use

\[
\int_0^\infty\frac{t^{\delta-1}}{1+t}\,dt
=\frac{\pi}{\sin(\pi\delta)},\qquad
\delta=\beta_{\rm spec}-m\in(0,1).
\]

Then

\[
\boxed{
G_E(s)=\sum_{n=0}^{m}(-1)^n\mu_n s^n
-\frac{\pi C}{\sin(\pi\beta_{\rm spec})}s^{\beta_{\rm spec}}
+o(s^{\beta_{\rm spec}}).
}
\]

The analytic coefficients depend on the entire measure. The leading nonanalytic coefficient depends on its threshold density. For `1<beta_spec<2`, the nonanalytic term is positive because `sin(pi beta_spec)<0`.

If `beta_spec=m` is an integer, only `mu_0,...,mu_(m-1)` are finite. The universal leading logarithm is

\[
\boxed{
G_E(s)=\sum_{n=0}^{m-1}(-1)^n\mu_n s^n
+(-1)^{m+1}C s^m\log(s/\lambda_{\rm ref})
+o\!\left(s^m\log(\lambda_{\rm ref}/s)\right),
}
\]

where `lambda_ref>0` carries the same units as `s`. Under only the stated `1+o(1)` density asymptotic, the remainder need not be `O(s^m)`: slowly varying corrections can produce weaker logarithms. A stronger density remainder, such as `O(lambda^epsilon)`, permits a finite convention-dependent analytic term at that order. The logarithmic coefficient itself is fixed.

## 3. A positive pole matched to two moments is a bound

Assume `mu_0,mu_1` are finite and positive. There is one positive single-pole function that matches the static value and first derivative:

\[
G_{1p}(s)=\frac{w}{\lambda_0+s},\qquad
\boxed{\lambda_0=\frac{\mu_0}{\mu_1},\quad
w=\frac{\mu_0^2}{\mu_1}.}
\]

Define a probability measure

\[
dP(\lambda)=\frac{d\rho(\lambda)}{\mu_0\lambda}.
\]

With `Y=1/lambda`, `E_P[Y]=mu_1/mu_0` and

\[
\frac{G_E(s)}{\mu_0}
=\mathbb E_P\left[\frac1{1+sY}\right].
\]

The function `1/(1+sY)` is strictly convex for `s>0,Y>0`. Jensen’s inequality proves

\[
\boxed{
G_E(s)\ge \frac{\mu_0^2}{\mu_0+s\mu_1}=G_{1p}(s),
\qquad s>0.
}
\]

Equality holds exactly when the contributing measure is supported at a single positive spectral value. Thus it is strict for a nontrivial continuum. This is an application of ordinary positive-measure/Stieltjes and Jensen methods, not a claimed new mathematical theorem.

For the normalized v294 trace-response shape

\[
\widehat\chi(s)=\frac{\mu_0-G_E(s)}s,
\]

the inequality reverses:

\[
0<\widehat\chi(s)\le
\frac{\mu_1}{1+s\mu_1/\mu_0}.
\]

Multiplying by the positive v294 prefactor `kappa_5^4/18` preserves this order. It does not establish a physical matter coupling.

## 4. Fully explicit counterexample

Choose dimensionless spectral units and

\[
d\rho(\lambda)=\lambda^{3/2}\mathbf 1_{[0,1]}(\lambda)d\lambda.
\]

This is a mathematical fixture, not a fitted or dimensioned v294 parameter point. Its threshold exponent corresponds to `a/k=1/2`, but its whole spectral density is not claimed to equal v294’s density.

Putting `lambda=x^2` gives

\[
G_E(s)=2\int_0^1\frac{x^4}{x^2+s}\,dx
=\boxed{\frac23-2s+2s^{3/2}\arctan(1/\sqrt{s})}.
\]

Its moments are

\[
\mu_0=\frac23,\qquad\mu_1=2,\qquad\mu_2=\infty.
\]

The exact one-pole fit is

\[
\boxed{\lambda_0=\frac13,\quad w=\frac29,\quad
G_{1p}(s)=\frac{2/9}{1/3+s}=\frac{2/3}{1+3s}.}
\]

Both have `G(0)=2/3` and `G'(0+)=-2`. Their next terms differ:

\[
G_E(s)=\frac23-2s+\pi s^{3/2}-2s^2+O(s^3),
\]

\[
G_{1p}(s)=\frac23-2s+6s^2+O(s^3).
\]

In particular, the continuum has a divergent second derivative at zero; the pole has a finite one. They cannot share the same exact derivative expansion.

The associated normalized trace responses are

\[
\widehat\chi_E(s)=2-2\sqrt{s}\arctan(1/\sqrt{s}),
\qquad
\widehat\chi_{1p}(s)=\frac2{1+3s}.
\]

Both tend to `2`, but the continuum’s leading static deficit is `pi sqrt(s)`, while the pole deficit is `6s`. This reproduces the v294 distinction between a finite static response and a potentially nonsmooth next derivative.

## 5. What the two moments cannot determine

For an interior point `x>0` where the density is regular, the stated continuation has

\[
\operatorname{Im}G_E(-x-i0)=\pi r(x).
\]

Thus, for the continuum fixture,

\[
\operatorname{Im}G_E(-x-i0)=\pi x^{3/2},\qquad0<x<1.
\]

The fitted pole instead has the distributional result

\[
\operatorname{Im}G_{1p}(-x-i0)=\pi\frac29\delta(x-1/3).
\]

Away from `x=1/3`, that pole limit is zero. Therefore exact matching of these two Euclidean moments does not identify the continuum discontinuity. If one writes `x=omega^2`, the delta transformation and any physical frequency-domain normalization must also be included; it cannot be treated as an ordinary constant density.

The sign shown follows the explicitly chosen `-i0` denominator. This is a **spectral discontinuity statement**. Calling it a measured absorptive susceptibility or a particle-production rate additionally requires a physical source coupling, a retarded prescription, state/vacuum specification, and normalization. Those ingredients have not been derived here. No real-time decay law is inferred from an imaginary-time/heat-semigroup decay.

## 6. Consequence for the TOE work

The v294 source supports a controlled statement about a declared scalar trace response. Because `beta_spec>1`, the static quantity and first Euclidean derivative can be finite despite a continuum reaching zero mass. If `1<beta_spec<=2`, replacing the response by a local Taylor expansion through second order fails at threshold; a nonanalytic term or explicit continuum must be retained. If `beta_spec>2`, more local coefficients exist, but the same finite-moment boundary eventually occurs.

The productive use is to preserve the derived scalar spectral asymptotic when constructing its EFT, and to test any proposed physical source coupling against that measure. This does not add the missing second source in v294 or complete the Standard Model/gravity crosswalk. N00AK’s finite Majorana-pole algebra is a separate object: complex-symmetric matrix residues do not satisfy the scalar positivity assumptions used in the Jensen bound.

## 7. Independent verification

The adjacent standard-library script `verify_continuum_bridge.py` performs:

- seven positive-axis comparisons between adaptive numerical quadrature and the independently evaluated exact toy formula;
- exact rational checks of both fitted moments;
- the `beta_spec=2` logarithmic fixture `G_E(s)=1/2-s+s^2 log((1+s)/s)`;
- 1,000 deterministic random positive-mixture tests of the Jensen bound;
- twelve finite-regulator spectral integrals, checked against the continuum cut limit and the pole’s off-atom zero limit.

All checks pass. The maximum exact-formula/quadrature discrepancy is `1.33e-15`. At regulator `eta=1e-4`, the largest continuum-limit discrepancy among the four test locations is `5.88e-4`; this is an explicitly finite-regulator convergence check, not a claim of an exact limiting numerical value. The analytic discontinuity and Jensen results are proved above.

Outputs: `verification_results.json`, `euclidean_comparison.csv`, and `spectral_comparison.csv`. No source files were changed, no experimental input was fabricated, and no release was published by this subtask.
