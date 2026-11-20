# Exact tensor UV residue helper

`tensor_engine.py` implements a general one-loop polynomial UV-pole functional
using exact SymPy arithmetic. This supplies a reproducible calculation tool; it
does not establish a new physical theory, interaction, or ultraviolet completion.

## Conventions

The default is Euclidean dimensional regularization, `d = 4 - 2 epsilon`,
denominators `D_i = (q - p_i)^2 + m_i^2`, and loop measure
`d^d q / (2 pi)^d`. Returned values are

\[
 R_E=16\pi^2\operatorname{Res}_{\epsilon=0} I_E.
\]

`masses2` takes squared masses. Exact integers, rational numbers and symbolic
expressions are accepted; floating-point inputs are rejected. Tensor components
have indices 0 through 3. Contracting four-dimensional indices gives the simple
UV residue, but does not retain finite evanescent contributions. Scalars in
dimensionally regularized scaleless integrals can have cancelling UV and IR
poles: this helper returns the UV polynomial, not the value of the full
scaleless integral. It is not an IR-pole detector or a finite-integral evaluator.

`MinkowskiUVIntegral` is a separate explicit wrapper for time-first metric
`(+---)` and denominators `(q-p_i)^2-m_i^2+i0`. It returns

\[
 R_M=(16\pi^2/i)\operatorname{Res}_{\epsilon=0} I_M,
\]

without vertex factors, propagator numerator factors, coupling constants or
diagram symmetry factors. It uses `p_E=(-i p_0,p_1,p_2,p_3)`,
`q_0=i q_E0`, and a factor `(-1)^N` from the denominators. This convention
gives a massive scalar tadpole `R_M=m^2`, and a triangle rank-two pole
`R_M^{mu nu}=g^{mu nu}/4`.

## General formula

For `N` denominators, introduce barycentric Feynman parameters,

\[
 P=\sum_i x_i p_i,\qquad
 \Delta=\sum_i x_i(m_i^2+p_i^2)-P^2.
\]

Let `S_(n,k)(P)` denote the sum over all **distinct, unnormalized** choices of
`k` disjoint metric pairs among `n` labels, with one factor of `P` for every
unpaired label. Shifting `q=l+P` and integrating the even powers of `l` gives

\[
 R_E^{\mu_1\cdots\mu_n}
 =\sum_{k=\max(0,N-2)}^{\lfloor n/2\rfloor}
 \frac{(-1)^{k+2-N}}{2^k(k+2-N)!}
 \int_{x_i\geq0}\!\!\delta(1-\sum_i x_i)\,d^Nx\;
 \Delta^{k+2-N} S_{n,k}(P).
\]

The pole coefficient comes from
`Gamma(N-k-2+epsilon)`. No extra `Gamma(N)` remains after cancellation
against the Feynman-parameter prefactor. Simplex monomials are integrated as

\[
 \int\delta(1-\sum_i x_i)\prod_i x_i^{a_i}d^Nx
 =\frac{\prod_i a_i!}{(N-1+\sum_i a_i)!}.
\]

The integral is exactly zero when the sum over `k` is empty. For example,
triangles of ranks zero and one and boxes below rank four have no UV pole.

For a three-denominator triangle define `S=sum_i p_i` and
`Q^{mu nu}=sum_i p_i^mu p_i^nu`. Useful checks are

\[
 R_E^{\mu\nu}=\delta^{\mu\nu}/4,\qquad
 R_E^{\mu\nu\rho}
 =\frac{1}{12}\bigl(\delta^{\mu\nu}S^\rho
 +\delta^{\mu\rho}S^\nu+\delta^{\nu\rho}S^\mu\bigr),
\]

\[
 R_E^{\mu\nu\rho\sigma}
 =\frac{1}{48}\sum_{6}\delta^{\mu\nu}
       (S^\rho S^\sigma+Q^{\rho\sigma})
 -\frac{4\sum_i m_i^2+3\operatorname{tr}Q-S^2}{96}
       \sum_3\delta^{\mu\nu}\delta^{\rho\sigma}.
\]

Here the first sum has six distinct choices of one metric pair; the second
has the three pairings of four indices. For a four-denominator box,

\[
 R_E^{\mu\nu\rho\sigma}
 =\frac{\delta^{\mu\nu}\delta^{\rho\sigma}
       +\delta^{\mu\rho}\delta^{\nu\sigma}
       +\delta^{\mu\sigma}\delta^{\nu\rho}}{24}.
\]

## API and verification

```python
import sympy as s
from tensor_engine import UVIntegral, MinkowskiUVIntegral

q = s.symbols('q0:4')
loop = UVIntegral([(0,0,0,0), (1,2,0,0), (2,0,1,0)], [0,0,3])
loop.component((0,0,1,1))       # one exact tensor component
loop.polynomial_residue(q[0]**2*q[1]**2, q)
loop.tensor(4)                  # independent symmetric components only
```

Run `python tensor_engine.py --self-test`. The producer suite contains 73
exact checks: simplex normalization; massive tadpoles and bubbles through rank
six; triangle and box coefficients; nonzero-mass denominator cancellation;
trace and momentum contractions; loop-routing shift covariance; denominator
permutations; and Minkowski signature and cancellation checks.

Eighteen cases also use an independent calculation with no Feynman
parameters: expand each propagator at large `|q|`, select the homogeneous term
with radial logarithmic degree, and average it on the unit three-sphere.
For a monomial `u_0^a0 ... u_3^a3`, odd exponents average to zero; even
exponents give `prod_i (a_i-1)!! / (2^h (h+1)!)`, where
`h=sum_i a_i/2`. This independently reproduces the tensor residue for
denominator counts one through five.

Three negative controls are required to fail: reversing the pole sign,
inserting a spurious factor of two, and normalizing the rank-four pairing sum.
These checks constrain normalization and signs; they do not replace a
complete model-specific diagram list, gauge-consistency analysis or comparison
with data.
