# TOE Tier-E v294 parent-lock, exact A1 zero resolvent, threshold law, inverse-moment ladder, and static trace-response theorem

## 1. Parent lock and scope

The canonical v293 parent used by this release is

```text
TOE_v293_COMPLETE.zip
SHA-256: 178b6378ec1169a1da57ab27d50f4eae5fa52f4e45b1a4db03b75b2ceb7bb567
bytes: 10111750
```

The archive passes ZIP integrity and safe-inventory checks.  v294 keeps the
minimal A1 action, the fixed-brane scalar quotient, the endpoint-complete
Hilbert form, and every v293 claim boundary unchanged.

v293 left two precise analytic questions open:

1. determine the boundary Weyl datum `G(0)`; and
2. determine the threshold law of the scalar spectral measure strongly enough
   to decide whether the inverse second moment is finite.

v294 closes both questions for the full locked A1 parameter family.  It does
not instantiate a numerical parameter point, add a second source action, or
construct an off-line parent matrix.

The central result is:

```text
G(0) has an exact positive hypergeometric family formula;
the static declared-trace response is finite;
the threshold density is proportional to lambda^(1+a/k);
the inverse moment I_n is finite exactly when n<2+a/k.
```

## 2. Locked A1 scalar system

The declared action remains

\[
\begin{aligned}
S_{\rm A1}={}&\frac1{2\kappa _5^2}
 \int_{M_+\cup M_-}\sqrt{-g}\,R
+\frac1{\kappa _5^2}\int_\Sigma\sqrt{-\gamma}(K_++K_-)\\
&-\int_{M_+\cup M_-}\sqrt{-g}
 \left[\frac12(\partial\Phi)^2+V(\Phi)\right]
-\int_\Sigma\sqrt{-\gamma}\,\lambda_b(\Phi),
\end{aligned}
\]

\[
W=W_0\sin(\Phi/v),\qquad
V=\frac12W_\Phi^2-\frac{2\kappa _5^2}{3}W^2,
\qquad
\lambda_b=2W+\frac{\mu_b}{2}(\Phi-\Phi_b)^2.
\]

Write

\[
a=\frac{W_0}{v^2}>0,
\qquad
k=\frac{\kappa _5^2W_0}{3}>0,
\qquad
x=a(u+u_0),
\]

so that on the half-line

\[
\bar\Phi'=\frac{W_0}{v}\operatorname{sech}x,
\qquad
A'=-k\tanh x.
\]

The endpoint parameters are

\[
x_0=au_0>0,
\qquad
\tau=\tanh x_0=\sin(\Phi_b/v)\in(0,1),
\]

\[
h=2A'(0^+)=-2k\tau<0,
\qquad
\beta=\frac2{\mu_b}>0,
\qquad
\mu_b>2a.
\]

The equation and energy forms inherited from v293 are

\[
m_s[f,g]=\int_0^\infty w\,\bar f g\,du
          +p_0\beta\,\bar f_0g_0,
\qquad
w=\frac1{\bar\Phi'^2},                                      \tag{2.1}
\]

\[
q_s[f,g]=\int_0^\infty p
 (\bar f'g'+4ka\,\bar f g)\,du
 -p_0h\,\bar f_0g_0,
\qquad
p=\frac{e^{2A}}{\bar\Phi'^2}.                               \tag{2.2}
\]

The positive self-adjoint scalar operator satisfies

\[
q_s[f,g]=m_s[f,A_sg],
\qquad
\sigma(A_s)=[0,\infty),
\qquad
\sigma_{\rm point}(A_s)=\varnothing.                        \tag{2.3}
\]

The endpoint evaluation vector `e_b` is defined by

\[
m_s[e_b,f]=f(0),
\]

and the boundary Weyl response is

\[
G(z)=\langle e_b,(A_s-z)^{-1}e_b\rangle_{m_s}
    =\int_0^\infty\frac{d\rho(\lambda)}{\lambda-z},
\qquad d\rho\ge0.                                           \tag{2.4}
\]

## 3. Exact subordinate zero-energy solution

Introduce the dimensionless quantities

\[
\gamma=\frac{k}{a}=\frac{\kappa _5^2v^2}{3}>0,
\qquad
\xi=\frac{1-\tanh x}{2},
\qquad
s=\xi(0)=\frac{1-\tau}{2}\in(0,1/2).                       \tag{3.1}
\]

The zero-energy radial equation reduces, after its endpoint factors are
removed, to

\[
\xi(1-\xi)F_{\xi\xi}
+\gamma(1-2\xi)F_\xi
-\frac{\gamma}{\xi(1-\xi)}F=0.                             \tag{3.2}
\]

Define

\[
H_\gamma(\xi)
={}_2F_1(2,2\gamma+1;\gamma+2;\xi),                         \tag{3.3}
\]

\[
Y_\gamma(\xi)=\xi(1-\xi)H_\gamma(\xi).                   \tag{3.4}
\]

Because all hypergeometric series coefficients are positive for
`gamma>0` and `0<=xi<1`, `H_gamma` and `Y_gamma` are positive on the
physical interval.  The unique subordinate zero-energy solution normalized
to unit endpoint trace is

\[
\boxed{
y_0(u)=\frac{Y_\gamma(\xi(u))}{Y_\gamma(s)},
\qquad y_0(0)=1.}                                            \tag{3.5}
\]

As `u` tends to infinity, `xi` is asymptotic to `exp[-2a(u+u0)]`, and

\[
y_0(u)\sim C_0e^{-2au}.                                     \tag{3.6}
\]

This is a subordinate threshold solution, not a normalizable eigenstate in
the equation Hilbert space.  Consequently zero remains in the continuous
spectrum; no bounded inverse `A_s^{-1}` on the full Hilbert space is claimed.

## 4. Exact zero-resolvent boundary datum

Set

\[
S_\gamma(\xi)=\partial_\xi\log H_\gamma(\xi).               \tag{4.1}
\]

The contiguous hypergeometric identity gives

\[
S_\gamma(\xi)
=\frac{2(2\gamma+1)}{\gamma+2}
 \frac{{}_2F_1(3,2\gamma+2;\gamma+3;\xi)}
      {{}_2F_1(2,2\gamma+1;\gamma+2;\xi)}.                 \tag{4.2}
\]

The zero-energy endpoint denominator is

\[
\boxed{
\Delta_0=-2p_0\left[(a+k)\tau
       +a\,s(1-s)S_\gamma(s)\right]<0.}                     \tag{4.3}
\]

Every term inside the square brackets is positive.  Therefore the locked
zero-energy branch has no endpoint zero or pole, and the zero-resolvent trace
is finite.  Direct evaluation yields

\[
\boxed{
G(0)=
\frac{W_0(1-\tau^2)}
{2\left[(1+\gamma)\tau
 +\dfrac{1-\tau^2}{4}S_\gamma(s)\right]}>0.}                \tag{4.4}
\]

Equation (4.4) is an exact family formula, not a unique number.  The locked
chain does not instantiate `W0`, `gamma`, `tau`, and `mu_b`.  A numerical
value must therefore wait for a separately declared parameter point with
units and provenance.

The formula is independent of `mu_b`.  This is not an omission: the endpoint
eigenparameter term is proportional to spectral parameter and therefore
drops out at zero energy.  `mu_b` re-enters the derivative of the resolvent
and the static trace response below.

For the exactly reducible check `gamma=1`,

\[
\boxed{G(0)=\frac{W_0}{2}(1-\tau).}                          \tag{4.5}
\]

## 5. Capacity identity and sharp elementary bounds

The trace functional is continuous in the energy-form norm.  Its exact
capacity is

\[
\boxed{
G(0)^{-1}=\inf_{f(0)=1}q_s[f,f],
\qquad
G(0)=\sup_{f\ne0}\frac{|f(0)|^2}{q_s[f,f]}.}                \tag{5.1}
\]

The minimizer is precisely the normalized subordinate profile `y0`.  Since
the energy form contains the positive endpoint term `-p0 h |f0|^2`,

\[
0<G(0)\le\frac1{-p_0h}.                                     \tag{5.2}
\]

Define

\[
A_*=\frac{W_0(1-\tau)}{1+\gamma},
\qquad
B_*=\frac{W_0(1-\tau^2)}{2(1+\gamma\tau)}.                 \tag{5.3}
\]

Then the hypergeometric logarithmic derivative obeys bounds equivalent to

\[
\boxed{
\min(A_*,B_*)\le G(0)\le\max(A_*,B_*).}                    \tag{5.4}
\]

Both inequalities are strict unless `gamma=1`, where the two comparison
values coincide with (4.5).  These are family bounds, not observational
constraints.

## 6. Finite inverse second moment and static declared-trace response

Let

\[
I_n=\int_0^\infty\lambda^{-n}\,d\rho(\lambda).              \tag{6.1}
\]

The first moment is the just-derived zero response:

\[
I_1=G(0)<\infty.                                             \tag{6.2}
\]

Differentiating the Weyl response at zero gives the exact identity

\[
\boxed{
G'(0)=G(0)^2\left[
 \int_0^\infty w(u)y_0(u)^2\,du+p_0\beta
 \right].}                                                   \tag{6.3}
\]

Since `w` grows as `exp(2au)` and `y0^2` decays as `exp(-4au)`, the
bulk integrand decays as `exp(-2au)`.  Thus

\[
\boxed{I_2=G'(0)<\infty.}                                   \tag{6.4}
\]

The exact bulk quadrature is

\[
\boxed{
\int_0^\infty wy_0^2\,du
=\frac1{8a^3v^2Y_\gamma(s)^2}
 \int_0^s H_\gamma(\xi)^2\,d\xi.}                         \tag{6.5}
\]

The declared mathematical probe in v293 couples as `-T F(0)`.  To avoid an
unproved matter crosswalk, v294 denotes its response by `chi_tr`:

\[
\boxed{
\chi_{\rm tr}(Q^2)
=\frac{\kappa _5^4}{18Q^2}
 [G(0)-G(-Q^2)],\qquad Q^2>0.}                               \tag{6.6}
\]

Its static limit now exists:

\[
\boxed{
\chi_{\rm tr}(0)
=\frac{\kappa _5^4}{18}I_2
=\frac{\kappa _5^4G(0)^2}{18}
 \left[
  \frac{\displaystyle\int_0^sH_\gamma(\xi)^2d\xi}
       {8a^3v^2Y_\gamma(s)^2}
  +p_0\beta
 \right]<\infty.}                                           \tag{6.7}
\]

Unlike `G(0)`, this static susceptibility depends on `mu_b` through
`beta=2/mu_b`.

## 7. Liouville tail and threshold spectral density

After the exact Liouville transformation to coordinate `zeta`, the scalar
operator has an inverse-square tail

\[
V_L(\zeta)=\frac{\nu^2-\tfrac14}{\zeta^2}
+O\!\left(\zeta^{-2-\epsilon}L(\zeta)\right),               \tag{7.1}
\]

where `L` is a harmless slowly varying logarithmic factor and

\[
\boxed{
\nu=1+\frac ak
=1+\frac{3}{\kappa _5^2v^2}>1,
\qquad
\epsilon=\min\!\left(1,\frac{2a}{k}\right)>0.}             \tag{7.2}
\]

The subordinate zero-energy transfer coefficient is nonzero because
`Delta0<0`.  Matching the normalized endpoint solution to the Bessel/Jost
basis gives

\[
\boxed{
\rho'(\lambda)=C_\rho\lambda^\nu[1+o(1)]
\quad(\lambda\downarrow0),\qquad C_\rho>0.}                 \tag{7.3}
\]

With

\[
d_\infty=
\frac{v\,k^{-(\nu-1/2)}}
     {4W_0\cosh x_0\,Y_\gamma(s)},                          \tag{7.4}
\]

the exact threshold constant in the locked normalization is

\[
\boxed{
C_\rho=
\frac{d_\infty^2}
{2^{2\nu-1}\Gamma(\nu)^2\Delta_0^2}.}                      \tag{7.5}
\]

The result is a threshold theorem for the locked scalar boundary measure.
It is not a particle-production rate, a cosmological power spectrum, or an
empirical fit.

## 8. Complete inverse-moment ladder

By (7.3), the lower endpoint of `I_n` behaves as

\[
\int_0 \lambda^{\nu-n}\,d\lambda.
\]

The upper endpoint is controlled by the finite total boundary spectral
weight

\[
\int_0^\infty d\rho(\lambda)=\frac1{p_0\beta}<\infty.       \tag{8.1}
\]

Therefore

\[
\boxed{
I_n<\infty
\quad\Longleftrightarrow\quad
n<\nu+1=2+\frac ak.}                                        \tag{8.2}
\]

At equality, the divergence is logarithmic.  In particular,

```text
I1: always finite,
I2: always finite,
I3: finite if a>k,
I3: logarithmically divergent if a=k,
I3: power divergent if a<k.
```

This ladder distinguishes the finite static response from the smoothness of
its next derivative.  It also prevents a formal Taylor series from being
used past the moment actually controlled by the A1 tail.

## 9. Exact response shape and endpoint asymptotics

Define the rescaling-invariant response

\[
\boxed{
U(s)=1-\frac{G(-s)}{G(0)},\qquad s\ge0.}                     \tag{9.1}
\]

Using the positive Stieltjes measure,

\[
U(s)=\frac{s}{G(0)}
 \int_0^\infty\frac{d\rho(\lambda)}{\lambda(\lambda+s)}.   \tag{9.2}
\]

Consequently

\[
U(0)=0,
\qquad
0<U(s)<1\ (s>0),
\qquad
U'(s)>0,
\qquad
U''(s)<0.                                                     \tag{9.3}
\]

More generally, the derivatives alternate with the Bernstein/Stieltjes sign
pattern.  The two endpoint laws are

\[
\boxed{
U(s)=\frac{I_2}{G(0)}s+o(s)
\quad(s\downarrow0),}                                       \tag{9.4}
\]

\[
\boxed{
1-U(s)\sim\frac1{p_0\beta\,G(0)\,s}
\quad(s\to\infty).}                                        \tag{9.5}
\]

The manner in which the trace susceptibility approaches its static value is
fixed by `alpha=a/k`:

\[
\chi_{\rm tr}(0)-\chi_{\rm tr}(s)=
\begin{cases}
\Theta(s^{a/k}),&0<a/k<1,\\
\Theta(s\log(1/s)),&a/k=1,\\
\Theta(s),&a/k>1,
\end{cases}
\qquad s\downarrow0.                                       \tag{9.6}
\]

The three cases are the same trichotomy as the `I3` moment gate in section 8.

## 10. Exact reducible fixture

The family point

```text
gamma=1,
tau=1/2,
W0=1,
v=1,
mu_b=4
```

is a dimensionless algebraic verification fixture, not a physical parameter
choice.  It gives

\[
G(0)=\frac14,
\qquad
\int_0^\infty wy_0^2du=\frac{37}{72},                       \tag{10.1}
\]

\[
I_2=\frac{85}{1152},
\qquad
\chi_{\rm tr}(0)=\frac{85}{2304}
\quad(\kappa_5^2=3,\ \text{implied by }\gamma=1,\ v=1).    \tag{10.2}
\]

These rational values are used only to detect sign, factor, endpoint, and
normalization regressions.

## 11. Operator-domain caution at threshold

The scalar operator remains positive with continuous spectrum beginning at
zero.  Therefore

```text
A_s^{-1} is not a bounded operator on the full equation Hilbert space.
```

The finite quantities in this release are vector-specific:

\[
\langle e_b,A_s^{-1}e_b\rangle=I_1<\infty,
\qquad
\langle e_b,A_s^{-2}e_b\rangle=I_2<\infty.                  \tag{11.1}
\]

They do not justify an unrestricted inverse, nor do they remove the
unbounded-threshold caution in the canonical map
`Fhat=(3sqrt(2)/kappa_5^2)A_s^(1/2)F`.

## 12. Why no parent graph coefficient follows

v294 closes a one-channel boundary response of the same scalar quotient
already normalized in v293.  It does not add a second independently prepared
source.  The following types remain distinct:

```text
F                  induced-metric scalar coordinate,
Fhat               canonical coordinate for the same scalar polarization,
P_F                phase-space cotangent dual of Fhat,
e_b                 endpoint-evaluation vector,
T                  declared mathematical boundary trace probe,
parent Pi          independently prepared parent source coordinate.
```

Neither a finite `G(0)` nor a finite `chi_tr(0)` promotes `P_F`, `e_b`, or
`T` to the missing parent `Pi`.  Thus the v292/v293 source-rank terminal
theorem survives unchanged:

\[
\boxed{
\text{the locked minimal-A1 data determine one source line, not a complete
two-coordinate parent source matrix}.}                       \tag{12.1}
\]

No `B`, `D_phys`, active-constraint inertia, Schur defect, or physical
`lambda=3` is exported.

## 13. Claim firewall

v294 certifies

```text
canonical v293 parent: locked,
exact subordinate zero-energy profile: derived,
exact family formula for G(0): derived and positive,
unique numerical G(0): absent without a parameter point,
trace capacity identity and analytic bounds: derived,
threshold density exponent: nu=1+a/k,
inverse-moment ladder: I_n finite iff n<2+a/k,
static declared-trace response: finite and exactly reduced to one quadrature,
complete two-parent-coordinate source map: absent,
physical B: absent,
physical D_phys: absent,
physical lambda=3: absent,
empirical or detection claim: absent.
```

The state remains

```text
PREREGISTERED_TEST_READY
certified lambda_total=0+0i
remaining target=3+0i
```

This is a mathematical consistency and spectral-analysis result inside the
declared A1 model.  It is not external confirmation of a theory of
everything.

## 14. Best v295 move

The next productive step is not another formal inversion.  Keep v294 locked
and pursue two separate lanes:

1. declare a fully dimensioned A1 parameter point, with independent
   provenance and uncertainty intervals, then evaluate (4.4), (6.7), and the
   response curve with interval arithmetic; and
2. if a two-coordinate parent graph is still desired, declare a new
   covariant Ward-consistent matter/source action that provides a genuinely
   independent second source row.

The numerical lane must not be used to manufacture the missing source row,
and the new-source lane must not retroactively alter the v294 scalar spectral
measure.
