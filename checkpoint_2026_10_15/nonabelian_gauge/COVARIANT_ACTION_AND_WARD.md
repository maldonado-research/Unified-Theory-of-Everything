# Matrix curvature, two-gauge completion, and the cubic action contact

15 October 2026. Independent literal-action Fourier audit of the selected
dimension-seven `bar L phi N_L` counterterm, where `phi=tilde H`.
The zero- and one-gauge pole coefficients are inputs from the prior checkpoint.
This note derives their non-Abelian local-action completion; it does not
integrate a loop or compute physical running.

## Conventions and matrix-valued local action

Use `(+---)`, `sigma^{mu nu}=i[gamma^mu,gamma^nu]/2`,
`D_mu=partial_mu+i A_mu`, and

\[
\mathcal A_\mu=g'y_LB_\mu1_2+gT^IW_\mu^I,
\qquad y_L=-\tfrac12,\qquad T^I=\tau^I/2.
\]

Both `L` and `phi` carry this fundamental connection. The row derivative is
`D_mu bar L=partial_mu bar L-i bar L A_mu`; `N` is neutral. Directly expanding
the commutator of two covariant derivatives gives

\[
[D_\mu,D_\nu]=i\mathcal F_{\mu\nu},\quad
\mathcal F_{\mu\nu}=\partial_\mu\mathcal A_\nu
-\partial_\nu\mathcal A_\mu+i[\mathcal A_\mu,\mathcal A_\nu],
\quad
\mathcal D_\nu\mathcal F^{\nu\mu}
=\partial_\nu\mathcal F^{\nu\mu}
+i[\mathcal A_\nu,\mathcal F^{\nu\mu}].
\]

For `U=exp(-i alpha)`, the transformations are
`phi -> U phi`, `bar L -> bar L U^{-1}`,
`A_mu -> U A_mu U^{-1}+i(partial_mu U)U^{-1}`,
and `F -> U F U^{-1}`. Infinitesimally,
`delta A_mu=partial_mu alpha+i[A_mu,alpha]`.
The plus sign in `i[A,A]` follows from the stated `D` convention.
In components, this convention gives
`F^K_mu nu=partial_mu W^K_nu-partial_nu W^K_mu
-g epsilon^{IJK}W^I_mu W^J_nu`, after factoring the overall `g`.

Let `F_flavor=C_(7,seed)Y_nu*`, retaining this flavor order. The candidate
counterterm action is

\[
\mathcal L_{\rm ct}=\frac{(F_{\rm flavor})_{\alpha i}}
{16\pi^2\epsilon}\left[O_J+
\frac14 i(D_\rho\bar L)\sigma^{\mu\nu}\gamma^\rho
\mathcal F_{\mu\nu}\phi N_L
+\frac13\bar L\gamma_\mu(\mathcal D_\nu
\mathcal F^{\nu\mu})\phi N_L\right]_{\alpha i}
+\mathrm{h.c.}
\]

The placement `(D bar L) F phi` matters. Weak matrices commute with the
gamma matrices but do not commute with each other. `O_J` uses the same six
placements as the prior literal-action audit, with full electroweak `D`:

\[
\begin{aligned}
O_J={}&\frac56\bar L\phi\,i\!\not\partial(-\partial^2)N_L
+\frac13(-D^2\bar L)\phi\,i\!\not\partial N_L
+\frac16\bar L(-D^2\phi)\,i\!\not\partial N_L\\
&-\frac16(\mathcal K\bar L)\phi(-\partial^2)N_L
+\frac12m_H^2\bar L\phi\,i\!\not\partial N_L
+\frac12m_H^2(\mathcal K\bar L)\phi N_L,
\qquad \mathcal K\bar L=-i(D_\mu\bar L)\gamma^\mu.
\end{aligned}
\]

These monomials reproduce the known one-B and one-W action kernels. They
are not a claim that every counterterm in the complete theory is listed.

## Functional two-gauge derivative

Take labelled outgoing modes
`A_j=X_j e_j exp(+i k_j.x)`, outgoing matter waves
`bar L~exp(+i r.x)`, `phi~exp(+i h.x)`, and incoming
`N~exp(-i p.x)`. Thus `p=r+h+sum k_j`. The `X_j` include the relevant
couplings; for a B mode `X_j=g'y_L 1_2`, and for `W^I`, `X_j=gT^I`.
Every kernel below is an **action kernel**. The corresponding vertex is
`+i` times the kernel and the common `F_flavor/(16 pi^2 epsilon)`.
The right projector `P_L` is implicit.

The known one-mode kernel is `K_1=X R_1`, where

\[
\begin{aligned}
R_1(r,h;k,e)={}&\left[-\frac{(2r+k)\cdot e}{3}
+\frac{(2h+k)\cdot e}{6}\right]\not p
+\left(\frac{p^2}{6}-\frac{m_H^2}{2}\right)\not e\\
&+\frac14[\not k,\not e]\not r
+\frac13[(k\cdot e)\not k-k^2\not e].
\end{aligned}
\]

Directly differentiating the two-mode terms of the matrix action gives

\[
\boxed{\begin{aligned}
K_2={}&\frac12\{X_1,X_2\}(e_1\cdot e_2)\not p\\
&+\frac14\left\{[X_1,X_2][\not e_1,\not e_2]\not r
-X_2X_1[\not k_1,\not e_1]\not e_2
-X_1X_2[\not k_2,\not e_2]\not e_1\right\}\\
&-\frac13[X_1,X_2]\not V,
\end{aligned}}
\]

\[
V=[(k_1+2k_2)\cdot e_1]e_2
-[(2k_1+k_2)\cdot e_2]e_1
+(e_1\cdot e_2)(k_1-k_2).
\]

The first line is the sum of the quadratic connections in the row and
scalar `D^2` terms. The second line contains both the nonlinear curvature
with the differentiated row and the two assignments of the row connection
to the linear curvature. The last line includes both `partial F_quad` and
`i[A,F_linear]`; dropping either produces the wrong coefficients of `V`.
Each derivative with respect to the two distinct mode variables is already
taken. An additional identical-field factor is not present.

Equivalently, write

\[
K_2=\frac12\{X_1,X_2\}R_2
+\frac12[X_1,X_2]C_2,
\]

\[
R_2=(e_1\cdot e_2)\not p
-\frac14\big([\not k_1,\not e_1]\not e_2
+[\not k_2,\not e_2]\not e_1\big),
\]

\[
C_2=\frac14\big([\not k_1,\not e_1]\not e_2
-[\not k_2,\not e_2]\not e_1\big)
+\frac12[\not e_1,\not e_2]\not r-\frac23\not V.
\]

`R_2` is symmetric and `C_2` antisymmetric under exchange of the two entire
gauge labels. Their products with the anticommutator and commutator are
therefore Bose symmetric. Pure B and mixed BW have `[X_1,X_2]=0`; they
reduce to the already fixed Abelian Lorentz polynomial times `X_1 X_2`.
For two W fields,
`{T^I,T^J}=delta^{IJ}1_2/2` and
`[T^I,T^J]=i epsilon^{IJK}T^K`.
The antisymmetric Lorentz polynomial is therefore essential when `I!=J`.

## Non-Abelian Ward recursion

At two gauge fields the functional variation gives

\[
\boxed{\begin{aligned}
K_2(r,h;k_1,k_2;k_1,e_2)={}&
X_2X_1R_1(r,h+k_1;k_2,e_2)\\
&-X_1X_2R_1(r+k_1,h;k_2,e_2)\\
&+[X_1,X_2]R_1(r,h;k_1+k_2,e_2).
\end{aligned}}
\]

The left side means replacing the first polarization by its momentum.
The last term arises from `delta A=i[A,alpha]` in the lower one-gauge
action. Its sign is positive with the outgoing Fourier convention above.
Both shifted terms and the merged-gauge term have the same incoming `p`.
Deleting the last term is an Abelian recursion and fails for noncommuting
generators. Exchanging the gauge labels gives the other identity.

If two local two-gauge degree-one kernels agree below two gauges and
obey this recursion, their difference is transverse in both independent
gauge momenta. Every monomial of such a difference contains both momenta
and has degree at least two. The degree-one difference therefore vanishes.
This is a conditional local uniqueness statement, not a second loop
calculation. It requires generic independent off-shell momenta, locality,
and a common complete lower-point sector.

## Cubic contact derived from the same action

`O_J` is at most quadratic in gauge fields. At three gauges only the row
connection times nonlinear curvature and the adjoint connection times
nonlinear curvature remain. Direct differentiation gives the constant
momentum-degree-zero action kernel

\[
\boxed{\begin{aligned}
K_3={}&-\frac14\sum_{i<j;\,\ell\notin\{i,j\}}
X_\ell[X_i,X_j][\not e_i,\not e_j]\not e_\ell\\
&-\frac13\sum_{i<j;\,\ell\notin\{i,j\}}
[X_\ell,[X_i,X_j]]
\left[(e_\ell\cdot e_i)\not e_j
-(e_\ell\cdot e_j)\not e_i\right].
\end{aligned}}
\]

There are three unordered pairs in each sum. These sums already include
all assignments of the three distinct external gauge fields. `K_3` is
Bose symmetric. Three B modes and two B modes with one W give zero.
For one B mode and two W modes, label the B mode by `ell`; the only
surviving term is

\[
-\frac14 X_B[X_{W_1},X_{W_2}]
[\not e_{W_1},\not e_{W_2}]\not e_B.
\]

The three-gauge Ward identity reads

\[
\begin{aligned}
K_3(k_1,e_2,e_3)={}&K_2(r,h+k_1;2,3)X_1
-X_1K_2(r+k_1,h;2,3)\\
&+K_2(r,h;[X_1,X_2],k_1+k_2,e_2;X_3,k_3,e_3)\\
&+K_2(r,h;X_2,k_2,e_2;[X_1,X_3],k_1+k_3,e_3).
\end{aligned}
\]

The notation in the last two lines specifies the generator, momentum,
and polarization of each remaining mode. The identity holds for each
choice of contracted leg. The Jacobi identity is needed when reducing
the nested commutators. A degree-zero local three-gauge ambiguity that
agrees through two gauges is separately transverse and vanishes.

The displayed action has no fourth-gauge contact. That is an action
statement; any assertion about higher-point loop ultraviolet finiteness
requires the corresponding power-counting and completeness assumptions.
In particular, none of these local statements determines a finite
three-gauge amplitude.

## Independent exact verification

Run `python3 -I -S -B -O verify_covariant_action.py` in this directory.
The standalone checker uses only the standard library, exact rational
real and imaginary parts, explicit Dirac/Pauli matrices, and polynomial
coefficients of distinct labelled Fourier modes. It imports no project
code and no loop residue. The covariant row/scalar derivatives, curvature,
and adjoint derivative are applied as polynomial differential operations
before comparison with the boxed simplified kernels.

Coverage includes every two-gauge momentum/polarization coefficient
(256 coefficient tests, using `E12 E21=E11` and `E21 E12=E22` to separate
the two ordered matrix words), all 64 cubic polarization coefficients for
non-orthogonal weak generators, and 12 additional rational off-shell
kinematic points. The latter test all six cubic label permutations,
both two-gauge Ward contractions, all three cubic Ward contractions,
Jacobi, and Abelian/mixed reductions. Deliberately wrong curvature and
adjoint signs, reversed row/curvature ordering, an omitted antisymmetric
sector, and an omitted cubic are rejected at noncommuting generic points.
The adjacent JSON file records the completed count and all negative
controls. These checks establish the stated local-action identities;
they must not be described as independent loop integrations, a full
anomalous-dimension calculation, or a new observable.

Completed result: **642 exact checks passed**, including **56 rejected
negative controls**. The test command succeeds with Python assertions
disabled (`-O`); failures use explicit exceptions.
