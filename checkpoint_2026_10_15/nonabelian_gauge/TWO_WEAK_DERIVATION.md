# Two weak fields and mixed BW: evaluated ordered one-loop pole

15 October 2026. Continuation of the literal seed action at source commit
`5c77bc5243b32da7b97a0ab86ea8571cb1cff5c0`. This is the coefficient linear
in `F=C_(7,seed)Y_nu*` and quadratic in electroweak gauge couplings in the
specified one-loop 1PI Green function. It is not a physical beta function,
finite amplitude, EOM reduction, or claim of unification or novelty.

## Result

Write `phi=tilde H`, `chi=L^c`, `T^I=tau^I/2`, `y_L=-1/2`,
`K=k_1+k_2`, and `p=r+h+K`. All momenta and polarizations are initially
generic and off shell, with outgoing `bar L,phi,A_1,A_2` waves and incoming
`N`. Let `E_i=slash epsilon_i`. For two generators `X_1,X_2`, stripping
their gauge couplings, the actual integrated matter-loop residue is

\[
\boxed{\mathscr R(X_1,X_2)
 =\frac12\{X_1,X_2\}R_2
 +\frac12[X_1,X_2]C_2,}
\]

\[
R_2=(\varepsilon_1\!\cdot\varepsilon_2)\not p
-\frac14\big([\not k_1,E_1]E_2+[\not k_2,E_2]E_1\big),
\]

\[
\boxed{\begin{aligned}
C_2={}&\frac14\big([\not k_1,E_1]E_2-[\not k_2,E_2]E_1\big)
 +\frac12[E_1,E_2]\not r\\
&+\frac23\big[((2k_1+k_2)\!\cdot\varepsilon_2)E_1
 -((k_1+2k_2)\!\cdot\varepsilon_1)E_2
 +(\varepsilon_1\!\cdot\varepsilon_2)(\not k_2-\not k_1)\big].
\end{aligned}}
\]

`R_2` is the previously evaluated Abelian kernel, while `C_2` is the new
ordered noncommuting combination calculated here. It is antisymmetric
under simultaneous exchange of legs 1 and 2; consequently the full
kernel is Bose symmetric. This statement includes arbitrary off-shell
polarizations, not just transverse gauge modes.

The full poles, with their literal action sign, are

\[
i\mathcal M_{W^IW^J}^{\rm UV}
=-\frac{i Fg^2}{16\pi^2\epsilon}\,
\bar u_L^{\,a}\mathscr R(T^I,T^J)_{ac}P_Lu_N\phi_c,
\]
\[
i\mathcal M_{BW^I}^{\rm UV}
=-\frac{i Fg'g}{16\pi^2\epsilon}\,
\bar u_L^{\,a}(y_LT^I)_{ac}R_2P_Lu_N\phi_c.
\]

Here the first polarization/momentum belongs to B in the mixed formula.
The external scalar tensor in the original `H*` variables is obtained by
right multiplication by epsilon. The matrix `F` acts in flavor space in
the displayed order; weak-matrix ordering is retained independently of flavor ordering.
The corresponding local counterterm action has the positive coefficient.

For SU(2), the generator identities give the explicit form

\[
\mathscr R(T^I,T^J)
=\frac{\delta^{IJ}}4\mathbf1_2 R_2
 +\frac{i\epsilon^{IJK}}2T^K C_2.
\]

For mixed BW, `[y_L 1,T^I]=0` and the coefficient is `y_LT^I R_2`.
This result includes both allowed assignments of the distinct gauge legs;
it carries no further factor of two or one half.

## Literal weak contractions and phases

The seed tensors in the phi basis are

\[
K^d_{ab;cd}=\delta_{ac}\delta_{bd},\qquad
K^x_{ab;cd}=\delta_{ad}\delta_{bc}.
\]

The source Yukawa contraction is `delta_ef`. If `C_be` is the weak matrix
on the conjugate fermion line and `S_df` the scalar matrix, summation of
the actual indices gives

\[
\sum_{bdef}K^x_{ab;cd}C_{be}S_{df}\delta_{ef}
=(SC^T)_{ac},
\quad
\sum_{bdef}K^d_{ab;cd}C_{be}S_{df}\delta_{ef}
=\delta_{ac}\operatorname{tr}(SC^T).
\]

For each gauge emission the kinetic matrices are `-X_i^T` on `chi` and
`X_i` on `phi`. The emitted one-gauge vertices both have the phase `-i`,
with the respective coupling and Lorentz numerator. The scalar kinetic
seagull vertex is

\[
+i g_1g_2\{X_1,X_2\}(\varepsilon_1\cdot\varepsilon_2).
\]

Order `ij` on one line means `i` is nearer Yukawa and `j` nearer the seed.
The source chain has phase `-1`. Each ordinary insertion adds
`i(-i)=+1`, while the extra scalar propagator and seagull add `i(i)=-1`.
After factoring the common `-Fg_1g_2`, the crossed factors are therefore:

| Diagram | Fermion matrix C | Scalar matrix S | Crossed factor, including relative phase |
|---|---|---|---|
| FF^ij | `(-X_j^T)(-X_i^T)` | `1` | `X_i X_j` |
| SS^ij | `1` | `X_j X_i` | `X_j X_i` |
| FS^ij | `-X_i^T` | `X_j` | `-X_j X_i` |
| SG | `1` | `{X_1,X_2}` | `-{X_1,X_2}(e_1.e_2)` |

The seagull Lorentz factor is kept outside `R_SG` below. The direct factor
for every row is the trace of the displayed crossed factor times `1_2`.
Those direct terms are UV finite individually at two gauge fields.
Tracelessness of a single generator is not used to remove them here.

## Complete seven-diagram routing

Define `D_a=(q-a)^2-m_H^2+i0`, `t=p-K`, and let the unit crossed seed be
`A_f^x(q)=(q-h-f)^2`. The unit direct seed is `A^d=t^2`; its trace factor
already supplies the weak contraction and is not to be multiplied by a
second trace weight. For `(i,j)=(1,2),(2,1)` the seven unweighted crossed
integrals are

\[
I_{FF}^{ij}=\int_q\frac{(q-h-K)^2(\not q-\not K)E_j
 (\not q-\not k_i)E_i\not q}
 {q^2(q-k_i)^2(q-K)^2D_p},
\]
\[
I_{SS}^{ij}=\int_q\frac{(q-h)^2\not q
 [(2p-2q-k_i)\cdot\varepsilon_i]
 [(2p-2q-2k_i-k_j)\cdot\varepsilon_j]}
 {q^2D_pD_{p-k_i}D_t},
\]
\[
I_{FS}^{ij}=\int_q\frac{(q-h-k_i)^2(\not q-\not k_i)E_i\not q
 [(2p-2q-k_j)\cdot\varepsilon_j]}
 {q^2(q-k_i)^2D_pD_{p-k_j}},\qquad
I_{SG}=\int_q\frac{(q-h)^2\not q}{q^2D_pD_t}.
\]

The direct integrals replace the crossed squared numerator by `t^2`.
Boxes then have loop degree at most three over four quadratic
denominators, and the seagull has loop degree one over three. Both are
UV convergent by at least one power. The direct sector's nonzero weak
traces therefore produce no simple UV pole. The crossed residue has
momentum degree one and is independent of `m_H^2`: a mass insertion
lowers the UV degree by two.

## Evaluating the commutator combination

The following primitive moments and individual Lorentz integrals are
reused from `two_gauge_leg/TWO_B_DERIVATION.md`; they are explicitly
attributed inputs/methodology, not a claim of a new independent
integration of the old Abelian problem. They follow from the
Feynman-simplex volume and first moment, in metric `(+---)`, stripping
`i/(16 pi^2 epsilon)`:

\[
\begin{array}{ll}
N=3:& \operatorname{Res}q^\mu q^\nu=g^{\mu\nu}/4,\qquad
 \operatorname{Res}q^\mu q^\nu q^\rho=
 \operatorname{sym}(gS)^{\mu\nu\rho}/12,\\
N=4:& \operatorname{Res}q^\mu q^\nu q^\rho q^\sigma=
 \operatorname{sym}(gg)^{\mu\nu\rho\sigma}/24,\qquad
 \operatorname{Res}q^5=\operatorname{sym}(ggS)/96,
\end{array}
\]

where `S` is the sum of denominator shifts. Lower ranks are UV finite
in these integrals. Only crossed numerator degrees four/five for boxes
and two/three for the seagull survive. Cancelling the seed's `q^2`
term against the unshifted fermion denominator also gives the following
individual residues. For `w_F=p+k_i-h`,

\[
R_{FF}^{ij}=\frac12(\not K E_jE_i+E_i\not k_iE_j)
-\frac16(\not w_F E_jE_i+E_i\not w_FE_j+E_jE_i\not w_F).
\]

For `u_i=2p-k_i`, `v_j=2p-2k_i-k_j`, `w_S=3p-2k_i-k_j-h`,

\[
R_{SS}^{ij}=\frac13[(e_i\cdot e_j)\not w_S
 +(w_S\cdot e_i)E_j+(w_S\cdot e_j)E_i]
-\frac12[(u_i\cdot e_i)E_j+(v_j\cdot e_j)E_i].
\]

For `u_j=2p-k_j`, `w_M=2p-k_j-h`,

\[
R_{FS}^{ij}=\frac12\not k_iE_iE_j
-\frac16(E_jE_i\not w_M+\not w_ME_iE_j)
+\left(\frac{w_M\cdot e_j}{3}-\frac{u_j\cdot e_j}{2}\right)E_i,
\qquad R_{SG}=\frac12(\not p+\not r).
\]

The literal weak table gives the two ordered coefficients

\[
U_{12}=R_{FF}^{12}+R_{SS}^{21}-R_{FS}^{21}
 -(e_1\cdot e_2)R_{SG},\quad
U_{21}=R_{FF}^{21}+R_{SS}^{12}-R_{FS}^{12}
 -(e_1\cdot e_2)R_{SG}.
\]

Thus the new integral to evaluate is explicitly

\[
C_2=U_{12}-U_{21}
=R_{FF}^{12}-R_{FF}^{21}-R_{SS}^{12}+R_{SS}^{21}
 +R_{FS}^{12}-R_{FS}^{21}.
\]

Substituting the integrated expressions above and reducing the ordered
Clifford products yields the boxed `C_2`. This evaluation precedes the
action comparison. It retains the part invisible in the Abelian sum.
The symmetric combination `U12+U21=R2` reproduces the previous result.
For an additional reproducible derivation, the checker expands each
routed numerator and applies the primitive tensor functional before
forming the new weak sum. It verifies the closed expression on every
coefficient allowed by UV power counting, not only sampled momenta.

## Explicit local action comparison and nonlinear sign

Introduce the coupling-valued matrix connection and curvature

\[
\mathcal A_\mu=g'y_L\mathbf1 B_\mu+gT^I W_\mu^I,
\quad D_\mu=\partial_\mu+i\mathcal A_\mu,
\quad \mathcal F_{\mu\nu}=\partial_\mu\mathcal A_\nu
 -\partial_\nu\mathcal A_\mu+i[\mathcal A_\mu,\mathcal A_\nu].
\]

This sign follows from `[D_mu,D_nu]=i mathcal F_mu nu`. Consequently
the component weak field strength is
`W^I_mu nu=partial_mu W^I_nu-partial_nu W^I_mu-g epsilon^IJK W^J_mu W^K_nu`.
The adjoint derivative is `D_nu mathcal F=partial_nu mathcal F
+i[mathcal A_nu,mathcal F]`. The positive cancelling action is

\[
\mathcal L_{ct}=\frac F{16\pi^2\epsilon}\left[
O_J+\frac i4(D_\rho\bar L)\sigma^{\mu\nu}\gamma^\rho
 \mathcal F_{\mu\nu}\phi N_L
+\frac13\bar L\gamma_\mu(D_\nu\mathcal F^{\nu\mu})\phi N_L
\right]+\mathrm{h.c.},
\]

where the six-placement operator from the original action audit is

\[
\begin{aligned}
O_J={}&\frac56\bar L\phi A(-\partial^2)N_L
+\frac13(-D^2\bar L)\phi A N_L
+\frac16\bar L(-D^2\phi)A N_L\\
&-\frac16(\mathcal K\bar L)\phi(-\partial^2)N_L
+\frac12m_H^2\bar L\phi A N_L
+\frac12m_H^2(\mathcal K\bar L)\phi N_L,\\
A={}&i\!\not\partial,\qquad
\mathcal K\bar L=-i(D_\mu\bar L)\gamma^\mu.
\end{aligned}
\]

Its quadratic Fourier kernels, with couplings stripped, are explicitly

\[
K_J=\tfrac12\{X_1,X_2\}(e_1\cdot e_2)\not p,
\]
\[
K_\sigma=\tfrac14\big([X_1,X_2][E_1,E_2]\not r
 -X_2X_1[\not k_1,E_1]E_2-X_1X_2[\not k_2,E_2]E_1\big),
\]
\[
K_{D\mathcal F}=\tfrac13[X_1,X_2]\big[
 ((2k_1+k_2)\cdot e_2)E_1-((k_1+2k_2)\cdot e_1)E_2
 +(e_1\cdot e_2)(\not k_2-\not k_1)\big].
\]

`K_sigma` includes both the row-derivative connection and nonlinear
curvature. `K_DF` includes both the derivative of the nonlinear curvature
and its adjoint connection. Their sum is precisely the integrated
`mathscr R`. An Abelian field-strength replacement would omit these
commutator terms and fails for distinct weak generators.

## 1PI topology completeness and Ward identity

The ungauged seed/Yukawa graph has two vertices and two internal matter
lines, hence one loop. In the literal action the composite differentiated
by the seed is neutral, so it supplies no gauge contact; the source
Yukawa has no derivative contact. At order two gauge couplings, two
ordinary matter-gauge insertions produce FF(2), SS(2), FS(2); one scalar
seagull produces SG(1). These exhaust the seven 1PI possibilities.

At the same order a matter one-gauge attachment joined to a pure weak
triple-gauge vertex can produce two external W fields. Its single
internal gauge propagator is a bridge: cutting it disconnects both
external W fields from the matter loop. It is therefore 1PR and is not
part of this declared 1PI kernel. A quartic pure-gauge vertex alone
cannot connect to this gauge-free seed/Yukawa graph at the selected
coupling order. Adding any further gauge attachment needed to avoid a
bridge exceeds the two-coupling monomial or adds another loop. There
are no ghost couplings to the seed/Yukawa matter vertices. External
matter-leg radiation is likewise 1PR. This is a valence argument for the
selected monomial; it is not an enumeration of the complete EFT.

Let `R_1(p,h;k,e)` be the previously integrated one-gauge Lorentz kernel,
where `r=p-h-k`. The new local kernel obeys the nonabelian 1PI Ward
identity

\[
\begin{aligned}
\mathscr R(X_1,X_2)\big|_{e_1=k_1}
={}&X_2X_1 R_1(p,h+k_1;k_2,e_2)
-X_1X_2 R_1(p,h;k_2,e_2)\\
&+[X_1,X_2]R_1(p,h;k_1+k_2,e_2).
\end{aligned}
\]

The last term is the nonabelian lower-gauge kernel, not a new two-gauge
loop topology. The old Abelian Ward difference is recovered when the
generators commute. The mass-dependent pieces of the three one-gauge
terms cancel. The second-leg identity follows by Bose symmetry. This
checker tests the integrated local Ward relation only; it does not claim
a new unintegrated Ward derivation.

## Reproduction, provenance, and limits

Run `python3 -I -S -B -O verify_two_weak.py --json`. Only the standard
library is used; the script imports no repository module. Gaussian weak
matrices, Clifford coefficients, and all comparisons use exact rational
arithmetic and explicit exceptions under `-O`.

The default run passes **7168 exact checks** over 256 coefficient cases
and 12 generic rational points. Of these, 2304 are explicitly labeled
inherited Lorentz/Abelian reproductions, 4768 are other positive checks
for this ordered calculation, and 96 are mutation rejections. The
machine-readable breakdown and two explicit sets of diagram residues
are saved in `TWO_WEAK_CHECKS.json`.

The Lorentz polynomial/moment engine and compact seven individual
residues are adapted from `two_gauge_leg/verify_two_b.py` at the stated
source commit. The script header records that dependence, and the result
counts label their reproduction as inherited. New tests cover all weak
indices for all sixteen ordered pairs of the three SU(2) generators and
hypercharge; all 256 independent linear-momentum/bilinear-polarization
coefficients of the new commutator; all weak and mixed kernels; generic
rational points; Bose symmetry; the nonabelian Ward recursion; and
rejections of missing conjugate transpose, reversed FF color order,
missing commutator, wrong nonlinear curvature sign, missing curvature
pieces, missing scalar seagull, and a half-size mixed seagull.

The 256 coefficient substitutions are evaluations of the already
defined UV polynomial. They do not set original exceptional integrals
to scaleless zero, separate IR poles by special kinematics, or determine
finite terms. Direct-seed and scalar-mass statements are consequences of
UV degree; the zero residue cannot remeasure the direct seed weight.
Four-dimensional Clifford identities suffice for the simple pole only.

The action passes this further selected two-gauge 1PI comparison.
Nothing here establishes complete operator mixing, a physical running
coefficient, finite photon amplitude, decay rate, or ultraviolet
completion. The independent audit is a separate artifact; this producer
checker must not be called independent of the earlier Lorentz engine.
