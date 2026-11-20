# Action-level field redefinition for the selected residual packet

15 October 2026. Derived against candidate `89f5867453bc7ed706d249dc82b393245b28d7ea`.
The inputs are the literal action and conventions of
`checkpoint_2026_10_15/operator_map/ACTION_SIGN.md` and the exact residual
certificate (5) of `electroweak_source/SELECTED_PACKET_REDUCTION.md`.
This note supplies the field substitutions implementing that certificate.
It is an action identity at first order in the selected coefficient, not a
new loop integration or a complete physical counterterm vector.

## Statement and orientation of the substitution

Write `kappa=1/(16 pi^2 epsilon)` and regard `kappa F` as first order.
For a formal derivation any real bookkeeping parameter may replace `kappa`.
All new fields occur on the right side of the substitutions below:

\[
\begin{split}
 L_{\alpha,\mathrm{old}}^a&=L_\alpha^a-\kappa\ell_\alpha^a,\\
 (\phi^a)^*_{\mathrm{old}}&=(\phi^a)^*-\kappa t_a,
 &\phi^a_{\mathrm{old}}&=\phi^a-\kappa t_a^*,\\
 B_{\mu,\mathrm{old}}&=B_\mu+\frac{\kappa g_y}{6}
       (j_\mu+j_\mu^\dagger),\\
 N_{Li,\mathrm{old}}&=N_{Li}-\kappa(x_i+y_i^c),
 &N_{Ri,\mathrm{old}}&=N_{Ri}-\kappa(y_i+x_i^c).
\end{split}                                                    \tag{6}
\]

The adjoint lepton substitution is the Hermitian conjugate of the first.
The explicit local coefficients are given below. No inverse mass occurs.
The map has a local perturbative inverse through this order, obtained by
reversing the shift signs and evaluating the coefficients on the old fields.
Equivalently, `N_old=N-kappa(w+w^c)` with `w=x+y`. This preserves `N=N^c`.
For `phi=epsilon H*`, `epsilon^12=+1` and `epsilon^2=-1`, the scalar
substitution is `H_old=H+kappa epsilon t^T`.

Let `S_0` be the **gauge-invariant renormalizable** part of the reference
action, including any explicitly retained renormalizable spectator
interactions, and let `S_ct=kappa integral(P_F+h.c.)`. The original
dimension-seven seed is excluded from `S_0`. Use precisely the sources
of this `S_0` in (5). Substitution gives

\[
 S_0[\mathrm{old}]+S_{\rm ct}[\mathrm{old}]
 =S_0[\mathrm{new}]+\kappa\int(P_{\rm red}+P_{\rm red}^\dagger)
 +O((\kappa F)^2),                                             \tag{7}
\]

modulo boundary terms and the Bianchi identity already used in (5).
`P_red` is exactly (1) of the input note. Sources are not discarded. For
complex `F`, the notation `O((kappa F)^2)` includes mixed `F,F*` products.
The direction `old=new-kappa coefficient` is essential: its inverse would
add the residual instead of removing it. Equation (7) is exact at first
order for the displayed `S_0+S_ct`. In the full EFT, varying a further
operator with coefficient `C_d=O(Lambda^{-(d-4)})` also produces
`O(kappa F C_d)=O(kappa Lambda^{-(d-1)})`, which is a separate inverse-scale
remainder, not generally a two-loop or quadratic-in-kappa remainder.
In particular varying the original seed gives `kappa F C_seed` at
`Lambda^-6`. These terms are outside this packet through `Lambda^-3`
under the fixed retained-mass counting. They can instead be retained by
using their full variational sources consistently in both sides. This
does not remove separate `C_6 C_5=O(Lambda^-3)` effects of other field
redefinitions.

## Action derivatives and their signs

Use `D_mu=partial_mu+i g_y B_mu` on `L,phi`, its opposite connection on
their adjoints, `A=i gamma^mu partial_mu`,
`K bar L=-i(D_mu bar L)gamma^mu`, and metric `(+---)`. The weak component
index is `a=1,2`; active and sterile indices are `alpha` and `i,j`.
Repeated indices in explicitly written sums have precisely that order.
Only the declared U(1) connection is retained; no SU(2) completion is claimed.
`M` is real positive diagonal as in the reference; the formulas also work
for constant real symmetric `M` with the ordered matrix products shown.

The source notation in the certificate has the following variational
meaning:

\[
\begin{split}
 e_{\bar L,\alpha}^a&=K\bar L_\alpha^a-\bar S_{L,\alpha}^a,
 &E_{L,\alpha}^a&=i\!\not D L_\alpha^a-S_{L,\alpha}^a,\\
 e_\phi^a&=-D^2\phi^a-(m_H^2+2\lambda\rho)\phi^a+S_\phi^a,
 &e_B^\mu&=\partial_\nu B^{\nu\mu}-g'J_Y^\mu,\\
 e_{Ri}&=A N_{Li}-M_{ij}N_{Rj}-S_{Ri},
 &e_{Li}&=A N_{Ri}-M_{ij}N_{Lj}-S^N_{Li}.
\end{split}
\]

In particular `e_phi=delta S_0/delta phi*`, with the sign shown, not its
negative. For the truncated Yukawa action, the explicit sources are

\[
\begin{split}
 S_{L,\alpha}^a&=Y_{\alpha j}\phi^a N_{Rj},
 &\bar S_{L,\alpha}^a&=\bar N_{Rj}Y_{\alpha j}^*(\phi^a)^*,\\
 S_\phi^a&=-\bar N_{Rj}Y_{\alpha j}^* L_\alpha^a,
 &S^N_{Li}&=Y_{\alpha i}^*(\phi^a)^*L_\alpha^a,
 &S_{Ri}&=(S^N_{Li})^c.
\end{split}                                                    \tag{8}
\]

Here `bar N_R=overline(N_R)`, and the sums in (8) are understood.
The minus in `S_phi` follows by differentiating
`-bar N_R Y^dagger phi^dagger L`. For additional retained renormalizable interactions,
replace these by the total variational sources. Higher-operator sources
can also be included, but then their inverse-scale order and the matching
terms in the action must be kept consistently as described after (7). Set `lambda=0` for the
literal potential in `ACTION_SIGN.md`. The gauge current includes every
retained charged field and the scalar connection terms.

For compatible Grassmann variations the integrated first variation is

\[
 \delta S_0=\int\!\left[
 e_{\bar L,\alpha}^a\delta L_\alpha^a
 +\delta\bar L_\alpha^a E_{L,\alpha}^a
 +\delta(\phi^a)^* e_\phi^a+(e_\phi^a)^*\delta\phi^a
 +e_B^\mu\delta B_\mu
 +\overline{\delta N_i}(e_{Ri}+e_{Li})\right].                 \tag{9}
\]

The last term has **no** extra `1/2`. Variation of both slots of the
Majorana quadratic action, with integration by parts in its kinetic
piece, cancels its action factor `1/2`. The Yukawa source and its charge
conjugate then give the Majorana equation `E_N=e_R+e_L`, with `E_N^c=E_N`.
The two terms for each complex boson or lepton are genuine conjugates;
one must not append another duplicate conjugate to (9).

## Explicit coefficients for the light fields

\[
\boxed{\ell_\alpha^a=\sum_i F_{\alpha i}\left[
 \frac{i}{3}\gamma^\mu(D_\mu\phi^a)A N_{Li}
 +\frac16\phi^a A^2N_{Li}
 +\frac{m_H^2}{2}\phi^a N_{Li}\right]}                         \tag{10}
\]

\[
\boxed{t_a=\frac16\sum_{\alpha i}F_{\alpha i}
 \bar L_\alpha^a A N_{Li},\qquad
 j_\mu=\sum_{\alpha i a}F_{\alpha i}
 \bar L_\alpha^a\phi^a\gamma_\mu N_{Li}.}                     \tag{11}
\]

`ell` is left chiral and transforms as `L`. `t_a` transforms as `phi*`;
its complex conjugate is the conjugate of the entire Grassmann-even
bilinear, including reversed fermion order. Thus, explicitly,
`t_a^*=(1/6) sum F* overline(A N_Li) L_alpha^a`.
The gauge shift is real. Its positive sign comes from the negative
`-g_y j_mu e_B^mu/6` term in (5).

## Fully expanded sterile coefficients

The following are columns, so that (6) is an actual field shift rather
than a formal instruction to delete equations of motion. Let
`c=2m_H^2/3+lambda rho/3`. The left- and right-chiral columns are

\[
\begin{split}
 x_i={}&\frac56\sum_{\alpha a}(F M^2)^*_{\alpha i}
                  (\phi^a)^* L_\alpha^a
       +c\sum_{\alpha a}F_{\alpha i}^*(\phi^a)^*L_\alpha^a\\
 &-\frac56\sum_{\alpha a}F_{\alpha i}^*\left[
       (\phi^a)^*D^2L_\alpha^a
       +2(D_\mu\phi^a)^*D^\mu L_\alpha^a
       +(D^2\phi^a)^*L_\alpha^a\right]\\
 &+\frac{i}{6}\sum_{\alpha a}F_{\alpha i}^*
             (\phi^a)^*\gamma^\mu D_\mu S_{L,\alpha}^a
  -\frac{i}{6}\sum_{\alpha a}F_{\alpha i}^*
             (D_\mu\phi^a)^*\gamma^\mu S_{L,\alpha}^a\\
 &-\frac16\sum_{\alpha a}F_{\alpha i}^*(S_\phi^a)^*L_\alpha^a
  -\frac{g_y}{12}\sum_{\alpha a}F_{\alpha i}^*
             (\phi^a)^*\sigma^{\mu\nu}L_\alpha^a B_{\mu\nu},\\[2mm]
 y_i={}&\frac{5i}{6}\sum_{\alpha a}(F M)^*_{\alpha i}\gamma^\mu
       \left[(\phi^a)^*D_\mu L_\alpha^a
                +(D_\mu\phi^a)^*L_\alpha^a\right]
  -\frac16\sum_{\alpha a}(F M)^*_{\alpha i}
                 (\phi^a)^*S_{L,\alpha}^a.
\end{split}                                                    \tag{12}
\]

`x_i=P_L x_i` and `y_i=P_R y_i`. Charge conjugation in (6) conjugates all
coefficients and fields in (12), with `psi^c=C bar psi^T` and the fixed
`C=i gamma^2 gamma^0`; it introduces no inverse mass or extra flavor
transpose. Formula (12), followed by that explicitly defined operation,
fixes every component of `N_old`. In particular `(FM)*`, not a reordered
`M F*`, appears. The constant-coefficient derivative does not act on `F`
or `M`; it does act on both factors in each neutral composite.

For (8) all expressions in (10)-(12) are local polynomials in the
reference fields and their derivatives. Leaving total sources displayed
allows additional spectator interactions without silently dropping their
induced operators. For example the two source-derivative terms in `x`
are explicitly `+i F* phi* gamma D(Y phi N_R)/6`
and `-i F*(D phi)* gamma Y phi N_R/6` in the truncated theory.

## Integrated-adjoint proof and Majorana completion

Define neutral rows

\[
 W_i=\sum_{\alpha a}F_{\alpha i}\bar L_\alpha^a\phi^a,
 \qquad V_i=\sum_{\alpha a}F_{\alpha i}\bar S_{L,\alpha}^a\phi^a,
 \qquad K_0 U=-i(\partial_\mu U)\gamma^\mu.
\]

These are row spinors; no adjoint is hidden in the letters `W,V`.
The two integrated derivative identities required here are

\[
 \int(K e_{\bar L})\,\phi A N_L
 =\int e_{\bar L}\,i\!\not D(\phi A N_L),
 \qquad \int U A E=\int(K_0 U)E.                              \tag{13}
\]

The charge of the first composite is `y_L`; the second composite is
neutral. The right side of the first identity expands to
`e_barL[i gamma^mu(D_mu phi) A N_L+phi A^2 N_L]`.
Consequently the original `-phi A^2 N_L/6` becomes the `+1/6` in (10).
For the sterile terms the integrated row coefficients are

\[
\begin{split}
 \bar x_i={}&\frac56[(W M^2)_i+K_0^2W_i]+cW_i
 +\frac13\sum_{\alpha a}F_{\alpha i}(K\bar S_{L,\alpha}^a)\phi^a
 -\frac16\sum_{\alpha a}F_{\alpha i}\bar L_\alpha^a S_\phi^a
 -\frac16K_0V_i-\frac{g_y}{12}W_i\sigma^{\mu\nu}B_{\mu\nu},\\
 \bar y_i={}&\frac56(K_0W)_j M_{ji}-\frac16V_jM_{ji}.
\end{split}                                                    \tag{14}
\]

Here `bar x=x^dagger gamma^0`, `bar y=y^dagger gamma^0`, and (12) is
exactly the adjoint of (14). In particular

\[
 K_0^2(\bar L\phi)=-\left[(D^2\bar L)\phi
        +2(D_\mu\bar L)D^\mu\phi+\bar L D^2\phi\right],
\]

and combining the two source derivatives in (14) gives
`-i F(D_mu bar S_L)gamma^mu phi/6`
`+i F bar S_L gamma^mu D_mu phi/6`. This records the derivative-adjoint
signs before any on-shell substitution.

Equation (5), after (13), is therefore precisely

\[
 P_F-P_{\rm red}=\sum_{\alpha a}e_{\bar L,\alpha}^a\ell_\alpha^a
       +t_a e_\phi^a-\frac{g_y}{6}j_\mu e_B^\mu
       +\bar x_i e_{Ri}+\bar y_i e_{Li}.                       \tag{15}
\]

There are no unlisted residuals in (15). To handle its conjugate sterile
channel, use the Grassmann scalar identity
`bar E w=overline(w^c) E` for a Majorana `E`. With `w=x+y`, chirality gives

\[
 \bar x e_R+\bar y e_L+\mathrm{h.c.}
 =\overline{(w+w^c)}\,(e_R+e_L).                              \tag{16}
\]

Equations (9), (15), and (16) prove (7), including the Majorana factor.
Using `-(w+w^c)/2` would cancel only half this residual packet. Treating
`N_L,N_R` as unrelated Dirac variables would not preserve the reference
Majorana constraint.

The literal gauge-fixed action additionally contains
`S_gf=-(1/2xi) integral(partial.B)^2`. Its variation is
`delta S_gf=-(1/xi) integral(partial.B)partial.delta B`, which is nonzero
for (6) in general. Equation (7) claims equivalence for the gauge-invariant
action, as does the input reduction. Applying the substitution to a
full gauge-fixed generating functional requires transforming its gauge
fixing and associated off-shell conventions as well.

## Quadratic terms on a neutral scalar background

At the unbroken origin, with all fields zero, each shift above starts at
field degree two or higher. The selected packet therefore induces no
quadratic normalization change there. This statement does not extend to
a nonzero scalar background.

Take `H=(0,u)`, `phi=(u,0)`, `u=v/sqrt(2)` real and constant, `B=0`,
with zero fermion background and the renormalizable sources (8). Let
`nu_L=L^1`, `c_v=2m_H^2/3+lambda u^2/3`. The linear fermion coefficients are

\[
\begin{split}
 \ell_\alpha^{1,(1)}&=u\sum_iF_{\alpha i}
                 (A^2/6+m_H^2/2)N_{Li},
 &\ell_\alpha^{2,(1)}&=0,\\
 x_i^{(1)}&=u\sum_\alpha\left[
     \frac56(FM^2)^*_{\alpha i}\nu_{L\alpha}
     +F_{\alpha i}^*(5A^2/6+c_v)\nu_{L\alpha}\right]
     +\frac{u^2}{6}\sum_{\alpha j}F_{\alpha i}^*Y_{\alpha j}A N_{Rj},\\
 y_i^{(1)}&=\frac{5u}{6}\sum_\alpha(FM)^*_{\alpha i}A\nu_{L\alpha}
     -\frac{u^2}{6}\sum_{\alpha j}(FM)^*_{\alpha i}Y_{\alpha j}N_{Rj}.
\end{split}                                                    \tag{17}
\]

These supply both active-sterile and sterile flavor transformations;
`A^2=-Box` produces momentum-dependent quadratic terms. The scalar and
gauge shifts have no linear part on this background, but they are
bilinear in fermions. If the scalar background is off shell, the scalar
tadpole `e_phi^(0)=-(m_H^2+2lambda u^2)phi` multiplied by its quadratic
shift contributes another fermion bilinear. A pure Jacobian transport of
the quadratic fermion kernel is valid only after accounting for this
term, or at a stationary background where `m_H^2+2lambda u^2=0`.
With `lambda=0`, an arbitrary nonzero `u` is not a stationary background
of the literal truncated scalar potential.

As a check on why normalization is needed, the original unreduced
selected packet has, on that constant background,

\[
 P_F^{(2)}=u\,F_{\alpha i}\bar\nu_{L\alpha}
            A(A^2+m_H^2)N_{Li},                              \tag{18}
\]

modulo integration by parts. Thus it changes the neutral kinetic kernel
already before the substitution. Equation (6) transports that change
into the representative (1), including its correlated mass and source
terms; it is not a declaration that every remaining kinetic matrix is
canonical. A further choice of fermion coordinates, all other retained
operators, and the complete parameter counterterms are needed for a
physical coefficient vector. The retained Majorana mass, source order,
and flavor order in this note do not authorize a beta function, pole
shift, or observable from this selected packet alone.

## Reproduction and scope of the checks

`verify_field_redefinition.py` is a standalone exact checker with no
repository imports. Its tests address the integrated derivative
adjoints, the explicit covectors, the Majorana scalar identity and half
factor, source signs, and the constant-background linearization. They
are action-algebra controls, not another integration of the source loops.
Run with `OMP_NUM_THREADS=1 python -I -B -O verify_field_redefinition.py`.
The optimized run with Python 3.12.14 and SymPy 1.14.0 passes 192 explicit
checks and rejects fourteen deliberate sign, source, chirality, factor,
flavor-slot, and background omissions. The saved result is
`FIELD_REDEFINITION_RESULTS.json`. Six arbitrary-jet trials use two active
and three sterile flavors with two weak components; the raw sterile
residual Fourier reconstruction supplies six additional trials. These
counts describe implementation controls within one derivation.
