# Literal-action sign of the one-hypercharge-leg pole

15 October 2026. An independent action/Fourier audit of the specified seed
sector. The previously evaluated triangle residues are inputs to this note;
they are not a complete anomalous dimension or observable.

## Result

Fix the conventions below, set `F=C_(7,seed) Y_nu*`, and keep the public ordered
residue `R=R_H-R_L=G_J+D`. The full poles and their cancelling local action are

\[
i\mathcal M_0^{\rm UV}=-\frac{i}{16\pi^2\epsilon}
 F\,\bar u_L\not J P_Lu_N,
\qquad
i\mathcal M_B^{\rm UV}=-\frac{i g'y_L}{16\pi^2\epsilon}
 F\,\bar u_L R^\nu P_Lu_N\,\varepsilon_\nu^*,
\]

\[
\boxed{\mathcal L_{\rm ct}=\frac{F_{\alpha i}}{16\pi^2\epsilon}
 \left[(O_J)_{\alpha i}+g'y_L\left(\frac14Q_\sigma+
 \frac13Q_{\partial B}\right)_{\alpha i}\right]+\mathrm{h.c.}}
\]

The external electroweak epsilon tensor in the amplitude is understood.
`y_L=-1/2` throughout. This sign is obtained from actual kinetic vertices,
charge conjugation and propagator phases, before imposing a Ward identity.

There is a necessary correction to the earlier charge language. For the
displayed `bar L ... N_L` channel, the open internal fermion is `L^c` and the
internal scalar is `tilde H`. Their charges are **`-y_L` and `y_L`**, respectively.
Consequently their physical charge-weighted sum is

\[
\operatorname{Res}[(-y_L)I_L+y_L I_H]=y_L(R_H-R_L),
\]

The earlier algebraic combination `y_L I_L+y_H I_H=y_H(I_H-I_L)`, with
`y_H=-y_L`, used the charges of the unconjugated fields. It is not the
physical weighting of this source open chain when `y_L` denotes the SM lepton
charge. The polynomial calculation of `R`, `G_J` and `D` does not change.
The previous action-phase caveat was concealing this missing charge-conjugate
identification; an arbitrary new phase is neither needed nor available once
the action and the zero-gauge calibration have been fixed.

## Complete conventions for this graph sector

Use metric `(+---)`, `d=4-2epsilon`, `S=T exp(i integral L_int)`, and Feynman
propagators with `+i0`. Define

\[
D_\mu\psi=(\partial_\mu+i g'y_\psi B_\mu)\psi,
\quad D_\mu\bar\psi=\partial_\mu\bar\psi-i g'y_\psi\bar\psi B_\mu,
\quad B_{\mu\nu}=\partial_\mu B_\nu-\partial_\nu B_\mu,
\]

\[
\gamma_5=i\gamma^0\gamma^1\gamma^2\gamma^3,
\quad P_{L,R}=(1\mp\gamma_5)/2,
\quad\sigma^{\mu\nu}=\frac{i}{2}[\gamma^\mu,\gamma^\nu],
\quad C=i\gamma^2\gamma^0.
\]

Thus `C gamma^{mu T} C^{-1}=-gamma^mu`, `C^T=-C`, and
`psi^c=C bar psi^T`, `bar psi^c=-psi^T C^{-1}`. The gauge transformation is
`psi -> exp(-i g'y_psi alpha) psi`, `B -> B+partial alpha`. No convention for
an anomalous closed chiral trace is needed in this open-chain calculation.

Let `L` be a left-chiral doublet, `H` a complex scalar doublet,
`tilde H=epsilon H*`, `epsilon^{12}=+1`, and `N=N^c` a Majorana singlet with
real positive diagonal mass `M`. The relevant literal reference action is

\[
\begin{aligned}
\mathcal L_{\rm ref}={}&-\frac14B_{\mu\nu}B^{\mu\nu}
-\frac1{2\xi}(\partial_\mu B^\mu)^2
 +\bar L i\!\not D L+(D_\mu H)^\dagger D^\mu H-m_H^2H^\dagger H\\
&+\frac12\bar N(i\!\not\partial-M)N
-\left(\bar L_\alpha(Y_\nu)_{\alpha i}\widetilde H N_{Ri}
 +\mathrm{h.c.}\right)
 +\left[\frac12(C_{7,\rm seed})_{\alpha\beta}
 \bar X_\alpha(-\Box)X_\beta^c+\mathrm{h.c.}\right],\\
X_\alpha={}&\widetilde H^\dagger L_\alpha,
\qquad X_\beta^c=\widetilde H^T L_\beta^c,
\qquad C_{7,\rm seed}=C_{7,\rm seed}^T.
\end{aligned}
\]

Here `y_H=+1/2`, `y_L=y_tildeH=-1/2`, `y_N=0`; `m_H^2` is the signed
parameter in `q^2-m_H^2+i0`. Other interactions may be added as spectators,
but none enters the selected coefficient linear in `C_(7,seed)`, `Y_nu*`
and one `g'`. This literal action specifies every vertex used here. No
sterile propagator occurs. The derivative acts on the entire neutral `X^c`,
so the seed has no gauge-field contact vertex.

For the displayed external channel the incoming `N` wave is `exp(-ip.x)`;
the outgoing `bar L`, `tilde H` and `B` waves are respectively `exp(+ir.x)`,
`exp(+ih.x)` and `exp(+ik.x)`. Thus `p=r+h+k`. Internal momentum `q` flows
from the Yukawa vertex to the seed in the charge-conjugate fermion line.
All external momenta are off shell until after the UV residue is defined.
The ordinary loop measure has UV normalization `i/(16 pi^2 epsilon)`.

The seed uses literal `-Box`; its action coefficient equals the public
positive-momentum vertex coefficient. With instead
`+(1/2) barX W1 Box X^c`, the translation is `C_(7,seed)=-W1`. This changes
`F` to `-W1 Y_nu*` in **both** poles and counterterms. It does not change
their relative sign.

The nondifferentiated normalization is anchored directly to the same
primary source: [Zhang, Eqs. (2)–(3)](https://arxiv.org/html/2405.18017v3#S2.E2)
uses `+(1/2) C5 O5+h.c.` and
`O5=bar L tildeH tildeH^T L^c=barX X^c`. Our literal seed inserts `-Box`
in precisely that monomial and orientation. Its action normalization is
therefore specified independently of the later Ward test.

## Seed normalization and the source chain

Differentiating the literal symmetric bilinear gives

\[
V_7=i(C_{7,\rm seed})^{\alpha\gamma}
\left[\epsilon^{ac}\epsilon^{bd}(r_a+h_c)^2
+\epsilon^{ad}\epsilon^{bc}(r_a+h_d)^2\right]P_R.
\]

The factor `1/2` is removed by exchanging the two fermion slots. The
Grassmann minus is cancelled by `C^T=-C`. The two identical scalar
derivatives generate the direct and crossed epsilon pairings; this is not
an extra overall factor of two. Momentum conservation equates the squared
momenta of complementary composites, so the two possible placements
inherited from the fermion-slot exchange agree. Contracting the Yukawa
epsilon yields weights two for the direct pairing and one for the crossed
pairing. The checker differentiates all four slot assignments explicitly.

For a charge-conjugated scalar bilinear,

\[
\bar N P_L L=\overline{L^c}P_LN^c.
\]

There is no additional minus: writing the right side with
`bar L^c=-L^T C^{-1}` produces a minus that is cancelled when the Grassmann
fields are returned to the original order. Since `N^c=N`, the Hermitian
conjugate Yukawa term can be written, without adding a second copy, as

\[
-\bar N P_L Y_\nu^\dagger\widetilde H^\dagger L
=-\overline{L^c}_\gamma(Y_\nu^*)_{\gamma i}
 \widetilde H^*N_{Li}.
\]

Its vertex in the source orientation is `-i epsilon^{bd}Y_nu* P_L`.
This gives the source ordering

`bar u_L [i C seed P_R] [i slash q/q^2] [-i Y_nu* P_L] u_N [i/D_p]`.

The product of its four factors of `i` is

\[
(i)(i)(-i)(i)=-1.
\]

Therefore the unintegrated zero-gauge graph is `-F T(p,h)`, precisely the
public calibration. At dimension five, replacing the seed polynomial by
one gives `-3 F integral slash q/[q^2 D_p]` and hence
`-i(3/2)F slash p/(16 pi^2 epsilon)`. This reproduces the explicit source
chain and pole in [Zhang, Appendix A, Eq. (35)](https://arxiv.org/html/2405.18017v3#A1.E35).
The derivative calculation itself is not supplied by that source.

## The actual gauge charges and ordered vertices

Write `chi=L^c` and `phi=tilde H`. Charge conjugation of the kinetic action,
using integration by parts for its derivative term, gives

\[
\bar L i\!\not D_{y_L}L
=\bar\chi i\!\not D_{-y_L}\chi
\quad\hbox{up to a total derivative}.
\]

Consequently the internal `chi` line has charge `-y_L`; `phi` has charge
`y_L`. Their propagators and gauge vertices follow from the above action:

| Line | Propagator | Emitted-B vertex |
|---|---|---|
| `chi=L^c`, charge `-y_L` | `i slash q/(q^2+i0)` with right-chiral restriction | `-i g'(-y_L) gamma^nu P_R` |
| `phi=tilde H`, charge `y_L` | `i/(s^2-m_H^2+i0)` | `-i g'y_L(s+s')^nu`, `s-s'=k` |

The scalar sign follows directly from
`i g'y_L[(partial phi^dagger)phi-phi^dagger(partial phi)]B` in its kinetic
Lagrangian. Incoming `phi` gives `-is`, outgoing `phi^dagger` gives `+is'`,
so the action polynomial is `-g'y_L(s+s')`, and the vertex is its product
with `i`. With `s=p-q`, `s'=p-q-k`, this gives the declared numerator
`(2p-2q-k)^nu` with the displayed charge and sign.

There is no concealed extra reversal of `q`. The propagator identity is

\[
C[S_L(-q)]^TC^{-1}=S_R(q).
\]

More decisively, conjugating the entire ordinary charged line gives

\[
C\{S_L(-q)V_{y_L}^\nu S_L(-q+k)\}^TC^{-1}
=S_R(q-k)V_{-y_L}^\nu S_R(q).
\]

Transposition reverses the product, while
`C gamma^{nu T} C^{-1}=-gamma^nu` reverses the vector charge. The resulting
open chain is exactly `(slash q-slash k) gamma^nu slash q`, with `q-k`
adjacent to the seed and `q` adjacent to the Yukawa vertex. Applying both
this conjugation and an additional charge or momentum sign would count
the same transformation twice.

The extra gauge insertion replaces one propagator by two and adds the
gauge vertex. Its extra scalar phase is `i(-i)=+1` for either line. Thus
the two unintegrated graphs are

\[
i\mathcal M_L=-g'F(-y_L)I_L,\qquad
i\mathcal M_H=-g'F y_L I_H,
\]

and their sum is `-g'y_L F(I_H-I_L)`. Multiplying by the pole factor `i`
proves the full one-B sign stated at the start. No Ward normalization was
used in this step, and no closed-fermion-loop minus occurs.

## Counterterm action and its one-B kernels

Retain the six derivative placements defining the public `O_J`:

\[
\begin{aligned}
O_J={}&\frac56\bar L\widetilde H A(-\partial^2)N_L
+\frac13(-D^2\bar L)\widetilde H A N_L
+\frac16\bar L(-D^2\widetilde H)A N_L\\
&-\frac16(\mathcal K\bar L)\widetilde H(-\partial^2)N_L
+\frac12m_H^2\bar L\widetilde H A N_L
+\frac12m_H^2(\mathcal K\bar L)\widetilde H N_L,\\
A={}&i\!\not\partial,\qquad
\mathcal K\bar L=-i(D_\mu\bar L)\gamma^\mu.
\end{aligned}
\]

Expanding its row and scalar derivatives, without assigning a Ward
normalization afterward, gives its one-B action kernel `g'y_L G_J^nu`,

\[
G_J^\nu=\left[-\frac{(2r+k)^\nu}{3}
+\frac{(2h+k)^\nu}{6}\right]\not p
+\left(\frac{p^2}{6}-\frac{m_H^2}{2}\right)\gamma^\nu.
\]

In particular `-D^2 barL` gives `-g'y_L(2r+k)^nu`, `-D^2 tildeH` gives
`+g'y_L(2h+k)^nu`, and `K barL` gives `-g'y_L gamma^nu`. These opposite
row/scalar signs are fixed by their covariant derivatives.

Define the gauge-invariant local monomials, with the shown gamma ordering,

\[
Q_\sigma=i(D_\rho\bar L)\sigma^{\mu\nu}\gamma^\rho
 \widetilde H N_L B_{\mu\nu},\qquad
Q_{\partial B}=\bar L\gamma_\mu\widetilde H N_L
 \partial_\nu B^{\nu\mu}.
\]

Their stripped one-B action kernels are respectively

\[
[\not k,\gamma^\nu]\not r,
\qquad k^\nu\not k-k^2\gamma^\nu.
\]

For example `sigma^{mu nu} B_{mu nu}` on the outgoing gauge mode gives
`-[slash k,slash epsilon]`; combining the `+ir_rho` row derivative and
the explicit `i` in `Q_sigma` fixes its positive displayed kernel. Hence
`g'y_L(Q_sigma/4+Q_partialB/3)` has kernel `g'y_L D`, with

\[
D^\nu=\frac14[\not k,\gamma^\nu]\not r
+\frac13(k^\nu\not k-k^2\gamma^\nu).
\]

The mixed-field action coefficient has no additional `1/2`. Its vertex
is `+i` times its action kernel and cancels the negative loop poles at zero
and one gauge field simultaneously. Flavor order remains
`(C_(7,seed)Y_nu*)_{alpha i}`; matrices are not commuted. The Hermitian
conjugate supplies the conjugate channel and must not be counted again
as a duplicate of the displayed channel.

## Two-gauge-field content of this chosen lift

The covariant derivative in `Q_sigma` includes a two-B contact. Together
with the quadratic connections in `O_L` and `O_H`, the action displayed
above predicts the following stripped two-B kernel, before multiplying
by `F/(16 pi^2 epsilon)` and the vertex factor `i`:

\[
(g'y_L)^2\left[(\varepsilon_1\!\cdot\varepsilon_2)\not p
-\frac14\left([\not k_1,\not\varepsilon_1]\not\varepsilon_2
+[\not k_2,\not\varepsilon_2]\not\varepsilon_1\right)\right]P_L.
\]

Here `p=r+h+k_1+k_2`. Both assignments of the two identical gauge fields
are included; the expression is Bose symmetric. This is the Fourier
content of the **chosen covariant action**, not a reported evaluation of
the two-B loop diagrams. That extra Green function remains a useful
independent check of its embedding in a larger counterterm calculation.

## Reproduction, negative controls, and limits

Run `python3 -I -S -B -O verify_action_sign.py` in this directory. It uses
only the standard library and exact Gaussian rational matrices, imports
no repository code, and reports 891 explicit checks over 32 generic
off-shell rational momentum choices. The checks include the literal
seed's four slot assignments, kinetic-action scalar vertex, the full
charge-conjugated ordered fermion chain, zero- and one-B action kernels,
the full phase/counterterm cancellation, and the two-B content of the
chosen lift. Wrong internal charges, an opposite counterterm sign,
omitting `D`, and reversing the fermion chain must all fail. The reversed
chain still passes the longitudinal contraction, showing why that
contraction cannot replace the charge-conjugation calculation.

The two already audited residue formulas `R_L,R_H` are explicit inputs to
this checker; these 891 checks must not be described as a third integration
of the two triangles. It tests the newly resolved action map and its
exact consequences. The charge-conjugation and Fourier identities are
stronger evidence for the relative sign than comparing two separately
normalized Ward polynomials.

The selected pole and its local cancelling action are now fixed in this
literal convention. No EOM reduction, subtraction-scheme finite term,
complete Green basis, two-B loop match, physical running coefficient, or
observable follows from this sign audit. In particular the singlet seed
does not authorize replacing the sourced light- and gauge-field EOM by
free equations when reducing the new field-strength operators.
