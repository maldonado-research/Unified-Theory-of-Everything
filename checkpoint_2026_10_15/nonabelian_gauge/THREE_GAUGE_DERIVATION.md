# Three external gauge fields: explicit 33-diagram UV residue

15 October 2026. A one-loop calculation for the single derivative seed and
conjugate Yukawa sector fixed by the literal action in
`checkpoint_2026_10_15/operator_map/ACTION_SIGN.md`. The loop integration
below is new; the lower-order action and vertex conventions are inputs.

## Result and exact scope

Write `F=C_(7,seed) Y_nu*`, and let `X_i` be the generator including its
coupling for gauge mode `i`. The `X_i` can be arbitrary noncommuting
matrices. For the electroweak theory they are `g' y_L I` for hypercharge
and `g_2 T^a`, with `T^a=tau^a/2`, for weak isospin. Define `E_i=slash(epsilon_i)` and
`d_ij=epsilon_i.epsilon_j`, with metric `(+---)`.

The sum of the 33 selected one-loop 1PI matter-line diagrams is

\[
\begin{split}
R_3={}&-\frac14\sum_{i<j;\,k\ne i,j}
 X_k[X_i,X_j][E_i,E_j]E_k\\
&-\frac13\sum_{i<j;\,k\ne i,j}
 [X_k,[X_i,X_j]](d_{ki}E_j-d_{kj}E_i).
\end{split}
\]

Each sum contains three terms. The full pole, in the fixed source-chain
convention, is

\[
i\mathcal M^{\rm UV}_{3A}
=-\frac{iF}{16\pi^2\epsilon}\bar u_L R_3P_Lu_N.
\]

The local cancelling action therefore has Fourier kernel
`+F R_3/(16 pi^2 epsilon)`. The external electroweak source contraction
and flavor labels are suppressed in this equation. The result agrees
with the independently expanded cubic kernel of the chosen covariant
action. It is a restricted off-shell UV statement; no finite amplitude,
physical beta function, or observable has been calculated here.

All `X_i` commute in the pure hypercharge specialization, so its cubic
residue vanishes. For a BWW specialization, the nested commutator term
vanishes and the surviving term is

\[
R_{BWW}=-\frac14X_B[X_{W_1},X_{W_2}]
 [E_{W_1},E_{W_2}]E_B.
\]

The WWW and BWW residues can be nonzero. No external momentum or scalar
mass appears in `R_3`; this follows independently from its UV degree zero.

An independent charged-basis projection of the actual loop sum is also
checked. With `X_A=e diag(0,-1)` and `X_±=g_2 T_±/sqrt(2)`, its upper
neutral-doublet component gives

\[
(R_{A W^+W^-})_{11}=e g_2^2\left\{
\frac18[E_A,E_-]E_+
+\frac16\left[d_{A-}E_++d_{A+}E_--2d_{+-}E_A\right]\right\}.
\]

The checker evaluates the loop color words directly at `e=1`,
`g_2=sqrt(2)` and proves all 64 Lorentz coefficients of this homogeneous
coupling-order identity; no floating-point square roots are needed.
This nonzero mixed local contact is compatible with a vanishing pure
neutral-photon pole. It is not a finite photon amplitude or a scattering
prediction; the remaining amplitude ingredients have not been supplied.

## Topology count, ordering, and physical vertex signs

The neutral seed has no explicit gauge contact. A fermion kinetic vertex
emits one gauge field. The scalar kinetic action emits either one gauge
field or a seagull pair; it has no cubic scalar-gauge contact.

With three single vertices, the number of fermion-line vertices is
`n_F=0,1,2,3`. Each distribution gives exactly six diagrams:

| Distribution | Choice and ordering | Number |
| --- | --- | --- |
| FFF | three fermion orderings | 6 |
| FFS | choose scalar label, order the two fermion labels | 6 |
| FSS | choose fermion label, order the two scalar labels | 6 |
| SSS | three scalar orderings | 6 |
| scalar seagull + scalar single | choose seagull pair, order pair/single along scalar line | 6 |
| scalar seagull + fermion single | choose seagull pair | 3 |

The total is 33. Line order `abc` always means Yukawa to seed. The open
fermion spinor chain is consequently read in the reverse order,
`Q E_c Q E_b Q E_a Q`, with `Q=slash(q)` at leading large loop momentum.

The scalar has generator `X_i`, while the charge-conjugate fermion has
generator `-X_i^T`. If the fermion and scalar emissions, read Yukawa to
seed, are `F=(f_1,...,f_r)` and `S=(s_1,...,s_t)`, respectively, their
crossed-seed color contraction is

\[
S_{\rm line} C_{\rm line}^T
=(-1)^rX_{s_t}\cdots X_{s_1}X_{f_1}\cdots X_{f_r}.
\]

This is a direct index contraction of the crossed source, not a choice
made by fitting the action. Its two-field specializations are
`FFij=Xi Xj`, `SSij=Xj Xi`, and `FSij=-Xj Xi`, agreeing with the separately
audited two-weak-field source convention.

The zero-field source phase is `(i)(i)(-i)(i)=-1`. An ordinary insertion
adds an extra propagator and vertex with phase `i(-i)=+1`. A scalar
seagull is `+i {Xi,Xj} d_ij`; its extra scalar propagator gives relative
phase `i(i)=-1`. The anticommutator already includes the two functional
gauge-field assignments. There is no further factor of two or global
symmetry factor. Its Abelian limit is the familiar relative coefficient
`-2 d_ij` after the common generator product is extracted.

For pair `(i,j)` and remaining label `k`, the physical color factors of
the three seagull topologies are

\[
-X_k\{X_i,X_j\},\quad-\{X_i,X_j\}X_k,
\quad+\{X_i,X_j\}X_k,
\]

for scalar pair-then-single, scalar single-then-pair, and fermion single,
respectively. The final plus combines the seagull phase with the
charge-conjugate fermion sign.

External Yang-Mills tree attachments to a matter-loop gauge insertion
are 1PR and are excluded. Internal gauge bridges connecting two distinct
matter points require higher gauge-coupling order. The displayed list is
the primitive matter-line 1PI content at the stated order, and does not
stand for the complete amplitude with external-leg or gauge-tree
radiation.

## UV extraction with an infrared separator

External momenta are initially generic and off shell. The crossed seed is
`q^2+O(q)`, while the direct seed is a polynomial in external momenta and
contains no loop momentum. Ordinary three-gauge graphs have five
denominators. Their leading numerator is degree six, giving UV degree
zero; every external-momentum correction lowers the degree and is UV
finite. Graphs with a seagull and one single vertex have four
denominators and numerator degree four, again degree zero. Thus only the
crossed seed's `q^2` term can contribute. Direct-seed terms are
individually UV finite here. The seed's earlier `2:1` pairing weight is
therefore not independently retested by this Green function.

Introduce a common positive auxiliary squared mass `M_IR^2` in all
leading denominators. Replacing a shifted physical propagator by its
unshifted massive leading term changes these degree-zero integrands by
UV-convergent terms. This replacement defines an IR-separated UV
functional and never sets an original massless integral to zero as
scaleless. Its radial normalization follows from

\[
\int\frac{d^dq}{(2\pi)^d}
\frac{(q^2)^{n-2}}{(q^2-M_{\rm IR}^2+i0)^n}
=\frac{i}{(4\pi)^{2-\epsilon}}
 (M_{\rm IR}^2)^{-\epsilon}
\frac{\Gamma(n-\epsilon)\Gamma(\epsilon)}
 {\Gamma(2-\epsilon)\Gamma(n)}
=\frac{i}{16\pi^2\epsilon}+O(1).
\]

Here `d=4-2epsilon`. This formula applies to `n=3,4,5` used below.
Mass-dependent changes affect finite terms only. Angular reduction uses

\[
\langle q^\mu q^\nu\rangle/(q^2)=g^{\mu\nu}/4,
\qquad
\langle q^\mu q^\nu q^\rho q^\sigma\rangle/(q^2)^2
=\operatorname{sym}(gg)/24.
\]

Four-dimensional angular and Clifford identities suffice for these
simple UV poles; `O(epsilon)` corrections change finite terms. The
script also retains the seed `q^2` and applies rank-six moments for
ordinary graphs, or rank-four moments for seagull graphs. That gives
exactly the same residues as cancelling `q^2` first. With the auxiliary
mass, the difference introduced by cancellation contains an extra
`M_IR^2` and is UV finite, so it cannot alter the pole.

## Individual integrated primitive residues

Strip the radial `i/(16 pi^2 epsilon)`, all color factors, and the
seagull's relative vertex phase. Let `(a,b,c)` be any permutation of
`(1,2,3)`. The ordinary residues are

\[
\begin{aligned}
R_{FFF}^{abc}&=\tfrac13(d_{ab}E_c-2d_{ac}E_b+d_{bc}E_a),\\
R_{FFS}^{ab|c}&=\tfrac12E_bE_aE_c
-\tfrac13(d_{ab}E_c+d_{ac}E_b-2d_{bc}E_a),\\
R_{FSS}^{a|bc}&=\tfrac13(d_{ab}E_c+d_{ac}E_b-2d_{bc}E_a),\\
R_{SSS}^{abc}&=-\tfrac13(d_{ab}E_c+d_{ac}E_b+d_{bc}E_a).
\end{aligned}
\]

These follow by applying the angular moments to, respectively,

\[
\begin{aligned}
&\langle Q E_c Q E_b Q E_a Q\rangle,\\
&-2\langle Q E_b Q E_a Q(q\cdot\varepsilon_c)\rangle,\\
&4\langle Q E_a Q(q\cdot\varepsilon_b)(q\cdot\varepsilon_c)\rangle,\\
&-8\langle Q(q\cdot\varepsilon_a)(q\cdot\varepsilon_b)
 (q\cdot\varepsilon_c)\rangle,
\end{aligned}
\]

with the appropriate radial powers stripped. For example,
`Q E Q=2(q.epsilon)Q-q^2 E` reduces all four expressions to the two
displayed angular moments. Both scalar seagull orderings and the
fermion-single seagull graph have the same unweighted Lorentz residue,

\[
R_{SG}^{ij|k}=-\tfrac12 d_{ij}E_k.
\]

The nine weighted seagull diagrams therefore sum to
`+(1/2) sum_(i<j) d_ij X_k {Xi,Xj} E_k`. Substituting all four ordinary
residues with their ordered color factors and adding this sum gives the
boxed cubic result above. This uses the actually integrated diagrams;
the action is used as a comparison target afterward.

## Reproduction and verification

Run `python3 -I -S -B -O verify_three_gauge.py`.

The standard-library checker imports no repository modules. It reports
**5805 total controls** (5798 positive checks and seven mutation rejections), including 2508 comparisons of
individual integrated diagrams with the displayed analytic residues and
2508 comparisons of cancelled/uncancelled seed extraction. It exhausts
all `4^3=64` polarization basis coefficients, giving exact equality of
each formal noncommuting color-word coefficient with the independently
derived cubic action. These 64 coefficients prove the trilinear Lorentz
identity for arbitrary polarizations and arbitrary generator matrices;
12 generic rational polarization triples are additional checks. The
direct loop `A W^+ W^-` neutral projection and all six Bose permutations
are checked at every point.

Seven deliberate mutations are rejected: dropping seagulls, reversing
their phase, dropping one FFF ordering, assigning the same generators to
both internal lines, reversing scalar/fermion color transport, omitting
the curvature-connection term, and omitting the adjoint-connection term.
Detection counts are recorded in `THREE_GAUGE_RESULTS.json`; individual
special polarization triples can have zero residues, so no false claim
of detection at every point is made. Checks remain active under `-O`.

The independent action/Fourier/three-gauge Ward audit is a separate
calculation and does not reuse these loop integrals. The separately authored matrix-cubature integration audit agrees on all
2,112 individual diagram evaluations and all 64 full color/Clifford
tensors; see `INDEPENDENT_CUBATURE_AUDIT.md`. That audit and its reproducible
comparison harness are distinct from this producer calculation.

## Bounded external-gauge closure and remaining work

For a primitive matter-line graph with `n` external gauge fields,
superficial UV degree is `3-n`: every ordinary insertion lowers it by
one, and a seagull replaces two such insertions and lowers it by two.
Consequently every such graph with `n>=4` is individually UV finite.
Together with the already computed zero-, one-, and two-gauge residues,
the three-gauge result completes the potentially divergent
external-gauge tower for this restricted one-loop source sector.

This is not a full EFT operator-basis closure. Other external matter
sectors, interaction-source terms, further coupling orders, scheme
choices, field normalization, and physical running remain to be
calculated. The result supplies a stronger gauge-consistency check and
a reproducible nonabelian UV counterterm, without supplying a new
experimental prediction.
