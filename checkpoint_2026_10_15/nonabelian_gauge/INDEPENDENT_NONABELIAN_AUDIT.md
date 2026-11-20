# Independent two-gauge-leg weak and mixed-sector UV audit

15 October 2026. Scope: the literal seed/conjugate-neutrino-Yukawa one-loop
matter 1PI sector of candidate.5, commit `5c77bc5`, with action conventions
in `checkpoint_2026_10_15/operator_map/ACTION_SIGN.md`.

The independent two-B angular extractor was adapted in a scratch directory.
No producer or repository module is imported at runtime. The non-Abelian
weak-index contractions, the ordered-color decomposition, the curvature and
covariant-action expansion, and the non-Abelian Ward recursion were derived
independently before comparison with the parallel producer and action audit.
This is a restricted simple-UV-pole result, not a complete EFT basis,
physical anomalous dimension, finite amplitude, or observable.

## Result

Include gauge couplings in the two external generator matrices, `M1` and
`M2`. For two weak fields use `Mi=g T^Ii`; for a hypercharge field use
`Mi=g' yL 1`, with `yL=-1/2`. Set `p=r+h+k1+k2`, and write
`Ei=slash(ei)`, `Ki=slash(ki)`, `e12=e1.e2`, `C=[M1,M2]`. Then the seven
matter diagrams give

\[
\begin{aligned}
R_2={}&\frac12\{M_1,M_2\}e_{12}\not p
-\frac14\bigl(M_2M_1[\not k_1,E_1]E_2
              +M_1M_2[\not k_2,E_2]E_1\bigr)\\
&+[M_1,M_2]\left\{\frac14[E_1,E_2]\not r
+\frac13\left[((2k_1+k_2)\cdot e_2)E_1
-((k_1+2k_2)\cdot e_1)E_2
+e_{12}(\not k_2-\not k_1)\right]\right\}.
\end{aligned}
\]

In the fixed action, source-chain, and Fourier conventions,

\[
i\mathcal M_2^{\rm UV}=-\frac{iF}{16\pi^2\epsilon}
\bar u_L R_2 P_Lu_N\,\phi,
\qquad F=C_{(7,\mathrm{seed})}Y_\nu^*.
\]

Its local counterterm has the opposite amplitude. The formula restores the
external weak matrices explicitly and includes all polarizations and generic
off-shell momenta. It agrees with the *quadratic* gauge-field expansion of
`O_J + (i DbarL sigma gamma F)/4 + (barL gamma D F)/3`.

For mixed B/W, the matrices commute and the result becomes
`g'g yL T^I R_BB`. There is no additional factor two: the seven labelled
external-leg diagrams and the anticommutator seagull already perform the
required field derivatives. For unequal weak indices the commutator survives;
that information is absent from the Abelian calculation.

## Literal color contractions

Use the seed tensors `Kd_ab;cd=delta_ac delta_bd` and
`Kx_ab;cd=delta_ad delta_bc`, with the source tensor `delta_ef`. The scalar
line matrices are `Mi`, while the charge-conjugate fermion matrices are
`Ci=-Mi^T`. For any ordered fermion chain `Cchain` and scalar chain `Schain`,
the crossed contraction is, explicitly,

\[
\sum_{bdef}\delta_{ad}\delta_{bc}
 (C_{\rm chain})_{be}(S_{\rm chain})_{df}\delta_{ef}
=(S_{\rm chain}C_{\rm chain}^{T})_{ac}.
\]

If label 1 is nearest Yukawa and label 2 nearest seed on a single line,
the diagram factors are therefore:

| Diagram | Ordered matrix, including mixed charge or seagull phase |
|---|---|
| FF12 | `M1 M2` |
| FF21 | `M2 M1` |
| SS12 | `M2 M1` |
| SS21 | `M1 M2` |
| FS12 | `-M2 M1` |
| FS21 | `-M1 M2` |
| SG | `-{M1,M2}` |

The unweighted SG residue includes `e12`. Its kinetic vertex is
`+i{M1,M2} g_mu_nu`, so its extra propagator supplies the relative minus.
For Abelian identical charges the anticommutator is `2 M1 M2`, reproducing
the existing relative `-2` seagull weight. Ordinary insertions have phase
`i(-i)=+1`; the zero-gauge source chain remains `-1`.

The transpose in the fermion representation reverses its matrix order
when contracted with the crossed tensor. Dropping the transpose fails for
`T2`, whose entries are imaginary. This audit sums every internal index
over exact Gaussian rationals for all nine WW, three BW, three WB, and one
BB generator pair, before forming any loop sum.

The direct seed pairing is proportional to `t^2`, with `t=p-k1-k2`; its
matrix traces need not vanish for two weak generators. Nevertheless all its
seven two-gauge diagrams are UV finite individually: the box numerators
have degree at most three rather than the minimum four for the logarithm,
and the triangle seagull has degree one rather than the minimum two.
Consequently this calculation does not revalidate the seed's earlier 2:1
pairing weight.

## Actual matter-loop integration and ordered coefficients

The scalar and fermion momenta and the seven routed numerators are those
listed in `TWO_B_DERIVATION.md`; the weak matrices do not change these
propagator denominators. No new matter topology is induced by replacing the
external Abelian labels with matrices. This is a statement about literal
kinetic vertices, not a replacement rule for an already-summed answer.

Write the independently angular-integrated, unweighted residues as
`FF12, FF21, SS12, SS21, FS12, FS21, SG`. Keeping the generator order through
integration gives

\[
R_2=M_1M_2 U_{12}+M_2M_1 U_{21},
\quad U_{12}=FF12+SS21-FS21-SG,
\quad U_{21}=FF21+SS12-FS12-SG.
\]

Each diagram is obtained by retaining its top two homogeneous numerator
powers in the large loop momentum q. For four denominators, the UV radial
logarithm is `A[N4]+A[2 q.sum(shifts) N5]`; for the seagull triangle it is
`A[N2]+A[2 q.sum(shifts) N3]`. The normalized d=4 angular moments use exact
Wick pairings, with denominators `4`, `4*6`, and `4*6*8` through rank six.
An infrared separator is implicit. This extracts a UV coefficient without
identifying an unseparated scaleless integral with a nonzero integral.

The q-independent parts of the crossed seed, and all scalar-mass terms in
the propagator expansion, lie below this pole's required degree. Thus the
residue is mass independent, jointly linear in `p,h,k1,k2`, and bilinear
in `e1,e2`. Four-dimensional Clifford reduction suffices for the simple
pole; this does not fix any finite or evanescent terms.

The resulting independent ordered coefficients agree with the target
obtained from the action:

\[
U_{12}=\frac12e_{12}\not p
-\frac14[\not k_2,E_2]E_1+S,
\quad
U_{21}=\frac12e_{12}\not p
-\frac14[\not k_1,E_1]E_2-S,
\]

where

\[
S=\frac14[E_1,E_2]\not r
+\frac13\left[((2k_1+k_2)\cdot e_2)E_1
-((k_1+2k_2)\cdot e_1)E_2
+e_{12}(\not k_2-\not k_1)\right].
\]

In particular `(U12+U21)` is the known two-B polynomial, while `(U12-U21)`
contains independent commutator information. The checker integrates all
seven actual routed numerators, compares each to its analytic angular
reduction, and only then compares these free ordered-color coefficients
with the covariant action.

## Curvature sign and quadratic covariant expansion

With matrix connection `A_mu` and `Dplus=partial+iA`,

\[
[D_\mu,D_\nu]=i\mathcal F_{\mu\nu},\qquad
\mathcal F_{\mu\nu}=\partial_\mu A_\nu-\partial_\nu A_\mu
+i[A_\mu,A_\nu].
\]

For `A=g T^I W^I`, the component nonlinear term is
`-g epsilon^IJK W_mu^J W_nu^K`. Importing the opposite curvature sign while
keeping this covariant derivative would violate the commutator identity.
The adjoint derivative is `D_adjoint F=partial F+i[A,F]`.

The six-term `O_J` gives `e12 slash(p) {M1,M2}/2` at two gauge fields.
The row connection of `i DbarL sigma gamma F /4` gives the two ordered
`-[k_i,E_i] E_j/4` terms. Its nonlinear curvature gives
`+[M1,M2][E1,E2] slash(r)/4`.

The derivative of the nonlinear curvature in `D_nu F^{nu mu}/3` gives

\[
\frac{[M_1,M_2]}3
 [((k_1+k_2)\cdot e_2)E_1-((k_1+k_2)\cdot e_1)E_2],
\]

and the connection `i[A_nu,F_linear^{nu mu}]/3` gives

\[
\frac{[M_1,M_2]}3
 [(k_1\cdot e_2)E_1-(k_2\cdot e_1)E_2
 +e_{12}(\not k_2-\not k_1)].
\]

Adding these independently obtained terms produces the displayed action
kernel. Each nonlinear structure can be removed separately in the checker;
each removal is detected on every generic noncommuting gauge pair.

## Non-Abelian Ward recursion and 1PI scope

Let `R1(p,h;k,e)` be the audited one-gauge polynomial with its generator
and coupling stripped, and `p` kept fixed. Contracting the first leg gives

\[
\begin{aligned}
R_2(e_1\to k_1)={}&M_2M_1 R_1(p,h+k_1;k_2,e_2)\\
&-M_1M_2 R_1(p,h;k_2,e_2)
+[M_1,M_2]R_1(p,h;k_1+k_2,e_2).
\end{aligned}
\]

The last term is essential: it comes from the transformation of the
remaining gauge field in the gauge-invariant 1PI effective action. The
first term has row momentum `r`, the second `r+k1`, and the combined-gauge
kernel again has row momentum `r`. The swapped identity holds for leg 2.
The three scalar-mass terms cancel in this combination. For commuting
B/W generators the last term vanishes, and the Abelian descent is recovered.

A tree triple-gauge vertex attached through a gauge propagator to the
one-W loop is *one-particle reducible*. Cutting its internal gauge line
disconnects the tree from the matter loop. Its nonlocal propagator does not
belong in the local two-W 1PI counterterm extracted here. The lower-one-W
commutator in the 1PI Ward identity is required without adding that graph.
Such a tree graph belongs in a connected scattering calculation and its
subdivergence is addressed by the lower-point counterterm there. No
claim about the full connected amplitude is made here.

## Standalone reproduction and limits

Run `python3 -I -S -B -O verify_two_gauge_angular.py`. Failures are explicit
exceptions and remain active under optimization. `INDEPENDENT_TWO_GAUGE_CHECKS.json`
records the completed checks and exact W1/W2 and B/W2 witnesses.

The completed default run passes **4,285 exact checks**:

- 144 literal internal-index color-contraction checks.
- 2,080 checks over 16 generic rational momentum/polarization samples,
  including both gauge-leg Ward recursions at three signed scalar masses.
- 2,048 exhaustive polynomial checks: seven individual diagram reductions
  plus both free ordered-color coefficients at each of 256 basis points.
- 13 mutation-detection assertions. Every one of the 96 generic
  noncommuting WW cases rejects omitted nonlinear curvature, reversed
  curvature sign, omitted adjoint connection, omitted lower-Ward kernel,
  and reversed FF color ordering. Wrong fermion transpose is detected in
  112 gauge/sample cases, omitted seagull in 150, and reversed mixed
  charge sign in all 256 gauge/sample cases.

The complete 256-point basis consists of one of the 16 momentum coordinates
of `(p,h,k1,k2)` and each of the 16 polarization-basis pairs. Because the
local residue has the degree established above, checking both free
ordered-color coefficients at every point proves the polynomial identity
for arbitrary external generator matrices at degree two. This is not a
complete enumeration of EFT operators.

This establishes an additional multi-W and mixed-B/W test of the chosen
covariant action within the explicitly selected matter loop. It does not
supply a complete counterterm set, an EOM-reduced operator basis, physical
running, finite photon matching, a decay rate, or a new measured prediction.
