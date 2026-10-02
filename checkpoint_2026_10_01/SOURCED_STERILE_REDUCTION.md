# Sourced sterile equations and a Schur-invariance control

1 October 2026. Independent AI-assisted algebra audit of the public routed
dimension-seven seed checkpoint. This note supplies a bounded operator identity
and an exact field-redefinition control. It does not calculate a new loop pole,
an anomalous dimension, a threshold correction, or a physical cancellation of
the full dimension-seven sector.

## Result and scope

For a declared two-component action, a cubic derivative on a surviving sterile
field has three sourced descendants in addition to the familiar mass-cubed
Yukawa term. The mass-cubed Yukawa term and one Weinberg contact have opposite
variations in the leading Schur combination. This cancellation holds for an
arbitrary rectangular coefficient matrix. Two derivative contacts remain in
this reduction. Their presence prevents promotion of this algebraic
cancellation to a statement about the complete physical seed amplitude.

The public routed checkpoint gives the stripped polynomial

\[
\not J_{\rm pole}=\left(\frac56p^2+\frac13r^2+\frac16h^2
+\frac12m_H^2\right)\not p
+\left(\frac12m_H^2-\frac16p^2\right)\not r.
\]

This note analyzes the sterile cubic derivative that can represent the
\(p^2\not p\) component, in an explicitly defined chiral convention. It does
not assert that this component alone is the complete gauge-covariant pole.
The four other momentum/mass dependences, the Higgs-mass multiple of the
sterile first derivative, and the completion with external gauge fields must
be included in a complete reduction. In particular, the number \(5/6\) below
is never treated as an inferred physical beta-function coefficient.

## Conventions and sourced equations

Use metric \((+---)\), Fourier convention \(e^{-ipx}\),
\(\sigma^\mu=(1,\boldsymbol\sigma)\), and
\(\bar\sigma^\mu=(1,-\boldsymbol\sigma)\). Let

\[
D=i\sigma^\mu\partial_\mu,\qquad
\bar D=i\bar\sigma^\mu\partial_\mu,\qquad
D\bar D=\bar D D=-\Box.
\]

Let \(n_a\) be left-handed sterile Weyl fields, and
\(\mathcal J_\alpha=\epsilon_{ij}L^i_\alpha H^j\) the left-handed gauge-singlet
composite. The light flavor index is \(\alpha\), and the sterile index is \(a\).
Spinor contractions are implicit. The transpose in a holomorphic bilinear is
a flavor transpose, not a Hermitian conjugate. For anticommuting Weyl fields
the contracted bilinear is symmetric under exchange of its two complete field
slots.

Declare

\[
\mathcal L\supset i n^\dagger\bar\sigma^\mu\partial_\mu n
-\frac12(n^T M n+\mathrm{h.c.})
-(\mathcal J^T Y n+\mathrm{h.c.})
+\frac12(\mathcal J^T C_5\mathcal J+\mathrm{h.c.}),
\]

where \(M=M^T=M^*>0\) is diagonal in a Takagi basis, \(Y\) is light-by-sterile,
and \(C_5=C_5^T\). Define the independent equation residuals

\[
E=D n^\dagger-Mn-Y^T\mathcal J,\qquad
\bar E=\bar D n-Mn^\dagger-Y^\dagger\mathcal J^\dagger.
\]

They vanish on the leading sourced equations of motion. This is a leading-EOM
reduction at first order in the higher-dimensional insertion; additional
interactions involving \(n\) would add sources. Standard Model gauge and Higgs
interactions do not change this sterile equation, but they enter the subsequent
light-field reductions. \(D\) on \(n\) or on the whole composite \(\mathcal J\)
is an ordinary derivative because both are gauge singlets.

## An exact off-shell identity

Applying \(D\) to the conjugate equation gives

\[
-\Box n=M^2n+M Y^T\mathcal J+Y^\dagger D\mathcal J^\dagger
+M E+D\bar E.
\]

Then apply \(-\Box\) to the definition of \(E\), substitute this last
identity, and use that the constant flavor matrices commute with derivatives:

\[
\boxed{D(-\Box)n^\dagger=M^3n+M^2Y^T\mathcal J
+M Y^\dagger D\mathcal J^\dagger+Y^T(-\Box)\mathcal J
+(-\Box+M^2)E+M D\bar E.}
\]

This is an exact differential identity for arbitrary off-shell fields when
\(E,\bar E\) have the definitions above. Discarding the last two terms is an
EOM reduction, which is legitimate at the specified order only with the
induced field redefinition and its correlated interactions. Setting the
sources to zero is a different operation. These are four-dimensional classical
operator identities. Translating a dimensionally regulated pole also requires
declared evanescent-operator conventions wherever those become relevant.

For an arbitrary light-by-sterile matrix \(F\), with mass dimension \(-3\),
define the mixed local operator

\[
\mathcal O_p=\mathcal J^T F D(-\Box)n^\dagger.
\]

Modulo the displayed leading EOM terms,

\[
\boxed{\mathcal O_p\simeq
\mathcal J^T F M^3n
+\mathcal J^T F M^2Y^T\mathcal J
+\mathcal J^T F M Y^\dagger D\mathcal J^\dagger
+\mathcal J^T F Y^T(-\Box)\mathcal J.}
\]

All four terms have operator dimension seven before multiplication by \(F\):
\([n]=3/2\), \([\mathcal J]=5/2\), \([D]=[M]=1\), and \([Y]=0\).
The third term is a derivative contact between a composite and its conjugate;
the fourth is a derivative contact between two unbarred composites. Both
must be retained until their own operator reduction and matching are done.
The equation does not state that either contact is separately an observable.

## The correlated dimension-five cancellation

Add \(t\mathcal O_p+\mathrm{h.c.}\), where \(t\) is a bookkeeping scalar.
In the action conventions above, the first two reduced terms give

\[
\delta Y=-t F M^3,\qquad
\delta C_5=t(FM^2Y^T+YM^2F^T).
\]

The two terms in \(\delta C_5\) arise from the symmetric Weyl bilinear and the
declared \(1/2\) normalization. No assumption that \(F\) is symmetric is used;
it need not even be square. At this algebraic order \(\delta M=0\). Thus

\[
\delta(YM^{-1}Y^T)
=-t(FM^2Y^T+YM^2F^T),\qquad
\boxed{\delta(C_5+YM^{-1}Y^T)=0.}
\]

Choosing \(F=\kappa W_1Y^*\) recovers the structure of the historical
conditional cancellation for any \(\kappa\). This does not derive that choice
of \(F\), its sign, a pole-to-running conversion, or the coefficient
\(\kappa\). Dropping the Weinberg contact makes this variation nonzero.
Using a Hermitian conjugate in place of the Majorana transpose generally also
fails. The last two derivative contacts in the boxed reduction remain to be
carried through the complete analysis.

## An exact finite field-redefinition control

An independent way to understand the cancellation is the finite local shift
\(n=n'+A\mathcal J\), where \(A\) has sterile-by-light shape and dimension
\(-1\). The holomorphic mass/contact part of the action changes exactly to

\[
Y'=Y+A^TM,\qquad
C_5'=C_5-YA-A^TY^T-A^TMA,\qquad M'=M.
\]

Direct multiplication establishes

\[
\boxed{C_5'+Y'M^{-1}Y'^T=C_5+YM^{-1}Y^T.}
\]

The quadratic \(-A^TMA\) term is necessary for finite shifts. The full action
also acquires mixed kinetic and derivative contact terms, including the
kinetic contribution quadratic in \(A\). Omitting them would not implement
the full field redefinition. At first order, taking
\(A^T=-tFM^2\) reproduces the two algebraic shifts above; this observation
alone does not reproduce all of \(\mathcal O_p\).

More generally, for a quadratic block kernel
\(\begin{pmatrix}C&B\\B^T&K\end{pmatrix}\), elimination of the second block gives
\(C-BK^{-1}B^T\). A triangular field shift changes
\(B\) to \(B+A^TK\) and \(C\) to
\(C+BA+A^TB^T+A^TKA\), leaving that Schur complement unchanged whenever the
inverse exists. The same proof applies to differential kernels with the
transpose interpreted as the integration-by-parts transpose; in momentum
space this includes momentum reversal. Applying it to the full Majorana
kinetic kernel requires retaining both chiral blocks. A zero of this
basis-change diagnostic does not establish a zero physical loop response.

## What the gauge-free polynomial does not determine

Replacing light derivatives by covariant derivatives can provide a local
gauge-covariant representative for the six displayed momentum/mass terms.
That construction is not unique from the three-point polynomial alone.
For example, the gauge-invariant operator

\[
\mathcal Q_B=\mathcal J^T G\sigma^{\mu\nu}D n^\dagger B_{\mu\nu}
+\mathrm{h.c.},\qquad [G]=-3,
\]

has field/operator dimension seven and vanishes identically when the external
hypercharge gauge field is set to zero. Here
\(\sigma^{\mu\nu}=\tfrac{i}{4}(\sigma^\mu\bar\sigma^\nu-
\sigma^\nu\bar\sigma^\mu)\). An arbitrary coefficient of this term is invisible
to the gauge-free three-point matching. Its single-gauge-boson vertex is
transverse because it contains \(B_{\mu\nu}\), so the longitudinal Ward identity
does not determine that coefficient either. Applying the sourced sterile EOM
converts it to a mass-dependent transition-dipole term and a contact with two
composites; it does not in general vanish. This example demonstrates a kernel
of the restriction to zero external gauge fields; it is not a claim to have
classified the full independent operator basis or computed \(G\).

## Chiral and normalization limits of the comparison

The explicit identity uses the undotted channel
\(\mathcal J^T D n^\dagger\). Its conjugate supplies the dotted channel
\(\mathcal J^\dagger\bar D n\) relevant to a four-component chain of the form
\(\bar u_L\not p P_Lu_N\), with the sterile Majorana spinor assembled from
\(n,n^\dagger\). Relative to the four-component convention
\(-\bar L Y_\nu\widetilde H N_R+\mathrm{h.c.}\), the undotted Yukawa matrix used
here is \(Y=Y_\nu^*\), with any common sign from the \(\epsilon\) convention
translated consistently. The coefficient \(F\) is likewise conjugated when
passing to the conjugate chiral channel. This observation identifies chirality and the need for both
conjugate sectors; it does not fix the action-to-amplitude sign or identify
\(F\) numerically with the checkpoint's pole coefficient. That final map must
respect its Fourier convention, external-leg orientation, and declared
Majorana normalization. No extra factor of two is inserted into the mixed
operator, which contains distinct composite and sterile fields. The factors
of \(1/2\) reside in the symmetric Majorana and composite contact terms as
shown explicitly.

## Independent computational controls

`verify_sourced_sterile_identity.py` uses only the Python standard library and
exact Gaussian-rational arithmetic, with genuinely complex entries in its
light-by-sterile flavor matrices. It performs 640 checks across 48 deterministic
trials plus 16 Pauli-matrix anticommutator checks. The trials test:

- The exact off-shell identity with arbitrary nonzero \(E,\bar E\).
- Exact solutions of both sourced equations at nonexceptional momenta.
- Failure of the free cubic substitution and of omission of each source.
- Correlated Schur cancellation for rectangular, complex \(F,Y\).
- Failure when the Weinberg contact is dropped or a required transpose is
  replaced by a Hermitian conjugate.
- Exact finite triangular-shift invariance and failure if its quadratic
  contact is omitted.

Run:

```sh
python3 -I -S -B -O verify_sourced_sterile_identity.py
```

The optimized run passed and its output is saved in
`sourced_sterile_identity_results.json`. The analytic derivation establishes
the identity under the declared assumptions. These arithmetic checks are
implementation controls, not proof by sampling, independent human peer review,
or a calculation of the unknown loop mixing vector. The note and code use only
the public checkpoint and the displayed action.

A separate AI-assisted reviewer independently checked the complex-transpose
Schur variation and reran the 640 optimized controls. That review found no
algebraic blocker and requested the explicit Yukawa-conjugation convention
stated above. It is distinct from external human peer review.

## Next required calculation

Match a complete gauge-covariant Green basis to the seed amplitude, including
the light-derivative structures. Apply the sourced light and sterile equations
consistently to all its terms. Carry the resulting contact, kinetic, field and
parameter shifts through the same matching combination before extracting
running or an observable. The identities here provide a precise control that
this fuller calculation must pass; they do not predict its result.
