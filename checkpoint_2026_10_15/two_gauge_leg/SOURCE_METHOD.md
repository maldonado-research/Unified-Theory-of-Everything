# Two-hypercharge-leg source and method audit

15 October 2026. This targeted note addresses a named gap after the zero- and
one-B calculations: the canonical scalar two-B vertex, the two-B Ward identity,
and whether an independent local degree-one two-B term could evade that identity.
It derives the relevant formulas from the same literal action. It is not a new
literature survey, an operator-novelty claim, or an evaluation of a three-B loop.

## 1. Primary sources checked

Gitte Elgaard-Clausen and Michael Trott, *On expansions in neutrino effective
field theory*, [arXiv:1703.04415v3](https://arxiv.org/abs/1703.04415v3),
[Eqs. (2)–(3)](https://arxiv.org/html/1703.04415v3#S2.EGx1), provides the
canonical scalar/fermion kinetic terms and explicitly declares
`D_mu = partial_mu + i g1 y B_mu + ...`. These agree with the project's D-plus
convention. The source's Eq. (6) also supplies the charge-conjugate sterile
Yukawa channels. Its kinetic action, rather than a convention-dependent table
of photon vertices, is the normalization checked here. The one- and two-B
vertices below are our direct differentiations of that action.

Ansgar Denner, Stefan Dittmaier and Georg Weiglein, *Application of the
Background-Field Method to the electroweak Standard Model*,
[arXiv:hep-ph/9410338v2](https://arxiv.org/abs/hep-ph/9410338v2),
[primary PDF](https://arxiv.org/pdf/hep-ph/9410338v2), Section 3,
Eq. (23), states the functional background-gauge invariance of the effective
action. Differentiation gives its vertex Ward identities; Eq. (36) is an
explicit Abelian fermion example. Section 4 immediately before Eq. (43)
explains that counterterms must satisfy the same Ward identities as bare and
renormalized vertices. These are method references: their broken-electroweak
vertex conventions are not imported into the present unbroken-hypercharge
kernel. Our charge and Fourier signs are derived below.

Only these two specific method/action sources were checked. The existing
same-day broader arXiv/INSPIRE survey was not repeated. In particular, the
polynomial uniqueness proof in Section 4 is supplied here and is not attributed
to either paper as a new theorem from those sources.

## 2. Canonical scalar vertices, including the seagull

Set `a = g' y_L`, `phi = tilde H`, `y_phi = y_L = -1/2`, and
`chi = L^c`, `y_chi = -y_L`. Use `D_mu phi = (partial_mu+i a B_mu)phi`,
metric `(+---)`, incoming scalar `exp(-i s.x)`, and outgoing conjugate scalar
`exp(+i s'.x)`. For emitted B momenta, `s-s' = sum k_i`. The canonical term is

    (D_mu phi)^dagger D^mu phi
      = partial_mu phi^dagger partial^mu phi
        + i a [(partial_mu phi^dagger)phi
               - phi^dagger(partial_mu phi)] B^mu
        + a^2 phi^dagger phi B_mu B^mu.

It gives the action kernels and vertices

    K_phi,B^nu(s,s') = -a(s+s')^nu,
    V_phi,B^nu       = -i a(s+s')^nu,
    K_phi,BB^mu nu   = +2 a^2 g^mu nu,
    V_phi,BB^mu nu   = +2 i a^2 g^mu nu.

The two in the last line is the functional-derivative factor for the two
identical B fields. It is already present in the vertex; do not multiply an
additional factor two when using that seagull in a labelled two-B amplitude.
The sign of the quadratic term is positive in the canonical scalar action.
A loop prefactor additionally includes all propagator and vertex phases;
copying the one-B scalar charge sign into a putative two-B seagull sign is
incorrect. Charge conjugation reverses the single-B fermion charge but not
its square. The canonical fermion kinetic action has a single-B vertex and
no fermion two-B seagull.

A useful scalar identity, directly from the displayed action kernel, is

    k_nu K_phi,B^nu(s,s-k)
      = -a [s^2-(s-k)^2].

The scalar mass cancels from this difference. At two B fields the canonical
seagull is required along with the two orderings of sequential scalar
insertions. For the selected bubble, one must also include both fermion-line
orderings and the two distributions with one B on each internal line. The
seed is a derivative of the neutral composite and has no B contact; the
nonderivative Yukawa has none either. These observations determine the
necessary topology classes, rather than evaluating their integrals.

## 3. Local Ward recursion with this Fourier convention

Let `K_n` be the local action kernel of `bar L phi N_L` with n labelled
outgoing B fields. All matter/gauge momenta are off shell; N is incoming and

    p = r+h+k_1+...+k_n.

The transformations are

    delta bar L = +i a alpha bar L,
    delta phi   = -i a alpha phi,
    delta N     = 0,
    delta B_mu  = partial_mu alpha.

Varying the action and taking the coefficient of an outgoing gauge parameter
with momentum k1 gives, for two B fields,

    k1_mu K_2^mu nu(r,h;k1,k2)
      = a [K_1^nu(r,h+k1;k2)-K_1^nu(r+k1,h;k2)].

The right-hand side includes the shifts of the charged external fields. A
single charged amplitude is consequently not separately transverse. An
Abelian background Ward identity applies to the gauge-invariant effective
action. The usual linear Abelian gauge-fixing term contributes only a pure
quadratic gauge term; differentiating in the mixed matter channel removes
it. No gauge propagator or closed anomalous chiral trace occurs in the
selected open-chain loop; no non-Abelian Slavnov–Taylor identity is being
assumed for it.

The same formula with 1 and 2 exchanged must hold. In stripping common
factors, take

    K_1 = a R_1,
    K_2 = a^2 R_2,
    R_1.epsilon = [-(2r+k).epsilon/3+(2h+k).epsilon/6] slash p
        +(p^2/6-m_H^2/2) slash epsilon
        +[slash k,slash epsilon] slash r/4
        +[(k.epsilon)slash k-k^2 slash epsilon]/3.

The previously chosen covariant action predicts

    R_2(epsilon1,epsilon2)
      = (epsilon1.epsilon2) slash p
        - ([slash k1,slash epsilon1] slash epsilon2
           +[slash k2,slash epsilon2] slash epsilon1)/4,

with the common right projector P_L understood. Contracting its first B leg
gives exactly

    R_2(k1,epsilon2)
      = (k1.epsilon2)slash p
        -[slash k2,slash epsilon2] slash k1/4
      = R_1(r,h+k1;k2,epsilon2)
        -R_1(r+k1,h;k2,epsilon2).

In this subtraction p remains fixed: both shifted one-B kernels have total
momentum p. The opposite leg follows by Bose exchange. This is an independent
check of the action's predicted kernel against the lower-point Ward recursion.
It must not be reported as an independent integration of the two-B diagrams.

## 4. A precise uniqueness result for the two-B local pole

Assume two candidate kernels have all the following properties:

1. They are local polynomials in independent off-shell momenta r,h,k1,k2,
   after eliminating p by momentum conservation. Coefficients may include
   the metric, epsilon tensors, constant gamma matrices, flavor matrices
   and masses, but no inverse momentum invariants or distributions.
2. They agree in the zero- and one-B sectors and obey the same linear
   Abelian Ward recursion, with the same field charges and conventions.
3. Their two-B terms have total momentum degree at most one, as appropriate
   to the selected dimension-seven L-H-N-B-B UV counterterm. There are no
   extra background spurions or higher-dimensional coefficients being
   reassigned to this sector.

Their difference `T^mu nu` then satisfies

    k1_mu T^mu nu = 0,       k2_nu T^mu nu = 0

as polynomial identities, not merely between on-shell polarizations.
Expand T by degree in k1. The term in `k1_mu T^mu nu` of degree one in k1
is `k1_mu T^mu nu|_(k1=0)`. Its coefficients must each vanish; hence every
component of T vanishes at k1=0. Every monomial of T therefore contains a
component of k1. The same reasoning with k2 shows that every monomial also
contains a component of k2. Equivalently, in the component polynomial ring,

    T belongs to (k1) intersection (k2) = (k1)(k2),

because these two sets of variables are disjoint. Every nonzero such
polynomial has momentum degree at least two. A degree-at-most-one T must
therefore vanish identically.

The proof works componentwise for the open spinor chain. Restricting to the
fixed odd-Clifford/chirality sector cannot create an exception. An epsilon
or dual-field-strength tensor supplies constant coefficients and does not
remove either required momentum; Bianchi and Schouten identities may
eliminate structures but cannot produce a lower-degree nonzero polynomial.
For example, an epsilon expression linear in k1 can be transverse in k1,
but cannot also be transverse in an independently variable k2. The proof
even allows constant tensors more general than Lorentz-invariant ones;
absence of fixed vectors is therefore more restrictive than it needs.

This does not say that two-B kernels of arbitrary momentum degree have no
new gauge-invariant structures. Products of two Abelian field strengths
start at two momenta and can evade lower-point matching when their field
content, chirality and canonical dimension permit them. They cannot be
added as a local degree-one ambiguity in the stated sector. Nor does it
say that the one-B kernel was fixed by its Ward identity: a homogeneous
one-B term can already have one momentum, and the lower-point transverse
remainder had to be calculated independently.

The independence and locality conditions matter. Setting k2=-k1 before
checking the identity loses an independent variable; restricting to on-shell
external states loses polynomial constraints. A nonlocal double-transverse
projector can contain inverse dot products and evade the degree argument.
Those cases are excluded from an off-shell local UV pole. If only a subset
of the lower-point Green basis has been matched, the corresponding higher-
point difference need not obey the homogeneous identity used in this proof.

Consequently a complete polynomial one-B match plus the correct Ward
identities fixes the degree-one two-B completion uniquely in this sector.
A separate two-B loop integration remains a valuable check of the diagrams,
phases, seagull normalization, and completeness assumptions; it is not
replaced by calling the uniqueness argument a loop calculation.

## 5. Three-B consequence, explicitly limited

The displayed chosen action contains at most two B fields: its charged
D^2 terms contain at most two connections, its Q_sigma term contains one
charged D and one Abelian field strength, and Q_partialB has one field
strength and no charged derivative. Therefore its three-B action kernel
vanishes. The Ward recursion checks this result: R_2 depends on r,h only
through p, so its row-shift and scalar-shift values are equal and their
difference is zero.

In the same dimension-seven sector a local three-B UV term has momentum
degree zero. If another completion agrees through two B fields, its
three-B difference is separately transverse and must vanish by the same
polynomial argument (already one independent gauge contraction kills a
constant). This is a conditional action/Ward consequence. No three-B loop
integral has been evaluated here. It is not a statement that the full
finite three-B amplitude, higher-dimensional local operators, or other
coupling/field sectors vanish.

## Retrieval and integrity record

All four requests below returned HTTP 200 on 15 October 2026 using inherited
managed HTTPS proxy and CA trust, with TLS verification enabled:

| Request | Response bytes | SHA256 |
| --- | ---: | --- |
| `https://arxiv.org/html/1703.04415v3` | 525092 | `b8e598a392bc19c2bbec98dff54487f62ceb9e891598d1e1db7044c9bbe3f910` |
| `https://export.arxiv.org/api/query?id_list=hep-ph/9410338` | 2729 | `79190160172eb8ba7f7d4da136ab1d1cb080cf9d2b38f0567491aec084178071` |
| `https://arxiv.org/pdf/hep-ph/9410338` | 341244 | `1ebacda7dd31d563b9af4bea30c5c517fbd6a996f12652eba5bdea636b1ddca5` |
| `https://arxiv.org/pdf/hep-ph/9410338v2` | 341244 | `1ebacda7dd31d563b9af4bea30c5c517fbd6a996f12652eba5bdea636b1ddca5` |

The metadata and PDF front page both identify the checked Ward source as
version 2. The unversioned and versioned PDF bytes agree. Downloaded papers
and full-text extractions are research-cache material and are not public
release files; this newly authored equation/method note is suitable for
review and publication. The source coverage is bounded to the named action
and Ward-method questions; no exhaustive-web, novelty, physical-running,
observable, or unification claim follows.
