# Primary comparisons for derivative-seed matter loops

23 October 2026, Pacific time. Input candidate: `f1f9a892e15545c70e4fa05827ea7fefc5a438eb`.

This review concerns a retained-sterile theory with one dimension-seven derivative seed and two neutrino-Yukawa interactions. The useful external comparisons address mass-dependent mixing, literal graph classification, and the additional work between an off-shell pole and physical running. None supplies the present seed's loop coefficient.

## Bounded coverage and retrieval limits

Three successful Firecrawl paper-index queries requested ten records each and returned **7, 4 and 10 records: 21 ranked appearances, 18 distinct papers and three duplicate appearances**. Each used the inclusive created/updated window 1 January 2025 through 23 October 2026. Two metadata inspections and three paper-passage requests returned eleven passages. Other search records received an abstract-level relevance screen.

A fourth, ordinary web-search attempt and a primary-URL scrape were rejected because the connected Firecrawl web/page account had insufficient credits; neither request ran, and neither was retried with changed parameters. The paper-index metadata/passage operations succeeded. Supported managed HTTPS retrieved two primary HTML bodies directly, both HTTP 200. Three previously saved primary bodies were reused and their hashes rechecked. Only the named sections of these five bodies were reviewed. The separate [provenance record](LITERATURE_PROVENANCE.json) records exact requests, counts, body versions and hashes; full texts, abstracts, raw outputs and personalized billing links are excluded from the public payload. This is a bounded literature assessment, not an exhaustive web search or external novelty determination.

## Mass-dependent mixing can populate lower operator dimensions

Renato M. Fonseca, Pablo Olgoso and José Santiago, [Renormalization of general Effective Field Theories: Formalism and renormalization of bosonic operators, 2501.13185v1](https://arxiv.org/html/2501.13185v1), newly retrieved body **22 January 2025**; metadata index update 24 January 2025.

[Section 3](https://arxiv.org/html/2501.13185v1#S3) explicitly reduces a general Green basis to a physical basis after canonical normalization. Its Eqs. (72)-(77) retain mass factors when redundant higher-dimensional interactions shift renormalizable masses, scalar couplings and Yukawas. These are coefficient-coordinate relations in a theory truncated at dimension six. They demonstrate why operator dimension alone does not separate the counterterm problem into independent blocks when dimensionful parameters are retained.

A more direct one-loop mixing benchmark is Jason Aebischer, Luigi C. Bresciani and Nudžeim Selimović, [Anomalous Dimension of a General Effective Gauge Theory II: Fermionic Sector, 2512.16890v1](https://arxiv.org/html/2512.16890v1), reused body **18 December 2025**. [Section 4.3.2, Eqs. (153)-(155)](https://arxiv.org/html/2512.16890v1#S4.SS3.SSS2) displays scalar-quartic running induced by fermionic higher-dimensional operators with explicit fermion masses. In particular, its dimension-five two-fermion/two-scalar insertion appears together with a mass and Yukawas, and its dimension-six two-fermion/three-scalar insertion can appear with two mass powers. The surrounding text also requires additional derivative-operator contributions when assembling the full quartic running.

**Applicability limit:** the first paper truncates its general action at mass dimension six. The second explicitly works through inverse-scale order `1/Lambda^2`, with `L=L4+L5+L6` in its Section 2. A single dimension-seven seed is order `1/Lambda^3`. These references therefore do not compute the present insertion, even when the target operator has dimension four or six. Their group, Weyl, coefficient and permutation conventions also require conversion before any numerical comparison. No template coefficient is copied into this checkpoint.

The elementary engineering-dimension check relevant here is

```text
[C7]=-3, [Y]=0, [M]=1,
[four-lepton operator]=6,
[two-derivative four-scalar operator]=6,
[four-scalar potential operator]=4.
```

Thus coefficient structures of the forms `C7 Y^2 M` for the first two targets and `C7 Y^2 M^3` for the last target are dimensionally allowed. The actual complex conjugations, noncommuting mass order, weak contractions, signs and coefficients must come from the literal diagram calculation. This is dimensional bookkeeping, not a proof that any particular pole survives. A mass-independent subtraction scheme does not make every UV counterterm independent of the masses retained in the action.

## A recent diagram benchmark separates topology from coefficient matching

Fabian Esser, Lukáš Gráf and Chandan Hati, [Cartography of LNV dim-9 SMEFT: Implications for Radiative Neutrino Masses and 0νββ, 2602.18395v1](https://arxiv.org/html/2602.18395v1), newly retrieved body **20 February 2026**; metadata index update 23 February 2026.

[Section 3.1](https://arxiv.org/html/2602.18395v1#S3.SS1) makes a useful methodological distinction: diagram classification with a specified external particle content does not, by itself, distinguish different internal weak or color contractions. Those operator coefficients must be determined by explicit matching. Accordingly, a count of allowed seed/Yukawa slot contractions is not an amplitude symmetry factor, and a field-content sector is not a unique physical operator coefficient. The derivative placement and full weak/flavor exchanges still need calculation.

[Section 2](https://arxiv.org/html/2602.18395v1#S2) explains matching a UV model onto an EFT near its heavy mass scale and subsequent RG evolution. [Section 4.1](https://arxiv.org/html/2602.18395v1#S4.SS1) performs dimension-five/seven/nine matching for an explicit heavy model, then converts an arbitrary Green basis using IBP, EoM and Fierz relations. That example assumes a common heavy mass; the stated hierarchical alternative requires staged matching.

**Use here:** this newer paper supports tracking correlated lower-dimensional operators before phenomenology and checking internal contractions after topology enumeration. Its heavy-model threshold matching is a different operation from extracting local UV poles with sterile fields still retained. Its predictions for neutrinoless double beta decay, neutrino masses and charged-lepton flavor violation are not predictions of the present packet. No matching software, benchmark model, fit or experimental bound from it is applied here.

The screened [D7RGESolver paper, 2505.06499](https://arxiv.org/abs/2505.06499), similarly concerns dimension-five/seven **pure-SMEFT** running. Its abstract-level screen does not provide a retained-sterile seed RGE; the software was not executed.

## From a cancelling pole to a scale derivative

Fonseca et al. [Section 4.1, Eqs. (89)-(93)](https://arxiv.org/html/2501.13185v1#S4.SS1) first canonically normalize and reduce the divergent action. It then writes the bare coefficient as

```text
a_i,bare = mu^((n_i-2) epsilon) (a_i+delta a_i),
```

where `n_i` counts fields and `d=4-2 epsilon`. The displayed counterterm cancels the divergent action coefficient; the one-loop scale relation is derived from the bare expansion. This independently benchmarks the need to distinguish an inverse-kernel pole, an action pole, a parameter pole and its beta function.

For the present primitive packets, the corresponding bare-mu powers and the Yukawa/seed homogeneity must be derived in their own coefficient convention. The reference's compact relation is stated after its full normalization/reduction steps; it cannot turn a single unreduced source or diagram residue into a complete physical running coefficient. In particular, no proper-vertex zero for one chosen external configuration implies that the other matter sectors vanish.

## Redundant coefficients and quantum sources remain relevant

Di Zhang, [Two-loop Renormalization Group Equations in the νSMEFT, 2504.00792v2](https://arxiv.org/html/2504.00792v2), reused body **5 May 2025**, supplies a retained-sterile reference through dimension five. [Section 3, Eqs. (13)-(15)](https://arxiv.org/html/2504.00792v2#S3) performs explicit field changes and states that the reduction formulas concern **bare** couplings before expansion into renormalized quantities and counterterms. Its separately ordered Majorana mass products are fixed by the field convention. [Section 6](https://arxiv.org/html/2504.00792v2#S6) retains derivatives with respect to redundant coefficients when deriving RG functions before restricting their coordinate values. The scope is different from the dimension-seven insertion, but the bookkeeping requirement remains pertinent.

Kevin Falls, [Equivalence of effective actions, 2504.17851v3](https://arxiv.org/html/2504.17851v3), reused body **22 August 2025**, [Section I, Eqs. (2)-(5)](https://arxiv.org/html/2504.17851v3#S1), distinguishes a microscopic nonlinear field redefinition from the transformation of the mean field that parametrizes a full quantum Legendre effective action. Beyond tree level, the latter generally includes quantum corrections. Thus a correct classical sourced action substitution does not, by itself, determine the complete transformed 1PI functional or composite-source renormalization. The served body also has a separate later dateline; the recorded arXiv body version date is retained here.

## UV/IR and four-fermion scheme controls

Zhang [Section 4, Eqs. (20)-(21)](https://arxiv.org/html/2504.00792v2#S4) discusses the spurious infrared divergences introduced by setting masses and momenta to zero and an auxiliary-mass propagator rearrangement. The physical Majorana matrix and an auxiliary IR mass have distinct roles. A scaleless integral vanishing in dimensional regularization does not establish that its UV residue vanishes. The present simple-pole calculation should retain off-shell kinematics or a controlled IR rearrangement and demonstrate regulator cancellation where used.

Fonseca et al. [Section 2.1](https://arxiv.org/html/2501.13185v1#S2.SS1) explains that four-fermion structures related in four dimensions by Fierz identities can produce evanescent shifts relevant to finite matching. Consequently an exact four-dimensional Grassmann/Fierz action check is valuable, but it does not fix a dimensional/evanescent prescription for finite amplitudes or later-loop running.

These comparisons support the calculation and its stated limits. Complete mixed matter-sector residues, their source-preserving physical reduction, other interactions, finite threshold matching and an observable remain separate requirements. No primary source in this bounded screen establishes the present seed's full physical beta function, external novelty, empirical confirmation or unification.
