# Literature benchmark for sourced action completion and off-shell Hessians

Accessed 16 October 2026. This is a bounded methodological review. It does not supply a selected-packet UV coefficient, a physical beta function, a finite amplitude or evidence of unification. The action algebra and its checks are separate artifacts in this checkpoint.

## Coverage and primary evidence

Four Firecrawl paper-index queries requested k=10 each and returned **6, 10, 10 and 10 records: 36 ranked appearances, 32 distinct papers and four duplicate appearances**. Three queries used the inclusive created/updated window 2025-01-01 through 2026-10-16. A fourth deliberately used 2018-01-01 through 2026-10-16 to recover foundational field-redefinition work. These are ranked, bounded searches, not an exhaustive web or novelty review.

Four selected papers received metadata inspections and fresh primary HTML retrievals, all HTTP 200. Three passage requests returned **3, 4 and 4 passages**. We read the specific body sections cited below and checked the retrieved byte hashes. A previously retrieved sterile-basis PDF was reused and its hash rechecked. No complete paper or supplement was exhaustively reviewed. Full bodies, verbatim abstracts and raw provider outputs stay in scratch storage; public provenance retains source identifiers, abstract hashes/character counts, our relevance summaries and retrieval receipts.

## A direct 2025 benchmark for background-source corrections

Timothy Cohen, Xu-Xiang Li and Zhengkang Zhang, [Geometric Building Blocks of Effective Field Theory Amplitudes, 2509.20449](https://arxiv.org/html/2509.20449), inspected body **v1, 24 September 2025**. The index update is 26 September 2025.

[Section 2.3, Eq. (28)](https://arxiv.org/html/2509.20449#S2.E28) distinguishes two conditions: the potential gradient vanishes at the physical vacuum, while inverse-propagator contributions disappear after the relevant vacuum, external-state and on-shell operations. These conditions cannot be silently imposed on an arbitrary off-shell background.

[Section 2.5, Eqs. (47)-(49)](https://arxiv.org/html/2509.20449#S2.E47) provides the most useful explicit benchmark for this continuation. Reexpressing the background source in nonlinear fluctuation coordinates generates additional quadratic and higher terms proportional to the potential gradient. Eq. (48) keeps the quadratic contribution `-Gamma^i_jk U_,i xi^j xi^k/2`, which restores the ordinary off-shell propagator. Eq. (49) gives the corresponding higher-vertex gradient terms. Dropping these terms before taking the stationary limit changes the off-shell object being calculated.

The paper treats scalar EFT geometry and tree-level amplitude building blocks. It does not calculate the retained-sterile packet or justify importing its scalar coefficients. Here it corroborates the need to transform background sources and retain gradient-dependent terms.

### The elementary Hessian identity to check

For real commuting scalar coordinates Q, a field map `Q^a=f^a(q)` and an action A, the ordinary chain rule gives

    E_a = delta A / delta Q^a,
    H_ab = delta^2 A / (delta Q^a delta Q^b),

    H'_ij = f^a_,i H_ab f^b_,j + E_a f^a_,ij.

For field theories, indices include spacetime and the products include the corresponding integrals and derivative kernels. The displayed form is a commuting-scalar check; it is not a prescription for replacing ordered Grassmann variations by ordinary matrix derivatives.

The second term vanishes at a stationary background, or for an affine map. In a sourced action with a source linear in the original Q, `E_a` includes that external source. A source chosen to hold a background fixed therefore also has to be transformed.

Crucially, for `f=q+eta h(q)`, the omitted term contains `eta E_a h^a_,ij`: it is already **linear in eta**. Working only to first order in the field-redefinition coefficient does not remove this Hessian correction. For a simple negative control, `A(Q)=t Q+m^2 Q^2/2` and `Q=q+kappa q^2/2` give `H'(0)=m^2+kappa t`, whereas a Jacobian-only calculation gives `m^2`. The stationary limit `t=0` and affine limit `kappa=0` recover the simpler result. This is a standard chain-rule consequence, not new mathematics.

## Tadpoles must accompany propagation away from a minimum

Timothy Cohen, Matthew Forslund and Andreas Helset, [Field Redefinitions Can Be Nonlocal, 2412.12247](https://arxiv.org/html/2412.12247), inspected body **v2, 3 October 2025**; first posted in December 2024, with an index update on 7 October 2025.

[Section 2.3.2](https://arxiv.org/html/2412.12247#S2.SS3.SSS2), especially Eqs. (2.44)-(2.50), explicitly shifts a scalar theory away from its minimum. New quadratic terms and tadpole insertions must be kept together: the resummed propagator agrees with the original one, and tadpole insertions restore a cubic interaction that appears absent in the shifted Lagrangian. [Section 3.1](https://arxiv.org/html/2412.12247#S3.SS1), Eqs. (3.1)-(3.2), states that correlation functions agree when sources are compared consistently; changing the source coupling computes a different quantity.

**Use here:** retain the tadpole/background-gradient contribution in an off-shell normalization or Hessian test. This methodological use does not require adopting nonlocal field maps or resumming beyond the packet's declared perturbative order.

## Classical field substitution and a quantum mean-field map are distinct

Kevin Falls, [Equivalence of effective actions, 2504.17851](https://arxiv.org/html/2504.17851), inspected body **v3, 22 August 2025**; index update 25 August 2025. The served body also contains a separate 2026 dateline; we do not treat that as a new arXiv version date.

[Section I, Eqs. (2)-(5)](https://arxiv.org/html/2504.17851#S1.E3) distinguishes a microscopic nonlinear field map from the corresponding map of mean fields. In particular, Eq. (3) states that ordinary Legendre effective actions are not generally related beyond tree level by simply substituting the classical microscopic map. Eq. (5) instead uses a quantum-corrected effective transformation, whose existence is studied in the paper.

[Section II, Eqs. (6)-(12)](https://arxiv.org/html/2504.17851#S2.E9) specifies the different source couplings and their Legendre transforms. The same underlying functional integral can generate correlations of different composite operators, with different effective-action arguments. Thus `mean(f(Q))` is not generally `f(mean(Q))` for a nonlinear f.

**Use here:** the elementary action/Hessian chain-rule check above is valid for a declared pulled-back functional. It must not be advertised as having computed the full transformed quantum 1PI functional, its composite-source renormalization or physical running. No Vilkovisky-DeWitt construction or quantum mean-field map has been implemented in this checkpoint.

## Source-retaining EoMs and EFT order

Juan Carlos Criado and Manuel Perez-Victoria, [Field redefinitions in effective theories at higher orders, 1811.09413](https://arxiv.org/html/1811.09413), inspected body **v2, 20 September 2019**; first posted November 2018, with index update 23 September 2019. This known foundational reference was inspected directly in addition to the bounded search results.

[Section 2, Eq. (2)](https://arxiv.org/html/1811.09413#S2.E2) transforms the action, source coupling and functional Jacobian together. [Section 3, Eq. (10)](https://arxiv.org/html/1811.09413#S3.E10) expands the actual field substitution and shows why EoM replacement reproduces only the linear term in the field change. Section 5 discusses how a change at one EFT order induces changes at higher orders, which are not recovered by simply applying a more accurate EoM to an isolated operator.

**Use here:** substitute the chosen fields into every term that contributes at the retained EFT and loop orders, then expand consistently. Keep interaction sources such as masses, Yukawa terms and the scalar potential. These interaction sources in Euler-Lagrange expressions are distinct from the external probe sources used to define correlation functions; both need consistent bookkeeping.

## Retained-sterile scalar/derivative inventory

Yi Liao and Xiao-Dong Ma, [Operators up to Dimension Seven in Standard Model Effective Field Theory Extended with Sterile Neutrinos, 1612.04527](https://arxiv.org/abs/1612.04527), reused PDF **v1, 14 December 2016**. The index's later update is not substituted for this body version.

Section 3 gives useful target classes: Eq. (14) has the two sterile/lepton/Higgs derivative-current operators; Eq. (15) has the charged-lepton/sterile two-Higgs two-derivative class; Eq. (18) has the sterile bilinear times `(H-dagger H)^2`; and Eq. (19) has two sterile-bilinear/Higgs derivative structures, including `(N^C N)(D_mu H)-dagger D^mu H`. They are references for organizing surviving terms, not a coefficient map for the current packet.

Convert `phi=tilde H`, sterile chirality/conjugation, epsilon contractions, derivative adjoints, flavor order, generator normalization and gauge-sign conventions before identifying a packet term with one of these representatives. EoM ellipses in a published basis relation are not permission to delete the explicit interaction sources in the present action. The selected tensor/current sectors previously checked also do not close this scalar/derivative sector or the full dimension-seven EFT.

The practical result of this review is an independent benchmark for the gradient-dependent off-shell correction and a clear separation between a sourced classical action identity, a quantum effective-action transformation and a physical prediction. No source in this screen supplies the missing packet-wide beta function or finite amplitude.
