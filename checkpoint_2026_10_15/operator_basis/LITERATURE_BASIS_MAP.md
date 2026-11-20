# Primary-source comparison for the selected dimension-seven basis sectors

Accessed 15 October 2026. This note compares operator definitions, flavor identities and field content. It does not calculate a counterterm, an anomalous dimension, a finite amplitude or an observable. The packet-specific algebra and tests are separate artifacts in this checkpoint.

## Bounded coverage

Three Firecrawl paper-index queries requested k=10 each. Two used the inclusive created/updated window 2025-01-01 through 2026-10-15; a third deliberately used 2000-01-01 through 2026-10-15 to recover Majorana operator-classification references. They returned **10, 6 and 10 records: 26 ranked appearances, 22 distinct papers and four duplicate appearances**. Exact queries, ranks, identifiers, scores and our relevance summaries are in `LITERATURE_PROVENANCE.json`.

The returned records received abstract-level screening. The comparison below instead rests on four primary bodies: three already retrieved in the preceding continuation and one newly fetched older benchmark. Their exact bytes were checked against SHA256 receipts. One new canonical metadata inspection was performed for the known older benchmark. No new passage-index requests or broad site crawls were needed. No selected 2025-2026 search hit supplied a more direct replacement for these basis definitions.

This is a focused follow-up, not an exhaustive search or a claim that any classification is new. Full bodies, verbatim abstracts and raw provider outputs stay in scratch storage. Public provenance contains identifiers, abstract hashes and character counts, receipts and our own screening summaries.

## Active-only two-lepton / two-Higgs / one-strength operators

Di Zhang, [Renormalization Group Equations for the SMEFT Operators up to Dimension Seven, 2306.03008](https://arxiv.org/html/2306.03008). The inspected body identifies **v2, 16 October 2023**; the paper-index update date is 18 October 2023.

[Table 1](https://arxiv.org/html/2306.03008#S2.T1) defines, with flavor labels restored,

    O_lHB^(alpha beta)
      = epsilon^(ab) epsilon^(de)
        (l_alpha,L^a C sigma_mn l_beta,L^d) H^b H^e B^mn,

    O_lHW^(alpha beta)
      = epsilon^(ab) (epsilon tau^I)^(de)
        (l_alpha,L^a C sigma_mn l_beta,L^d) H^b H^e W^(I,mn).

Here `tau^I` denotes the source's Pauli matrices, written `sigma^I` there, to distinguish them from its Lorentz tensor `sigma_mn`. The table has **no explicit gauge coupling** multiplying these monomials. Its independent flavor counts are

| Source monomial | Flavor count for n active flavors | n=1 | n=3 |
| --- | ---: | ---: | ---: |
| O_lHB | n(n-1)/2 | 0 | 3 |
| O_lHW | n^2 | 1 | 9 |

The relation `O_lHB^(alpha beta) + O_lHB^(beta alpha) = 0` appears explicitly in [Section 2.1, Eq. (5)](https://arxiv.org/html/2306.03008#S2.E5). [Table 2](https://arxiv.org/html/2306.03008#S2.T2) keeps the antisymmetric B representative and an unrestricted ordered W coefficient. These are counts for the displayed chiral monomials and their flavor coefficients; the Hermitian-conjugate completion is not a second independently assigned set of coefficients.

This agrees structurally with the present tensor inventory. In the conjugated-scalar convention `phi=tilde H`, the B weak tensor `phi^a phi^b` is symmetric, leaving a flavor-antisymmetric coefficient. The ordered W tensor

    P^(I,ab) = (T^I phi)^a phi^b,       T^I = tau^I/2,

has weak-symmetric and weak-antisymmetric parts. The former pairs with flavor-antisymmetric coefficients and the latter with flavor-symmetric coefficients. Their dimensions add to

    n(n-1)/2 + n(n+1)/2 = n^2.

Thus one conventional ordered W monomial can encode both flavor sectors. Calling them two projected sectors does not double the source's `n^2` count.

### A separate check on W flavor projection

Nicole F. Bell et al., [Model Independent Bounds on Magnetic Moments of Majorana Neutrinos, hep-ph/0606248](https://arxiv.org/abs/hep-ph/0606248). The newly retrieved PDF identifies **v2, 10 October 2006**; the index update of 11 November 2009 is distinct.

Section II, Eqs. (7)-(8), uses B and W operators with **explicit g-prime and g**, respectively, and a Pauli insertion in the W operator. Eq. (10) defines

    O_W^(+/-)_(alpha beta)
      = [O_W_(alpha beta) +/- O_W_(beta alpha)]/2.

The text states that the ordered W monomial is neither flavor symmetric nor flavor antisymmetric. Appendix A, Eq. (A4), relates the second displayed W contraction to the flavor transpose of the first; the discussion of Eq. (A5) identifies another frequently used contraction with `2 O_W^-`. These explicit projection factors matter when comparing representative coefficients.

Section II following Eqs. (12)-(13) also states that the flavor-symmetric `O_W^+` does not contribute to the tree-level neutrino magnetic interaction, whereas the B operator and `O_W^-` do. Therefore the existence of a symmetric W flavor sector is compatible with the absence of a diagonal Majorana magnetic moment in this sector. It does not turn a selected multi-gauge UV contact into a physical magnetic-moment prediction. We have not reproduced this paper's loop matching, naturalness estimates or phenomenological bounds.

### Normalization and orientation conditions

The agreement asserted here is structural and concerns flavor dimensions and projection identities. A coefficient equality requires additional conversion:

- **Scalar and fermion orientation:** the source has two unbarred lepton fields and two H fields. The present representative uses barred leptons, charge-conjugated leptons and two `phi=tilde H` fields. Choose the appropriate Hermitian conjugate and retain the resulting flavor order, epsilon contractions and spinor order.
- **Generator normalization:** the source's Pauli insertion is `tau^I=2 T^I`. The present W operator uses T, so a factor of two must be accounted for before comparing coefficients. This condition alone does not fix the final sign or flavor transpose.
- **Gauge convention:** 2306.03008 Section 2 defines `D=partial-i g1 Y B-i g2 T^I W^I`. The present packet uses all-Dplus. Bell instead explicitly writes weak curvature with nonlinear term `-g epsilon^(abc) W^b W^c`, and uses `epsilon=-i tau2`. A full field/strength conversion, not an isolated sign replacement, is required.
- **Couplings and Lorentz conventions:** the modern table has no explicit g factors, while Bell's Eqs. (7)-(8) include them. Fix the Lorentz-sigma normalization, charge-conjugation convention, Hermitian completion and coefficient summation convention before using any numerical map.

No source coefficient or beta-function entry has been imported into the selected representative `-F Y_nu^T/12` from these comparisons. That coefficient remains an input established by the packet calculation; it is not established by matching the literature counts.

## Retained-sterile four-fermion / Higgs current classes

The present current descendants have the schematic form

    (bar L gamma_mu N_L) phi (bar f gamma^mu f),

with singlet or weak-triplet contractions where allowed. With `N_L=N_R^c`, these must be compared with the appropriate conjugated orientation in sterile bases written with N_R. The vector-current representation is a useful representative, not automatically the published scalar-Fierz basis monomial.

Yi Liao and Xiao-Dong Ma, [Operators up to Dimension Seven in Standard Model Effective Field Theory Extended with Sterile Neutrinos, 1612.04527](https://arxiv.org/abs/1612.04527), inspected PDF **v1, 14 December 2016**, gives the relevant scalar-Fierz representatives in **Section 3, Eq. (25)**. Hao-Lin Li et al., [Operator Bases in Effective Field Theories with Sterile Neutrinos: d <= 9, 2105.09329](https://arxiv.org/abs/2105.09329), inspected PDF **v1, 19 May 2021**, provides a separately organized inventory in **Section 4.3**. The following is a field-content correspondence, not an equality of coefficients:

| Spectator current f | Liao-Ma Eq. (25) family | Li et al. Section 4.3 |
| --- | --- | --- |
| L | O_LNLLH | Eq. (4.84), l^2 lbar N H with explicit flavor projectors |
| Q | O_QNQLH1 and O_QNQLH2 | Eq. (4.78), two l N q qbar H contractions |
| e_R | O_eLNeH | Eq. (4.85), ebar e l N H |
| u_R | O_uLNuH | Eq. (4.81), l N ubar u H |
| d_R | O_dLNdH | Eq. (4.83), dbar d l N H |

Liao-Ma explains that its vector bilinears are transformed into scalar ones using Fierz identities. Its Table 4 tracks the two quark-doublet contractions and the repeated-lepton family. The Li et al. Young projectors make the repeated-field flavor decomposition explicit. In particular, for the repeated-L sector, weak-current contractions and flavor permutations must be reduced together; counting each visible singlet/triplet current as a separate arbitrary flavor tensor can overcount. A flavor-diagonal or universal gauge current is also a restricted coefficient slice of these general flavor inventories.

These references establish suitable target classes for the surviving source terms. A complete map still requires an explicit Fierz conversion with weak and color contractions, field ordering, sterile conjugation, flavor permutations and Hermitian completion fixed. This note does not classify every retained-sterile dimension-seven sector: sources with additional sterile fields, derivatives or different external matter content require their own treatment. Neither the active tensor counts nor the current-class correspondence proves closure of the full EFT or supplies its physical running.
