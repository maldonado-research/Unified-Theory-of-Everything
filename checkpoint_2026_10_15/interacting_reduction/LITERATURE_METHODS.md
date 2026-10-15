# Bounded literature follow-up for interacting source reduction

Accessed 15 October 2026. This is a methodological follow-up to candidate.6, not an exhaustive web review, a novelty assessment, or a calculation of the selected packet's physical anomalous dimension. Exact query/result provenance, screening dispositions, body sections, version dates and retrieval receipts are in `LITERATURE_PROVENANCE.json`.

## Coverage and what was actually read

Four complementary Firecrawl paper-index searches requested k=10 each and returned **10, 4, 6 and 10 records: 30 ranked appearances, 23 distinct canonical papers and seven duplicate appearances**. Two searches used the inclusive created/updated window 2025-01-01 to 2026-10-15. Two deliberately widened that window to 2024 and 2020 to recover basis and mixed-order reduction work that recent-only searches miss. An index update date is not a preprint version date.

Five selected papers received canonical metadata inspection; four received passage requests. Those requests returned **0, 4, 3 and 4 passages**. The index lacked full text for 2105.09329. Direct managed HTTPS supplied two primary HTML bodies and two primary PDFs, with an additional abstract page used to check the 2105.09329 version history. Its HTML-body endpoint returned 404 before the successful PDF retrieval. The already saved 2504.00792v2 HTML was reused and its SHA256 rechecked. The selected sections below were read; no whole paper or supplement was exhaustively reviewed, and no external symbolic package or published RGE solver was run.

The five primary bodies are the evidence for the methodological discussion. Other search hits received an abstract-level relevance screen only. The 23 papers are this round's distinct search records, not 23 newly discovered papers; some appeared in earlier repository searches. Full texts, verbatim abstract metadata, provider passage dumps and raw provider outputs remain in scratch storage rather than being redistributed. Public provenance retains abstract SHA256 hashes and character counts, source identification, query/rank/date metadata and our own relevance summaries.

## Corrected retained-sterile basis: retain the generated sources

Yi Liao and Xiao-Dong Ma, [Operators up to dimension seven in standard model effective field theory extended with sterile neutrinos, 1612.04527](https://arxiv.org/abs/1612.04527). The retrieved PDF identifies **v1, 14 December 2016**; the index separately reports an update on 19 July 2017. This source was reached as a cited antecedent of the inspected 2407.16751 paper, not as an additional search result.

Section 3, Eq. (17), selects four operators containing a sterile/lepton current, a **derivative on the scalar**, and either B or W strength or its dual. Table 3 and Eqs. (32)-(34) explain why some earlier operators containing a dual strength and a **derivative on the lepton** are redundant. Their reduction generates an additional neutrino-Yukawa-weighted scalar/sterile tensor operator involving B strength. Thus a derivative moved through a chiral identity or an interacting equation of motion does not justify deleting the interaction source.

The paper's abstract and introduction report seven redundant dimension-six operators, sixteen redundant dimension-seven operators and two missing dimension-seven operators in the older classifications it revisits. These are the authors' comparisons in their conventions. We have not recomputed their entire basis, and none of these counts is a closure certificate for the selected packet here. In particular, the redundancy of the derivative-on-lepton dual-strength representatives does **not** eliminate every dual-strength scalar-derivative representative.

Hao-Lin Li et al., [Operator Bases in Effective Field Theories with Sterile Neutrinos: d <= 9, 2105.09329](https://arxiv.org/abs/2105.09329), supplies a later independent classification. The retrieved PDF identifies **v1, 19 May 2021**, consistent with the checked arXiv submission history; the index update of 24 November 2021 is separate. Sections 2-3 explain amplitude/operator correspondence, integration-by-parts and identical-field flavor projection. Section 4.3 lists the dimension-seven sterile operators. These supply a target inventory against which to map the packet's surviving sources, including conjugation and flavor symmetry. They do not provide its counterterm coefficients.

**Use here:** map each surviving weak-current, scalar-source, tensor and multi-fermion term to an explicit chosen retained-sterile basis. Preserve field order, flavor order and Hermitian completion. Translate the source's epsilon, chirality, dual-strength, gauge and charge-conjugation conventions before comparing coefficients.

## Mixed EFT orders: a dimension-six change can affect dimension seven

Di Zhang, [Renormalization Group Equations for the SMEFT Operators up to Dimension Seven, 2306.03008](https://arxiv.org/html/2306.03008), retrieved body **v2, 16 October 2023**; the index update is 18 October 2023.

[Section 2](https://arxiv.org/html/2306.03008#S2), especially the EoMs through order Lambda^-1 and [Table 3](https://arxiv.org/html/2306.03008#S2.T3), shows why redundant dimension-six operators removed with Higgs/lepton equations can also contribute to physical dimension-seven coefficients when the Weinberg operator is present. [Appendix A](https://arxiv.org/html/2306.03008#A1), Eqs. (41)-(42) and the following reduction relations, explicitly retains these mixed-order terms. The conclusion states that the corresponding dimension-six counterterms are required along with dimension-seven counterterms.

A schematic functional expansion explains the bookkeeping. For

    S = S0 + Lambda^-1 S1 + Lambda^-2 S2 + Lambda^-3 S3 + ...,
    Phi -> Phi + Lambda^-2 f2,

the change through this order contains

    Delta S = Lambda^-2 (delta S0 / delta Phi) f2
            + Lambda^-3 (delta S1 / delta Phi) f2 + O(Lambda^-4).

Here the functional variation includes all fields and adjoints with their actual ordering; this schematic line is not a commuting rule for fermion fields. Removing a dimension-six redundant term against S0 leaves a possible dimension-seven term from S1. A complete power-counting ledger must track EFT order and loop order independently. The expansion is a standard field-substitution identity, not newly proposed mathematics.

**Use here:** retain the normalization/source changes acting on any present lower-dimensional operators. The optional dimension-six-normalization times dimension-five benchmark in the previous checkpoint is an example of this bookkeeping. This paper is **SMEFT without a dynamical sterile field**; its coefficient formulas cannot be copied into the present retained-sterile theory. Threshold matching would be a separate calculation if N were integrated out.

## Bare reduction and RG projection are different steps

Di Zhang, [Two-loop Renormalization Group Equations in the nuSMEFT, 2504.00792v2](https://arxiv.org/html/2504.00792v2), **5 May 2025**, is the previously identified source benchmark, rechecked rather than presented as a new discovery.

[Section 3, Eqs. (13)-(15)](https://arxiv.org/html/2504.00792v2#S3) explicitly shifts the sterile and lepton fields with mass, derivative and composite Higgs/lepton sources. It states that the resulting reduction relations are for **bare couplings**; the renormalized coupling and counterterm maps must be obtained by expansion at the chosen loop order. [Section 6](https://arxiv.org/html/2504.00792v2#S6), particularly the discussion following Eq. (51), demonstrates that setting redundant coefficients to zero before differentiating the renormalization constants can spoil pole-consistency relations. Differentiate while retaining that dependence, then project. Flavor-dependent Yukawa RG functions also depend on allowed field rotations, whereas flavor invariants remove this coordinate ambiguity.

**Use here:** keep an explicit sourced bare coefficient map, perform canonical normalization, expand the map consistently in poles and loop order, and only then extract the RG equations in a declared flavor convention. This source covers nuSMEFT through **dimension five**, at two loops; it is a workflow benchmark, not the selected dimension-seven result.

## A lower-dimensional running benchmark and its limit

Marco Ardu and Xabier Marcano, [Completing the one-loop nuSMEFT Renormalization Group Evolution, 2407.16751](https://arxiv.org/html/2407.16751), retrieved body **v2, 5 August 2024**; the index update is 6 August 2024.

[Table 1](https://arxiv.org/html/2407.16751#S2.T1) specifies dimension-six sterile operators and Hermitian-conjugation counting. [Section 3](https://arxiv.org/html/2407.16751#S3), Eqs. (5)-(7), describes the Yukawa-dependent one-loop matrix and the projection of redundant counterterms into that physical basis. The calculation includes mixed gauge/Yukawa and Higgs/Yukawa contributions; the authors combine it with existing gauge contributions at dimension six.

**Use here:** a future comparison for the induced dimension-six sector, after converting conventions and flavor order. Its one-loop linear EoM projection is not a replacement for the mixed-order bookkeeping above, and its matrix does not determine our dimension-seven running. No individual entry from that matrix has been numerically or symbolically reproduced here.

## Result of this screen

The actionable advance from this literature work is a stricter reduction protocol: preserve the interaction terms generated by identities; map them to corrected sterile bases; keep induced lower-order products; distinguish a bare field substitution from RG differentiation; and check a final flavor-invariant quantity. The explicit packet algebra and its numerical or symbolic verification are separate artifacts in this checkpoint.

No paper found in this bounded screen establishes the selected packet's full physical beta function, a finite photon amplitude, a unique non-Abelian transverse completion or a unification result. The external-gauge hierarchy already checked in candidate.6 does not close the additional external-matter sectors. Finite matching would also need a declared dimensional/evanescent convention beyond four-dimensional identities used for selected simple poles.
