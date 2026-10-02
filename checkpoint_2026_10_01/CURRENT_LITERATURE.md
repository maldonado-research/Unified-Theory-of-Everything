# Focused primary-source audit for the October 2026 operator checkpoint

Research cutoff: **1 October 2026**. Access and rechecking: **1 October 2026 Pacific**. This is a bounded literature review of sourced equations of motion (EOM), sterile-neutrino effective theories, dimension-seven bases, sequential matching and operator mixing. It is not an exhaustive web search, a novelty determination, or an independent reproduction of published loop calculations.

## Outcome and consequence for this checkpoint

The primary literature supports the methodological choices made in the checkpoint: keep the sterile and light-field interaction sources; distinguish an off-shell local representative from a complete independent Green's basis; declare the entire effective action and power counting before reduction; and carry the resulting correlated coefficients to a physical matching combination. None of the sources inspected establishes the checkpoint's dimension-seven seed coefficient as a physical beta function, or establishes a unified theory.

A particularly direct benchmark is **arXiv:2405.18017v3, Appendix A, Eqs. (35)–(38)**. Its dimension-five Weinberg-insertion loop needs a redundant derivative counterterm. Applying the sourced sterile EOM produces **both** a Yukawa-shaped term and a Weinberg-shaped contact. This verifies that the same structural issue addressed here already occurs in the published lower-dimensional calculation. It does not fix the higher-derivative seed normalization.

No contradiction with the checkpoint's stated normalizations was established by this review. No coefficient-level translation to a full published dimension-seven sterile basis was completed. In particular, the reference defines its derivative operator without the explicit `i` used in the present checkpoint; its redundant coefficient therefore contains an `i`. These conventions must be translated before comparing signs. Its `3/2` in the dimension-five pole cannot be substituted for the derivative seed's `5/6` projection.

## Established sources rechecked from primary records and relevant text

### 1. Sequential seesaw matching and EOM sources

Gitte Elgaard-Clausen and Michael Trott, **On expansions in neutrino effective field theory**. First submitted 13 March 2017; version inspected **v3, 17 June 2018**. [Record](https://arxiv.org/abs/1703.04415), [versioned full text](https://arxiv.org/html/1703.04415v3), [journal](https://doi.org/10.1007/JHEP11(2017)088).

Checked: metadata, abstract, Section 2.1 and relevant matching passages in Sections 3–4. Section 2.1 explicitly states that sequential decoupling must use EOM containing the heavy neutrinos that remain in the spectrum. Matching is at tree level through dimension seven. This is a direct reason to retain remaining sterile sources between thresholds and a benchmark for the eventual tree matching packet. It supplies neither the missing one-loop gauge completion nor the current seed's finite threshold.

### 2. A corrected sterile-neutrino operator basis

Yi Liao and Xiao-Dong Ma, **Operators up to Dimension Seven in Standard Model Effective Field Theory Extended with Sterile Neutrinos**. **v1, 14 December 2016**, the sole version listed. [Record](https://arxiv.org/abs/1612.04527), [full text](https://arxiv.org/html/1612.04527v1), [journal](https://doi.org/10.1103/PhysRevD.96.015012).

Checked: metadata, abstract and relevant basis-construction discussion. The paper uses EOM, integration by parts, Fierz identities and group identities and checks counting by Hilbert series. Relative to earlier sterile-neutrino literature it removes redundant operators and adds missing dimension-seven operators. Its operator counts depend on flavor and Hermitian-conjugation conventions.

Implication: start a completeness audit with this corrected sterile-neutrino basis, then build whatever redundant off-shell set the calculation actually requires. The six checkpoint monomials have not been shown independent under integration by parts and therefore should remain a **local representative**, as currently labeled. A gauge-field-strength structure invisible to a three-point restriction is compatible with this need for more external-leg sectors; the present review has not mapped that structure to a specific minimal-basis coefficient.

For comparison, the earlier Subhaditya Bhattacharya and José Wudka paper, **Dimension Seven Operators in Standard Model with Right handed Neutrinos**, [arXiv:1505.05264v2](https://arxiv.org/abs/1505.05264v2), first submitted 20 May 2015 and v2 dated 21 September 2016, was checked at metadata/abstract level only. It should not be treated as an independently complete corrected basis without the later relations and relevant erratum.

### 3. Dimension-five threshold running with sterile fields retained

Di Zhang, **Threshold Effects on the Massless Neutrino in the Canonical Seesaw Mechanism**. First submitted 28 May 2024; **v3, 1 April 2025**. [Record](https://arxiv.org/abs/2405.18017), [versioned full text](https://arxiv.org/html/2405.18017v3).

Checked: metadata, abstract, dimension-five Green's basis in Eq. (4), sequential-decoupling discussion and Appendix A, Eqs. (35)–(40). The appendix explicitly gives the redundant derivative counterterm and the associated Yukawa and Weinberg contributions after sourced EOM reduction. The symmetric contact involves both flavor orientations; suppressing one orientation loses part of the coefficient.

The paper's one-loop rank result concerns the minimal/canonical type-I seesaw and the stated running/matching assumptions. It does not show that arbitrary extra operators, finite thresholds or higher loops preserve a massless state. Use the rank statement within that domain, and use the appendix as a convention-sensitive lower-order check, not as a derivation of the present dimension-seven result.

### 4. SMEFT Green's and physical bases through order Λ^-3

Di Zhang, **Renormalization Group Equations for the SMEFT Operators up to Dimension Seven**. First submitted 5 June 2023; **v2, 16 October 2023**. [Record](https://arxiv.org/abs/2306.03008), [full text](https://arxiv.org/html/2306.03008v2).

Checked: metadata, abstract and the discussion of basis reduction and EOM in Section 2. The paper explicitly includes dimension-five corrections to EOM where redundant dimension-six structures feed the dimension-seven result. It distinguishes the Green's basis, independent under integration by parts and algebraic identities but redundant under field redefinitions, from a physical basis.

Implication: the October checkpoint's renormalizable source examples are useful identities but do not yet constitute reduction of a complete sequential EFT containing lower effective operators. A declared power counting must determine which lower-dimensional source corrections contribute. This paper is SMEFT after the sterile fields are removed; its basis cannot simply replace a theory that still propagates sterile fields.

Di Zhang, **Revisiting Renormalization Group Equations of the SMEFT Dimension-Seven Operators**. First submitted 17 October 2023; **v2, 31 January 2024**. [Record](https://arxiv.org/abs/2310.11055), [full text](https://arxiv.org/html/2310.11055v2).

Checked: metadata, abstract and the calculation strategy. It obtains the same-dimension running through off-shell diagrams, Green-basis counterterms, reduction and physical-basis running, with partial checks using other tools and gauges. This is an eventual implementation benchmark; no ancillary code was executed and no anomalous-dimension entry was independently reproduced here.

### 5. Limits of EOM substitution and finite field changes

Juan Carlos Criado and Manuel Pérez-Victoria, **Field redefinitions in effective theories at higher orders**. First submitted 23 November 2018; **v2, 20 September 2019**. [Record](https://arxiv.org/abs/1811.09413), [full text](https://arxiv.org/html/1811.09413v2).

Checked: metadata, abstract and relevant passages in Sections 3–5 on the difference between EOM substitution and higher-order field transformations. EOM substitution at first order does not license dropping the higher powers generated by a finite change of variables. This supports retaining the quadratic contact in the checkpoint's exact triangular-change/Schur test; it does not elevate that test to a full EFT loop calculation.

Abdurrahman Barzinji, Michael Trott and Anagha Vasudevan, **Equations of Motion for the Standard Model Effective Field Theory: Theory and Applications**. First submitted 17 June 2018; **v2, 7 December 2018**. [Record](https://arxiv.org/abs/1806.06354), [full text](https://arxiv.org/html/1806.06354v2).

Checked: metadata, abstract and relevant discussion/conclusion. Lower-dimensional contact operators modify EOM and can affect higher-dimensional matching coefficients. This provides a second methodological source for the reference-action requirement; it is not a substitute for computing the full field redefinition at the required perturbative order.

## Previously cited 2026 leads revalidated

These are recent relative to the older foundations, but were already present in the September checkpoint. They are not newly discovered October results.

- **Carla Biggio, Marta Fuentes Zamoro, Xu Li, Luca Merlo and Luca Ottonello**, *How to Identify a Majoron: Effective Field Theories of Spontaneous Lepton Number Breaking*, **arXiv:2608.11522v1, 12 August 2026**. [Record](https://arxiv.org/abs/2608.11522), [full text](https://arxiv.org/html/2608.11522v1). Checked metadata/abstract, matching-order discussion and concluding scope statement. It compares removal of a heavy radial mode and seesaw fields in two orders. The conclusion explicitly states that matching was performed at tree level and that one-loop work is needed for consistent running. The model assumes a complex symmetry-breaking scalar and a surviving angular Majoron. Its correlated observables must not be imported into a real-scalar or explicit-Majorana-mass model without deriving the different matching.
- **J. de Vries, S. Fajfer, L. P. S. Leal, O. Sumensari and R. Zukanovich Funchal**, *Neutrinoless Double-Beta Decays from Operator Mixing*, **arXiv:2608.07657v1, 7 August 2026**. [Record](https://arxiv.org/abs/2608.07657), [full text](https://arxiv.org/html/2608.07657v1). Checked metadata/abstract and the discussion of logarithmic evolution in Section III and conclusions. A coefficient with no direct decay contribution can contribute through mixing; iterated one-loop mixing can make a double logarithm the first nonzero contribution in a flavor channel. This is distinct from having computed the full two-loop anomalous dimension. A viable observable test requires the correlated coefficient vector, electroweak/hadronic matching and stated nuclear inputs, not a universal single-operator bound.
- **Weiyi Deng, Chengcheng Han, Wuzhou Yin and Tong Ju**, *Right-Handed Neutrino Production by an Axion-like Inflaton: Implications for Leptogenesis*, **arXiv:2607.13592v1, 15 July 2026**. [Record](https://arxiv.org/abs/2607.13592), [full text](https://arxiv.org/html/2607.13592v1). Metadata and abstract were rechecked; the full HTML was retrieved but detailed equations were not re-audited this session. Derivative axion production, Pauli blocking and repeated production are a separate cosmological application. They do not repair a missing local matching coefficient and are not adopted as this project's underlying model.

## Search after the September checkpoint, including failures

A successful [arXiv title search for neutrino](https://arxiv.org/search/?query=neutrino&searchtype=title&abstracts=show&order=-announced_date_first&size=200) was screened for matching, sterile fields, dimension-seven bases and operator mixing. The [recent hep-ph list](https://arxiv.org/list/hep-ph/recent?skip=0&show=2000) returned 247 entries announced 28 September–2 October 2026; entries after the 1 October cutoff were excluded from the assessment. A title search misses relevant work without the search word in its title; this is a material coverage limitation.

One subsequent result inspected at record/abstract level was **Yan Shao and Zhen-hua Zhao**, *Purely flavon driven leptogenesis for exactly degenerate Dirac or Majorana neutrino mass matrix in the type-I seesaw model*, [arXiv:2609.35299](https://arxiv.org/abs/2609.35299), first submitted **28 September 2026**. It adds flavon interactions or heavy vectorlike mediators to generate decay asymmetries in specially degenerate textures. Those additional interactions are not in the present reference action. Its abstract does not provide the missing dimension-seven Green-basis reduction; no extension of the project is justified from it alone.

No subsequent result that resolves the seed's full gauge completion or physical beta coefficient was verified in the successful searches. This does **not** establish that no such publication exists.

Failed requests were recorded rather than interpreted as an absence of papers: the export.arxiv.org API ID-list request timed out; arxiv.org/api queries for sterile/dimension, neutrino/Green and date-restricted seesaw returned HTTP 429; several multiword search pages and all-field searches for dimension-seven, nuSMEFT, sterile and field redefinitions timed out; the September monthly hep-ph listing returned HTTP 404. The successful abstract/full-text routes above were used instead. This review used primary arXiv pages; it did not use Google Drive, Consensus or Firecrawl and claims no coverage from them.

## Concrete next benchmark

Implement the dimension-five calculation and its two sourced descendants from arXiv:2405.18017v3 in the project's own explicit action conventions as a calibration problem. Then declare the retained-field dimension-seven EFT, identify the relevant corrected basis and all lower-order insertions, and match the additional gauge and multi-field sectors. Only after complete counterterm and finite-threshold bookkeeping should a physical running or matching claim be tested. Literature agreement on the method does not certify the desired coefficient.

Only newly authored notes and links belong in any public update. Retrieved paper HTML and extracted text remain in a separate local literature cache and are not publication artifacts.
