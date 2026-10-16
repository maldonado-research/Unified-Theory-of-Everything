# Bounded literature follow-up: field redefinitions and canonical normalization

Accessed 2026-10-15. This is a targeted methodological screen, not an exhaustive web review, a novelty assessment, or evidence for unification. Machine-readable query and result provenance is in `LITERATURE_PROVENANCE.json`.

## Coverage and retrieval limits

Four Firecrawl paper-index queries requested 12 results each within the provider's inclusive created/updated window 2025-01-01 to 2026-10-15. They returned **12, 12, 4 and 7 ranked records: 35 appearances, 27 distinct canonical papers and 8 duplicate appearances**. Four papers received metadata inspection and in-body passage review. Five passage requests returned 18 passages in total.

Four Firecrawl arXiv-abstract scrape requests failed because the connected account reported insufficient credits. A fallback using ordinary managed HTTPS then successfully retrieved all four selected primary arXiv HTML documents with HTTP 200. The local bytes were checked against SHA256 receipts, and the specific sections below were read to corroborate the methodological claims. No complete paper or supplement was exhaustively reviewed. Full HTML stays in scratch storage and is not part of the publication artifacts.

The date filter is creation **or update**, so a 2024 identifier is retained transparently in Q4; it is not presented as first posted in 2025.

## Primary-source receipt cross-check

All four saved HTML byte hashes match the retrieval receipts. Version dates below are read from each served HTML header, independently of index update dates.

| Primary URL | UTC retrieval | Observed version | Bytes | SHA256 |
| --- | --- | --- | --- | --- |
| [2504.00792v2](https://arxiv.org/html/2504.00792v2) | 2026-10-15T21:01:35.160444+00:00 | 2504.00792v2, 05 May 2025 | 1506975 | `27a02c94440ffb3936bdadb89efd50185353b3b85d37e05d9ff26d70378a02d2` |
| [2512.15866](https://arxiv.org/html/2512.15866) | 2026-10-15T21:01:35.161046+00:00 | 2512.15866v1, 17 Dec 2025 | 1378254 | `215a651053c6a95873bf3644396ce4af7a0db2056811a31c186a86dad8b5864a` |
| [2512.16890](https://arxiv.org/html/2512.16890) | 2026-10-15T21:01:35.806447+00:00 | 2512.16890v1, 18 Dec 2025 | 1380293 | `47d272ade8c65f712382eb1a0ecac96c397b0e4a30fb86055cb6ddc3ff37d466` |
| [2504.18537](https://arxiv.org/html/2504.18537) | 2026-10-15T21:01:35.913841+00:00 | 2504.18537v1, 25 Apr 2025 | 1003322 | `b22f8d43c18cbef65ab7cbdc81a703c4e7be13c527bceb054eb2613c5b451814` |

The receipts establish the exact source bytes inspected; they do not reproduce the papers' results or certify that every result applies to this model.

## Exact search protocol

- **Q1**, k=12; returned 12: Sterile neutrino effective field theory nuSMEFT one loop renormalization group equations field redefinitions canonical normalization dimension six
- **Q2**, k=12; returned 12: Effective field theory derivative operators interacting equations of motion field redefinitions canonical normalization anomalous dimensions physical running
- **Q3**, k=12; returned 4: Right handed neutrino effective field theory dipole photon gauge invariant one loop mixing hypercharge weak isospin
- **Q4**, k=12; returned 7: Field redefinitions effective theories higher orders equations of motion renormalization derivative fermion source operators

All queries used from=2025-01-01, to=2026-10-15 and date accessed 2026-10-15. The tools supplied canonical paper IDs but no per-query search UUID. Per-query ranks, scores and IDs are preserved in JSON.

## Selected methodological benchmarks

### 1. Di Zhang, Two-loop Renormalization Group Equations in the νSMEFT

[arXiv:2504.00792v2](https://arxiv.org/abs/2504.00792v2); retrieved body header states 05 May 2025. Metadata reports first posting 1 April 2025 and index update 7 May 2025; these are distinct date fields. Scope: νSMEFT through dimension five at two loops.

[Section 3](https://arxiv.org/html/2504.00792v2#S3), Eqs. (13)-(14), gives explicit right-handed-neutrino and lepton-doublet shifts. The sterile shift contains a mass-dependent linear term, a derivative of the charge-conjugated sterile field, and a Higgs/lepton source. The doublet shifts by a Higgs/sterile composite. The induced relations affect the sterile Majorana mass, neutrino Yukawa, Weinberg coefficient and Higgs/sterile coefficient. Those are **bare-coupling** relations; renormalized coefficients and counterterms follow by order-by-order expansion.

[Section 6](https://arxiv.org/html/2504.00792v2#S6), especially Eq. (51) and the discussion through Eq. (57), warns that setting redundant coefficients to zero before differentiating renormalization constants loses contributions and spoils pole-consistency relations. Restore the full dependence, differentiate, then project onto vanishing redundant coefficients. Anti-Hermitian flavor rotations can change flavor-dependent Yukawa RG functions. Flavor invariants avoid that flavor-coordinate ambiguity.

**Use here:** the derivative/source field shifts and the order of RG projection are directly applicable as methodology. The paper does **not** supply the complete physical coefficient map or anomalous dimension for the present dimension-seven packet.

### 2. Fonseca, Olgoso and Santiago, Renormalization of general Effective Field Theories: Renormalization of fermionic operators

[arXiv:2512.15866](https://arxiv.org/abs/2512.15866): live HTML identifies **v1, 17 December 2025**; the paper index separately reports an update date of 19 December 2025.

[Section 3.1](https://arxiv.org/html/2512.15866#S3.SS1), Eqs. (13)-(15), follows a physical UV boundary condition into off-shell divergences with nontrivial kinetic terms and redundant operators, then canonical normalization and physical-basis reduction, then bare/renormalized coupling splitting. Some **bosonic** Green-basis divergences from [the companion, arXiv:2501.13185](https://arxiv.org/abs/2501.13185), are required even for fermionic physical beta functions.

The authors report computations with Matchmakereft and group-theory/model checks. This review did not execute those packages, obtain all supplements, or specialize their templates.

**Use here:** an independent workflow benchmark for source retention, field normalization and reduction. A fermion-only packet may not contain all inputs required for a physical beta function.

### 3. Aebischer, Bresciani and Selimovic, Anomalous Dimension of a General Effective Gauge Theory II: Fermionic Sector

[arXiv:2512.16890](https://arxiv.org/abs/2512.16890): live HTML identifies **v1, 18 December 2025**. The index separately reports an update date of 16 September 2026 and [JHEP09(2026)159](https://doi.org/10.1007/JHEP09(2026)159); the DOI page itself was not retrieved. The index update date does not establish a newer inspected arXiv version.

The introduction describes a physical basis through dimension six for arbitrary gauge groups and scalar/fermion/vector content, retaining all Lagrangian parameters through order 1/Λ² and mixing across dimensions. [Section 3](https://arxiv.org/html/2512.16890#S3), Eq. (41), extracts ultraviolet anomalous dimensions from on-shell form factors **after subtracting infrared anomalous dimensions**. Eq. (42) describes the multi-insertion expansion. Appendix A.1 covers collinear terms, while soft logarithms are removed separately.

**Use here:** a possible independent on-shell cross-check after conventions and representations are mapped. A vanishing selected on-shell/neutral UV structure cannot be read as a finite photon amplitude, and an on-shell loop pole must be separated from infrared singularities.

### 4. Assi, Helset, Pagès and Shen, Renormalizing Two-Fermion Operators in the SMEFT via Supergeometry

[arXiv:2504.18537](https://arxiv.org/abs/2504.18537): live HTML identifies **v1, 25 April 2025**; the paper index separately reports an update date of 28 April 2025.

The formalism includes mixed boson-fermion one-loop graphs and manifests covariance under **nonderivative** field redefinitions. [Section 5](https://arxiv.org/html/2504.18537#S5) reduces Green-basis outputs using explicit operator relations and Matchete as cross-checks. [Section 6](https://arxiv.org/html/2504.18537#S6) explicitly identifies higher-derivative operators and derivative field redefinitions as a remaining challenge for the standard geometric picture.

“Supergeometry” refers to commuting and anticommuting field-space coordinates; it does not require a physically supersymmetric theory.

**Use here:** organize mixed-loop structure and independently check reductions where its hypotheses apply. It does not automatically eliminate the retained derivative-source packet or establish a geometric unification of forces.

## Consequences for the next calculation

1. Expand explicit, gauge-covariant field shifts through the full interacting leading-order action, including Yukawa sources and Majorana masses. Free equations of motion alone omit induced terms.
2. Carry canonical normalization into every coupling and operator coefficient and state the EFT/loop order retained. Use multiple flavors or complex parameters to detect mistakes hidden by commuting one-flavor limits.
3. Retain redundant-coefficient dependence through RG differentiation and only then project. Record anti-Hermitian flavor freedom and use flavor invariants for any claim of physical running.
4. A complete physical coefficient needs all source sectors required by the reduction, including bosonic counterterms where applicable. Specialize a published general template before transferring numerical coefficients.
5. For an on-shell check, subtract infrared terms and include external-state normalization and any matching/input conventions needed by the observable.
6. Finite matching and an actual observable remain necessary before a finite photon amplitude, physical dipole prediction or experimental bound can be claimed.

These are established consistency requirements and implementation targets, not a new mathematical theory or a demonstrated breakthrough.

## Complete abstract-screen inventory

The four papers above received indexed passage review plus targeted corroboration in primary HTML. Other entries are screened from indexed titles/abstracts; their reported results have not been independently reproduced. No abstracts, full paper passages or HTML copies are included here.

| Paper | Query | Relevance and disposition |
| --- | --- | --- |
| [arxiv:2504.00792](https://arxiv.org/abs/2504.00792) | Q1, Q2, Q3, Q4 | Two-loop nuSMEFT through dimension five; derivative redundant coefficients, canonical/field redefinition shifts and flavor dependence are directly relevant; not the complete physical coefficient map for the present dimension-seven packet. |
| [arxiv:2501.17134](https://arxiv.org/abs/2501.17134) | Q1 | General one-loop bosonic physical-operator RGEs; useful only after matching conventions and field content; no fermionic completion established from this abstract. |
| [arxiv:2604.21492](https://arxiv.org/abs/2604.21492) | Q1 | Finite-temperature cosmology applications; not needed to establish current operator reduction or physical running. |
| [arxiv:2510.24843](https://arxiv.org/abs/2510.24843) | Q1 | Dimension-six nuSMEFT leptogenesis application; relevance comes after an independently normalized physical operator basis. |
| [arxiv:2601.19974](https://arxiv.org/abs/2601.19974) | Q1 | Two-loop dimension-six SMEFT and evanescent scheme treatment; ordinary SMEFT field content is not the retained-sterile theory. |
| [arxiv:2601.15901](https://arxiv.org/abs/2601.15901) | Q1 | Threshold/input-scheme/tadpole corrections warn that measured broken-phase parameters require finite matching; no immediate sterile-basis result. |
| [arxiv:2510.08682](https://arxiv.org/abs/2510.08682) | Q1, Q2 | Baryon-number violation and NLL scheme cancellation; sector differs from neutral sterile-photon packet. |
| [arxiv:2602.05461](https://arxiv.org/abs/2602.05461) | Q1 | Broad dimension-eight SMEFT thesis; not a direct field-redefinition coefficient benchmark. |
| [arxiv:2512.08827](https://arxiv.org/abs/2512.08827) | Q1, Q2 | Two-loop general bosonic EFT; omits fermionic sector needed here. |
| [arxiv:2608.07657](https://arxiv.org/abs/2608.07657) | Q1 | Neutrinoless double-beta mixing effects; cannot transfer bounds without operator/flavor map. |
| [arxiv:2512.15866](https://arxiv.org/abs/2512.15866) | Q1, Q4 | General one-loop fermionic EFT through dimension six; passages explicitly require nontrivial kinetic terms, full Green basis, bosonic poles, canonical normalization, and physical-basis shifts. |
| [arxiv:2502.14303](https://arxiv.org/abs/2502.14303) | Q1 | Sterile EFT but baryon-number violating nucleon decays; no direct selected packet calculation. |
| [arxiv:2511.15609](https://arxiv.org/abs/2511.15609) | Q2 | Truncation and field-redefinition uncertainty; reinforces power-counting caution, screened only at abstract level. |
| [arxiv:2512.04064](https://arxiv.org/abs/2512.04064) | Q2 | Alternate two-loop extraction via compactified 5D matching; methodological candidate only, not applied. |
| [arxiv:2511.16740](https://arxiv.org/abs/2511.16740) | Q2 | Scalar phi-four multiloop EFT/CFT; not a fermion derivative-source reduction. |
| [arxiv:2502.14030](https://arxiv.org/abs/2502.14030) | Q2 | Bosonic counterpart of selected on-shell fermionic template; may be necessary to assemble a complete map. |
| [arxiv:2512.16890](https://arxiv.org/abs/2512.16890) | Q2, Q4 | On-shell physical one-loop arbitrary gauge EFT through dimension six; passages identify soft/collinear IR subtraction, single/double insertion mixing and cross-dimension templates. |
| [arxiv:2608.19316](https://arxiv.org/abs/2608.19316) | Q2, Q4 | Two-loop current-current four-fermion templates; only abstract screened, no universality claim independently checked. |
| [arxiv:2607.18166](https://arxiv.org/abs/2607.18166) | Q2 | Noninjective scalar field redefinitions; distinct global issue from local perturbatively invertible EFT redefinitions. |
| [arxiv:2604.14741](https://arxiv.org/abs/2604.14741) | Q2 | OPE operator-mixing computation; abstract-level alternative, not applied. |
| [arxiv:2512.15927](https://arxiv.org/abs/2512.15927) | Q2 | Einstein-Maxwell shockwaves do not establish sterile-neutrino normalization or physical photon dipole. |
| [arxiv:2501.13150](https://arxiv.org/abs/2501.13150) | Q3 | Electromagnetic leptogenesis with heavy/light sterile dimension-five dipoles; field content and gauge-invariant operator differ from an SM-doublet dimension-six dipole. |
| [arxiv:2508.08903](https://arxiv.org/abs/2508.08903) | Q3 | W measurements constrain low-mass sterile interactions; requires independent observable matching first. |
| [arxiv:2609.28415](https://arxiv.org/abs/2609.28415) | Q3 | Dark U(1) kinetic mixing model; do not equate dark-photon mixing with a neutrino electromagnetic dipole. |
| [arxiv:2504.18537](https://arxiv.org/abs/2504.18537) | Q4 | Covariant mixed boson-fermion one-loop formulas and SMEFT dimension-eight examples; covariance explicitly nonderivative and higher derivatives remain a stated challenge. |
| [arxiv:2409.08984](https://arxiv.org/abs/2409.08984) | Q4 | Higher-derivative SUSY EFT order reduction; different field content. First posted 2024 despite provider date filter allowing created/updated hits; metadata update not independently verified. |
| [arxiv:2501.13185](https://arxiv.org/abs/2501.13185) | Q4 | Green/physical general EFT bases and reduction shifts through dimension six; selected fermionic paper states some of its bosonic poles are required. |

## Publication boundary

This screen supports a more careful calculation and identifies reusable methods. It supplies no completed anomalous dimension for the current packet, finite photon amplitude, new unification mechanism, global novelty claim, or empirical confirmation.

