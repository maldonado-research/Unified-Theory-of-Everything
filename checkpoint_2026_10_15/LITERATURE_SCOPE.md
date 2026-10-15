# Bounded primary-source audit for the hypercharge-leg checkpoint

Research cutoff and access: **15 October 2026**. This is a focused search for
sterile-neutrino EFT renormalization, dimension-seven operator mixing, and
one-loop gauge matching. It is not an exhaustive web survey or a novelty
assessment. Earlier requests returned proxy HTTP 403; later normal HTTPS
requests to arXiv and INSPIRE returned HTTP 200. The earlier no-fresh-access
assessment is superseded by the successful retrievals described below.

## Outcome

Complete dimension-seven **SMEFT** running and numerical solvers already
exist. Sterile-neutrino dimension-five/six running and implementations also
exist. These are different scopes: a calculation with a surviving sterile
external field cannot be identified directly with a pure-SMEFT anomalous
dimension after sterile removal.

The inspected material did not provide a coefficient-level reproduction of
this checkpoint's restricted off-shell `C_(7,seed) Y_nu^* g′` four-point UV
residue. No full convention/basis translation was completed, so the search
establishes neither novelty nor independent literature validation of that
residue. The defensible project result remains the calculated local
polynomial and its computational checks, with its stated restricted scope.

## Successful searches and coverage

All queries below returned HTTP 200 on 15 October 2026. Counts are the
services' reported totals; screened counts are returned records, with overlap
between searches. Screening used titles and relevant abstracts, not every
paper's full text. Broad searches returned unrelated work and are not a
completeness certificate.

| Service and exact query | Total | Screened |
|---|---:|---:|
| arXiv: `all:sterile AND all:"dimension-seven"`, latest submitted first | 3 | 3 |
| INSPIRE: `sterile neutrino dimension seven`, most recent first | 6 | 5 |
| arXiv: `(all:sterile OR all:"right-handed" OR all:nuSMEFT OR all:NRSMEFT) AND (all:renormalization OR all:anomalous)`, latest submitted first | 568 | 40 |
| INSPIRE: `("dimension-seven" OR "dimension seven" OR "dimension-7") AND ("one-loop" OR renormalization)`, most recent first | 97 | 40 |
| INSPIRE: `(nuSMEFT OR NRSMEFT OR "sterile neutrino") AND ("renormalization group" OR "anomalous dimension")`, most recent first | 24 | 24 |
| INSPIRE: `(title:renormalization OR title:renormalisation OR title:anomalous) AND (title:neutrino OR title:sterile)`, most recent first | 208 | 50 |

A separate successful arXiv ID lookup returned all five foundational records:
`2405.18017v3`, `2306.03008v2`, `2310.11055v2`, `1703.04415v3`, and
`1612.04527v1`. This revalidated their titles and latest version timestamps
returned by the API. Four additional full HTML texts were retrieved
successfully and selected sections inspected, as specified below. Full
responses and equation extractions remain in a local cache outside the
public repository. No downloaded paper or private archive excerpt is a
publication artifact.

## Calibration checked against fresh primary equations

**Di Zhang, *Threshold Effects on the Massless Neutrino in the Canonical
Seesaw Mechanism*, [arXiv:2405.18017v3](https://arxiv.org/abs/2405.18017v3)**,
latest API timestamp 1 April 2025. Fresh HTML Appendix A,
[Eqs. (35)–(38)](https://arxiv.org/html/2405.18017v3#A1.E35), confirms the
`-3/2` dimension-five momentum pole and the redundant derivative counterterm.
Its sourced EOM produces **both** a Yukawa descendant and a symmetric Weinberg
descendant, including both flavor orientations. Its derivative operator has
no explicit `i`; action normalization and derivative conventions must be
translated before comparing coefficients. The dimension-five `3/2` does not
supply the derivative-insertion result or the transverse four-point term.

The other four foundational records were freshly checked at metadata level
here; their equation-level discussion is inherited from the public
[1 October ledger](../checkpoint_2026_10_01/CURRENT_LITERATURE.md).
Elgaard-Clausen and Trott, *On expansions in neutrino effective field theory*,
[1703.04415v3](https://arxiv.org/abs/1703.04415v3), addresses retained-field
sources and sequential tree matching. Liao and Ma,
[1612.04527v1](https://arxiv.org/abs/1612.04527v1), supplies the corrected
sterile-neutrino operator basis. Zhang's
[2306.03008v2](https://arxiv.org/abs/2306.03008v2) and
[2310.11055v2](https://arxiv.org/abs/2310.11055v2) address pure-SMEFT
Green/physical reduction and dimension-seven running. None of those metadata
checks constitutes an independent reproduction of a loop calculation.

## Relevant returned papers and actual reading

**Yi Liao, Xiao-Dong Ma, Hao-Lin Wang and Xiang Zhao, *RGE solver for the
complete dim-7 SMEFT interactions and its application to 0νββ decay*,
[arXiv:2505.06499v2](https://arxiv.org/html/2505.06499v2), first submitted
9 May 2025; [JHEP 08 (2025) 138](https://doi.org/10.1007/JHEP08(2025)138).**
Fresh full HTML inspected: introduction, Section 2, and Appendix A's dipole
RGEs. Section 2 explicitly notes a corrected charged-lepton Yukawa mixing
term into `O_LHB,LHW` relative to earlier work; numerical comparisons should
use corrected expressions. Appendix A,
[Eqs. (15)–(16)](https://arxiv.org/html/2505.06499v2#A1.EGx14), gives
`C_LHW` and `C_LHB` running, including mixing from derivative operators.
This is a useful future check **after** matching to pure SMEFT; it does not
contain the retained sterile external leg of the present four-point
function. Section 2 also states that multiple effective-operator insertions
are neglected in the implementation. No solver code was run here.

**Jason Aebischer, Tejhas Kapoor and Jacky Kumar, *A package for
renormalization group running in the SMEFT with sterile neutrinos*,
[arXiv:2411.07220v1](https://arxiv.org/html/2411.07220v1), first submitted
11 November 2024; [EPJC publication](https://doi.org/10.1140/epjc/s10052-025-14127-8).**
Fresh full HTML inspected: Sections 2.1, 3, and conclusion. Its inventory is
dimension five/six, and its `wilson` implementation is therefore not a
ready-made solution of the present dimension-seven problem. The inspected
version explicitly omits Majorana-mass-dependent dipole-to-Weinberg mixing
and notes an upcoming inclusion of missing Yukawa terms. Those version
limits preclude treating a package name as evidence of complete running at
the requested order. Its reference to Ardu and Marcano,
[2407.16751](https://arxiv.org/abs/2407.16751), is a further dimension-six
renormalization lead; that paper was not separately inspected here.

**Lukáš Gráf, Chandan Hati, Ana Martín-Galán and Oliver Scholer, *Flavorful
Lepton-Number-Violating SMEFT: Loop Mixing, 0νββ and Rare Meson Decays*,
[arXiv:2609.28382v1](https://arxiv.org/html/2609.28382v1), first submitted
23 September 2026.** This is the newest directly relevant dimension-seven
mixing result returned by the targeted INSPIRE query. Fresh full HTML
inspected: Sections 1–2 and 8, plus the stated low-energy-input limitation
in Section 5. It studies pure-SMEFT flavor mixing and complementary nuclear
and meson observables, with a scalar-leptoquark UV example. It does not
supply a surviving-sterile seed coefficient. Its single-operator/flavor
assumptions and hadronic inputs cannot be transferred to a correlated
coefficient packet without a new matching/observable analysis.

**Di Zhang, *Renormalization Group Equations for the Dimension-7 SMEFT
Operators*, [arXiv:2409.02622v1](https://arxiv.org/html/2409.02622v1), first
submitted 4 September 2024; [conference contribution](https://doi.org/10.22323/1.476.0776).**
Fresh full HTML inspected: Sections 1–2 and summary. This summarizes the
established Green-basis/physical-basis strategy and pure-SMEFT one-loop
running through order Lambda^-3. It reinforces the need for lower-dimensional
EOM corrections; it is not an independent retained-sterile calculation.

Additional returned leads checked at record/abstract level only include
Chala and Titov, *One-loop matching in the SMEFT extended with a sterile
neutrino*, [2001.07732](https://arxiv.org/abs/2001.07732), and de Vries et al.,
*Neutrinoless Double-Beta Decays from Operator Mixing*,
[2608.07657](https://arxiv.org/abs/2608.07657). These do not supply an audited
coefficient comparison in this session.

## Consequence for the next calculation

A longitudinal Ward contraction leaves transverse local terms undetermined.
The new four-point residue should next be translated into an independent
retained-sterile operator basis, including its correlated sourced descendants.
Only after the remaining counterterms, field/parameter renormalization and
threshold matching are included can physical running or an observable be
tested. Existing dipole/mixing literature constrains that program; finding a
nonzero transverse residue is not by itself discovery of a new operator,
physical law, unified theory, or experimentally established effect.

No search or primary-text request in this bounded recovery batch failed.
Earlier HTTP 403 failures remain part of the access history. Neither
successful retrieval nor an absence of a matching result in the screened
subset establishes comprehensive literature coverage or novelty.
