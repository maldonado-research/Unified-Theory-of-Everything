# Literature comparison for derivative transport and field-content projection

Accessed 16 October 2026. This note checks operator classes, source retention, mass conventions and the limits of helicity arguments. It does not calculate the packet's UV coefficient, a physical anomalous dimension or a finite amplitude.

## Bounded coverage

Four Firecrawl paper-index queries requested k=10 each and returned **9, 5, 2 and 10 records: 26 ranked appearances, 24 distinct papers and two duplicate appearances**. Three queries used the inclusive created/updated window 2025-01-01 through 2026-10-16; a fourth deliberately widened the window to 2016 to recover sterile operator classifications. The query strings, ranks, identifiers and our screening summaries are preserved in `LITERATURE_PROVENANCE.json`.

Two newly selected papers received metadata inspections, passage requests returning four passages each, and fresh primary HTML retrievals, both HTTP 200. Four previously retrieved primary bodies were rechecked for the sections below, with their saved byte hashes verified. Only the specified sections were reviewed; complete papers, supplements and external software were not exhaustively checked. Other search hits received an abstract-level screen. Raw papers, verbatim abstracts, passage dumps and provider results remain in scratch storage; public provenance keeps metadata, hashes and our own summaries.

## The sterile dimension-seven derivative/strength class

Yi Liao and Xiao-Dong Ma, [Operators up to Dimension Seven in Standard Model Effective Field Theory Extended with Sterile Neutrinos, 1612.04527](https://arxiv.org/abs/1612.04527), inspected PDF **v1, 14 December 2016**.

Section 3, Eq. (17), gives four representatives with a sterile/lepton vector bilinear, a derivative on the Higgs, and a B or W strength or its dual. Eqs. (32)-(34) then show explicitly that moving between dual-strength representatives with a derivative on the lepton involves gamma identities **and an interacting EoM contribution**. Their example generates a Yukawa-weighted scalar/sterile tensor-strength operator. Thus derivative placement and a dual-strength identity cannot be treated independently of the interaction sources. The fact that some derivative-on-lepton representatives are redundant does not eliminate the scalar-derivative classes in Eq. (17).

Hao-Lin Li et al., [Operator Bases in Effective Field Theories with Sterile Neutrinos: d <= 9, 2105.09329](https://arxiv.org/abs/2105.09329), inspected PDF **v1, 19 May 2021**, provides a separately organized check. Section 4.3, Eqs. (4.60)-(4.63), lists the chiral-strength/lepton/sterile/Higgs-derivative types. Section 2.2 explains amplitude/operator correspondence; Eqs. (4.13)-(4.14) display chirality-dependent strength contractions with two-component sigma matrices.

**Use here:** the calculated image of `QVh` and `QEh` can be compared structurally with this derivative/strength sector. The packet's explicit projector

    X_+ = (X + i dual(X))/2

is defined and checked in the packet's own metric, epsilon, gamma and chirality conventions. We do not transfer an `F_L` or `F_R` name, a dual sign or a numerical coefficient from a source. A full coefficient map would additionally need `phi=tilde H`, weak-index orientation, Hermitian conjugation, Pauli versus `T=tau/2`, gauge-field sign conventions and derivative adjoints converted explicitly.

The literature supplies a target class, not the selected row-transport formula. That formula must differentiate every covariant factor before applying the declared lepton, sterile and gauge EoMs. Its current/source descendants and ordered residuals are part of the action identity. A boundary term or Bianchi identity must be declared separately from an EoM residual. Using free equations or deleting the gauge-current contribution would change the identity being tested.

## A zero-sterile slice is not a threshold EFT

Wouter Dekens, Jordy de Vries, Kaori Fuyuto, Emanuele Mereghetti and Guanghui Zhou, [Sterile neutrinos and neutrinoless double beta decay in effective field theory, 2002.07182](https://arxiv.org/html/2002.07182), newly retrieved body **v2, 14 June 2020**; the index update is 15 July 2020.

[Section 2, Eq. (1)](https://arxiv.org/html/2002.07182#S2.E1) explicitly retains dynamical right-handed sterile fields and a general complex symmetric mass matrix. Its operator tables select terms relevant to its double-beta-decay application and are not used here as a complete derivative-sector classification.

Actual elimination of heavy neutrinos is a separate operation in [Section 2.3.1, Eqs. (40)-(41)](https://arxiv.org/html/2002.07182#S2.SS3.SSS1): it solves a heavy-field equation using an inverse mass and a specified expansion in momentum over that mass. That is a threshold approximation, not the act of retaining monomials with zero explicit sterile fields.

The present `pi_0(P_red)` instead means the polynomial field-content slice with no sterile fields or their derivatives. The sterile field remains in the theory. No inverse mass expansion, path integration over N or new low-energy matching coefficient follows from taking this slice.

### Why projection and an EoM field map can fail to commute

The following is an elementary explanation, not a coefficient taken from a paper. A Yukawa interaction contains

    L_Y = -bar L Y_nu phi N_R + h.c.

A permitted first-order field shift `delta N_R=eta chi_R[L,H]`, with no N in the composite chi, generates

    pi_0(delta L_Y) = -eta bar L Y_nu phi chi_R + h.c.

even though `pi_0(L_Y)=0`. Applying the slice first therefore loses information needed for the subsequent variation. Interacting EoMs likewise connect monomials of different sterile degree. A zero-N projection is consequently representative dependent under such action-coordinate changes.

Juan Carlos Criado and Manuel Perez-Victoria, [Field redefinitions in effective theories at higher orders, 1811.09413](https://arxiv.org/html/1811.09413), rechecked body **v2, 20 September 2019**, supplies the general methodological foundation: [Section 2, Eq. (2)](https://arxiv.org/html/1811.09413#S2.E2) transforms sources as well as fields, and [Section 3, Eq. (10)](https://arxiv.org/html/1811.09413#S3.E10) expands the actual field substitution. The full action, rather than a prematurely restricted field-content slice, determines the induced terms at the retained order.

**Use here:** first evaluate the sources and complete the selected field map, then take and label the zero-N slice of that chosen representative. It is not an invariant projection of the EoM equivalence class. Its completeness is relative to the declared reference interactions: a neutrino-Yukawa-only reference and a reference with Standard Model spectator Yukawas are different inputs that must be identified explicitly.

## Complex symmetric masses require chirality-specific products

Dekens et al. [Section 2.2.1, Eqs. (16)-(18)](https://arxiv.org/html/2002.07182#S2.E17) diagonalizes its symmetric neutrino mass matrix by **unitary congruence**, `U^T M_nu U=m`, rather than unitary similarity. The field vector and dagger placements are explicitly fixed there; they should not be copied into a different left/right convention without conversion.

Di Zhang, [Two-loop Renormalization Group Equations in the nuSMEFT, 2504.00792v2](https://arxiv.org/html/2504.00792v2), rechecked body **5 May 2025**, gives another practical ordering benchmark. [Section 3, Eqs. (13)-(14)](https://arxiv.org/html/2504.00792v2#S3) contains sterile/lepton source field shifts and mass corrections involving separately ordered products such as `M M-dagger G` and `G M-dagger M`. These are bare reduction formulas through dimension five, not the present dimension-seven transport coefficients.

For the explicitly declared convention

    A N_L = M N_R + S_R,
    A N_R = M-dagger N_L + S_L,
    A = i slash(partial),        M^T=M,        partial M=0,

direct iteration gives

    A^2 N_L = M M-dagger N_L + M S_L + A S_R.

The product acting on N_R is instead `M-dagger M`. The two agree in a real diagonal convention but cannot generally be interchanged. For example, with `N_R=V N_R'`, `N_L=V* N_L'`, one has

    M' = V^T M V,                F' = F V*,
    M' M'-dagger = V^T (M M-dagger) V*,
    F' M' M'-dagger = F (M M-dagger) V*.

These are covariance checks for the declared constant-mass convention. They explain what the ordered `R/T/U` transport must respect; they do not determine that transport's coefficients. They require neither mass inversion nor a nonzero mass eigenvalue and therefore should not be replaced by a heavy-sterile expansion.

## A current helicity theorem is a check with additional premises

Luigi C. Bresciani and Nudzeim Selimovic, [All-Order Helicity Selection Rules in Effective Field Theories, 2609.20914](https://arxiv.org/html/2609.20914), newly retrieved body **v1, 17 September 2026**; index update 21 September 2026.

The setup preceding Eqs. (1)-(2) uses a **four-dimensional massless theory**, a common operator dimension and an **irreducible physical basis**. Operator weights are obtained from minimal on-shell configurations, including their particle count and total helicity. [Eq. (11)](https://arxiv.org/html/2609.20914#S0.E11) gives a higher-loop zero subject to an operator-length relation and the stated weight inequalities. The theorem's scheme statement also requires finite changes that preserve the length selection rule.

The discussion preceding [Eq. (6)](https://arxiv.org/html/2609.20914#S0.E6) identifies exceptional amplitudes from non-holomorphic Yukawa sectors and explains how the weight bound is relaxed. Its example uses the product of Standard Model up- and down-type Yukawas. The mere appearance of a complex-conjugated Yukawa matrix is not, by itself, the paper's definition of this exception.

**Use here:** after a complete physical-basis calculation, helicity and length can provide useful independent checks on proposed anomalous-dimension entries. The present retained-mass, source-preserving off-shell transport has not established all those premises. A chiral curvature image or a raw `5/3 K_S` derivative coefficient therefore supplies neither a physical renormalization coefficient nor an all-order non-renormalization result. No beta-function entry is imported from this paper.

This review supports keeping the derivative/strength image, source transport, complex-mass ordering and field-content slice as separately defined checks. It establishes no novelty, complete-EFT closure, physical running or finite photon prediction.
