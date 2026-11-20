# Primary benchmarks for pure neutrino-Yukawa one-loop poles

15 October 2026, Pacific time. Retrieval receipts and the existing checkpoint directory use 16 October UTC. This comparison concerns the renormalizable retained-sterile theory and an external-field contribution to a dimension-seven counterterm. It does not establish the full dimension-seven anomalous dimension or a new physical effect.

## Bounded coverage

Four Firecrawl paper-index queries requested k=10 each and returned **6, 5, 5 and 10 records: 26 ranked appearances, 20 distinct papers and six duplicate appearances**. Three queries used the inclusive created/updated window 2025-01-01 through 2026-10-16; the fourth deliberately used 2000-01-01 onward to recover established seesaw running references. This is a bounded screen, not an exhaustive review or a novelty assessment.

Two papers received new metadata inspections and primary HTML retrievals, both HTTP 200. A previously retrieved 2025 primary body was reused and its saved hash checked. Three passage requests supplement the direct body review; their counts and exact questions are in `LITERATURE_PROVENANCE.json`. Only the named sections were reviewed. Other returned records received an abstract-level relevance screen. Full papers, abstracts and raw provider outputs remain in scratch storage; public provenance keeps identifiers, hashes, our summaries and request receipts.

## A directly comparable unbroken-theory benchmark

Di Zhang, [Two-loop Renormalization Group Equations in the nuSMEFT, 2504.00792v2](https://arxiv.org/html/2504.00792v2), rechecked body **5 May 2025**. Despite the title's two-loop emphasis, **Appendix A, Eq. (58)** explicitly collects the one-loop counterterms needed here.

The reference [Section 2](https://arxiv.org/html/2504.00792v2#S2) uses `d=4-2 epsilon`, the neutrino interaction `-bar L Y tilde H N_R+h.c.`, the symmetric Majorana mass in `-bar N_R^c M N_R/2+h.c.`, and scalar potential `+m^2 H-dagger H+lambda(H-dagger H)^2`. The bare fields use square-root wave-function factors. Thus a tabulated `delta Z` is the kinetic renormalization matrix, not the `delta Z/2` coefficient multiplying a single external field.

Define, in the present orientation with Y having active rows and sterile columns,

    H_L=Y Y-dagger,      H=Y-dagger Y,
    T=Tr(H),            S=Tr(H M-dagger M),
    kappa=1/(16 pi^2 epsilon).

Extract only the pure neutrino-Yukawa monomials from [Appendix A](https://arxiv.org/html/2504.00792v2#A1): omit gauge, other Yukawa and scalar-quartic contributions and set all dimension-five and redundant-operator coefficients to zero. The reference then gives

| Quantity | Pure-Y pole |
| --- | --- |
| delta Z_L | `-kappa H_L/2` |
| delta Z_NL | `-kappa H^T` |
| delta Z_NR, by conjugation | `-kappa H` |
| delta Z_phi, for phi=tilde H | `-kappa T` |
| additive scalar parameter delta m^2 | `kappa(m^2 T-2 S)` |
| additive Majorana parameter delta M | `kappa(H^T M+M H)/2` |
| additive Yukawa parameter delta Y | `kappa(3 Y H/4+T Y/2)` |

The scalar result is written additively: the source tabulates `m^2 delta Z_m2`, whose product remains meaningful without dividing by m^2. Complex mass and flavor orders must be preserved; in particular `S=Tr(H M-dagger M)` follows from the displayed Majorana/Yukawa orientation.

### Scalar amplitude versus scalar mass parameter

The scalar mass counterterm in the action includes the field factor as well as the mass-parameter counterterm:

    delta m^2 + m^2 delta Z_phi = -2 kappa S.

For the declared kinetic/mass inverse propagator `p^2-m^2`, its two-point counterterm vertex is therefore

    i[-kappa T p^2+2 kappa S].

It cancels a loop UV term `+i kappa[T p^2-2 S]`. This last translation is our algebraic comparison in the declared amplitude convention, not a separate reproduction of every self-energy diagram in the reference. The sign of a potential parameter, an inverse-propagator correction and the associated Feynman amplitude must not be equated without this conversion.

### Proper Yukawa vertex versus parameter counterterm

Combining the parameter and external-field terms gives

    delta Y + delta Z_L Y/2 + Y delta Z_NR/2
            + delta Z_phi Y/2 = 0.

Thus the published counterterms are consistent with **zero net proper pure-Y Yukawa-vertex UV pole**. A nonzero Yukawa parameter counterterm is still required because the fields renormalize. This comparison does not assert that every finite vertex form factor vanishes, nor replace the direct check of all allowed Majorana contractions. Mass-insertion graphs and their UV degree must be considered in that diagram calculation.

## An independent established beta-function check

Stefan Antusch, Joern Kersten, Manfred Lindner, Michael Ratz and Michael Andreas Schmidt, [Running Neutrino Mass Parameters in See-Saw Scenarios, hep-ph/0501272](https://arxiv.org/html/hep-ph/0501272), newly retrieved body **v3, 21 March 2005**. The index update of 9 July 2009 is a separate date.

[Section 2, Eq. (2)](https://arxiv.org/html/hep-ph/0501272#S2.E2) writes its Yukawa matrix with sterile rows, in `-bar N_R Y_source L tilde H-dagger`, and its Majorana term as `-bar N_R M_source N_R^c/2+h.c.`. Relative to the present action these identify `Y_source=Y-dagger` and `M_source=M-dagger`. Its **Appendix D.1, Eqs. (D.33b)-(D.33c)** then reproduces, after this conversion and extraction of the pure-Y terms,

    16 pi^2 beta_Y = 3 Y H/2 + T Y,
    16 pi^2 beta_M = H^T M + M H.

The same results follow from the pole residues above using the declared `d=4-2 epsilon` bare-coupling expansion. For the scalar parameter that expansion gives

    16 pi^2 beta_m2 = 2 m^2 T - 4 S.

These are standard renormalizable pure-Y contributions. Gauge, other Yukawa, scalar and higher-dimensional contributions are outside this extraction. The threshold matching and numerical running software in the older paper were not used or executed.

## A newer full-model renormalization reference, with different coordinates

Jihong Huang and Shun Zhou, [One-loop Renormalization of the Type-I Seesaw Model in the Modified Minimal-subtraction Scheme, 2507.21691](https://arxiv.org/html/2507.21691), newly retrieved body **v2, 4 January 2026**; index update 6 January 2026.

[Section 3, Eq. (3.3)](https://arxiv.org/html/2507.21691#S3.Ex10) uses field factors `1+delta Z/2`. The following discussion specifies `d=4-2 epsilon` and `Delta=1/epsilon-gamma_E+ln(4 pi)`; modified MS subtracts that UV combination. Section 3 also declares the Fleischer-Jegerlehner tadpole prescription. Section 4 treats physical masses and mixing matrices in the broken-phase mass basis, including the altered Majorana-neutrino self-energies and Higgs running.

This is a useful newer benchmark for complete one-loop bookkeeping, Majorana contractions and the separation of field and parameter counterterms. Its scalar potential is `-m_source^2 rho+lambda_source rho^2/4`, whereas the direct benchmark above uses `+m^2 rho+lambda rho^2`. Therefore neither its physical Higgs-mass beta function nor a mass-basis self-energy is copied into the present unbroken-basis calculation. A comparison would need the potential, vacuum/tadpole prescription, mass/mixing rotations and finite subtraction conventions converted first. No such full broken-to-unbroken conversion or finite matching has been performed here.

## UV poles must remain separated from infrared zeros

Zhang [Section 4](https://arxiv.org/html/2504.00792v2#S4), especially Eqs. (20)-(21), explains that setting masses and external momenta to zero can introduce spurious IR divergences, and gives an auxiliary-mass rearrangement with the necessary local counterterms. A dimensionally regulated scaleless integral vanishing does not by itself show that its UV residue vanishes.

For the present one-loop checks, keep a suitable off-shell scale or an explicitly controlled IR rearrangement while extracting the pole. An auxiliary regulator mass is distinct from the retained Majorana matrix M. In particular, the physical trace S in the scalar pole cannot be removed merely to simplify the IR prescription. The use of a four-dimensional chiral numerator identity for a simple pole does not independently fix finite or evanescent-scheme terms.

## What external-field renormalization fixes for a C7 coefficient

For a coefficient matrix C with two active-lepton and two scalar external fields, inserting only the checked field factors gives the compensating residue

    r_C,field = (H_L C+C H_L^T)/4 + T C.

For this term, quadratic in Y, the declared bare-coupling expansion yields

    (16 pi^2 beta_C)_field = 2 r_C,field.

This factor follows from the EFT and loop-order bookkeeping for that residue; it is not a rule that an arbitrary bare pole can always be doubled without checking its coupling powers and engineering-dimension factors. The matrix order and transpose follow the declared two-lepton coefficient convention.

This is only the **external-field contribution**. Operator-specific proper C7 vertex graphs, other counterterms, mixing, field-map transport and the chosen physical basis can add terms. A vanishing proper renormalizable Yukawa vertex does not establish vanishing proper C7 vertices. None of the primary references above supplies the full C7 beta function for the selected packet, and this field contribution is not presented as that result.
