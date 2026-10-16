# Electroweak normalization and retained-sterile research update

Checked 15 October 2026. This is a bounded extension of the existing October
15 source ledger. It checks the weak-generator normalization needed for a
one-W continuation, the sourced weak EOM, and selected newer renormalization
and dipole studies. It does not claim an exhaustive search, novel operators,
a physical anomalous dimension, or an observational constraint on this
project's seed coefficient. Downloaded papers and response bodies remain in
a cache outside the public repository.

## 1. Separate Pauli matrices, generators, and operator normalizations

Liao and Ma, *Operators up to Dimension Seven in Standard Model Effective
Field Theory Extended with Sterile Neutrinos*,
[1612.04527v1, Eq. (17)](https://arxiv.org/html/1612.04527v1#S3.EGx10),
contains the established four retained-sterile operators

    O_NLB1 = epsilon_ij (N C gamma^mu L^i) (D^nu H^j) B_mu nu,
    O_NLB2 = epsilon_ij (N C gamma^mu L^i) (D^nu H^j) dual(B)_mu nu,
    O_NLW1 = (epsilon tau^I)_ij (N C gamma^mu L^i) (D^nu H^j) W^I_mu nu,
    O_NLW2 = (epsilon tau^I)_ij (N C gamma^mu L^i) (D^nu H^j) dual(W)^I_mu nu.

The displayed definitions contain no explicit gauge coupling, i, or mixed
bilinear factor 1/2. Both dual operators belong to this chosen independent
basis. The redundant dual operators in comparison Table 3 have derivatives
on other fields; they are not the same representatives. These are inherited
source facts, rechecked here rather than new results of this search.

The normalization of the weak matrix must be declared separately. The paper
does not explicitly spell out a Pauli-versus-generator equation in its body.
The standard Warsaw convention, used as the project's explicit dictionary,
is directly stated by Grzadkowski et al.,
[1008.4884v3, Section 2, Eq. (2.2)](https://arxiv.org/pdf/1008.4884v3):

    tau^I are Pauli matrices; T^I = tau^I/2,
    D_mu = partial_mu + i g' Y B_mu + i g T^I W_mu^I + ... .

Consequently, if O_W[Pauli] uses epsilon tau and O_W[T] uses epsilon T,

    O_W[Pauli] = 2 O_W[T],
    C_W[T] = 2 C_W[Pauli]

for the same action. This factor applies separately to the ordinary and
dual operators and is not supplied by the gauge coupling. Index orientation,
charge conjugation, dual sign and Hermitian conjugation still require their
own translation. A one-W loop must use the actual generator T, not infer its
normalization from an operator label.

## 2. The weak-current descendants are fixed by the action

The same Warsaw source's
[Eq. (5.1)](https://arxiv.org/pdf/1008.4884v3) gives, in its D-plus convention,

    (D^rho W_rho mu)^I = (g/2) [H^dagger i<->D_mu^I H
                            + barL gamma_mu tau^I L
                            + barQ gamma_mu tau^I Q].

Its Eq. (2.3) defines the Higgs bidirectional derivative with the Pauli
matrix. Equivalently, using T=tau/2 throughout,

    (D_rho W^rho mu)^I = +g J_W^{I,mu},
    J_W^{I,mu} = i H^dagger T^I D^mu H
                  - i (D^mu H)^dagger T^I H
                  + sum_flavors [barL gamma^mu T^I L
                                  + barQ gamma^mu T^I Q].

The quark current already includes its color sum. Right-handed SM singlets
and a gauge-singlet sterile add no weak-current term. This is the EOM of the
gauge-invariant classical action before gauge fixing; both D on the adjoint
field strength and D in the Higgs current carry their non-Abelian terms.

Zhang,
[2306.03008v2, Eq. (3)](https://arxiv.org/html/2306.03008v2#S2.EGx3),
provides a second check in D-minus conventions, with D^nu W_mu nu and Pauli
matrices sigma^I on the right. Changing the covariant-derivative convention
and exchanging the antisymmetric field-strength indices are distinct sign
operations. Its pure-SMEFT coefficient RGEs are not retained-sterile RGEs.

These equations determine which correlated contacts a weak-field divergence
operator generates. They do not determine the loop coefficient multiplying
that operator. Likewise a one-W calculation alone does not check the
non-Abelian multiple-W vertices of a completed action.

## 3. A newer retained-sterile renormalization benchmark sharpens the method

Di Zhang, *Two-loop Renormalization Group Equations in the nuSMEFT*,
[2504.00792v2](https://arxiv.org/html/2504.00792v2),
[JHEP 06 (2025) 106](https://doi.org/10.1007/JHEP06(2025)106), calculates
two-loop renormalizable and dimension-five running. It does **not** provide
the dimension-seven mixing required by the present calculation.

Its Section 3 is nevertheless directly relevant. The three redundant
dimension-five operators in
[Eq. (12)](https://arxiv.org/html/2504.00792v2#S3.EGx12) include a sterile
second-derivative bilinear and two mixed lepton-Higgs-sterile derivatives.
[Eq. (13)](https://arxiv.org/html/2504.00792v2#S3.EGx13) gives explicit field
redefinitions; [Eq. (14)](https://arxiv.org/html/2504.00792v2#S3.EGx14)
tracks their correlated effects on M_N, Y_nu, C_5 and C_HN. Eqs. (16)-(17)
demonstrate how different integrations by parts followed by naive EOM
replacement can obscure the flavor structure of the sterile derivative
operator. This is a reason to verify a fully indexed action-level field
redefinition, including canonical normalization, rather than project only
one selected tensor.

Section 6 and Eqs. (49)-(54) address apparent divergent two-loop RG functions
in flavor-dependent field/Yukawa coordinates. The paper restores the
dependence on redundant coefficients before differentiation, or considers
flavor invariants. This does not imply an inconsistency of the EFT or a new
physical divergence. It supplies a useful future diagnostic, not an extra
counterterm or beta function already calculated for this project.

The separate order-Lambda^-3 issue remains: dimension-five sources acting on
redundant dimension-six terms can contribute to dimension-seven physical
coefficients. Zhang 2306.03008v2 Section 2, Eq. (3), Table 3 and Appendix A
establish that point for pure SMEFT. A retained-sterile calculation must
derive its own complete version, including the surviving sterile sources.

## 4. Recent phenomenology is useful only after matching

**Natascia Vignaroli,** *New avenues for the neutrino dipole
portal exploration at the energy frontier*,
[2507.01130v2](https://arxiv.org/html/2507.01130v2), studies dimension-six
electroweak dipoles. Its
[Eq. (1)](https://arxiv.org/html/2507.01130v2#S2.EGx1) includes explicit
g',g and defines its tau as a generator, sigma/2. Its
[Eqs. (2)-(3)](https://arxiv.org/html/2507.01130v2#S2.EGx2) give correlated
photon/Z/W dipoles and multiple-boson contacts. Section 4.1 explains that
omitting the latter spoils electroweak gauge invariance in the studied
processes. Its collider sensitivities are projections under its mass,
mixing, flavor and EFT-validity assumptions, not measured exclusions of
this project's seed.

For a separate dimension-six action convention with Pauli matrices and no
explicit couplings,

    L_6 = c_B barL sigma^mu nu N tildeH B_mu nu
           + c_W barL sigma^mu nu N tau^I tildeH W_mu nu^I + h.c.,

setting tildeH=(v/sqrt(2),0)^T and using e=g s_theta=g' c_theta gives

    d_gamma = v/sqrt(2) (c_theta c_B + s_theta c_W),
    d_Z     = v/sqrt(2) (-s_theta c_B + c_theta c_W),
    d_W     = v c_W.

The source's coefficients are c_B=g' C_B and c_W=g C_W/2, which reproduces
its Eq. (3). The dictionary is derived using the stated standard
electroweak relation and the displayed operator normalizations.
The two neutral combinations rotate together, and cancelling d_gamma does
not generally remove d_Z or the charged-current channel. This dictionary
does not yet match any computed dimension-seven seed coefficient to c_B,c_W.

**M. F. Mustamin and M. Demirci,** *Constraints on active-sterile neutrino transition magnetic
moments from low-energy electronic recoils at direct detection experiments*,
[2511.00701v2](https://arxiv.org/html/2511.00701v2), uses
[Eq. (3)](https://arxiv.org/html/2511.00701v2#S2.EGx3)
L_int=(mu/2) bar(nu_L) sigma N F+h.c. Its normalization therefore differs
by a factor of two from an action written with d bar(nu_L) sigma N F.
Section IV.2 and Table 2 report solar-neutrino upscattering limits conditional
on sterile mass and flavor assumptions. The low-mass XENONnT
flavor-independent limit is 0.95e-11 mu_B at 90% CL with two degrees of
freedom; the manuscript explicitly notes that its analysis does not include
the full detector response simulation. Such a number is not a bound on the
unmatched dimension-seven seed.

**Tousib Ahmed, Ayan Chattaraj, Anirban Majumdar, Newton Nath and Rahul Srivastava,**
*Precision Tests of SM and new physics with the COHERENT
Ge-mini and TEXONO data*,
[2608.10104v1](https://arxiv.org/html/2608.10104v1), Section II.4,
[Eq. (31)](https://arxiv.org/html/2608.10104v1#S2.E31), gives
mass-dependent dipole-portal nuclear and electron cross sections; Sections
III-IV specify detector inputs, nuisance parameters and constraints.
These are prospective benchmark ingredients for a declared observable
analysis after matching, not an independent validation of the UV pole.

Finally, *Enhanced active-sterile neutrino polarizability at the intensity
frontier*, [2512.07691v2](https://arxiv.org/html/2512.07691v2), has a
potentially misleading dimension label for the present purpose. Its Section
I, Eq. (1) is a dimension-seven **nuLEFT** two-photon operator. The same
paragraph states that its electroweak-invariant parent is dimension eight.
It is a different sector from the one-field-strength dimension-seven
nuSMEFT operators above, and from the two-B contact required by covariant
completion of the existing seed calculation. It supplies no justification
for relabeling that contact as a new polarizability prediction.

## Search and retrieval coverage

Four new INSPIRE queries restricted by the service's `date 2025->2026`
filter returned 4, 36, 14 and 9 records in total, with respectively 4, 25, 14
and 9 screened records: 52 entries, 46 distinct records. Titles and returned
abstracts were screened; the dipole query was capped at its 25 newest
results. The date filter can select a recent publication of an older
preprint, so arXiv first-submission and update timestamps were checked
separately. No coverage beyond the 15 October 2026 research date is claimed.

The four queries target retained-sterile dimension seven, retained-sterile
dipoles, retained-sterile renormalization and electroweak dipole observables.
They are an explicit bounded extension, not a repetition represented as a
new comprehensive survey. Eight arXiv IDs were looked up, seven versioned
HTML texts were retrieved, and the Warsaw versioned PDF and its metadata
were fetched. The four search requests and ten primary metadata/text
requests all returned HTTP 200 through the inherited managed HTTPS proxy
with TLS verification enabled. URLs, safe query strings, reported totals,
byte counts and SHA256 digests appear in SOURCE_REQUEST_RECEIPTS.json.
All response bodies remain cache material outside the public repository.

The inspected subset supplies established operators, normalization checks,
a sharper field-redefinition benchmark and observable leads. It does not
establish absence of a matching calculation elsewhere, novelty, or complete
dimension-seven retained-sterile running.
