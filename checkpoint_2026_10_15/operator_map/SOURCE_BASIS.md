# Primary-source constraints on the retained-sterile hypercharge map

Checked 15 October 2026. This is a targeted equation-level extension of the
public October 15 literature ledger, addressing one named gap: which existing
retained-sterile dimension-seven operators and sourced gauge EOM constrain the
new transverse polynomial. It is not a novelty search or a complete EFT matching
calculation. This note is newly authored; downloaded paper full texts remain
outside the repository.

## 1. The relevant operators are already in the corrected basis

Yi Liao and Xiao-Dong Ma, *Operators up to Dimension Seven in Standard Model
Effective Field Theory Extended with Sterile Neutrinos*,
[arXiv:1612.04527v1](https://arxiv.org/abs/1612.04527v1).
Fresh HTML was retrieved successfully and Section 3, item (4),
[Eq. (17)](https://arxiv.org/html/1612.04527v1#S3.EGx10), its comparison
Table 3, and Eqs. (32)–(34) were read.

With the source's right-handed sterile field N, its four independent
N-psi-H-D-X operators are

    O_NLB1 = epsilon_ij (N C gamma^mu L^i) (D^nu H^j) B_mu nu,
    O_NLB2 = epsilon_ij (N C gamma^mu L^i) (D^nu H^j) dual(B)_mu nu,
    O_NLW1 = (epsilon tau^I)_ij (N C gamma^mu L^i) (D^nu H^j) W^I_mu nu,
    O_NLW2 = (epsilon tau^I)_ij (N C gamma^mu L^i) (D^nu H^j) dual(W)^I_mu nu.

Here `N C gamma L` is the source's charge-conjugate bilinear notation, and the
dual includes the conventional factor 1/2. Eq. (17) has **no explicit i and no
gauge coupling** multiplying the operators. The displayed definitions also have
no factor 1/2 for the mixed N-L bilinear. Do not carry the seed's symmetric
identical-flavor normalization onto these mixed bilinears without deriving the
translation. We use Pauli matrices tau^I and generators T^I=tau^I/2 in the
project's explicit dictionary. The source does not spell out that
normalization equation in its body; the direct Warsaw convention check
and coefficient translation are recorded in the subsequent
[electroweak source audit](../electroweak_source/ELECTROWEAK_SOURCE_CHECK.md).
The four operators are non-Hermitian; the conjugate channel must be retained
when constructing a real action. Source Section 3 states that its dimension-seven
operators carry two units of lepton number under its field conventions.

A local hypercharge-field-strength structure in the present calculation is
therefore compatible with an established operator family. Finding it does not
constitute discovery of a previously absent operator. Mapping a particular
open-chain polynomial to these coefficients still requires consistent field
orientation, charge conjugation, flavor, Fourier and action conventions.

There is an important Table 3 distinction. The older list contains operators
with a derivative on N or L, with and without dual field strengths. Its four
dual operators are marked redundant in that comparison. This does **not** say
that O_NLB2 or O_NLW2 in the new Eq. (17) are redundant: both new dual operators
have the derivative on H and belong to the independent chosen basis.

## 2. A published sourced identity shows what reduction must keep

Liao–Ma [Eq. (34)](https://arxiv.org/html/1612.04527v1#S3.EGx22) explicitly
reduces a dual-B operator with a derivative on L. With flavor order suppressed
as in their displayed calculation, define

    A_L = overline(N^C) gamma^mu (tilde(H)^dagger D^nu L) B_mu nu,
    A_Ldual = overline(N^C) gamma^mu (tilde(H)^dagger D^nu L) dual(B)_mu nu,
    T_NNB = (H^dagger H) (overline(N^C) sigma^mu nu N) B_mu nu.

Their final equality is

    2 A_Ldual = -2 i A_L + Y_N T_NNB,

using the **sourced** lepton EOM. Thus replacing a dual derivative-on-L
structure generates a sterile-pair dipole contact, rather than just deleting
one tensor. The source's sterile-pair dipole is the family O_NNB in
[Eq. (20)](https://arxiv.org/html/1612.04527v1#S3.EGx12).

The source uses sigma_mu nu = i[gamma_mu,gamma_nu]/2. Its Eqs. (32)–(33)
fix the epsilon/gamma5 translation, including

    sigma_mu nu P_L = (i/2) epsilon_mu nu rho sigma sigma^rho sigma P_L,

with the opposite sign on P_R. An overall dual-tensor sign must be checked
against these identities and the project's chosen epsilon convention; it must
not be inferred from an operator's name. Conjugating or transposing the
open-chain channel also changes the relevant matrix order. The source's
flavor-suppressed displayed Eq. (34) is not a substitute for a fully indexed
project coefficient identity.

## 3. The sourced hypercharge EOM and its sign are directly checkable

Di Zhang, *Renormalization Group Equations for the SMEFT Operators up to
Dimension Seven*,
[arXiv:2306.03008v2](https://arxiv.org/abs/2306.03008v2).
Fresh HTML was retrieved successfully and the renormalizable action, covariant
derivative, EOM and Green-basis discussion in Section 2 were read.
[Eq. (3)](https://arxiv.org/html/2306.03008v2#S2.EGx3) gives

    D^nu B_mu nu = g1 J_Y,mu,
    J_Y,mu = (1/2) H^dagger i<->D_mu H
               + sum_f Y(f) overline(f) gamma_mu f.

The sum runs over the SM chiral fermion multiplets and flavors. Their
hypercharges are Y(Q)=1/6, Y(u)=2/3, Y(d)=-1/3, Y(L)=-1/2, Y(e)=-1.
A gauge-singlet sterile contributes no hypercharge-current term. A charged
fermion vector current has its own representation sums; do not insert an
extra color or weak factor into a current that already sums those components.

Zhang declares **D_mu = partial_mu - i g1 Y B_mu - ...**. Accordingly,
raising indices and exchanging the antisymmetric field-strength indices gives

    partial_mu B^mu nu = -g1 J_Y^nu       [Zhang's D-minus convention].

For an action declared instead with D_mu = partial_mu + i g' Y B_mu + ...,
variation of the same kinetic terms gives

    partial_mu B^mu nu = +g' J_Y^nu       [D-plus convention].

The scalar term has weight **+1/2 within J_Y** in both definitions. A change
in covariant-derivative sign and an exchange of the field-strength indices
are separate sign changes. This is a direct check against mistaking the
source's D^nu B_mu nu for the project's partial_mu B^mu nu.

With Fourier exp(-ikx), the one-B kernel of partial_mu B^mu nu is

    -k^2 delta^nu_alpha + k^nu k_alpha.

It therefore supplies precisely the k^nu slash(k)-k^2 gamma^nu structure
in the project's transverse polynomial. Reducing it using the sourced gauge
EOM generates the correlated Higgs-current and all charged-fermion-current
contacts. If its prefactor already carries one g', those contacts acquire one
additional g'. Those contacts can occupy different external-field and formal
coupling sectors from the four-point loop residue; that does not justify
silently removing them from the reduced effective action.

The quoted EOM is a gauge-invariant classical-action identity before adding
gauge-fixing terms. Zhang's paper concerns SMEFT after sterile removal. Its
renormalizable hypercharge identity remains useful with neutral steriles
retained, but its pure-SMEFT anomalous dimensions do not automatically become
retained-sterile anomalous dimensions. At the linear dimension-seven order,
a dimension-seven redundant operator uses the renormalizable EOM; complete
order-Lambda^-3 reductions may additionally require dimension-five EOM
corrections acting on dimension-six redundant coefficients.

## 4. Retained sterile sources have a separate direct benchmark

Gitte Elgaard-Clausen and Michael Trott, *On expansions in neutrino effective
field theory*, [arXiv:1703.04415v3](https://arxiv.org/abs/1703.04415v3).
Fresh HTML was retrieved successfully and Section 2.1, Eqs. (7)–(9), was read.
The text explicitly retains sterile states that have not yet been integrated
out. Its lepton EOM contains the charged-lepton Yukawa, surviving sterile
Yukawa, and Weinberg-operator sources. Its sterile EOM contains both its mass
and its light-field Yukawa sources. These are required when integration by
parts moves a derivative onto a retained sterile field.

This paper declares D=partial+i g1 y B+... in
[Eq. (3)](https://arxiv.org/html/1703.04415v3#S2.EGx2), opposite Zhang's
D-minus convention. The difference is conventional, but mixing source
coefficient signs without translating it would be an error.

## 5. What the present polynomial can safely establish

The displayed project remainder is

    D^nu = [slash(k),gamma^nu] slash(r)/4
               + [k^nu slash(k)-k^2 gamma^nu]/3.

Writing V_a^nu = a^nu slash(k)-(a.k)gamma^nu and
E_a^nu = Gamma(k,nu,a), the Clifford identity gives

    D^nu = (V_r^nu + E_r^nu)/2 + V_k^nu/3.

This identifies a field-strength/fermion-derivative channel and a
derivative-of-field-strength channel. In a generic off-shell, homogeneous
quadratic, odd-Clifford one-B ansatz built only from r,h,k, there are five
transverse directions: V_r,V_h,V_k,E_r,E_h. The grade-one ansatz has 15
structures (six gamma^nu times scalar products and nine a^nu slash(b));
its contraction has rank 12. The three grade-three structures have
contraction rank one and transverse dimension two. This kinematic counting
is restricted to that stated tensor sector. It does not prove a complete
independent gauge-invariant Green basis or establish independence modulo
EOM.

The accompanying [OPERATOR_MAP.md](OPERATOR_MAP.md) supplies a scoped
translation to the independent H-derivative B/dual-B family while tracking
lepton-EOM, sterile-EOM and gauge-EOM contacts. The current note does not claim a fully indexed action-level
coefficient map, new physics, a physical beta function, or a new observable.
Its definite advance is replacing an inherited metadata-level basis reference
with primary operator definitions and an explicit sourced-reduction benchmark.

## Retrieval record and scope

The three versioned primary HTML requests above returned HTTP 200 through the
inherited managed HTTPS proxy with TLS verification enabled. Version-specific
response byte counts and SHA256 digests are recorded in
[SOURCE_REQUEST_RECEIPTS.json](SOURCE_REQUEST_RECEIPTS.json).
Downloaded full texts and local text extractions are research cache material,
not public release files. Only this newly authored note and primary-source
links are suitable for inclusion in the repository. No broad search queries
were repeated; the named gap was the retained-sterile basis and sourced
hypercharge EOM normalization. No novel operator or missing-literature claim
is made from this bounded check.
