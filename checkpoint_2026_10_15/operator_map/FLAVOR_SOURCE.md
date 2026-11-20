# Flavor support of the sourced dipole packet

15 October 2026. Conditional matrix identities for the declared source action,
not a calculation of a net observable dipole or physical beta function.

Let Y be the active-by-sterile neutrino Yukawa matrix, C=C^T the symmetric
active-by-active positive-momentum seed, and F=C Y*. Keep fixed retained
masses and the action conventions of ACTION_SIGN.md. Set H=Y Y†. Transpose
and Hermitian conjugation have distinct roles throughout.

The lepton-source route produces a sterile tensor with matrix Y†F. Its
antisymmetric part is zero because

    (Y† C Y*)^T = Y† C Y*.

The charged-lepton term of the same hypercharge source contraction has
H_scalar† tildeH_scalar=0 for a single commuting Higgs doublet. These are
two specific cancellations. They do not justify discarding total EFT
sources, or the sources from reducing a derivative on the sterile field.

The sterile-source route instead produces an active composite tensor with
matrix F Y^T=C H^T. The antisymmetric coefficient is

    K = (F Y^T-Y F^T)/2 = (C H^T-H C)/2,     K^T=-K.

Only this antisymmetric flavor part contracts the same-chirality composite
tensor bilinear. In general it does not vanish. For example, with
Y=diag(1,2) and C=[[0,1],[1,0]],

    K = [[0,3/2],[-3/2,0]].

This coefficient is one source descendant of the restricted insertion.
Other operators, direct diagrams and correlated counterterms can contribute
to the same physical structure. The displayed nonzero matrix is not an
invariant low-energy transition moment by itself.

## Exact kernel and rank within these flavor conventions

Diagonalize the Hermitian positive-semidefinite H as H=U diag(lambda_i) U†,
with lambda_i real and nonnegative, and use the unitary congruences

    C_hat=U† C U*,        K_hat=U† K U*.

Then

    K_hat_ij = (lambda_j-lambda_i) C_hat_ij/2.

Consequently the source tensor vanishes precisely when C_hat has support
only within equal-eigenvalue blocks of H. A null block and any degenerate
positive-eigenvalue block are both allowed. In particular H proportional
to the identity gives K=0 for every symmetric C. This is a difference-of-
eigenvalues statement; it must not be replaced by a sum-of-eigenvalues
source-closure identity from a different calculation.

If the distinct eigenspaces have multiplicities d_a and n active flavors,
the complex-linear map Sym_n(C) to Skew_n(C), C -> K, has

    dim(kernel) = sum_a d_a(d_a+1)/2,
    rank        = (n^2-sum_a d_a^2)/2.

Proof: in the diagonal basis each independent off-diagonal symmetric entry
is multiplied by a nonzero number exactly when its two eigenvalues differ;
diagonal entries and entries within a degenerate block are killed. Each
nonzero entry contributes one independent skew entry. The congruence is
invertible, so it preserves rank and kernel dimension. No inverse of Y,
H or a retained sterile mass was used.

Frobenius norms give an additional exact diagnostic:

    ||K||_F² = (1/4) sum_ij (lambda_j-lambda_i)² |C_hat_ij|².

If at least two distinct eigenvalues exist, the squared smallest and largest
nonzero spectral gaps bound this expression by one quarter of those gaps
squared times the norm squared of the entries joining unequal-eigenvalue
blocks. If every eigenvalue is equal, K is identically zero and the nonzero-
gap bound is inapplicable. These are standard finite-dimensional linear
algebra consequences applied to the present source packet; no novel
fundamental mathematics is claimed.

## Orders and limits

The loop factor multiplies C Y*. The new active contact has two Yukawa
factors and one hypercharge coupling. Gauge-EOM contacts have an additional
hypercharge coupling. Both are correlated descendants of the chosen
operator reduction even though an unreduced coupling-order projection can
omit them. Their absence from the original four-external-field integral is
not permission to call its reduced action source-complete.

The statements require the displayed seed symmetry and source conventions.
Arbitrary F can retain a sterile antisymmetric tensor as well. Additional
light doublets or higher-dimensional terms can change the source packet;
only terms beyond the explicitly retained inverse-scale order may be
discarded by power counting. No full operator basis, finite threshold,
renormalized observable or complete running is established here.
