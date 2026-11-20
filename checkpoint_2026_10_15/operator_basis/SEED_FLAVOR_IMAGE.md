# Flavor image of the selected symmetric seed

15 October 2026. This refines the [previous active-source spectral criterion](../operator_map/OPERATOR_MAP.md)
by using both weak tensor channels. It is exact finite-dimensional linear
algebra for the selected source, not a physical dipole bound or a new
fundamental mathematical law. Input: candidate
`34acf849172231282a0a6e7a07a81157b2d48799` and its sourced representative.

## Source, conventions and covariance

Hold the Yukawa matrix `Y` fixed. With `n` active flavors, any number of
retained sterile flavors and symmetric complex seed `C=C^T`, the source is

    F=C Y*,     H_Y=Y Y^dagger >= 0,
    K=F Y^T=C H_Y^T,
    S=(K+K^T)/2,     A=(K-K^T)/2.                          (F1)

`H_Y` is a flavor matrix, distinct from the Higgs field. In the
[active tensor basis](ACTIVE_TENSOR_BASIS.md), `A` supplies the B and
flavor-antisymmetric W source, while `S` supplies the flavor-symmetric W
source. The actual independent-pair coefficients also contain the declared
gauge couplings, `-1/12` representative coefficient and pair factors.
The statements below concern raw `S,A` in (F1), before those normalizations.

Under a unitary active-flavor change,

    C -> U C U^T, H_Y -> U H_Y U^dagger, K -> U K U^T.

Thus symmetry and antisymmetry use transpose, not Hermitian adjoint. An
arbitrary Hermitian matrix is insufficient for the results below: the
positive-semidefinite premise follows specifically from `H_Y=Y Y^dagger`.

## The symmetric weak channel fixes the antisymmetric source

Choose an eigenbasis `H_Y=diag(h_i)` with `h_i>=0`, and transform `C,K`
by the same congruence. Direct multiplication gives

    K_ij=C_ij h_j,
    S_ij=C_ij(h_i+h_j)/2,
    A_ij=C_ij(h_j-h_i)/2.                                (F2)

For each supported pair,

    A_ij=((h_j-h_i)/(h_i+h_j)) S_ij     if h_i+h_j>0,
    S_ij=A_ij=0                        if h_i=h_j=0.      (F3)

The earlier antisymmetric spectral-gap criterion is inherited. The further
statement is that the symmetric weak source retains every seed slot visible
to this packet and determines its antisymmetric source at fixed `H_Y`,
including rank-deficient and degenerate cases. A zero antisymmetric source
does not require a zero symmetric weak source. For one nonzero Yukawa
singular value and one active flavor, `S` can be nonzero while `A=0`.

In arbitrary unitary flavor coordinates, let `P_0` be the orthogonal
projector onto `ker H_Y`. The exact image of the map in (F1) is characterized by

    H_Y K=K^T H_Y^T,        K P_0^T=0.                   (F4)

Necessity follows directly from (F1). For sufficiency work in the eigenbasis.
The second condition sets every zero-eigenvalue column of `K` to zero.
For `h_j>0`, set `C_ij=K_ij/h_j`; the first condition makes these entries
symmetric on the positive subspace. Complete mixed zero/positive entries by
symmetry and choose the zero/zero block arbitrarily. This constructs a
symmetric `C` with `K=C H_Y^T`. Only nonzero Yukawa eigenvalues are divided
in this proof; no active or sterile mass inverse occurs.

The first equation of (F4) alone is insufficient when `H_Y` is singular.
For example `H_Y=diag(0,2)` and `K=diag(1,0)` obey it but cannot be a source
of the form (F1), since their zero-eigenvalue column is nonzero.

Equivalently, the supported `A` is the unique solution of

    H_Y A+A H_Y^T=S H_Y^T-H_Y S,
    P_0 S P_0^T=P_0 A P_0^T=0.                          (F5)

In the eigenbasis the Sylvester operator multiplies an `ij` entry by
`h_i+h_j`. Its only zero denominators occur on the excluded zero/zero block.
For an indefinite Hermitian matrix the premise fails: `diag(1,-1)` with an
off-diagonal symmetric seed gives `S=0` and `A!=0`.

## Exact source-image dimension

Let `r=rank H_Y`. The kernel in symmetric `C` consists precisely of its
arbitrary symmetric block on the `n-r` zero eigenspace. Therefore

    dim_C image(C -> K)
      = n(n+1)/2 - (n-r)(n-r+1)/2
      = r(2n-r+1)/2.                                    (F6)

The map `C -> S` has the same rank by positivity and (F2). In comparison,
the antisymmetric map has rank equal to the number of unordered pairs
with unequal `h_i`, as in the earlier spectral check. For three active
flavors the full selected source-image ranks are `0,3,5,6` for `r=0,1,2,3`.
These are complex coordinate counts at fixed `Y`, not numbers of particles
or observables. They are smaller than the unrestricted active tensor sector:
its B/W monomial space has dimensions `n(n-1)/2` and `n^2`, respectively.

## A sharp relative algebraic norm bound

For nonzero `H_Y`, define

    eta_H=(h_max-h_min)/(h_max+h_min),      0<=eta_H<=1.

Equation (F3) gives a sharp bound on the raw selected source projections:

    ||A||_F <= eta_H ||S||_F <= ||S||_F.                 (F7)

Each off-diagonal ratio is bounded by `eta_H`, while diagonal entries
contribute only to `S`. Both Frobenius norms are unchanged by unitary
congruence, so (F7) is independent of the chosen unitary flavor basis.
A symmetric seed supported only on an off-diagonal pair of extremal
eigenvectors saturates the first inequality when `h_max>h_min`. When the
extrema coincide, `A=0` and equality holds trivially. For singular nonzero `H_Y`,
`eta_H=1`; for positive definite `H_Y`, `eta_H<1`. If `H_Y` is proportional
to the identity, `A=0`. If `H_Y=0`, both sources vanish and set `eta_H=0`.

This is a matching-consistency check for this selected source in the
declared representative. It is not an inequality between measured magnetic
moments, normalized operator coefficient vectors, or full amplitudes. Other
operators, field changes and other loop sectors can contribute additional
coefficients. Unitary flavor covariance does not make these source entries
invariant under arbitrary EFT field redefinitions.

## Reproduction

Run `python -I -B -O verify_seed_flavor_image.py` in the pinned SymPy
environment. It imports no repository code and uses explicit failures.
Exact linear maps and their constraint-space ranks are checked for 45
spectra with one through five active flavors, including zero, repeated and
distinct eigenvalues. Complex rational unitary rotations test the required
transposes, the projector condition and norm invariance. Additional controls
verify source ordering with an actual rectangular complex Yukawa matrix,
supported-channel reconstruction, sharpness, and the failure of extending
the result to an indefinite Hermitian matrix. The adjacent JSON separates
positive controls from rejected mistakes and records every tested spectrum.

The general statements follow from the displayed proof; the finite checks
test their implementation and edge cases. No new loop integration, full
physical coefficient map, RG function, finite matching or external novelty
claim is made here.
