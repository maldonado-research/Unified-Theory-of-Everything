# Independent algebraic active-tensor basis

15 October 2026. Derived independently from the candidate.7 matrix source,
before reading the candidate.8 producer's basis or checker. The task is the
complete **algebraic nonderivative monomial sector** with two same-chirality
active leptons, two identical commuting `phi=tilde H` scalars, and one
hypercharge or SU(2) field strength. This is not a complete EFT quotient under
derivative equations of motion, integration by parts, or field redefinitions.

For `n` active flavors the answer is `n(n-1)/2` independent hypercharge
operators and `n^2` independent SU(2) operators, per chiral orientation.
Their coefficients are complex; appending Hermitian conjugates makes the
action real and does not double this complex operator count. At `n=3` this
sector has **3 B operators and 9 W operators**.

## Explicit operators and algebraic completeness

Use `D=partial+i A`, `T^I=tau^I/2`, `y_L=-1/2`,
`sigma^{mu nu}=i[gamma^mu,gamma^nu]/2`, and `X_alpha=phi^dagger L_alpha`.
Let

    B_alpha beta = bar X_alpha sigma^{mu nu} X_beta^c B_mu nu,
    W_alpha beta = (bar L_alpha T^I phi)
                            sigma^{mu nu} X_beta^c W^I_mu nu.       (1)

Neither coupling constants nor flavor coefficients are included in these
operator definitions. The field strengths retain the full nonlinear
`D=partial+i A` convention. The two W channels below are named for their
**flavor** exchange symmetry:

    W_A,alpha beta=(W_alpha beta-W_beta alpha)/2,
    W_S,alpha beta=(W_alpha beta+W_beta alpha)/2.                    (2)

The same-chirality spin block of `sigma^{mu nu} C` is symmetric. Exchanging
the two Grassmann lepton slots therefore changes the sign of their combined
flavor/weak coefficient. The commuting scalar pair is symmetric in its two
weak slots. The allowed weak tensors are

    P_B = phi phi^T,
    P_W^I=(T^I phi)phi^T,
    P_W,S^I=(P_W^I+(P_W^I)^T)/2,
    P_W,A^I=(P_W^I-(P_W^I)^T)/2.                                (3)

The flavor-antisymmetric W operator uses `P_W,S`; the flavor-symmetric W
operator uses `P_W,A`. This opposite pairing is essential. In particular,
`B_beta alpha=-B_alpha beta`, while `W_alpha alpha` can be nonzero.

Completeness follows before choosing a scalar background. The weak scalar
pair is the triplet `Sym^2(2)=3`; the two dual lepton weak slots are
`2* tensor 2*=1_A+3_S`. In the B channel, only `3_S tensor 3` contains
a singlet, with multiplicity one. In the W channel one needs a triplet to
contract with the field strength: `1_A tensor 3` and `3_S tensor 3` each
contain one triplet. These give precisely the two tensors in (3).
Contractions written using epsilon tensors are included in this count.

For the Lorentz slots, two same-chirality Weyl spinors have a symmetric
spin-one component and an antisymmetric spin-zero component. A field
strength contracts the spin-one component uniquely. Its opposite duality
does not yield an additional same-chirality invariant. Explicit chiral
Clifford matrices in the checker give rank three for the six antisymmetric
`mu nu` spin matrices. The spin-one intertwiner equations have nullity one.
Thus adding an operator with a dual field strength does not double the
complex count in this fixed-chirality sector.

An independent normalized basis is consequently

    B_alpha beta       for alpha<beta,
    W_A,alpha beta     for alpha<beta,
    W_S,alpha beta     for alpha<=beta.                            (4)

The counts are respectively `n_A=n(n-1)/2`, `n_A`, and
`n_S=n(n+1)/2`. For `n=1,2,3,4`, the total counts are `1,5,12,22`.
This construction proves linear independence and spanning under the stated
algebraic identities. It does not assert independence after adjoining a
larger derivative/EOM sector: a curvature can itself be represented by a
covariant-derivative commutator, and the previously established source and
field-redefinition relations can transport these coefficients.

## Map of the selected source and independent-pair normalization

Write `K=F Y_nu^T`, `K_A=(K-K^T)/2`, and `K_S=(K+K^T)/2`, keeping
flavor order. The existing source action is

    P_source=-1/12 sum_alpha beta K_alpha beta
                      [g'y_L B_alpha beta+g W_alpha beta].         (5)

The `-1/12` belongs to the combined selected packet. It must not be replaced
by the isolated transverse-packet coefficient `-1/4`. There is no extra
Majorana half factor in this mixed source action.

With the basis normalization in (2)-(4), (5) becomes

    P_source = sum_alpha<beta [
          -g'y_L/6 (K_A)_alpha beta B_alpha beta
          -g/6 (K_A)_alpha beta W_A,alpha beta
          -g/6 (K_S)_alpha beta W_S,alpha beta]
          -g/12 sum_alpha K_alpha alpha W_S,alpha alpha.           (6)

The off-diagonal factors of two arise from summing both ordered flavor
pairs. There is no such doubling for a diagonal symmetric operator. Using
unhalved definitions of `W_A/W_S` changes these coefficients, so the basis
normalization must always accompany the map.

For nonzero gauge couplings and arbitrary `K`, the source map has rank
`n^2` inside the `n^2+n_A` dimensional tensor sector. Its B and W_A
coefficients are correlated by the same `K_A`. Thus even arbitrary `K`
does not independently populate every B and W operator. Setting a coupling
to zero changes the rank of this particular source map, not the algebraic
operator count.

For the selected symmetric seed, `F=C_seed Y_nu^*`, `C_seed=C_seed^T`,
and `H_Y=Y_nu Y_nu^dagger`,

    K=C_seed H_Y^T,
    K_A=(C_seed H_Y^T-H_Y C_seed)/2,
    K_S=(C_seed H_Y^T+H_Y C_seed)/2.                              (7)

This further constrains the source image and does not change the basis
count. In a basis diagonalizing `H_Y` with nonnegative eigenvalues `h_i`,
the symmetric channel has entries `C'_ij(h_i+h_j)/2`, while the
antisymmetric channel has `C'_ij(h_j-h_i)/2`. Where the sum is nonzero,
the former determines the latter; if both eigenvalues vanish, both
entries vanish. This is a source-coordinate statement, not a relation
among independent operators or a physical Wilson-coefficient prediction.

## Neutral and charged component projections

At `phi=(u,0)` with real `u`, define
`N_alpha beta=bar nu_Lalpha sigma^{mu nu}(nu_Lbeta)^c` and
`C_alpha beta=bar e_Lalpha sigma^{mu nu}(nu_Lbeta)^c`, leaving their
two Lorentz indices implicit. The former obeys `N_beta alpha=-N_alpha beta`.
The latter exchanges two different weak components and is not flavor
antisymmetric. Define the full field-strength combination
`W^-_mu nu=(W^1_mu nu+i W^2_mu nu)/sqrt(2)`. Then

    B_alpha beta=u^2 N_alpha beta B,
    W_alpha beta=u^2[N_alpha beta W^3/2
                                      +C_alpha beta W^-/sqrt(2)],
    W_A,alpha beta=u^2[N_alpha beta W^3/2
                         +(C_alpha beta-C_beta alpha)W^-/(2sqrt(2))],
    W_S,alpha beta=u^2(C_alpha beta+C_beta alpha)W^-/(2sqrt(2)).    (8)

In particular `W_S,alpha alpha=u^2 C_alpha alpha W^-/sqrt(2)` survives
for one flavor. All flavor-symmetric W operators are invisible to this
neutral-component projection. A neutral-only count therefore misses an
entire independent channel. The sign in `W^-` follows from the lower-left
entry of `T^I W^I`; it is not inferred from a convention for a free curl.
The nonlinear terms remain part of the field strength.

Equations (6) and (8) reproduce the previous matrix-source projection.
They do not perform canonical normalization, broken-phase matching, or
calculation of an observable.

## Independent exact verification

Run with SymPy 1.14.0:

```sh
python -I -B -O independent_tensor_basis.py > INDEPENDENT_TENSOR_BASIS_RESULTS.json
```

The standalone core imports no repository code or producer result. It
solves all three infinitesimal SU(2) invariance conditions on the entire
weak coefficient spaces: dimension 12 and rank 11 for B, dimension 36
and rank 34 for W. It separately verifies the exchange eigenvalues and
that the displayed tensors span the resulting nullspaces. The Lorentz
spin-one coefficient space has dimension 12 and constraint rank 11.

For one through four flavors it builds full sparse polynomials in exterior
Grassmann generators, the unfixed commuting scalar components, and all six
two-form components. Exact coefficient ranks establish independence of (4).
Further tests reconstruct arbitrary complex nonsymmetric `K` using the
independent-pair normalizations, check source-image ranks, and evaluate
neutral and charged projections. Negative controls reject a missing weak
channel, doubling ordinary/dual structures, omitted diagonal W tensors,
lost off-diagonal factors, wrong charged phase, reordered flavor matrices,
and the wrong inherited packet coefficient.

The optimized receipt records **80 exact checks** and **25 rejected
negative controls**. These certify algebraic representation and component
calculations only. No loop integral, beta function, complete retained-
sterile EFT basis, derivative-EOM quotient, or physical prediction is
claimed.
