# Independent algebraic sectors and explicit selected-source maps

15 October 2026. Continuation of candidate
`34acf849172231282a0a6e7a07a81157b2d48799`, prepared with AI assistance.
This checkpoint supplies independent operator coordinates for two specified
nonderivative field-content sectors of the previously reduced action.
Completeness within these algebraic sectors is distinct from closure of
the full retained-sterile EFT under derivative equations of motion and RG.
No additional loop is integrated here.

## Active and retained-mass tensor sectors

The [active tensor derivation](ACTIVE_TENSOR_BASIS.md) fixes the full
four-dimensional Lorentz/weak invariant spaces for two active leptons,
two conjugated Higgs scalars and one electroweak strength. Grassmann exchange
gives `n(n-1)/2` independent complex B coordinates and `n^2` W coordinates
for `n` active flavors. The W space contains both antisymmetric and symmetric
flavor channels. With three active flavors the counts are **3 B plus 9 W**.
Dual strengths add no independent complex chiral coordinate.

The [independent tensor audit](INDEPENDENT_TENSOR_BASIS.md) separately solves
the full invariant equations, expands scalar and Grassmann polynomials,
and checks ordered-to-independent flavor-pair factors. A distinct retained
mass sector has `n n_N` coordinates per strength and keeps the ordered
coefficient `F M`, with no mixed-operator half factor. The tables map the
explicit `-1/12` curvature/source contribution in the chosen representative;
they do not turn it into a physical dipole prediction.

## Four-fermion current sectors

The [current basis](CURRENT_OPERATOR_BASIS.md) adds the Lorentz Fierz and
identical-field quotient missing from the earlier weak-index contraction.
For two identical barred lepton slots the exact relation is

    Oex^L_(alpha beta gamma i)=O0^L_(beta alpha gamma i).

The independent ordered lepton basis has dimension `n_L^3 n_N` and the
selected gauge-current contribution has coefficient

    C^L_(alpha beta gamma i)
       =(g^2-g'^2)F_(alpha i)delta_beta,gamma/24
           -g^2 F_(beta i)delta_alpha,gamma/12.

The flavor permutation matters. One active flavor gives the nonzero
coefficient `-(g^2+g'^2)F/24`. Quark doublets have two weak contractions;
right-handed singlets have one. The note gives their exact source maps,
color normalization and conversion to explicitly defined scalar representatives.
The [separate current audit](INDEPENDENT_CURRENT_BASIS.md) reconstructs each
formal flavor/coupling coefficient and gives explicit nonsingular component
minors, with constructive witnesses that work for arbitrary flavor counts.

## Symmetric-seed image and consistency relation

The [flavor-image proof](SEED_FLAVOR_IMAGE.md) imposes the actual source
`K=C_seed H_Y^T`, `H_Y=Y_nu Y_nu^dagger>=0`. In the `H_Y` eigenbasis,
`S=(K+K^T)/2` and `A=(K-K^T)/2` obey

    A_ij=((h_j-h_i)/(h_i+h_j))S_ij    when h_i+h_j>0,
    S_ij=A_ij=0                     when h_i=h_j=0.

At fixed `H_Y` the symmetric weak source determines the antisymmetric
source on its support. The exact image has complex dimension
`r(2n-r+1)/2`, where `r=rank H_Y`, and is characterized by global
transpose/projector constraints. A sharp relative Frobenius bound holds
for these raw selected-source projections. It is an algebraic consistency
check, not a magnetic-moment bound, an experimental prediction, or a
constraint on the full EFT coefficient vector.

## Reproduction and next scope

Run the five standalone checks in the pinned root environment:

```sh
python -I -S -B -O verify_active_tensor_basis.py
python -I -B -O independent_tensor_basis.py
python -I -S -B -O verify_current_operator_basis.py
python -I -S -B -O independent_current_basis.py
python -I -B -O verify_seed_flavor_image.py
```

The independent tensor audit and seed-image checker use pinned SymPy; the
other three use the standard library. Every checker uses explicit failures
and imports no project calculation. Receipts distinguish invariant-space
ranks, full polynomial expansions, component witnesses and negative controls.
They are internal AI-assisted derivation/implementation checks, not human
peer review or independent experimental evidence.

The [literature comparison](LITERATURE_BASIS_MAP.md) records three focused
searches, 22 distinct paper records and four checked primary bodies. The
operator classes and flavor counts agree with existing classifications;
this is not discovery of previously unknown operators. Actual coefficient
conversion requires the declared chirality, scalar epsilon, gauge, Pauli
normalization, conjugation and flavor conventions. Full texts and raw
abstracts are not redistributed.

The coefficient tables here transport the displayed curvature-source and
gauge-current terms only. Remaining scalar/source derivatives, other
sterile sectors and field-coordinate changes can contribute to these same
operator classes after further reductions. The source-reducing and local
normalizing maps must be composed on the full action. Completing that map,
including lower-order insertions, remaining loop/field/parameter counterterms
and a dimensional/evanescent prescription, precedes physical running and
finite matching. No complete unification, physical dipole, empirical
validation, new fundamental mathematics or external novelty is claimed.
The GitHub candidate remains for review; Zenodo still identifies its
unchanged earlier baseline.
