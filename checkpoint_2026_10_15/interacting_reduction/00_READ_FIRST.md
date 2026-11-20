# Interacting electroweak reduction and a local normalization convention

15 October 2026. Provisional continuation of candidate
`d822e9f4d70ad5b455436f40abc3c2b45909bbef`, prepared with AI assistance.
The input is the already matched zero-/one-/two-/three-gauge local action;
this continuation evaluates its source-retaining field changes and selected
interacting normalization. It is not another loop integration.

The [matrix reduction](ELECTROWEAK_REDUCTION.md) carries the full
`SU(2)_L x U(1)_Y` curvature into the combined packet. It retains the scalar,
lepton, sterile and gauge-current sources and gives the real indexed weak-field
shifts required to implement the reduction in the reference action. The
explicit curvature coefficient remains `-1/12` in this representative. Its
matrix ordering, source descendants and field changes belong together; the
isolated coefficient is not a physical dipole prediction.

The [independent action and Grassmann audit](INDEPENDENT_EW_AUDIT.md) checks
the non-Abelian derivative and variation signs. It also resolves a distinction
from the Abelian source. If `K=F Y_nu^T` and
`P_ab=(T^I phi)_a phi_b`, the relevant exchange projection is

    (K tensor P - (K tensor P)^T)/2
       = K_antisymmetric tensor P_symmetric
         + K_symmetric tensor P_antisymmetric.

Consequently a flavor-symmetric weak source can survive, including a
one-flavor charged component. The flavor-only antisymmetrization valid for
the hypercharge tensor cannot be imposed on the weak tensor. This is a
source term in the selected representative, not a measured transition moment.

The [current descendants](CURRENT_DESCENDANTS.md) explicitly contract the
weak generators into ordered fermion currents and scalar derivative currents.
On a constant neutral scalar background the selected real field changes have
`delta A=0`, nonzero `delta Z`, and nonzero charged `delta W` directions.
The zero photon-coordinate shift does not remove the higher-point photon
UV contact calculated in the prior checkpoint.

The [interacting normalization](INTERACTING_CANONICAL.md) extends the
previous neutral quadratic normalization by a local covariant differential
operator. Direct substitution of the actual six-placement packet separates
its normalization anticommutator from an explicit interaction remainder.
Derivatives act on scalar fields and on entire neutral composite sources;
the induced mass and source terms remain part of the action. The selected
normalization has at most one gauge connection, so it leaves the matched
two-/three-gauge local kernels unchanged in this convention and order.
These are action-coordinate statements, not finite scattering results.

The source reduction and local normalization are two specified coordinate
maps of the original packet. Combining their coefficient tables requires
composing the field changes on the full action and its sources at the declared
order. The normalization calculation alone does not establish canonical
normalization of the source-reduced representative.

The [literature note](LITERATURE_METHODS.md) and
[provenance ledger](LITERATURE_PROVENANCE.json) record bounded searches and
primary-source comparisons. Source conventions and basis identities must be
converted before comparing coefficients; no unrelated coefficient is imported.

Run the four standalone checks in the pinned root Python environment:

```sh
python -I -B -O verify_electroweak_reduction.py
python -I -B -O independent_ew_audit.py
python -I -B -O verify_current_projection.py
python -I -B -O verify_interacting_canonical.py
```

The root runner executes the same stages on disposable copies. Individual
JSON receipts state the actual coverage and separate positive controls from
rejected mistakes. These are internal computational checks and independently
implemented AI-assisted derivations, not independent human peer review.

The matrix reduction passes 1,410 positive exact controls and rejects 20
mutations; the separate action/Grassmann audit passes 122 and rejects 81
negative controls. The weak-current projection passes 61 and rejects 20;
the interacting-normalization checker passes 85 and rejects 33. These counts
describe implementation controls, not counts of physical discoveries.

The selected source map and local normalization convention do not complete
the full retained-sterile operator basis. Other external-matter loop sectors,
lower-order insertions, the remaining field/parameter/contact counterterms,
and a subtraction/evanescent prescription are still required before physical
running, finite matching or an observable can be extracted. No new
fundamental mathematics, completed unification, empirical validation or
external novelty is established. The earlier Zenodo DOI continues to identify
only its unchanged baseline; this continuation is review-stage GitHub work.
