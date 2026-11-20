# Explicit field changes and neutral canonical normalization

15 October 2026. Provisional continuation prepared with AI assistance.
Input: candidate `89f5867453bc7ed706d249dc82b393245b28d7ea` and the
[selected packet's exact residual certificate](../electroweak_source/SELECTED_PACKET_REDUCTION.md).

This continuation implements the next action-level step for that selected
hypercharge packet. It does not evaluate another loop or complete physical
running. The earlier one-W and neutral-photon pole results keep their scope.

The [indexed field redefinition](FIELD_REDEFINITION.md) gives local shifts
of the lepton, Higgs, Majorana sterile and real hypercharge gauge fields.
Substitution in the declared gauge-invariant renormalizable action converts
the original selected packet into its source-retaining representative,
through first order in its coefficient and the declared inverse-scale order.
The change includes conjugate sources, differentiated composites and ordered
mass factors. The Majorana-real completion has no additional factor `1/2`.
The [separate action audit](INDEPENDENT_ACTION_AUDIT.md) derives the same
covectors and tests the Grassmann, chirality, reality and derivative signs.

The [neutral quadratic calculation](CANONICAL_NORMALIZATION.md) starts from
the actual unreduced six-placement packet. For constant real
`phi=(u,0)` and zero gauge background it gives

    delta L^(2) = kappa u F bar nu_L A(-Box+m_H^2)N_L + h.c.,
    Z_as(p^2) = u(p^2+m_H^2)F,     F=C_seed Y_nu*,

where `K=I+kappa Z` has that off-diagonal block and its Hermitian conjugate.
The local perturbative normalization `q=(I-kappa Z/2)q_c` induces the
correlated holomorphic mass and source changes

    delta Mcal = -(Z^T Mcal+Mcal Z)/2,
    delta G = -G Z/2.

The note declares the mass-conjugation convention and retains the
momentum-dependent mass/source terms. Independent pole kernels and
source-dressed propagator identities check the coordinate transformation
without requiring a mass inverse. The momentum expansion does not justify
interpreting extra roots of a truncated polynomial as new particles.

An off-stationary scalar background adds a further Hessian term from the
scalar tadpole times the nonlinear scalar shift. Its absence cannot be
assumed when comparing the original and reduced packets. The notes state
when a simple fermion congruence is sufficient and when that extra term
must be retained.

The [lower-order closure benchmark](LOWER_ORDER_CLOSURE.md) makes a separate
missing contribution executable. For declared optional dimension-six
sterile normalization `K=I+t^2 rho Z_s` and dimension-five scalar mass
`t rho C`, the same field change induces

    Delta D_7 = -(Z_s^T C+C Z_s)/2.

These optional coefficients are inputs to the benchmark, not newly computed
seed-loop coefficients. A separate [Grassmann review](INDEPENDENT_CLOSURE_REVIEW.md)
checks this nonzero order-`Lambda^-3` contact, including a singular mass.

Run the five standalone checkers in the pinned root SymPy environment:

```sh
python -I -B -O verify_field_redefinition.py
python -I -B -O independent_action_audit.py
python -I -B -O verify_canonical.py
python -I -B -O verify_lower_order_closure.py
python -I -B -O independent_grassmann_closure_check.py
```

The root runner executes them on disposable copies. Focused JSON receipts
describe exact implementation controls and rejected mistakes; they are not
experimental observations or independent human peer review.

The indexed-shift checker passes 192 exact controls and rejects 14 mutations;
the separate action audit passes 33 controls and rejects 21. The canonical
checker passes 82 and rejects 24, the lower-order benchmark has 124 controls
including 21 rejected mutations, and its separate component review has seven
controls. These counts distinguish the implementations and their scope;
they are not counts of new physical predictions.

The [literature follow-up](LITERATURE_METHODOLOGY.md) records four bounded
2025–2026 queries: 35 appearances, 27 distinct papers. Four selected primary
texts were retrieved and relevant sections checked after an initial
Firecrawl scrape-credit failure. Their versions, byte counts and hashes,
query results and coverage limits are in [the provenance record](LITERATURE_PROVENANCE.json).
The newer general EFT templates support the need for canonical normalization,
bosonic contributions and infrared subtraction. They do not supply the
missing dimension-seven coefficient of this calculation.

The next task is to assemble a complete retained-sterile coefficient basis,
lower-order insertions, remaining diagram and field/parameter counterterms,
and a subtraction/evanescent prescription. The full interacting
canonicalization, two-W/mixed-BW matching, physical running, finite matching
and an observable remain open. No mathematical novelty, completed
unification, finite photon prediction or new physical dipole coefficient is
claimed. This is review-stage GitHub progress; the earlier Zenodo DOI still
archives only its unchanged baseline.
