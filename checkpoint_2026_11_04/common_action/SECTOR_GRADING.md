# Sector, spurion and bare-mu grading of the common action inputs

4 November 2026 (Pacific). Input research commit `41f9b709c0fafdefdc9e1b990caa5421ee99f7dc`. This note supplies exact bookkeeping for the earlier source counterterm, independently integrated matter poles, known field inputs and their selected classical transports. It calculates no new loop, coefficient, physical operator quotient or full beta function.

## Fix the field and spurion conventions

Use the literal [earlier source action](../../checkpoint_2026_10_15/interacting_reduction/ELECTROWEAK_REDUCTION.md), [matter poles](../../checkpoint_2026_10_23/seed_yukawa/00_READ_FIRST.md), [seed field input](../../checkpoint_2026_10_16/yukawa_counterterms/SEED_FIELD_TRANSPORT.md) and [real matter map](../../checkpoint_2026_10_30/real_matter/REAL_MATTER_TRANSPORT.md). The retained sterile field remains in the theory. A mass factor multiplying a higher-dimensional coefficient can generate a lower-dimensional local counterterm without integrating out that field.

Assign formal global lepton charges

```text
q(L)=q(N_R)=+1, q(N_L)=-1, q(phi)=0,
q(Y)=q(lambda)=q(m_H^2)=0.
```

The reference Majorana action `-bar N_L M N_R/2+h.c.` requires `q(M)=-2`. With `X=phi^dagger L`, both `bar X C5 X^c/2` and `bar X C7(-Box)X^c/2` require `q(C5)=q(C7)=+2`. Conjugation reverses every charge. These assignments are a spurion classification of the literal action, not a claim that its fixed numerical parameters preserve a physical lepton-number symmetry.

Consequently `C7 Y* M Y^dagger` is neutral, whereas the earlier `F=C7 Y*` source has charge `+2`. The two packets belong to different field/spurion grades. The potential invariant must retain the ordered complex mass word `C7 Y* M M^dagger M Y^dagger`: its three mass factors have net charge `-2`, not `-6`. A literal `M^3` is valid only in the previously declared real diagonal specialization; it must not replace the ordered word for general complex masses.

## Canonical dimensions and explicit mu powers

In `d=4-2epsilon`, use scalar dimension `1-epsilon`, fermion dimension `3/2-epsilon`, and derivative dimension one. A zero-external-gauge field monomial with `n_s` scalar slots, `n_f` fermion slots and `n_D` derivatives therefore has

```text
D_4=n_s+3n_f/2+n_D,
D_d=D_4-epsilon(n_s+n_f),
[c]_d=4-D_4+epsilon(n_s+n_f-2).
```

The bare recipient coefficient convention is `c0=mu^(w_out epsilon)(c+delta c)`, where `w_out=n_s+n_f-2`. Derivatives affect `D_4` but do not contribute to `w_out`. Conjugated couplings have the same explicit mu weight as their unconjugated partners. In the declared input convention, `w(Y)=1`, `w(C5)=w(C7)=w(lambda)=2` and `w(M)=w(m_H^2)=0`.

The table grades a displayed coefficient orientation; Hermitian partners have the opposite lepton charge and identical dimensions and mu weights. `Y^n` denotes the total Yukawa degree, not an assertion that noncommuting matrix factors can be reordered. `M^3` in the mass-degree column counts three mass factors and denotes the ordered palindrome where appropriate.

| Selected input or transport | `(n_f,n_s,n_D)` | Field-monomial `D_4` | Coefficient mass dimension | Input monomial | `w_in` | `w_out` | Coefficient lepton charge |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Seed external-field input | `(2,2,2)` | 7 | -3 | `C7 Y^2` | 4 | 2 | +2 |
| Earlier derivative-three source | `(2,1,3)` | 7 | -3 | `C7 Y*` | 3 | 1 | +2 |
| Earlier mass-dressed source | `(2,1,1)` | 5 | -1 | `m_H^2 C7 Y*` | 3 | 1 | +2 |
| Primitive four-lepton contact | `(4,0,0)` | 6 | -2 | `C7 Y^2 M` | 4 | 2 | 0 |
| Primitive four-scalar derivative | `(0,4,2)` | 6 | -2 | `C7 Y^2 M` | 4 | 2 | 0 |
| Primitive mixed derivative | `(2,2,1)` | 6 | -2 | `C7 Y^2 M` | 4 | 2 | 0 |
| Primitive quartic potential | `(0,4,0)` | 4 | 0 | `C7 Y^2 (M M^dagger M)` | 4 | 2 | 0 |
| Radial mass-squared potential change | `(0,4,0)` | 4 | 0 | `m_H^2 C7 Y^2 M` | 4 | 2 | 0 |
| Induced neutrino Yukawa change | `(2,3,0)` | 6 | -2 | `C7 M Y^3` | 5 | 3 | 0 |
| Induced scalar sextic | `(0,6,0)` | 6 | -2 | `lambda C7 Y^2 M` | 6 | 4 | 0 |
| Conditional independent `C5` insertion | `(2,4,0)` | 7 | -3 | `C5 C7 Y^2 M` | 6 | 4 | +2 |
| Conditional `C5` source-map derivative part | `(2,3,2)` | 8 | -4 | `C5 F*`, `F=C7 Y*` | 5 | 3 | 0 |
| Conditional `C5` source-map mass-squared part | `(2,3,0)` | 6 | -2 | `m_H^2 C5 F*` | 5 | 3 | 0 |
| Renormalizable lepton field input | `(2,0,1)` | 4 | 0 | `Y^2` | 2 | 0 | 0 |
| Renormalizable scalar field input | `(0,2,2)` | 4 | 0 | `Y^2` | 2 | 0 | 0 |
| Majorana mass input | `(2,0,0)` | 3 | +1 | `M Y^2` | 2 | 0 | -2 |
| Yukawa parameter/field combination | `(2,1,0)` | 4 | 0 | `Y^3` | 3 | 1 | 0 |
| Quartic field-factor change | `(0,4,0)` | 4 | 0 | `lambda Y^2` | 4 | 2 | 0 |

The `m_H^2 F` source has a canonical dimension-five **field monomial**, even though it occurs in the same fixed-mass, seed-linear derivative expansion as the dimension-seven source. A mass-dressed field monomial must not be mislabeled by the dimension of the seed from which its coefficient originated. Likewise, the conditional `C5` terms are linear in each separately supplied coefficient; they are not generated by assuming a value of `C5` from the seed loop.

The common source normalizer also acts on `L5`. Its derivative-containing cross term is proportional to `rho C5 F* P`, with `P=-Box+m_H^2` and one retained sterile slot. The two pieces therefore have canonical field-monomial dimensions eight and six, respectively, while both have recipient mu weight three. The coefficient orientation `C5 F*` is neutral: `q(C5)=+2` and `q(F*)=-2`. The conjugate sterile/composite field slots are neutral together. The common-action note supplies the actual ordered tensors, literal half and Hermitian completion; grading alone does not calculate them.

Canonical field dimension is separate from inverse-scale counting. If `C7=O(Lambda^-3)`, an independent `C5=O(Lambda^-1)`, and retained masses are fixed in the expansion, **all three conditional `C5 C7` rows have inverse-scale degree four**. This includes the dimension-seven `C5 C7 M` matter term: its explicit retained mass raises the mass dimension of its coefficient without changing its fixed-mass inverse-scale order. They can be retained as a separately declared mixed-input extension, but they do not belong to a strict `Lambda^-3` truncation under these assignments. The receipt records `3 n_C7+n_C5` separately from field dimension, loop order and sterile-mass degree. A different mass-scaling assumption must be stated before changing that conclusion.

Independent probe sources are retained in the adjacent action calculation. They are not assigned recipient RG coefficients in this table: doing so would require their own declared `d`-dimensional source conventions. Gauge vertices also require expanding their elementary fields and explicit gauge couplings before applying a corresponding count; the table is restricted to the displayed zero-gauge field contents.

## Why the weight difference is two

For a local vertex `v` with `N_v` elementary scalar/fermion slots, its bare coupling weight is `w_v=N_v-2`. A connected graph with `E` external slots, `I` internal lines and `V` vertices obeys

```text
sum_v N_v=E+2I,     L=I-V+1,
sum_v w_v-w_out=(E+2I-2V)-(E-2)=2L.
```

This is a general topological and dimensional identity. At one loop the difference is two, despite the recipient exponents `0,1,2,3,4` in the table. Applying a one-loop field change to a reference interaction inherits the same homogeneous weight difference. Applying it to an existing one-loop pole instead starts at second loop order; its product must not be counted again in a one-loop assembly.

By contrast, applying a one-loop seed-linear field change to the **tree** seed remains at one loop, but is quadratic in `C7` and has fixed-mass inverse-scale degree six. Its exclusion from the present seed-linear action is an inverse-scale/seed-degree restriction, not a claim that every product involving a field change is two-loop order.

For an already signed **parameter** pole `delta c=r/(16pi^2 epsilon)`, differentiate the bare coefficient at fixed `c0`. The canonical contribution to the derivative of `r` is `-w_in epsilon r`, while the outer prefactor supplies `+w_out r/(16pi^2)`. Thus, in the specified unreduced coordinates,

```text
16pi^2 beta_c | signed residue r=(w_in-w_out)r.
```

The difference two does not fix `r`, turn a positive loop action into its negative cancelling action, close missing operator mixing or supply a physical beta function. The earlier source `+kappa P` was already a counterterm, while the primitive matter `+kappa Q` was a loop pole; the common cancelling action uses `+kappa P-kappa Q`. The action derivation must fix that sign before applying this weight conversion. Redundant-coordinate dependence must also be retained through the appropriate RG differentiation before physical projection.

Equal dimensions and weights are necessary but insufficient for adding coefficient tables. For example the four-lepton and four-scalar derivative terms both have `D_4=6`, coefficient dimension `-2` and recipient weight two, yet their field content differs. The original seed input and conditional `C5` transport both have `D_4=7` and lepton charge `+2`, yet they have different scalar counts and recipient weights. Flavor, weak and derivative tensors, source coordinates and the independent EOM/IBP quotient require further checks. In particular, consistent grading does not prove that a proper graph exists: the earlier seed-shaped proper four-point graph is absent at the declared Yukawa degree even though elementary power counting alone permits its grade.

## Reproduce

Run `python -I -B -O verify_sector_grading.py` in the root pinned environment. The standalone checker imports no repository calculations and uses explicit guards active under optimization. Its [receipt](SECTOR_GRADING_RESULTS.json) records the complete 18-row ledger, reference-action spurion checks, the general symbolic dimension and connected-graph identities, valence consistency examples and adverse controls for incorrect recipient weights, mass-spurion substitutions, field grades, inverse-scale truncation and loop-order counting. These are bookkeeping checks, not new loop integrations, independent operator classifications or independent human peer review.
