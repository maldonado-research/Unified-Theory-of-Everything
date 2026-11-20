# One-loop Yukawa inputs for the retained-sterile theory

15 October 2026 (Pacific; 16 October UTC). Provisional GitHub candidate.11, continuing input commit `0770c304c7e05f213dd796d3fd5a194107f38858`.

The preceding checkpoints matched selected gauge/seed poles and transported their source-dependent action. This continuation computes a different missing input: the actual one-loop renormalizable neutrino-Yukawa field and mass poles. Two independently authored integration routes agree. A separate proper-vertex analysis closes the pure-neutrino-Yukawa part of the renormalizable parameter running. Applying the field poles to the derivative seed then supplies its external-field counterterm and RG contribution.

These are controlled ingredients for the next physical calculation. The renormalizable running recovered here is an established benchmark, not a new fundamental law. The full dimension-seven beta function, finite photon amplitude, observable prediction and unification remain unestablished.

## Calculated results

Use `d=4-2 epsilon`, `kappa=1/(16 pi^2 epsilon)` and the literal interaction `-bar L Y phi N_R+h.c.`, with `phi=tilde H`. Define

```text
H_L=Y Y^dagger, H_R=Y^dagger Y,
T=Tr H_R, Q=Tr(H_R M^dagger M).
```

The scalar mass parameter is signed by `-m_H^2 phi^dagger phi`. Gauge interactions, other Yukawas, the scalar quartic and higher-dimensional insertions are absent in this loop sector. Constant complex symmetric `M` is retained; no sterile field is integrated out.

1. **Two-point UV integration.** [Feynman-parameter integration](YUKAWA_TWO_POINT.md) gives `delta Z_L=-kappa H_L/2`, `delta Z_NR=-kappa H_R` and `delta Z_phi=-kappa T`. The scalar inverse-kernel pole is `+i kappa(T p^2-2Q)`. The proper chirality-flipping fermion mass pole vanishes, while the Majorana mass parameter still runs through its bare-field factors. A separate [large-momentum/radial calculation](INDEPENDENT_TWO_POINT_AUDIT.md) reproduces these results with a common auxiliary IR mass, whose dependence cancels from the UV residues. Literal phases, the closed fermion sign, fixed versus summed weak components and the Majorana normalization are checked.
2. **Proper Yukawa vertex.** [Explicit cubic graph and chirality analysis](YUKAWA_VERTEX.md), with an [independent vertex audit](INDEPENDENT_VERTEX_AUDIT.md), finds no logarithmic UV pole in the proper renormalizable `bar L phi N_R` vertex at pure order `Y^3`. The permitted opposite-Majorana orientation requires a mass insertion and is UV finite. This statement concerns its UV pole, not its finite amplitude.
3. **Restricted renormalizable running.** Bare-field and parameter expansion, followed by the canonical-epsilon scale derivative, recovers

   ```text
   16 pi^2 beta_Y = (3/2)Y H_R + T Y,
   16 pi^2 beta_M = H_R^T M + M H_R,
   16 pi^2 beta_mH2 = 2m_H^2 T - 4Q.
   ```

   These are the one-loop pure-neutrino-Yukawa contributions in fixed flavor coordinates and minimal subtraction. They are not full Standard Model plus sterile-neutrino running. The action's scalar mass coefficient pole is `-2 kappa Q`; the parameter pole is instead `kappa(m_H^2 T-2Q)`.
4. **Derivative-seed external-field contribution.** For the original symmetric half-normalized seed, [literal four-field substitution](SEED_FIELD_TRANSPORT.md) gives

   ```text
   delta C7_parameter | fields
     = kappa[(H_L C7+C7 H_L^T)/4 + T C7],
   16 pi^2 beta_C7 | fields
     = (H_L C7+C7 H_L^T)/2 + 2T C7.
   ```

   Varying-scalar Fourier/Grassmann checks keep derivatives on the full composite. The factor two between this simple-pole residue and RG contribution follows from `Y0=mu^epsilon(...)` and `C7,0=mu^(2epsilon)(...)`. Proper seed vertices, source/operator mixing, gauge and other interaction sectors are still required. This is an external-field input, not the complete seed counterterm or physical anomalous dimension.

Four bounded literature searches screened twenty distinct records. Three inspected primary bodies support the counterterm/scale conventions, including a January 2026 revision with a different broken-phase coordinate system; the required conversion limits are retained.

The complex mass trace uses `M^dagger M` because it belongs to right-handed propagation; the preceding left-handed source recurrences use `M M^dagger`. Their different order is fixed by the chiral action and checked by complex-unitary covariance.

## Read and reproduce

| Material | Role |
| --- | --- |
| [Two-point calculation](YUKAWA_TWO_POINT.md) · [receipt](YUKAWA_TWO_POINT_RESULTS.json) | Exact parameter integrals, literal action poles, field versus parameter masses |
| [Independent two-point audit](INDEPENDENT_TWO_POINT_AUDIT.md) · [receipt](INDEPENDENT_TWO_POINT_RESULTS.json) | Asymptotic radial/angular integration with UV/IR separation |
| [Vertex and parameter map](YUKAWA_VERTEX.md) · [receipt](YUKAWA_VERTEX_RESULTS.json) | Exhaustive cubic topology/field flow and restricted beta benchmarks |
| [Independent vertex audit](INDEPENDENT_VERTEX_AUDIT.md) · [receipt](INDEPENDENT_VERTEX_RESULTS.json) | Separate chirality/topology and UV-degree controls |
| [Seed field transport](SEED_FIELD_TRANSPORT.md) · [receipt](SEED_FIELD_RESULTS.json) | Bare composite substitution and correctly weighted RG input |
| [Focused literature assessment](LITERATURE_YUKAWA_POLES.md) · [provenance](LITERATURE_PROVENANCE.json) | Bounded current searches, primary comparisons and applicability limits |

Run the five adjacent Python scripts in the root pinned environment. The producer two-point/vertex scripts and the independent vertex audit use the standard library; the independent two-point audit and seed checker use pinned SymPy. They import no project modules, use explicit failure guards and support isolated optimized execution. The root runner includes all five in both modes. The complete publication replay is [recorded separately](../../verification/YUKAWA_COUNTERTERM_REPLAY_LINUX.json); stage counts are implementation controls, not additional independent physical results.

## Next calculation and publication state

The next priority is proper one-loop diagrams with a derivative-seed insertion at the remaining Yukawa orders, with enough external matter fields to distinguish the complete local operator/source coefficients. Those poles must be combined with the new field/parameter inputs in the same subtraction, evanescent and field-coordinate convention, and only then reduced to a physical coefficient vector. Continue to retain lower-order insertion effects and differentiate redundant coefficients before projecting them away. Actual finite sterile-threshold matching and one declared observable come after that closure.

Earlier scientific files and receipts are preserved. AI-assisted derivations and independent implementations within this chat are not independent human peer review. This candidate remains in research PR #2, with website PR #58 following the exact research commit; both stay open for review. No merge, live deployment, GitHub release, new Zenodo version or DOI is made. DOI `10.5281/zenodo.22319172` continues to identify only the unchanged N00AK-r1 baseline.
