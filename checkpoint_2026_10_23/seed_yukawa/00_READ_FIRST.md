# Derivative-seed Yukawa matter loops and lower-dimensional poles

23 October 2026, Pacific. Provisional GitHub candidate.12; input commit `f1f9a892e15545c70e4fa05827ea7fefc5a438eb`.

This continuation evaluates new matter amplitudes containing the original derivative seed. It proves which proper graphs exist at one seed and two neutrino-Yukawa insertions, then integrates the four-lepton and four-Higgs sectors by two independent methods each. Retaining the Majorana mass allows a dimension-seven seed to generate dimension-six contacts and a scalar-potential pole. These are local UV counterterms in a specified retained-sterile theory; no particle is integrated out.

The calculations advance a restricted part of the operator mixing. They do not establish a complete physical running vector, a finite amplitude, a new fundamental law, external novelty or unification.

## What changed

1. **The proper graph classification is now explicit.** The [literal slot/bridge enumeration](SEED_TOPOLOGY.md) exhausts the one-loop sector with one original seed and total Yukawa degree two. The original two-lepton/two-Higgs proper vertex has no allowed graph. Its fixed, unreduced Green-coordinate coefficient is therefore supplied by the preceding field-pole calculation in this monomial. Connected sterile-source bridges remain and cannot be relabelled as external-field factors or discarded during a physical EOM reduction.
2. **New four-lepton pole.** The [parameter calculation](FOUR_LEPTON_POLE.md) and [independent auxiliary-mass/Gaussian audit](INDEPENDENT_FOUR_LEPTON.md) agree on a nonzero contact proportional to `C7` and `R=Y* M Y^dagger`. The anomalous sterile propagator supplies the retained mass, and the derivative acting on the entire scalar/lepton composite supplies the logarithmic loop numerator. Identical-flavor statistics do not erase it.
3. **New four-Higgs poles.** The [parameter integral](FOUR_SCALAR_POLE.md) and [independent angular/radial audit](INDEPENDENT_FOUR_SCALAR.md) agree on a mass-cubic potential pole and a mass-linear derivative interaction. Their complex flavor invariants and full conjugate derivative placement are retained. No scalar quartic or physical Higgs-mass prediction follows without specifying the remaining action and reduction.
4. **Counterterm and source transport.** A [normalized Grassmann-pair calculation](GREEN_MIXING_AND_SOURCES.md) gives a rank-at-most-two primitive flavor matrix per weak triplet. Explicit bare-mu differentiation supplies the projected MS inputs in declared coordinates. A local scalar phase change can remove the imaginary derivative placement only with its induced Yukawa and external-probe terms retained.

## Convention-fixed loop results

Use `d=4-2epsilon`, `kappa=1/(16 pi^2 epsilon)`, `phi=tilde H`, `P0=-Box` on the entire neutral composite, and

```text
L_seed=(1/2) C_alpha beta (bar L_alpha phi) P0(phi^T L_beta^c)+h.c.,
C=C^T,  L_Y=-bar L Y phi N_R+h.c.,
L_M=-(bar N_L M N_R+h.c.)/2, M=M^T.
```

For the primary four-lepton projection define `epsilon_spin=[[0,1],[-1,0]]` and

```text
P_alpha beta gamma delta
 =sum_(u,v)(bar l_alpha^u epsilon_spin bar l_beta^{v T})
           (l_gamma^{u T} epsilon_spin l_delta^v),
R=Y* M Y^dagger.
```

The four-lepton local loop action is

```text
L_loop,4L^UV=+(kappa/2) C_alpha beta R_gamma delta
                          P_alpha beta gamma delta+h.c.
```

Its cancelling action has the opposite sign. The vector-current Fierz translation is checked for four-dimensional physical spinors; it is not a completed dimensional/evanescent basis.

For the scalar result, set

```text
rho=phi^dagger phi,
K1=Tr(C Y* M Y^dagger),
K3=Tr(C Y* M M^dagger M Y^dagger),
O=(phi^dagger partial_mu phi)^2
       -rho partial_mu phi^dagger partial^mu phi.
```

Then

```text
L_loop,4phi^UV=kappa[K3 rho^2+K1 O]+h.c.
```

Again the cancelling action has the opposite sign. `K1` has dimension `-2` and `K3` dimension zero. The ordered mass palindrome is essential for complex `M`; both expressions transform covariantly without a mass inverse and vanish at `M=0`. The imaginary part of `K1` multiplies a distinct Hermitian derivative packet and must not be dropped.

These statements follow from actual UV integration, literal Wick/closed-loop signs, seed and expansion factors, weak components, Bose/Grassmann identities and complex-flavor covariance. Coefficient agreement between independently authored routes is distinct from independent human peer review.

## Read and reproduce

| Material | Role |
| --- | --- |
| [Topology and derivative placement](SEED_TOPOLOGY.md) · [receipt](SEED_TOPOLOGY_RESULTS.json) | Complete literal proper/reducible graph classification in the declared sector |
| [Four-lepton pole](FOUR_LEPTON_POLE.md) · [receipt](FOUR_LEPTON_RESULTS.json) | Actual simplex/radial integral and ordered local counterterm |
| [Independent four-lepton audit](INDEPENDENT_FOUR_LEPTON.md) · [receipt](INDEPENDENT_FOUR_LEPTON_RESULTS.json) | Independent UV extraction, Gaussian/action factors and component identities |
| [Four-Higgs pole](FOUR_SCALAR_POLE.md) · [receipt](FOUR_SCALAR_RESULTS.json) | Off-shell polynomial, local derivative/potential action and constant-background sign control |
| [Independent four-Higgs audit](INDEPENDENT_FOUR_SCALAR.md) · [receipt](INDEPENDENT_FOUR_SCALAR_RESULTS.json) | Auxiliary-mass angular/radial integral and separate action reconstruction |
| [Green mixing and sources](GREEN_MIXING_AND_SOURCES.md) · [receipt](GREEN_MIXING_RESULTS.json) | Pair normalization/image, primitive MS inputs and an explicit source-retaining phase change |
| [Literature methods and limits](LITERATURE_SEED_YUKAWA.md) · [provenance](LITERATURE_PROVENANCE.json) | Eighteen distinct records in three successful bounded queries; five primary bodies checked |

Run the six adjacent Python programs in the root pinned environment. The topology checker uses only the standard library; the other five use pinned SymPy. All use explicit guards, import no repository modules and support isolated optimized execution. The [full root replay](../../verification/SEED_YUKAWA_REPLAY_LINUX.json) includes all six; implementation control counts do not extend the scientific scope of a calculation. The literature assessment records the unsuccessful credit-blocked web retrieval attempts and the successful primary HTTPS checks; it is not an exhaustive web or novelty review.

## What remains and publication state

The immediate missing amplitude is the **mixed two-lepton/two-Higgs Majorana-mass triangle**. Its two derivative placements differ: internal legs in the same neutral composite lose their loop momentum from `P0` and are UV finite, while opposite-composite legs have a potentially divergent momentum-linear term. A degree bound is not that term's coefficient. Its integration is needed before the three allowed proper matter sectors can be collected.

The original three-point source counterterms and their interacting EOM/field maps must also be combined with the new poles before a physical basis is chosen. Other couplings, additional higher-operator/lower-order insertions, canonical and quantum-source transport, subtraction/evanescent conventions, finite matching and one declared observable remain open. A zero original proper four-point graph does not give a physical beta-function zero. The new rank statement applies only to the single primitive contact.

Earlier scientific files and receipts are preserved. Research PR #2 precedes the source-pinned website update in PR #58; both remain open for review. No merge, live deployment, tag/release, new Zenodo version or DOI is made. DOI `10.5281/zenodo.22319172` continues to identify only the unchanged N00AK-r1 baseline.
