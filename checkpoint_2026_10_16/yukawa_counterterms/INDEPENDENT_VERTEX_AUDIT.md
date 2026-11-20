# Independent proper-Yukawa-vertex and pole-to-RG audit

15 October 2026, Pacific client date; 16 October UTC directory date.
Input: `0770c304c7e05f213dd796d3fd5a194107f38858`.

This independent audit derives the proper pure-neutrino-Yukawa one-loop
`bar L phi N_R` UV vertex from the literal action, before using a known
beta function. It then combines that result with the separately completed
[`YUKAWA_TWO_POINT.md`](YUKAWA_TWO_POINT.md) pole inputs. The executable
imports no producer code or saved coefficients. This is AI-assisted internal
verification, not independent human peer review.

## Literal field flow and every cubic one-loop topology

Use the unbroken zero-scalar-background action

```
L_Y=−bar L_alpha^a Y_(alpha i)phi_a N_Ri+h.c.,
phi=tilde H, d=4−2epsilon, kappa=1/(16 pi^2 epsilon),
eta=(+---), C=i gamma^2 gamma^0, N=N^c.
```

Write the two vertices as `A=bar L phi N_R` and its conjugate `B`.
In independent left-Weyl fields, `A` has one barred active component,
one barred sterile component and `phi`; `B` has the unbarred components
and `phi^dagger`. There is no mixed-operator half at these unlike-field
Yukawa vertices. Majorana mass propagators can connect two sterile
components of the same orientation. They do not create a `phi phi`
propagator or an active charged-lepton Majorana propagator.

For three cubic vertices, a connected one-loop three-point graph has
three internal edges. A proper triangle must put exactly one external
leg at each vertex. The selected `bar L`, `phi` and `N_R` source legs
each occur at an `A` vertex, forcing `AAA`. After those external legs
are removed, the internal scalar pair has equal charge and the active
fermion pair is barred/barred. Both propagators are absent in this
reference vacuum. Thus no proper triangle exists for that selected
dimension-four source orientation.

This is checked by exhaustive field-slot enumeration, not just by a
picture: the script generates all eight assignments of `A/B` to three
vertices, all compatible external placements, all internal pairings,
and every vertex relabeling. It tests connectedness and removes each
edge in turn to identify bridges. It includes the anomalous sterile
mass propagators during this enumeration.

There are exactly **three connected selected-source graphs**, each with
a bridge and therefore one-particle reducible:

| Bridge | Loop subgraph | Weak factor between fixed external components |
|---|---|---|
| active lepton | active-field two-point insertion | `delta_ab` |
| scalar | scalar two-point insertion | `delta_ab` |
| sterile fermion | sterile-field two-point insertion | `2 delta_ab` |

The sterile bubble alone has a free internal weak-doublet sum. These
external-leg chains supply field-renormalization effects; they must not
be counted again as proper Yukawa-vertex counterterms. The receipt records
the actual edges, external assignments and bridge tests. A superficial
degree is assigned to proper graphs and loop subgraphs, rather than
incorrectly scaling an external bridge momentum as loop momentum.

## The opposite Majorana orientation is a separate UV-finite triangle

Allowing the external unbarred sterile component, rather than the selected
`N_R` source, permits one `AAB` proper triangle. Its scalar and active
lines obey charge flow, but the two remaining sterile components have
the same barred orientation. Their propagator must contain the Majorana
mass. Its weak contraction is `delta_ab`, with no free factor two.
Turning off anomalous propagators removes this triangle in the ordinary
Dirac comparison. The other connected graph for this Majorana orientation
has a bridge.

The ultraviolet statement does not require evaluating a finite form
factor. The two fermionic numerators that distinguish the possibilities
are checked as actual chiral matrices:

```
P_R(slash q+M)P_L slash q P_R=q^2 P_R,
P_R(slash q+M)P_R slash q P_L=M slash q P_L.              (V1)
```

The first is a nonzero logarithmic control. The second is the
mass-suppressed Majorana orientation: its kinetic numerator vanishes
between the two equal projectors. Multiplication on the right by `P_R`
also kills this second chain, consistently with the source distinction.
With generic shifted momenta the mass chain is
`M slash(q+s)P_L`; it vanishes at `M=0`.

Three propagator denominators scale as `q^6`. In four dimensions the
normal numerator gives radial behavior `dq/q`, whereas the one-mass
numerator gives `dq/q^2`, which is UV integrable. Additional mass or
external-momentum insertions only lower the degree. At one loop this
proper triangle has no smaller UV loop subgraph. Finite Majorana form
factors, including possible IR questions on special kinematics, are not
calculated or claimed to vanish.

## An independent real-component and large-momentum control

Splitting `phi=(h+i a)/sqrt(2)` must preserve the relative phase of its
real components. Stripping only their common `1/sqrt(2)`, the `A` endpoint
weights are `(1,i)` and the conjugate `B` weights `(1,−i)`. Hence

```
sum_real (A A)/2=0,       sum_real (A B)/2=1.
```

The first relation is the complex-scalar `phi phi` obstruction in a
different representation. A would-be ordinary chiral chain has the
nonzero first numerator in (V1), but its same-charge internal scalar
endpoints cancel between `h` and `a`. Extracting the external `phi`
source also retains the `A` coupling and rejects its conjugate `B`
coupling. These source and internal-charge controls are checked
separately; an allowed diagram is not manufactured by replacing the
complex scalar with a single real component.

The executable evaluates twelve exact off-shell momentum choices. It
constructs the complete numerator polynomials after `q -> t q`, retaining
the shifted constant terms and the sterile mass. The coefficient of
`t^−4`, the potentially logarithmic integrand term, follows from the
actual degree-two numerator divided by `(q^2)^3`. The unweighted ordinary
control is nonzero, its charge-weighted value is zero, and the anomalous
mass triangle has zero degree-two coefficient. Incorrectly conjugating a
same-charge scalar endpoint, omitting the imaginary scalar, or replacing
the mass projector by the kinetic projector gives a nonzero adverse
control.

This uses a Wick direction with nonzero `q^2` and generic external shifts.
A positive auxiliary mass can separate the leading UV expansion from
the IR: the ordinary control
`integral q_E^2/(q_E^2+mu_aux^2)^3` has normalized simple-pole residue one.
Gamma recursion reduces its residue to `(d/2)/Gamma(3) -> 1` at `d=4`.
No scaleless integral is assigned zero and then used to erase a UV pole.
The zero proper vertex follows from charge/chirality, not from such an
integration prescription or from fitting a published beta function.

## Independent parameter-pole and four-dimensional RG reconstruction

The completed two-point inputs, recomputed by other standalone checks,
have exactly the orientations used here:

```
H_L=Y Y^dagger, H_N=Y^dagger Y, T=tr H_N,
delta Z_L=−kappa H_L/2,
delta Z_NR=−kappa H_N,
delta Z_phi=−kappa T.
```

These are external inputs to this independent vertex executable, not
two-point integrations repeated by it. Their action convention is
`L0=(1+delta Z_L/2)L`,
`N_R0=(1+delta Z_NR/2)N_R` and
`phi0=(1+delta Z_phi/2)phi`. The vanishing proper vertex then fixes

```
delta Y=kappa r_Y,
r_Y=H_L Y/4+Y H_N/2+T Y/2
   =3Y H_N/4+T Y/2.                                   (V2)
```

The field-induced action coefficient is the negative of (V2). A proper
vertex zero therefore does not mean a zero Yukawa parameter counterterm.

For completeness the independent field-flow enumeration also checks
the sterile two-point orientations: an anomalous mass source would need
same-oriented scalar and active contractions and has no graph; the
ordinary kinetic source has its one bubble and weak multiplicity two.
There is no proper chirality-flipping sterile mass pole in this coupling
sector. With the mass action
`−[bar N_L M N_R+h.c.]/2`, cancellation of its bare-field contribution
therefore requires

```
delta M=kappa r_M,
r_M=(H_N^T M+M H_N)/2.                                 (V3)
```

The transpose is essential for a general complex symmetric retained
mass. A complex singular three-sterile test matrix is allowed; no mass
inverse is used. Conversion through the left-Weyl convention
`Mcal=M*`, `Z_NL=Z_NR*` gives the same result.

The pole-to-RG conversion is also derived. For
`Y0=mu^epsilon(Y+delta Y)`, the residue `r_Y` is homogeneous of Yukawa
degree three, while the bare prefactor has weight one. The canonical
`−epsilon Y` term differentiates every `Y` and conjugate `Y`, giving
`(3−1)r_Y=2r_Y`. For `M0=M+delta M` there is no explicit mu power;
`r_M` has Yukawa degree two, giving `2r_M`. Thus, in these fixed field
coordinates and this minimal-subtraction convention,

```
16 pi^2 beta_Y | pureY =3Y H_N/2+T Y,
16 pi^2 beta_M | pureY =H_N^T M+M H_N.                  (V4)
```

These formulas are checked with rectangular complex Yukawa matrices,
noncommuting flavor products, exact rescalings of the actual residues,
and the anomalous-dimension form
`gamma_L Y+Y gamma_NR+gamma_phi Y` and
`gamma_NR^T M+M gamma_NR`. Omitting the epsilon factor two, using the
wrong chiral sterile kinetic orientation, or replacing the Majorana
transpose by an adjoint fails.

The conditional labels in the machine receipt identify the dependence
of (V2)-(V4) on supplied two-point poles. Those inputs are now completed
and match the displayed values. This audit supplies an independent
proper-vertex and RG consistency check for the **pure-neutrino-Yukawa
renormalizable coupling sector**. It does not include gauge couplings,
other Yukawas, scalar quartic effects or higher-dimensional insertions.
Keeping a running Majorana mass diagonal would require an additional
flavor-coordinate rotation and its correlated Yukawa change.

## Reproduction and remaining limits

Run

```
python -I -S -B -O independent_vertex_audit.py
```

The standard-library checker uses exact Gaussian rationals, literal
graph slots, bridge tests, weak component sums, explicit chiral matrices
and large-momentum numerator polynomials. Its
[`INDEPENDENT_VERTEX_RESULTS.json`](INDEPENDENT_VERTEX_RESULTS.json)
records **184 positive checks and 52 adverse controls**, 236 total.
These are algebraic/asymptotic controls rather than 236 separate loop
integrations. Failures use explicit exceptions and no project modules
are imported.

The calculation does not determine finite triangle amplitudes, pole
masses, a full Standard Model beta function, or any dimension-seven
proper seed counterterm. In particular the zero renormalizable proper
vertex cannot be transferred to a derivative-seed insertion. Its proper
graphs, complete operator/source mixing and subtraction conventions
remain necessary before a physical higher-dimensional coefficient vector
can run or an observable can be predicted.
