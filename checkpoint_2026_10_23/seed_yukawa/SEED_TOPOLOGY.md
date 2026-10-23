# Proper derivative-seed Yukawa topology and its derivative placements

23 October 2026, Pacific time. Input commit:
`f1f9a892e15545c70e4fa05827ea7fefc5a438eb`.

This note exhausts the literal one-loop, one-original-seed, two-neutrino-Yukawa
vertex sector in the unbroken zero-scalar-background theory. The original
seed-shaped four-point proper vertex has no graph at this coupling order.
Other four-point field contents do have allowed Majorana-mass triangles.
The derivative placement distinguishes a UV-finite mixed contraction from
ones that may require local counterterms. Topological allowance and power
counting do not supply those pole coefficients; their actual integrations
are separate calculations.

Here “two Yukawa vertices” means total Yukawa degree two and includes all
`AA`, `AB`, `BA` and `BB` orientations below. It does not mean two identical
unconjugated Yukawa matrices.

## Action, slots and retained mass

Use the preceding [literal action](../../checkpoint_2026_10_15/operator_map/ACTION_SIGN.md)
and [renormalizable Yukawa conventions](../../checkpoint_2026_10_16/yukawa_counterterms/YUKAWA_VERTEX.md):

```text
d=4-2epsilon, eta=(+---), phi=tilde H,
P0=(i slash partial)^2=-Box,
S=(1/2)C_ab (bar L_a phi) P0(phi^T L_b^c),  C=C^T,
L_seed=S+S^dagger,
A=bar L Y phi N_R,  B=bar N_R Y^dagger phi^dagger L,
L_Y=-(A+B),  N_L=N_R^c.
```

Gauge vertices, other Yukawas, scalar quartic vertices and other higher
operator insertions are excluded. Scalar and Majorana masses remain in the
propagators; no sterile field is integrated out. The mass parameter is a
constant complex symmetric matrix. No mass inverse is required.

For graph enumeration the slots are

```text
S : (bar L_0, bar L_1, phi_0, phi_1),
A : (bar L, phi, nbar),
B : (L, phi^dagger, n),
n=N_L,  N_R=n^c.
```

The two neutral seed composites pair slots `(0,2)` and `(1,3)`. The second
is differentiated. These slot labels specify the literal monomial before
identical-field differentiation; they are not distinguishable particles.
The conjugate seed interchanges all field orientations.

The ordinary weak-charged propagators join `L` to `bar L` and `phi` to
`phi^dagger`. The sterile kinetic propagator joins `n` to `nbar`. The
anomalous sterile propagators may join `n n` or `nbar nbar`, with a mass
numerator instead of a slash numerator. Charged same-orientation scalar
or light-fermion propagators do not exist in this reference vacuum.

The left-Weyl holomorphic mass in this convention is `Mcal=M*` in
`-n^T epsilon Mcal n/2+h.c.`. Thus the leading large-momentum `n n`
propagator contains `Mcal^dagger=M`, using `M=M^T`; its conjugate contains
`M*`. This orientation follows from the displayed chiral action rather
than replacing a general complex mass by a real diagonal number.
The exact chiral control

```text
P_R slash q P_R=0,             P_R M P_R=M P_R
```

and its charge-conjugate endpoint version distinguish the anomalous edge
from the ordinary kinetic edge. The pole degree below counts this mass
numerator as degree zero in loop momentum.

## Exhaustive graph proof

There are ten literal slots. A connected one-loop graph with three vertices
has three internal lines, six internal endpoints and four external fields:

```text
I-V+1=3-3+1=1,    E_external=4+3+3-2I=4.
```

A one-particle-irreducible graph has no bridge. In a connected graph on
three vertices with three edges, every vertex of a bridgeless graph has
internal degree at least two. The sum of internal degrees is six, so
each degree is exactly two and the proper graph is the triangle. Generic
self-loops and parallel edges are included before this conclusion: the
checker finds fourteen connected structural candidates satisfying the
vertex capacities and precisely one proper structural candidate.
A self-loop cannot evade the bridge criterion here.

The seed has no sterile slot. Each Yukawa triangle vertex must therefore
connect to the seed through a charged slot. For `S`, both possible seed
fields require a conjugate charged slot: only `B` supplies `L` or
`phi^dagger`. Both triangle vertices are consequently `B`. Their remaining
internal slots can connect to one another only through the anomalous
`n n` propagator. A kinetic sterile edge would need opposite sterile
orientations. The conjugate seed analogously requires `AA` and `nbar nbar`.

This proves all of the following in the declared sector:

- The selected `S` proper triangle has two conjugate Yukawa vertices and
  one retained-Majorana mass numerator; the conjugate graph is also kept.
- No allowed proper graph has an external sterile field.
- No allowed proper graph has the original external field content
  `bar L bar L phi phi` or its conjugate.
- At `M=0` every allowed proper triangle in this sector disappears.

The executable independently generates every choice of six internal slots,
every perfect pairing of those slots, all four Yukawa orientations, both
seed orientations, every bridge test and all external field contents.
It finds the following counts per seed orientation:

| Seed | Yukawa orientations | Connected labelled slot contractions | Proper | Reducible |
|---|---|---:|---:|---:|
| `S` | `AA` | 0 | 0 | 0 |
| `S` | `AB` | 8 | 0 | 8 |
| `S` | `BA` | 8 | 0 | 8 |
| `S` | `BB` | 36 | 12 | 24 |
| `S^dagger` | `AA` | 36 | 12 | 24 |
| `S^dagger` | `AB` | 8 | 0 | 8 |
| `S^dagger` | `BA` | 8 | 0 | 8 |
| `S^dagger` | `BB` | 0 | 0 | 0 |

These are **labelled literal-slot contractions**, not physical diagram
multiplicities, graph symmetry factors or amplitude coefficients. The
`1/2` in the identical-fermion seed and `1/2!` from two Yukawa insertions
must be handled in an amplitude calculation rather than inferred from
this table. Exchanging Yukawa labels likewise does not create a second
physical contribution by fiat.

## Allowed external sectors and actual seed momentum

There are three external field sectors. Their local degree bounds depend
on the *full composite* momentum in `P0`, not on a rule that every seed
insertion supplies two loop powers.

Orient the triangle cycle `seed -> vertex 1 -> vertex 2 -> seed`.
The two incoming loop momenta at the seed have coefficients `-q` and
`+q`; all external seed slots have zero coefficient of `q`.
The differentiated composite momentum has loop coefficient equal to
**the sum** of those on `bar L_1` and `phi_1`. Momentum conservation makes
the other composite momentum its negative. Squaring either composite
therefore gives the same full polynomial.

| Seed internal charged lines | External field content | Labelled proper slots per oriented seed | Seed degree in `q` | Other numerator powers | Superficial UV degree |
|---|---|---:|---:|---:|---:|
| Two light-fermion lines | `phi phi phi^dagger phi^dagger` | 2 | 2 | 2 | 2 |
| One light and one scalar, in opposite seed composites | `bar L L phi phi^dagger` | 4 | 2 | 1 | 1 |
| One light and one scalar, in the same seed composite | `bar L L phi phi^dagger` | 4 | 0 | 1 | -1 |
| Two scalar lines | `bar L bar L L L` | 2 | 2 | 0 | 0 |

Every row also has one anomalous sterile mass numerator. With three
quadratic propagator denominators, the formula is

```text
omega=4 + (light-fermion slash powers) + (seed loop powers) - 2(3).
```

For the same-composite mixed case, both internal momenta occur in the same
neutral composite. Their loop parts cancel. The seed factor is a squared
external momentum, and the remaining mass/light numerator supplies only
one loop power. This graph is UV convergent. Applying `P0` to the lepton
alone would destroy that cancellation and produce an incorrect divergent
term. For the opposite-composite mixed case there is one internal leg in
each composite; the quadratic seed loop momentum remains and the bound is
`omega=1`. Odd leading integrands can remove the leading term, while a
local term linear in external momenta may still have a logarithmic pole.
The degree is not a calculation of that coefficient.

The four-lepton case has degree zero and can generate mass-weighted local
four-fermion terms. The four-scalar case has degree two and can in principle
generate two-derivative terms or scalar-mass-dependent potential terms.
Chiral, weak, Bose and flavor sums may further constrain their actual
residues. Those sums and loop integrations must be performed before any
specific generated coefficient is reported.

All weak contractions are explicitly joined from the seed, Yukawa and
propagator Kronecker deltas. Every weak component is fixed by external
indices; there is **no free weak-doublet trace** in these proper graphs.
At fixed external weak components the tensor is zero or one. The direct
and crossed identical-slot contractions remain distinct until the full
operator vertex is differentiated; neither fact justifies appending a
blanket weak factor two.

For the mixed row, the same-composite contraction carries the weak pairing
`(bar L phi)(phi^dagger L)`, whereas the opposite-composite contraction
carries `(bar L L)(phi^dagger phi)`. This describes weak indices only;
spinor derivatives and flavor orders are separate.

## The target proper zero and the surviving reducible source bridge

The absence of a proper original `bar L bar L phi phi` graph is exact at
one seed and total Yukawa degree two. It is independent of an IR regulator,
a scaleless integral prescription, four-dimensional Fierz identities or
an epsilon-suppressed Dirac reduction.

The target does have connected **reducible** graphs. For each labelled
`AB` or `BA` orientation the eight contractions split into:

```text
2 lepton bridges: external-lepton self-energy attachments,
2 scalar bridges: external-scalar self-energy attachments,
4 sterile bridges: earlier C7 Y* three-point source bubble,
                   joined to a Yukawa tree vertex by a sterile line.
```

The sterile-bridge contribution must not be relabelled as an ordinary
external-field factor. The earlier three-point pole and its local source
counterterm are already part of the reference calculation. Their sourced
EOM or field-coordinate transport can feed the seed-shaped coefficient
when a physical operator basis is chosen. Conversely, such a bridge is
not a proper four-point loop and must not be added to the proper vertex.

Using the separately completed
[field renormalization inputs](../../checkpoint_2026_10_16/yukawa_counterterms/SEED_FIELD_TRANSPORT.md),
the target's projected parameter pole in the **fixed original unreduced
Green coordinates**, at linear `C7` and total pure-Yukawa degree two, is
therefore supplied by its field factors:

```text
H_L=Y Y^dagger, T=Tr(Y^dagger Y),
delta C7_parameter | target Green = kappa r7,
r7=(H_L C7+C7 H_L^T)/4+T C7,
kappa=1/(16 pi^2 epsilon).
```

With the already checked `C7,0=mu^(2epsilon)(C7+delta C7)` and
`Y0=mu^epsilon(Y+delta Y)` bookkeeping, its four-dimensional projected
coefficient has the previously derived contribution

```text
16 pi^2 beta_C7 | target Green, pureY
  =(H_L C7+C7 H_L^T)/2+2T C7.
```

The new result in this note is the proper-topology certificate, not a new
integration of those field poles. This closes that specific unreduced
four-point counterterm projection. It does **not** close the generated
four-scalar, mixed or four-lepton sectors, the earlier three-point source
counterterms, their EOM/field map, other couplings or the physical Wilson
coefficient vector. It must not be advertised as a complete physical seed
beta function, finite matching result or observable prediction.

## Reproduction and UV/IR scope

Run `python -I -S -B -O verify_seed_topology.py`. The standard-library
checker uses exact literal combinatorics, weak-component equations,
varying-external-momentum seed polynomials, Gaussian-rational chiral
matrices and explicit failure guards. Its adjacent
[receipt](SEED_TOPOLOGY_RESULTS.json) records **960 positive comparisons
and 41 adverse controls**. Many are rows of exhaustive weak-component
systems, not separate loop integrations or physical results. It imports
no repository modules and has no removable assertions.

The proper target zero follows from absence of graphs. For the allowed
triangles use generic nonexceptional off-shell external momenta or a
controlled positive auxiliary IR mass during a separate pole integration.
A vanishing dimensionally regulated scaleless integral must not be used
to erase a possible UV pole. No amplitude phase, finite term, tensor pole
coefficient, human peer-review outcome or new fundamental physics is
claimed by this topology/degree check.
