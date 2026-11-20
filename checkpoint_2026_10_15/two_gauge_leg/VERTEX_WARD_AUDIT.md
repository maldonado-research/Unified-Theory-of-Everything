# Independent seven-diagram two-B vertex and topology audit

15 October 2026. This audit derives the unintegrated two-hypercharge-field
sector of the literal action in `operator_map/ACTION_SIGN.md`, starting from
its zero-gauge normalization `-F`, where `F=C_(7,seed) Y_nu*`. It does not
integrate the diagrams. The two independent loop producers perform that step.
The action and one-B conventions are inputs; no repository module is imported
by the standalone checker.

Use `D=partial+i g' y B`, metric `(+---)`, and outgoing `r,h,k1,k2`, so
`p=r+h+k1+k2`. The internal fermion is `chi=L^c` with charge `-yL`, while
the scalar is `phi=tilde H` with charge `yL`. The loop momentum `q` is the
fermion momentum adjacent to the Yukawa vertex and flows toward the seed.
All external momenta are generic and off shell.

## Seven kernels and their relative weights

Define

```
K = k1+k2,  t=p-K,  s=p-q,
E_i = slash(e_i),  D_a=(q-p+a)^2-m_H^2,
S_a = 2t^2+(q-h-a)^2.
```

All terms below carry the same rightmost `P_L` and the external epsilon
tensor inherited from the seed. The direct and crossed scalar contractions
have already been combined in `S_a`; they are not extra graph multiplicities.
The label `ij` on ordered emissions means that gauge leg `i` is nearest
the Yukawa vertex on the stated line. On a mixed graph, `i` is on the
fermion and `j` is on the scalar.

```
FF_ij = S_K slash(q-K) E_j slash(q-k_i) E_i slash(q)
        / [q^2 (q-k_i)^2 (q-K)^2 D_0],

HH_ij = S_0 slash(q) [(2s-k_i).e_i] [(2s-2k_i-k_j).e_j]
        / [q^2 D_0 D_{k_i} D_K],

FH_ij = S_{k_i} slash(q-k_i) E_i slash(q) [(2s-k_j).e_j]
        / [q^2 (q-k_i)^2 D_0 D_{k_j}],

SG = S_0 slash(q) (e_1.e_2) / [q^2 D_0 D_K].
```

The complete stripped integrand is

```
R2 = FF_12+FF_21 + HH_12+HH_21 - FH_12-FH_21 - 2 SG,
i M_BB = -F (g' yL)^2 integral R2 P_L.
```

There are seven graphs: two orderings on the fermion, two on the scalar,
the two mixed assignments, and one scalar contact. For each sequential
insertion the new propagator and gauge vertex give `i*(-i)=+1`; thus the
relative charge products are `(+1,+1,-1)` for `(FF,HH,FH)` after factoring
out `yL^2`. There is no closed-fermion-loop minus.

The canonical scalar action supplies

```
(D phi)^dagger D phi
 = (partial phi)^dagger partial phi
 + i g'yL [(partial phi^dagger)phi-phi^dagger(partial phi)] B
 + (g'yL)^2 phi^dagger phi B_mu B^mu.
```

Differentiating the last term with respect to the two labeled external
fields produces `2(g'yL)^2 g^{mu nu}`. Its vertex is `+2i(g'yL)^2 g^{mu nu}`.
Replacing one scalar propagator by two adds one factor `i`; hence the
relative contact factor is `i*(+2i)=-2`. The factor two is already the
identical-field functional derivative, and must not be canceled by another
`1/2!` or duplicated as a second seagull graph. This agrees with expanding
the quadratic action using `B=b1 e1+b2 e2`: its `b1*b2` coefficient is
`2 e1.e2`.

The seed differentiates the neutral composite as a whole, so it has no
one-B or two-B contact. The Yukawa interaction has no gauge contact and
the minimal fermion kinetic term has no two-B contact. Hypercharge has
no Abelian gauge self-vertex. At one loop and the specified powers
`C7 Y_nu* g'^2`, each gauge interaction therefore attaches to the original
two-edge seed/Yukawa loop. Radiation from an external leg is one-particle
reducible and is outside this 1PI vertex. An added gauge propagator or
other interaction would change the loop order or coupling order. These
statements establish completeness of this fixed graph sector, not of a
general effective-action calculation.

## Exact finite-integrand Ward identity

Let `R1(p,h;k,e;q)=H-L` denote the already specified one-B integrand, with
`q` adjacent to Yukawa:

```
H = [2(p-k)^2+(q-h)^2] slash(q) [(2p-2q-k).e]
    / [q^2 ((q-p)^2-m2) ((q-p+k)^2-m2)],

L = [2(p-k)^2+(q-h-k)^2] slash(q-k) slash(e) slash(q)
    / [q^2 (q-k)^2 ((q-p)^2-m2)].
```

The exact rational identity before integration is

```
R2(p,h;k1,k2;k1,e2;q)
 = R1(p,h+k1;k2,e2;q) - R1(p,h;k2,e2;q)
 + R1(p-k1,h;k2,e2;q) - R1(p-k1,h;k2,e2;q-k1).
```

It follows from

```
P_{k1} slash(k1) P_0 = P_{k1}-P_0,
Delta_{k1} [(2s-k1).k1] Delta_0 = Delta_{k1}-Delta_0,
P_a=slash(q-a)/(q-a)^2,  Delta_a=1/D_a.
```

For clarity, contraction of the two FF orderings telescopes to

```
S_K [P_K E2 P_{k1} - P_{k2} E2 P_0] Delta_0.
```

Write `A2(a)=(2s-2a-k2).e2`. Contraction of the two HH orderings is

```
S_0 P_0 [Delta_K A2(k1) Delta_{k1}
         -Delta_{k2} A2(0) Delta_0
         +2(k1.e2) Delta_K Delta_0].
```

The last term is canceled precisely by the contracted `-2 SG`. The
mixed graphs leave the endpoint terms that organize into the displayed
four one-B integrands. The final pair is a shift difference. It vanishes
only under a translation-invariant regulator such as dimensional
regularization; generally it is nonzero at any fixed `q`. Hence the exact
integrated identity and its UV-pole restriction are

```
k1.R2 = R1(p,h+k1;k2)-R1(p,h;k2).
```

Exchanging labels supplies the second Ward identity. These identities
retain the full seed polynomial and all masses. They do not infer the
transverse loop residue from a longitudinal test.

The known action predicts the stripped two-B residue

```
V2 = (e1.e2) slash(p)
     -1/4 ([slash(k1),E1] E2 + [slash(k2),E2] E1).
```

Its contraction is `(k1.e2)slash(p)-[slash(k2),E2]slash(k1)/4`, exactly
the same difference of the known one-B local kernel. The checker tests
this consequence but labels it an action-kernel check, not an independent
integration of the seven diagrams.

## Bose symmetry and open-line ordering

Exchange of `(k1,e1)` and `(k2,e2)` swaps each pair `FF`, `HH`, `FH` and
leaves `SG` fixed. Thus both the unintegrated total and the chosen local
action are Bose symmetric with no loop shift.

The two-insertion fermion numerator follows from conjugating the ordinary
ordered line, using `C gamma^{mu T} C^{-1}=-gamma^mu` and `C^2=-1`.
Transposition reverses the entire product, each vector insertion changes
its charge sign, and the conjugated momenta are `(q-K,q-k1,q)` in the
displayed left-to-right source ordering. There is no additional momentum
reversal after this step. The checker verifies this product identity as
a Dirac matrix identity before imposing either Ward contraction. All
tests also apply after the rightmost chiral projector.

## No extra local two-B freedom invisible at zero and one B

This is a statement in the fixed field/chirality sector and at local
momentum degree at most one, before on-shell conditions or EOM quotienting.
Suppose two local lifts agree at zero and one B, and let their two-B
difference be `Delta^{mu nu}`. The off-shell Ward hierarchy makes it
homogeneously transverse in both gauge legs. Eliminate `p` using momentum
conservation and regard `r,h,k1,k2` as independent. At degree at most one,
the most general matrix polynomial can be written

```
Delta^{mu nu} = A^{mu nu}(r,h) + B^{mu nu rho} k1_rho
               + C^{mu nu rho} k2_rho,
```

where `A` contains all terms independent of both gauge momenta and is at
most linear in the other momenta; `B,C` are constant in momenta. The
first Ward identity sets `A=0` and `C=0` by the independent monomials
`k1`, `k1*r`, `k1*h`, and `k1*k2`, and requires `B` to be antisymmetric
in its first and third indices. The second Ward identity then sets
`B=0`, since the coefficients of every independent `k2_nu*k1_rho`
vanish. This argument is entrywise in spinor space; epsilon and chiral
Dirac structures introduce no exception. Bose symmetry is compatible
with the conclusion but is not needed to prove it.

Equivalently, a local difference invisible at lower gauge multiplicity
would need two field strengths and thus at least two powers of external
gauge momenta. That degree is unavailable here. A chirality-changing
mass-spurion sector, nonlocal form factors, on-shell-only Ward tests,
other external fields, evanescent finite effects, and a complete EFT
operator classification are not being asserted by this proof.

## Standalone reproducibility and limits

Run

```
python3 -I -S -B -O verify_two_b_vertices.py
```

The checker uses only Python's standard library and exact rational
arithmetic. Complex Dirac matrices are represented by independent
8 by 8 real matrices, rather than importing either loop producer's
Clifford implementation. Generic rational loop/external momenta and
signed scalar mass-squared values avoid singular denominators; singular
candidates are discarded before counting an accepted test point.
Failures use explicit exceptions and remain active under `-O`.

It checks the two exact rational Ward identities including shifts,
the scalar telescoping remainder and seagull factor, diagram-by-diagram
Bose exchange, the open-line charge-conjugation order, phases, the
rightmost chirality projection, and the local action's Ward/Bose
consequences. Actual mutated kernels omit, halve, or flip the seagull;
flip mixed charges; use the wrong FF seed; reverse the FF chain; drop
an ordered FF graph; or discard the shift difference prematurely.
Every mutation must fail at one or more tested generic points.

The output JSON records the realized point counts, individual check
counts, and the number of points detecting each mutation. These are
vertex/routing checks. Only the separate loop integrations can establish
the value of the new UV residue.
