# Independent action audit of the selected matter normalization

30 October 2026 (UTC). Input commit `5094584ebf402736b0bc5b597f6769f640178fd1`.
This audit was derived from the preceding primitive mixed and four-scalar loop
results before reading the new producer's implementation. It performs classical
first-order action algebra with those earlier loop poles; it evaluates no new
loop. Independent authorship and executable controls are internal AI-assisted
verification, not human peer review or empirical evidence.

## Conventions and independent derivation

Use the preceding [mixed-loop result](../../checkpoint_2026_10_26/mixed_yukawa/MIXED_SEED_POLE.md),
[four-scalar result](../../checkpoint_2026_10_23/seed_yukawa/FOUR_SCALAR_POLE.md)
and [phase map](../../checkpoint_2026_10_26/combined_matter/COMBINED_PHASE_TRANSPORT.md).
Write

```text
rho=phi^dagger phi, phi=tilde H, epsilon=i sigma^2,
K=C7 Y_nu* M Y_nu^dagger=G+i A,
G=G^dagger, A=A^dagger, g=Tr G, t=Tr A,
J_mu=i[phi^dagger D_mu phi-(D_mu phi)^dagger phi],
K3=Tr(C7 Y_nu* M M^dagger M Y_nu^dagger),
kappa=1/(16 pi^2 epsilon_UV).
```

The epsilon matrix and dimensional-regularization pole are distinct objects.
For one arbitrary derivative direction let `z=phi^dagger D phi`. Then
`z=(partial rho-iJ)/2`. Consequently the complete Hermitian four-scalar
derivative pole is exactly

```text
kappa[(g+i t)(z^2-rho |D phi|^2)+h.c.]
 =kappa[g((partial rho)^2-J^2)/2
          -2g rho |D phi|^2+t(partial rho)J].            (A1)
```

The same identity holds direction by direction with the declared Lorentz
contraction. It neither assumes a real primitive invariant nor discards its
conjugate derivative placement.

Define `F=G/2+iA` and `c=g+i t`. The complete selected field substitution
means **old fields equal new fields plus the displayed delta**:

```text
delta L=kappa rho F L,       delta bar L=kappa rho bar L F^dagger,
delta phi=kappa rho c phi,   delta phi^dagger=kappa rho c* phi^dagger.
```

Its lepton real part is nonunitary and its scalar real part is radial. With
`L_kin,L=(i/2)[bar L slash D L-(D bar L)slash L]`, direct independent
variation of all four factors gives

```text
delta L_kin,L
 =+(kappa rho/2)bar L G i<->slash D L
   -kappa(partial_mu rho)bar L A gamma^mu L.             (A2)
```

In particular the two radial derivative-of-rho terms cancel inside the
symmetric kinetic expression; the phase terms add. Using an unsymmetrized
kinetic term instead requires retaining its total derivative consistently.
For the scalar kinetic action the product rule gives

```text
delta L_kin,phi
 =kappa[2g rho |D phi|^2+g(partial rho)^2-t(partial rho)J],
delta rho=2kappa g rho^2.                                (A3)
```

Adding (A2)-(A3) to the input poles removes precisely the selected norm-times-
kinetic terms and the two imaginary-current terms. The remaining derivative
terms from these primitive poles are

```text
kappa J_mu bar L G gamma^mu L/2
 +kappa g[3(partial rho)^2-J_mu J^mu]/2.                  (A4)
```

Equation (A4) is a local representative, not an assertion that the full scalar
metric or all interacting fields are canonically normalized. A complete
independent EOM/IBP basis and the loop's external-gauge completion are separate
tasks. Covariant product rules establish the displayed action identity, not
additional gauge-vertex loop residues.

## Potential, interactions and independent probes

For the reference action `L_pot=-V(rho)`, the radial map contributes
`-2kappa g rho^2 V'(rho)`. In the established convention
`V=m_H^2 rho+lambda rho^2`, combining only this contribution with the primitive
potential pole produces

```text
kappa[2(Re K3-g m_H^2)rho^2-4g lambda rho^3].             (A5)
```

The sextic term is retained. Neither coefficient in (A5) is the total physical
scalar coefficient or a full renormalizable counterterm. There may be other
couplings and diagrams. No positivity or nonzero sign is assumed for `g`.

The tilde convention requires `H=-epsilon phi*`, so the actual Higgs change is
`delta H=kappa rho c* H`: radial parts agree and phase parts have opposite
signs. For the reference interactions

```text
L_Y=-bar L Y_nu phi N_R-bar L Y_e H e_R
    -bar Q Y_u phi u_R-bar Q Y_d H d_R+h.c.,
```

the complete first-order variation, with singlet and quark fields unchanged,
is

```text
delta L_Y=-kappa rho bar L [F^dagger+c I]Y_nu phi N_R
          -kappa rho bar L [F^dagger+c*I]Y_e H e_R
          -kappa rho c  bar Q Y_u phi u_R
          -kappa rho c* bar Q Y_d H d_R+h.c.              (A6)
```

The flavor matrix acts before the Yukawa matrix. In general that ordering
cannot be reversed. For `L_probe=bar eta L+bar L eta+j^dagger phi+phi^dagger j`,
the independently held sources induce

```text
delta L_probe=kappa rho[bar eta F L+bar L F^dagger eta
                         +c j^dagger phi+c*phi^dagger j]. (A7)
```

These interactions are retained in the transformed action. Sterile kinetic
and Majorana terms have no direct variation because the sterile fields and
parameters are unchanged in this selected substitution. The separately
computed one-loop four-lepton action also stays unchanged at this order:
applying a one-loop field shift to an already one-loop term starts at
`kappa^2`. Independently supplied tree operators require their own variations.

## Conditional lower-order insertion

Suppose the reference action independently contains
`L5=bar X C5 X^c/2+h.c.`, where `X=phi^dagger L` and `C5=C5^T`.
Then

```text
delta X=kappa rho B X,
B=F+c*I=G/2+gI+i(A-tI),
D=B^dagger=G/2+gI-i(A-tI).
```

The conjugate composite changes as `delta bar X=kappa rho bar X D` and
`delta X^c=kappa rho D^T X^c`. Direct variation therefore gives

```text
delta L5=(kappa rho/2)bar X [D C5+C5 D^T]X^c+h.c.         (A8)
```

The transpose in the charge-conjugate slot and the literal half are essential.
`D C5+C5 D^T` is symmetric even when `D` is not Hermitian. Because `D` has
mass dimension minus two, this is a dimension-seven term conditional on the
independent dimension-five input. It is not a claim that the primitive loops
generated `C5`. Setting `G=g=0` recovers the earlier phase result with its
minus-imaginary sign.

## Independent executable and controls

Run `python -I -B -O independent_real_audit.py` in the root pinned environment.
The standalone [program](independent_real_audit.py) imports no project
calculations and uses explicit guards rather than removable assertions. Its
[receipt](INDEPENDENT_REAL_RESULTS.json) records **30 positive checks and 22
adverse or witness controls**.

The scalar part proves (A1)-(A3) as exact polynomials in the real components of
an arbitrary complex two-component scalar and its independent first jet. It
also varies a general cubic polynomial in `rho` to check the potential chain
rule. The recorded finite witness has `rho=10`, `partial rho=-8`, `J=4`,
`g=t=5` and residual `440`; both current and norm gradient are nonzero.

The fermion part represents fields and derivative jets by independent exterior
generators, uses noncommuting complex Hermitian flavor matrices, and varies
the literal symmetric kinetic action for all four chiral sigma matrices.
It checks its Hermiticity, the four Yukawa interactions, independent scalar
and lepton probes, and the actual odd-field `C5` monomial including its
Hermitian conjugate. Deliberate errors include half radial normalization,
wrong scalar-current sign, same Higgs phase, reversed Yukawa ordering,
deleted probe/radial sources, adjoint in place of transpose, omitted scalar
contribution to `C5`, and removal of the half. Nonzero witness guards are
included explicitly in the adverse-control count; these counts do not enlarge
the scientific scope.

The audit agrees with the selected classical representative (A4)-(A8). Earlier
three-point source terms, full retained-sterile operator reduction, other
insertions, quantum/composite-source transport, gauge-loop completion,
evanescent/subtraction definitions, complete physical running, finite matching
and one observable remain open. No new law, unification, empirical confirmation
or novelty claim follows from these action identities.
