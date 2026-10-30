# Real matter normalization with interactions retained

30 October 2026 (UTC). Input research commit
`5094584ebf402736b0bc5b597f6769f640178fd1`.
This is a first-order **classical action-coordinate calculation** using the
previous mixed and four-Higgs loop poles. It integrates no new loop and does
not complete a physical operator basis, quantum source map, running or matching.

## Conventions and chosen coordinate map

Use the [mixed pole](../../checkpoint_2026_10_26/mixed_yukawa/MIXED_SEED_POLE.md),
[four-scalar pole](../../checkpoint_2026_10_23/seed_yukawa/FOUR_SCALAR_POLE.md)
and [preceding phase map](../../checkpoint_2026_10_26/combined_matter/COMBINED_PHASE_TRANSPORT.md).
Write

```text
phi=tilde H, rho=phi^dagger phi, T=(D_mu phi)^dagger D^mu phi,
z_mu=phi^dagger D_mu phi,
J_mu=i[z_mu-z_mu*], u_mu=partial_mu rho=z_mu+z_mu*,
R=Y_nu* M Y_nu^dagger, K=C7 R,
G=(K+K^dagger)/2, A=(K-K^dagger)/(2i),
g=Tr G=Re Tr K, t=Tr A=Im Tr K,
K3=Tr(C7 Y_nu* M M^dagger M Y_nu^dagger),
kappa=1/(16 pi^2 epsilon).
```

`G,A` are Hermitian constant flavor matrices with mass dimension `-2`;
`g,t` are real. Lorentz contractions use `(+---)`. The reference scalar
action is `T-V(rho)`, with `V=m_H^2 rho+lambda rho^2`; the signed mass
parameter and this quartic normalization agree with the earlier action map.
Use the symmetrized lepton kinetic action
`bar L i<->slash D L/2`, which differs from `bar L i slash D L` by a boundary.

The real portions of the two known derivative poles are

```text
L_mixed,real^UV = (kappa/2)[-rho bar L G i<->slash D L
                           +J_mu bar L G gamma^mu L],
L_scalar,real^UV=kappa g[(u_mu u^mu-J_mu J^mu)/2-2rho T].   (1)
```

The covariant notation is an algebraic extension of the zero-external-gauge
action. It is not an independent calculation of external-gauge loop vertices.
The scalar identity follows from `z=(u-iJ)/2`, so no component of the
Hermitian-conjugate derivative placement has been dropped.

Define the direction of the substitution explicitly: **old fields equal new
fields plus the displayed delta**. First choose its real part,

```text
delta_r L=(kappa/2)rho G L, delta_r bar L=(kappa/2)rho bar L G,
delta_r phi=kappa g rho phi, delta_r phi^dagger=kappa g rho phi^dagger.
```

Direct product rules give

```text
delta_r L_kin,L=(kappa/2)rho bar L G i<->slash D L,
delta_r T=kappa g[2rho T+u_mu u^mu],
delta_r rho=2kappa g rho^2.                              (2)
```

Consequently this particular coordinate map removes the displayed
`rho bar L G i<->slash D L` and `rho T` terms, but leaves

```text
L_derivative,res^UV=(kappa/2)J_mu bar L G gamma^mu L
       +kappa g[3u_mu u^mu/2-J_mu J^mu/2]
 = (kappa/2)J_mu bar L G gamma^mu L
       +2kappa g[(partial rho)^2-|phi^dagger D phi|^2].    (3)
```

This is a specified local normalization step, not complete canonical
normalization on all scalar backgrounds. The residual derivative structures
must remain in any subsequent EOM/IBP or geometric analysis.

The previous phase map can be added at this order. Its composition with
(2) is

```text
delta L=kappa rho F L,      F=G/2+iA,
delta bar L=kappa rho bar L F^dagger,
delta phi=kappa rho b phi,  b=g+it=Tr K,
delta phi^dagger=kappa rho b* phi^dagger.                 (4)
```

For `epsilon=i sigma^2`, `phi=epsilon H*`, the inverse relation is
`H=-epsilon phi*`; hence `delta H=kappa rho b* H`.
The imaginary part of (4) cancels precisely the two selected phase-current
terms already calculated. The two steps commute through first order in
`kappa`; their cross-terms would be second order.

## Potential, Yukawa interactions and independent probes

The radial part now changes the potential, which the earlier pure phase
left invariant:

```text
delta L_potential=-2kappa g rho^2 V'(rho)
                =-2kappa g m_H^2 rho^2-4kappa g lambda rho^3.
```

Including the unchanged primitive mass-cubic pole, the displayed scalar
potential polynomial becomes

```text
2kappa[Re K3-g m_H^2]rho^2-4kappa g lambda rho^3.          (5)
```

This coefficient is in the chosen field coordinates; it is not a Higgs
mass, physical quartic, sextic prediction or total counterterm.

For the full declared reference Yukawa action

```text
L_Y=-bar L Y_nu phi N_R-bar L Y_e H e_R
    -bar Q Y_u phi u_R-bar Q Y_d H d_R+h.c.,
D=G/2+g I, B=A-t I,
```

the combined substitution (4), leaving quarks and singlet fermions fixed,
produces

```text
delta L_Y=-kappa rho bar L[D-iB]Y_nu phi N_R
          -kappa rho bar L[D-i(A+t I)]Y_e H e_R
          -kappa rho(g+it)bar Q Y_u phi u_R
          -kappa rho(g-it)bar Q Y_d H d_R+h.c.            (6)
```

In particular, even the one-active-flavor neutrino phase cancellation of
the previous checkpoint does not erase the radial Yukawa variation. Setting
`A=t=0` does not remove (6) when `G` is nonzero.

The independent probe action
`L_probe=bar eta L+bar L eta+j^dagger phi+phi^dagger j`
also changes:

```text
delta L_probe=kappa rho[bar eta F L+bar L F^dagger eta
                      +b j^dagger phi+b* phi^dagger j]. (7)
```

No Yukawa or probe source in (6)-(7) is silently set to zero. Other
interactions and composite sources in a fuller reference action must also
be transformed.

## Conditional lower-order insertion

If an **independently supplied** symmetric dimension-five interaction
is present,

```text
L5=(1/2)bar X C5 X^c+h.c., X=phi^dagger L, C5=C5^T,
```

then `delta X=kappa rho phi^dagger U L`, where
`U=F+b*I=D+iB`. Its barred and charge-conjugated slots give

```text
delta L5=(kappa rho/2)bar X [
     D C5+C5 D^T-i(B C5+C5 B^T)]X^c+h.c.                 (8)
```

The transpose on the second slot is essential. The matrix in brackets
is symmetric. Its coefficient has mass dimension `-3`, so this is a
dimension-seven contribution conditional on the independent `C5`.
It is not a coefficient newly generated by the previous seed loop.

## Selected primitive matter ledger and limits

At the declared order, the primitive four-point family is represented by
the residual derivatives (3), potential (5), full source variations
(6)-(7), and the previously integrated four-lepton contact

```text
L_4L^UV=(kappa/2)C7_alpha beta R_gamma delta
               P_alpha beta gamma delta+h.c.,
P_alpha beta gamma delta=sum_(u,v)
 (bar l_alpha^u epsilon_spin bar l_beta^{vT})
 (l_gamma^{uT} epsilon_spin l_delta^v).
```

The contact is unchanged by (4) at this order: substituting the
one-loop field shift into a one-loop pole contributes only at `kappa^2`.
Substitution into the original `C7` seed is quadratic in that seed and is
outside the declared one-seed sector. Equation (8) is added only when
that independent lower-order input belongs to the reference action.
The cancelling action has the opposite of the collected pole in the
same coordinates, with the field-counterterm convention fixed consistently.

This ledger collects the displayed **primitive proper four-point matter
inputs** under one classical field substitution. The earlier three-point
source packet, its reducible sterile bridges, and renormalizable field,
mass and coupling counterterm inputs remain separate. Additional
interactions, lower-order insertions, a complete retained-sterile EOM/IBP
and evanescent basis, full quantum 1PI/composite-source transport,
finite matching, physical running and an observable remain unresolved.
No sterile degree of freedom is integrated out. Neither a residual
coordinate coefficient nor the algebraic removal of a kinetic term is
evidence of unification, new fundamental mathematics or empirical confirmation.

Run `python -I -B -O verify_real_matter.py` in the root pinned symbolic
environment. The standalone checker uses exact scalar first jets,
literal exterior Grassmann lepton kinetics, full weak-doublet Yukawa
and probe slots, and the dimension-five insertion. Its adjacent receipt
states the actual positive and deliberately incorrect controls; it imports
no project calculations and uses explicit failures active under optimization.
