# Curvature of the selected scalar derivative packet

30 October 2026 (UTC). This is a first-order geometric diagnostic of the previously calculated real four-scalar pole and its radial field substitution. It evaluates no new loop. The geometry is the geometry of scalar field coordinates in a selected local action, not spacetime geometry or an empirical gravitational result. The divergent pole coefficient is treated as a formal perturbative parameter; this note does not construct a finite physical sigma model or assert positivity.

## Literal real coordinates and radial substitution

Use the zero-external-gauge scalar conventions of [the primitive four-scalar calculation](../../checkpoint_2026_10_23/seed_yukawa/FOUR_SCALAR_POLE.md). Define

```text
phi=(q0+i q1, q2+i q3)/sqrt(2),
rho=phi^dagger phi=q^Tq/2,
S=(partial_mu phi)^dagger partial^mu phi=(partial_mu q)^T(partial^mu q)/2,
J_mu=i[phi^dagger partial_mu phi-(partial_mu phi)^dagger phi],
g=Re Tr(C7 Y_nu* M Y_nu^dagger), h=kappa g.
```

Let `J0=diag([[0,-1],[1,0]],[[0,-1],[1,0]])`; this constant real matrix obeys `J0^T=-J0`, `J0^2=-I`. The current is `J_mu=q^T J0 partial_mu q`. The symbol `J0` is a scalar-coordinate complex structure, distinct from the current `J_mu`.

The real scalar derivative pole is

```text
Delta L_old=h[(partial rho)^2/2-J_mu J^mu/2-2rho S].       (1)
```

For `L_kin=(1/2)g_ij(q) partial_mu q_i partial^mu q_j`, its metric correction is

```text
g_old=I+h H_old,
H_old=q q^T-(J0 q)(J0 q)^T-(q^Tq)I.                     (2)
```

Under the actual coordinate substitution

```text
q_old=q+h xi(q), xi(q)=(q^Tq)q/2,
```

the pullback, with its full Jacobian retained, gives

```text
g_new=I+h H_new+O(h^2),
H_new=3q q^T-(J0 q)(J0 q)^T,
H_new-H_old=partial xi+(partial xi)^T.
```

Equivalently, `delta phi=h rho phi` changes the reference scalar kinetic term by `h[(partial rho)^2+2rho S]`. Thus the remaining derivative interaction is

```text
Delta L_new=(h/2)[3(partial rho)^2-J_mu J^mu].             (3)
```

The substitution also changes the potential, interactions and sources. Those changes must accompany this kinetic diagnostic; they are retained in the adjacent [matter transport calculation](REAL_MATTER_TRANSPORT.md). Equation (3) does not assert canonical normalization of the full interacting theory.

## Curvature calculation and obstruction

Fix the convention

```text
R^i_(jkl)=partial_k Gamma^i_(lj)-partial_l Gamma^i_(kj)
           +Gamma^i_(km) Gamma^m_(lj)-Gamma^i_(lm) Gamma^m_(kj),
R_ijkl=g_im R^m_(jkl).
```

The metric corrections (2)-(3) are quadratic in `q`. Both metrics equal `I` and have vanishing first derivatives at `q=0`. Their curvature there, to first order in `h`, is therefore determined by their second derivatives:

```text
R_ijkl(0)=(h/2)[H_il,jk+H_jk,il-H_ik,jl-H_jl,ik].         (4)
```

Differentiating either metric gives the same tensor,

```text
R_ijkl(0)=h{3(delta_ik delta_jl-delta_il delta_jk)
            +2 J0_ij J0_kl+J0_ik J0_jl-J0_il J0_jk}
          +O(h^2).                                      (5)
```

This follows directly from `partial_k partial_l(q_i q_j)=delta_ik delta_jl+delta_il delta_jk` and the analogous identity for `(J0q)_i(J0q)_j`. The metric difference in (2)-(3) is a coordinate Lie derivative of the flat metric, whose first-order Riemann tensor vanishes. No sampled field configurations are needed for these identities.

For arbitrary tangent vectors `X,Y`, contraction of (5) gives

```text
R(X,Y,X,Y)=3h[X^2 Y^2-(X dot Y)^2+(X^T J0 Y)^2]+O(h^2).
```

For unit orthogonal vectors the sectional curvature is consequently

```text
K(X,Y)=3h[1+(X^T J0 Y)^2]+O(h^2).                        (6)
```

A complex line, for example the `q0,q1` plane, has sectional curvature `6h`; a totally real plane, for example `q0,q2`, has `3h`. Because `J0` is orthogonal, `|X^T J0 Y|<=1`. The two cases illustrate the exact direction dependence rather than selecting special scalar backgrounds. With the declared sign convention, `Ricci_ij(0)=12h delta_ij` and `R_scalar(0)=48h`.

If `g` is nonzero, this first-order curvature is nonzero. An invertible local scalar reparametrization with no derivatives cannot turn this selected kinetic metric into a constant flat metric through this order: curvature transforms as a tensor. The radial substitution removes the displayed `rho S` coefficient while leaving the nonflat derivative packet (3). If `g=0`, the entire real scalar packet discussed here vanishes. The sign of `g` is unrestricted; none of these identities is a positivity statement.

This is a restricted obstruction to flattening the selected scalar kinetic metric. It does not obstruct source-preserving EOM or derivative-dependent field redefinitions into other operators, prove a complete independent EFT basis, or show that other loop and operator sectors cannot change the total curvature. Fermions, retained sterile fields, lower-order operator insertions, probe sources, gauge completion, subtraction conventions and finite matching must still be included in a full physical calculation.

## Reproduction

Run `python -I -B -O verify_scalar_geometry.py` in the root pinned environment. The standalone checker imports no project calculations and uses explicit guards that remain active under optimization. Its [receipt](SCALAR_GEOMETRY_RESULTS.json) records **29 positive check groups and eight adverse controls rejected**. The groups cover literal complex scalar first jets, the full radial pullback, all 256 curvature entries from both metrics, curvature symmetries and Bianchi identities, arbitrary-plane polynomial contraction and the two sectional curvatures. Component equality counts in the receipt describe implementation coverage, not independent physical validations or evidence for the hypothesis.
