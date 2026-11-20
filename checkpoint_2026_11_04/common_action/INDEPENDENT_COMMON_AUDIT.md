# Independent audit of the common action and neutral source map

4 November 2026 (Pacific). Input research commit `41f9b709c0fafdefdc9e1b990caa5421ee99f7dc`.

This audit was derived from the preceding counterterm, matter-pole and bare-field notes before comparing the new producer's implementation. It checks a first-order classical action convention and its explicit neutral quadratic projection. It evaluates no loop and does not supply a physical running vector, finite matching or a vacuum prediction. Independent AI-assisted derivation and computation are not independent human peer review.

## Signs and the chosen map

The input [electroweak packet](../../checkpoint_2026_10_15/interacting_reduction/ELECTROWEAK_REDUCTION.md) is already a **cancelling action**, `+kappa P`. In contrast, the [matter input](../../checkpoint_2026_10_30/real_matter/REAL_MATTER_TRANSPORT.md) is a **loop pole**, `+kappa Q`. The common cancelling action therefore contains `+kappa P-kappa Q`, together with the separately supplied renormalizable action counterterms. The sign of the matter normalizer reverses when it is used on this cancelling action.

Let `vC` be the specified differential normalizer of the original three-point packet and `w` the positive-loop matter map of the preceding checkpoint. The selected common substitution is

```text
old fields = common fields + kappa(vC-w),
bare fields = common fields + kappa(z+vC-w),
```

where `z` denotes the standard renormalizable bare-field factors. Products of a one-loop counterterm and a one-loop substitution start at second loop order. The variation of an independently present **tree** operator is a different contribution and cannot be dropped on that argument.

The source-reducing direction `vS` is an alternative map of the original packet. To use its table, the bridge is `vC-vS-w`; it is not `vC-w` applied again on top of `vS`. The full earlier bridge is supplied by the [action-map calculation](../../checkpoint_2026_10_16/action_completion/FIELD_MAP_COMPOSITION.md). This audit reconstructs its neutral active-field component and does not independently re-expand every interacting `vS` term.

## Explicit neutral action

Set `phi=(u,0)` with real constant `u`, zero gauge background and zero fermion background. The scalar need not be stationary. Use left-Weyl fields `q=(nu_L,N_L)` and

```text
F=C7 Y*, P=p^2+m_H^2, H=u P F,
K=C7 Y* M Y^dagger=G+iA,
G=G^dagger, A=A^dagger, g=Tr G, t=Tr A, b=g+it,
Z=[0 H; H^dagger 0],
D_s=Z+diag(u^2 G,0),
R_z=diag(-H_L/2,-H_R^T),
H_L=Y Y^dagger, H_R=Y^dagger Y.
```

The original six placements give the coefficient `5/6+1/3-1/6=1` multiplying `p^2` and `1/2+1/2=1` multiplying `m_H^2`. There is no additional Majorana half on the off-diagonal block `H`. The positive matter loop has active kinetic coefficient `-u^2 G`; its cancelling action has `+u^2 G`. The total selected neutral kinetic coefficient before the new substitution is consequently

```text
K_neutral=I+kappa(D_s+R_z).
```

The fermion and background substitutions from the common map are

```text
q_old=(I+kappa J)q,
J=-D_s/2-i u^2 diag(A,0),
delta phi=-kappa u^3 b (1,0).
```

Since `J+J^dagger=-D_s`, direct kinetic substitution gives

```text
K_common=I+kappa R_z.
```

The standard renormalizable counterterms remain explicit. This result does not say that the entire quadratic action has become canonically normalized. Removing `R_z` would be an additional field convention, with corresponding tree-operator and parameter transport.

For the holomorphic mass action `-q^T epsilon Mcal q/2+h.c.`, define

```text
Mcal=[0 B; B^T D], B=uY*, D=M*.
```

The scalar background shift contributes `delta B_scalar=-u^3 b*Y*`. Independent variation of both odd fermion slots therefore gives

```text
delta Mcal=J^T Mcal+Mcal J
          +[0 delta B_scalar; delta B_scalar^T 0].
```

Its explicit blocks are

```text
delta M_LL=-(H* B^T+B H^dagger)/2,
delta M_LN=-H*D/2
           -u^3[G^T/2+iA^T+(g-it)I]Y*,
delta M_NN=-(H^T B+B^T H)/2.
```

The transpose is not an adjoint. The scalar-background term is essential even though the original neutral calculation is a fermion Hessian. These are coordinate-induced coefficients in the selected action, not independently predicted masses. Complex symmetric and singular retained masses are allowed; no mass inverse enters this calculation.

For a fixed row source `S q+h.c.`, the new coefficient is `S+kappa S J`. The checker varies actual independent odd probe fields. Dropping this source transport would change response functions even if a spectral comparison happened to agree.

## A supplied dimension-five coefficient

If an independent symmetric `C5` is included as `L5=+bar X C5 X^c/2+h.c.`, its holomorphic term is **positive**, `+(u^2/2)nu_L^T epsilon C5*nu_L+h.c.`. The declared mass convention has a negative half, so its mass block is `M5=diag(-u^2 C5*,0)`. Applying exactly the same fields and background gives

```text
delta M5=J^T M5+M5 J
         +diag(+2u^4 b* C5*,0).
delta M5_LN=+(u^2/2) C5* H.
```

The positive mixed **mass block** in the second line enters the negative-half mass action. It therefore agrees with the negative mixed **interaction** contact of the producer's common action. Retaining only the radial matter change would miss it. The checker explicitly constructs the positive `L5` bilinear and rejects the opposite mass-block sign before varying its two slots. Its `P0` and `m_H^2` pieces have different canonical operator dimensions, as explained in the [common-action derivation](COMMON_ACTION_MAP.md). The `C5` input is independent, and these mixed-input variations are not new loop-generated coefficients. Under strict counting `C7=O(Lambda^-3)`, `C5=O(Lambda^-1)` with fixed masses, their product is outside a truncation at `Lambda^-3`; retaining them conditionally requires stating that different counting choice.

## Renormalizable action coefficients and the tree seed

The [pure-Yukawa input](../../checkpoint_2026_10_16/yukawa_counterterms/YUKAWA_TWO_POINT.md) fixes

```text
z_L=-H_L/4, z_NR=-H_R/2, z_phi=-T_Y/2,
T_Y=Tr H_R, Q_M=Tr(H_R M^dagger M).
```

The parameter residues

```text
rY=H_L Y/4+Y H_R/2+T_Y Y/2,
rM=(H_R^T M+M H_R)/2,
rm_H^2=m_H^2 T_Y-2Q_M
```

obey the action identities

```text
rY+z_L Y+Y z_NR+z_phi Y=0,
rM+z_NR^T M+M z_NR=0,
rm_H^2+2m_H^2 z_phi=-2Q_M.
```

Thus a nonzero parameter pole cannot simply be appended as another proper action vertex. In particular the proper Majorana mass coefficient is zero in this restricted input sector despite nonzero `rM`.

The standard field factors act on the original **tree** derivative seed as

```text
delta C7_fields=-kappa r7,
r7=(H_L C7+C7 H_L^T)/4+T_Y C7.
```

Their inverse would produce `+kappa r7`. This is first loop order and linear in the seed. It is not a product of two one-loop poles. The independently identified parameter subset `+kappa r7` cancels the displayed bare-field variation in that convention; neither term alone is the total physical seed counterterm. The audit constructs the literal positive tree-seed bilinear, checks its negative-half mass-block conversion, then varies both odd slots and both scalar factors on its constant background. It does not itself test nonconstant scalar derivatives in the seed; the earlier [seed-field program](../../checkpoint_2026_10_16/yukawa_counterterms/verify_seed_field_transport.py) checks those separately.

## Reproduction and limits

Run `python -I -B -O independent_common_audit.py` in the pinned root environment. The standalone program imports no project calculations and uses explicit failures active under optimization. Its adjacent [receipt](INDEPENDENT_COMMON_RESULTS.json) gives the final counts. The unrestricted scalar first-jet part checks the loop-versus-cancelling sign directly. The neutral part uses two complex flavor examples, including a singular complex symmetric mass, actual exterior fermion components, all four chiral spin matrices and independent probes.

Adverse controls include the wrong matter sign, a halved mixed block with an explicit nonzero sensitivity witness, omission of the scalar background shift or probe, incorrect mass adjoints, deletion of the conditional `C5` mixed contact, confusion of parameter and action mass poles, deletion of tree-seed inverse transport, and adding the source and normalizing maps instead of their bridge. Their number is implementation coverage rather than a count of independent physical validations.

This verifies a selected common classical convention and a specified projection. Full derivative/EOM and evanescent reduction, remaining couplings and amplitudes, quantum composite-source transport, complete physical running, finite matching and an observable remain open. A constant scalar background is not a sterile threshold integration, a vacuum condition or empirical evidence for unification.
