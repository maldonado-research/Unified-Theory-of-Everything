# External-field input to the derivative-seed counterterm

15 October 2026 (Pacific; 16 October UTC). Input commit: `0770c304c7e05f213dd796d3fd5a194107f38858`.

This note uses the independently integrated [pure-Yukawa two-point poles](YUKAWA_TWO_POINT.md). It derives their contribution to the coefficient of the original derivative seed, with the scalar derivatives retained. This is one counterterm and RG input. It does not supply the proper dimension-seven vertex, its operator mixing, a complete physical seed beta function or a finite observable.

## Conventions and literal bare-field substitution

Use `d=4-2 epsilon`, `kappa=1/(16 pi^2 epsilon)`, the unbroken theory and the same literal action as the two-point calculation. Define

```text
H_L = Y Y^dagger,  H_R = Y^dagger Y,  T = Tr H_R,
X = phi^dagger L,  phi = tilde H,
P0 = (i slash partial)^2 = -Box,
L_seed = (1/2) bar X_alpha C7_alpha beta P0 X_beta^c + h.c.,
C7 = C7^T.
```

The half is the coefficient convention of the seed with its identical fermionic slots. The loop-derived field poles in this coupling sector are

```text
delta Z_L = -kappa H_L/2,
delta Z_phi = -kappa T.
```

The constant field substitutions are `L0=(1+delta Z_L/2)L` and `phi0=(1+delta Z_phi/2)phi`. There are two lepton factors and two scalar factors in the seed. Applying them to the displayed literal action gives a coefficient shift

```text
Delta C7_fields = (delta Z_L C7 + C7 delta Z_L^T)/2
                 + delta Z_phi C7
               = -kappa r7,
r7 = (H_L C7 + C7 H_L^T)/4 + T C7.
```

Thus the parameter pole needed to cancel this field-induced action shift is

```text
delta C7_parameter | fields = +kappa r7.
```

This pole is a defined part of the parameter counterterm, not the sum of all one-loop seed contributions. The transpose on the second flavor slot is essential; replacing it by a dagger is wrong for a general complex Yukawa matrix. Symmetric `C7` gives symmetric `r7`. No extra Majorana half multiplies this coefficient.

Constant field matrices commute with `P0`, but the scalar field itself varies. The operator is `P0(phi^T L^c)`, including its scalar and mixed derivatives. Replacing it by `phi^T P0 L^c` would change the seed. The checker implements the whole product in a finite Fourier/exterior algebra, checks the unintegrated and integrated substitutions, and rejects that omission on nonconstant scalar backgrounds. Its two examples check the component transport and derivative implementation; the general matrix formula follows directly from the linear bare-field substitution above.

## Canonical-epsilon conversion to an RG input

Write the dimensionful renormalized seed coefficient with mass dimension `-3`, and use

```text
Y0 = mu^epsilon (Y + delta Y),
C7,0 = mu^(2 epsilon) (C7 + delta C7).
```

In `d` dimensions the displayed operator has dimension `7-4 epsilon`, so its coefficient has dimension `-3+2 epsilon`. These exponents fix the conversion; a simple pole is not directly its four-dimensional beta function.

Treat `Y` and its conjugate as independent formal coordinates while differentiating, then restore conjugation. The pole `r7` has canonical weight four: two Yukawa factors of weight one and one `C7` of weight two. Differentiating the bare coefficient at fixed `C7,0` therefore gives

```text
16 pi^2 beta_C7 | fields = (4-2) r7
                        = (H_L C7 + C7 H_L^T)/2 + 2 T C7.
```

Equivalently, with `gamma=mu d(ln Z)/dmu /2`,

```text
16 pi^2 gamma_L = H_L/2,  16 pi^2 gamma_phi = T,
beta_C7 | fields = gamma_L C7 + C7 gamma_L^T + 2 gamma_phi C7.
```

The last equation expresses the external-field contribution only. The [proper renormalizable Yukawa vertex](YUKAWA_VERTEX.md) is a different amplitude; its vanishing pure-Yukawa pole cannot be transferred to the higher-derivative seed vertex. Proper seed graphs and the complete retained-sterile operator/source map are still required before a physical coefficient vector can run.

## Relation to the preceding source projection

The preceding [active field-content projection](../derivative_projection/ACTIVE_SOURCE_PROJECTION.md) transported a selected sourced action and obtained a raw half-normalized seed-shaped projection `5 kappa K_S/3`. That was a representative-dependent action contribution with different specified inputs. The present `r7` follows from new renormalizable loop poles and the bare seed's own field factors. Neither selected quantity is a total `C7` counterterm. They must not be equated, summed into a physical result without a consistent complete action map, or used to reinstate the withdrawn all-odd physical-running claim.

## Executable checks

Run `python -I -B -O verify_seed_field_transport.py` in the pinned symbolic environment. The standalone program imports no repository modules and uses explicit failures rather than removable assertions. [Its receipt](SEED_FIELD_RESULTS.json) records 17 positive checks and 15 rejected mutations: varying-scalar substitutions, integrated action coefficients, flavor symmetry and covariance, scalar multiplicity, the seed half, the required flavor transpose and the epsilon-to-RG factor. This is algebra using independently computed loop inputs; it is not another loop integration or independent human peer review.
