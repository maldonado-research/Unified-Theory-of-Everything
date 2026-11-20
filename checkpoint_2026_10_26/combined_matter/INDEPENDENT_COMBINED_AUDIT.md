# Independent audit of the correlated scalar and lepton phase change

26 October 2026 (Pacific). Input commit: `9b1e18450c0c1e18fb6f917105368b05ca123df2`.

This is an independently implemented action-variation audit of the selected phase-current contributions from the preceding mixed and four-scalar poles. It imports no producer code, does not recompute either loop, and does not establish a complete physical or quantum effective-action map. All displayed statements are first order in the loop coefficient and linear in the original seed.

## Field variation and kinetic terms

Use `phi=tilde H`, `epsilon=i sigma^2=[[0,1],[-1,0]]`, so `phi=epsilon H*` and the inverse is **`H=-epsilon phi*`**. Write

```text
R=Y_nu* M Y_nu^dagger,  K=C7 R,
A=(K-K^dagger)/(2i),    t=Tr A,    rho=phi^dagger phi,
J_phi^mu=i[phi^dagger D^mu phi-(D^mu phi)^dagger phi].
```

`A` is Hermitian and `t` real. Its matrix entries have dimension minus two. The independently integrated zero-gauge loop actions contain the selected terms

```text
kappa (partial_mu rho) bar L A gamma^mu L
 + kappa t (partial_mu rho) J_phi^mu,
kappa=1/(16 pi^2 epsilon).
```

Take the simultaneous infinitesimal changes

```text
delta L=i kappa rho A L,    delta bar L=-i kappa rho bar L A,
delta phi=i kappa rho t phi, delta phi^dagger=-i kappa rho t phi^dagger.
```

The inverse-Higgs convention gives `delta H=-i kappa rho t H`. Since `rho` is a gauge singlet and `A,t` are constant in spacetime, literal product differentiation gives

```text
delta(bar L i slash D L)=-kappa (partial_mu rho) bar L A gamma^mu L,
delta[(D_mu phi)^dagger D^mu phi]=-kappa t (partial_mu rho)J_phi^mu,
delta rho=0.
```

These variations cancel the two displayed phase-current terms in the first-order local action. They do not cancel the other Hermitian-current, four-fermion, potential or source terms. The covariant product rule here is an action identity; zero-gauge loop results alone do not establish the corresponding external-gauge loop vertices. No statement about a complete quantum Jacobian or Legendre-source transformation follows from this classical action variation.

## Literal interaction and probe changes

For the reference renormalizable interactions, with weak/spin/color contractions suppressed,

```text
L_Y = -bar L Y_nu phi N_R -bar L Y_e H e_R
      -bar Q Y_u phi u_R -bar Q Y_d H d_R + h.c.,
```

apply the two field changes to the fields actually present, without substituting free equations. The first-order interaction variation is

```text
delta L_Y = +i kappa rho bar L (A-t I)Y_nu phi N_R
            +i kappa rho bar L (A+t I)Y_e H e_R
            -i kappa rho t bar Q Y_u phi u_R
            +i kappa rho t bar Q Y_d H d_R + h.c.       (1)
```

The opposite phases of `phi` and `H` are essential. Matrix products retain their displayed order. `N_R`, `e_R`, `Q`, `u_R` and `d_R` were not changed in this particular coordinate choice; possible further field changes constitute additional work. These terms are induced sources, not new independent loop integrations.

For an independent scalar probe `j` and lepton probe `eta`, define

```text
L_probe=j^dagger phi+phi^dagger j+bar eta L+bar L eta.
```

Their direct variation is

```text
delta L_probe=i kappa rho t(j^dagger phi-phi^dagger j)
             +i kappa rho(bar eta A L-bar L A eta).     (2)
```

These expressions retain the probes as independent sources. A full action can contain further interactions and lower-order operators, all of which must also be varied. Source-dependent terms cannot be discarded merely because the displayed kinetic phase currents cancel.

## An independently supplied lower-order insertion

If the reference action independently includes a symmetric dimension-five input

```text
L5=(1/2) bar X C5 X^c+h.c., X=phi^dagger L, C5=C5^T,
```

set `B=A-t I`. The actual composite variations are `delta X=i kappa rho B X`, `delta bar X=-i kappa rho bar X B` and `delta X^c=-i kappa rho B^T X^c`. Thus

```text
delta L5=-(i kappa rho/2) bar X(B C5+C5 B^T)X^c+h.c.    (3)
```

This is a dimension-seven contact induced by the field change acting on an **independently supplied** dimension-five input. Its presence is conditional on that input and its coefficient convention; it is not an independently calculated seed-loop amplitude. The checker expands the literal two-flavor, two-weak-component odd-field composites and rejects omission of the scalar phase and hence of `t I`.

## A rank-safe source condition

The neutrino term in (1) vanishes precisely when `(A-t I)Y_nu=0`. Let `P` be the orthogonal projector onto `im Y_nu`, `Q=1-P`, and `r=rank Y_nu`. Since

```text
R Q=Y_nu* M Y_nu^dagger Q=0,
Q K Q=Q K^dagger Q=0, hence Q A Q=0,
```

the source condition is equivalent to `A P=t P`. Hermiticity then gives `P A=t P`, and the cross blocks vanish; combined with `Q A Q=0`, this implies `A=t P`. Taking its trace gives `(r-1)t=0`:

- If `r>=2`, the condition requires `t=0` and `A=0`.
- If `r=1`, it is equivalent to `A=t P`, which can be nonzero.
- If `r=0`, `Y_nu=0` gives `R=K=A=0`.

This proof uses no inverse sterile mass and remains valid for singular Yukawa matrices. Under the active flavor rotation of the reference action, `P` and `A` transform by the same unitary similarity, so the condition is coordinate covariant. An explicit rank-one nonzero witness and rank-zero/rank-one/full-rank support checks accompany the proof. These are conditions on the selected induced neutrino source; they are not a physical cancellation theorem or an observable.

For one active flavor, automatically `A=t`, so `delta X=0` and the neutrino interaction and (3) are unchanged. The charged-lepton and quark terms in (1), and probes in (2), can still be nonzero. This one-flavor statement has no force for a general multiflavor theory.

## Reproduction and limits

Run `python -I -B -O independent_combined_audit.py` in the pinned SymPy environment. It uses independent exterior Grassmann generators for lepton, sterile, quark and probe fields; it directly varies the literal Yukawa, kinetic first jets, source probes and optional `C5` composite, with complex matrices and singular-Y witnesses. Four explicit left-Weyl sigma components check the lepton kinetic variation. The adjacent receipt records exact controls and rejected mutations. It imports no repository modules and uses explicit failures, rather than assertions removed by optimization.

The finite component checks exercise the analytic first-order identities above; they are not a new loop calculation or independent human peer review. Complete operator mixing, other interactions, source and gauge-vertex transport at the quantum level, a dimensional/evanescent prescription, finite threshold matching and an observable remain open.
