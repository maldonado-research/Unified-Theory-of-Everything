# Complex Majorana masses and covariant source recurrences

16 October 2026. Continuation of candidate
`4ed8560e019dfe74ed6c5178b0f990fe3ccb3008`. The original sourced packet
used a real diagonal retained Majorana mass. This note extends its
ordered chiral recurrences, integrated field-map rows and coefficient
transport to a constant complex symmetric mass. It does not integrate
out a sterile field or calculate a new loop or physical beta function.

## Opposite chiral equations contain the adjoint

Keep the original four-component mass orientation:

```text
L_mass=−(1/2)[bar N_L M N_R+bar N_R M^dagger N_L],
M=M^T, N_R=N_L^c, A=i slash partial.
```

The half is the Majorana action normalization; varying both constrained
slots removes it from the Euler equations. Including the actual
variational interaction sources gives

```text
e_R=A N_L−M N_R−S_R,
e_L=A N_R−M^dagger N_L−S_L^N.                          (M1)
```

For symmetric `M`, `M^dagger=M*`, so the two equations have the required
charge-conjugation relation when the sources do. The complex mass in
the holomorphic left-Weyl mass matrix would instead be `M*`, as declared
in the earlier normalization calculation. Those two conventions must not
be identified by keeping the same unconjugated matrix name.

`M` is a constant parameter. Hence `A` commutes with mass multiplication,
but matrix products need not commute with one another. Successive exact
substitutions in (M1) give

```text
R=M N_R+S_R,
T=M M^dagger N_L+M S_L^N+A S_R,
U=M M^dagger M N_R+M M^dagger S_R+M A S_L^N+A^2 S_R,

A N_L=R+e_R,
A^2 N_L=T+M e_L+A e_R,
A^3 N_L=U+M M^dagger e_R+M A e_L+A^2 e_R.              (M2)
```

Derivatives act on complete source composites. No source, residual,
mass inverse or on-shell condition is deleted. Since the sterile fields
and these composite sources are gauge singlets, `A^2=−partial^2` is
appropriate here. A field-dependent mass would introduce additional
derivative terms and is outside this calculation.

The selected sourced representative keeps its previous displayed form
with `R,T,U` replaced by (M2). In its off-shell residual replace
`M^2 e_R` by `M M^dagger e_R`. These are the only real-mass powers in
the sterile recurrences; the earlier Clifford, scalar, gauge-current
and covariant-IBP identities do not depend on a real mass. The other
ordered factors `M e_L`, `M S_L^N` and `F M` remain in their given order.
For real symmetric `M`, (M2) reduces to the original powers `M^2,M^3`.

## Integrated field-map rows

Use the earlier neutral row composites

```text
W_i=sum_alpha a F_alpha i bar L_alpha^a phi^a,
V_i=sum_alpha a F_alpha i bar S_Lalpha^a phi^a,
K0 W=−i(partial_mu W)gamma^mu,
c=2m_H^2/3+lambda rho/3.
```

After moving derivatives off the residuals, the sterile rows are

```text
bar x=(5/6)(W M M^dagger+K0^2 W)+c W
       +(1/3)[(K bar S_L)phi]_F−(1/6)[bar L S_phi]_F
       −K0 V/6−[bar L Fcal_mu nu phi sigma^{mu nu}]_F/12,
bar y=(5/6)(K0 W)M−V M/6.                              (M3)
```

Their adjoints define left-chiral `x` and right-chiral `y`. The
substitution remains `N_old=N−kappa(x+y+(x+y)^c)` with no extra half.
For example the adjoint of the first mass term in `bar x` is
`M M^dagger W^dagger`, rather than `M^2 W^dagger`. This follows from
the actual residual (M2); it is not an independent convention for the
mass ordering. The other field shifts retain their original form.

## Unitary flavor coordinates and coefficient projections

Let the old fields be `L=V L'`, `N_L=U N_L'`, with constant unitary
`V,U`, and keep the scalar convention fixed. Charge conjugation requires
`N_R=U* N_R'`. Direct substitution into the action fixes

```text
Y'=V^dagger Y U*,
M'=U^dagger M U*,
C_seed'=V^dagger C_seed V*,
F'=C_seed' Y'^*=V^dagger F U.                           (M4)
```

In particular the Majorana mass has a congruence, not a unitary
similarity. It remains symmetric. The sources and residuals transform as

```text
S_R'=U^dagger S_R,       (S_L^N)'=U^T S_L^N,
e_R'=U^dagger e_R,       e_L'=U^T e_L,
R'=U^dagger R, T'=U^dagger T, U_source'=U^dagger U_source.
```

The last line uses `U_source` for the third source composite to
distinguish it from the unitary matrix. Contracting with `F'` then
transforms the active output with `V^dagger`. The integrated rows in
(M3) transform with `U` and `U*`, respectively; their adjoints and
Majorana completion therefore transform consistently too.

The adjacent [active projection](ACTIVE_SOURCE_PROJECTION.md) and the
retained tensor-source sectors involve the ordered matrices

```text
K=F Y^T,
K_M=F M M^dagger Y^T,
G=F M Y^dagger,
H_N=Y^dagger F.                                        (M5)
```

Under (M4), `K,K_M` transform as `V^dagger K V*` and
`V^dagger K_M V*`; their symmetric scalar/Box projections have the
same congruence. `G` transforms as `V^dagger G V`, so its Hermitian
kinetic projection is also covariant. The retained tensor matrix obeys
`F'M'=V^dagger(F M)U*`. Finally
`H_N'=U^T H_N U`; it is symmetric for the actual symmetric seed
`F=C_seed Y*`. These laws also check the source symmetries used in the
[curvature transport](CURVATURE_DERIVATIVE_BASIS.md).

The real-basis `F M^2 Y^T` in the active projection is thus generalized
to `F M M^dagger Y^T`. Replacing every mass by an ordinary complex
matrix while leaving `M^2,M^3` unchanged violates these laws. These are
coefficient-coordinate checks within the selected action, not physical
predictions or a proof that every EFT coefficient is known.

## Reproduction and scope

Run `python -I -B -O verify_complex_mass_covariance.py` with root pinned
SymPy. The standalone exact checker uses three active and two sterile
flavors, three generic off-shell momenta, two complex rational unitary
sterile changes, a complex active change, and generic, singular and zero
complex symmetric masses. It constructs arbitrary chiral fields and
independent nonzero interaction sources, forms the literal residuals,
and tests all three recurrences, their chiral covariance and the actual
`Y^T X^c,Y^dagger X` source orientations. No mass is inverted.

Checks of the integrated rows use independent derivative/source row
coefficients; they test their transformation law rather than redoing
every graded action variation. Deliberate wrong adjoints, real-mass
powers and a Majorana similarity in place of congruence must fail.
Adjacent JSON records the exact controls. Failures remain active under
optimization and no repository calculations are imported.

This algebra permits a general complex symmetric retained mass in the
selected source formulas. It supplies no quantum functional Jacobian,
subtraction prescription, finite matching, physical running or observable.
