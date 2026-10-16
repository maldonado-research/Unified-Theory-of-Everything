# Charged-lepton and quark Yukawa transport

16 October 2026. A classical first-order completion of the spectator terms
left symbolic in the previous field maps. The input is candidate
`5f9889b76661faa155c681a5ecfac5d036398bc8`. This calculation adds no loop,
physical Wilson coefficient, or new operator-classification claim.

## Declared interactions and scalar convention

Add the Standard Model spectator Yukawa interactions to the reference action:

```text
L_spec=−bar L_alpha (Y_e)_alpha p H e_Rp
       −bar Q_r (Y_d)_r s H d_Rs
       −bar Q_r (Y_u)_r s phi u_Rs+h.c.,
phi=epsilon H*, epsilon=[[0,1],[-1,0]], H=−epsilon phi*.
```

Weak indices contract between each barred doublet and its scalar. Quark
color contracts within each monomial; summing color does not multiply a
coefficient by an additional factor of three. Yukawa matrices are arbitrary
complex matrices in the displayed order. The underlying quantum theory
still needs the other retained interactions and counterterms.

Write `A_N=i slash partial`, `F=C_seed Y_nu*`. The already specified
source-reducing substitution has

```text
delta_S L_alpha^a=−ell_alpha^a,
ell_alpha^a=sum_i F_alpha i [
   (i/3)gamma^mu(D_mu phi)^a A_N N_Li
   +(1/6)phi^a A_N^2 N_Li+(m_H^2/2)phi^a N_Li],
delta_S phi^a=−t_a*, delta_S H^a=epsilon_ab t_b,
t_a=(1/6)sum_alpha i F_alpha i bar L_alpha^a A_N N_Li.
```

Here `delta` denotes the coefficient of the common real first-order
parameter `kappa`. The even bilinear is conjugated in reverse fermion order:
`t_a*=(1/6)sum F_alpha i* overline(A_N N_Li)L_alpha^a`.
All right sides use new fields. No additional mixed-operator or Majorana
half is inserted. Quarks and charged singlets have zero field shifts.
The sterile substitution does not enter `L_spec` directly.

## Direct first-order action substitution

Expanding the displayed action in both its scalar and doublet slots gives

```text
delta_S L_spec =
   +overline(ell_alpha)(Y_e)_alpha p H e_Rp
   −bar L_alpha^a (Y_e)_alpha p epsilon_ab t_b e_Rp
   −bar Q_r^a (Y_d)_r s epsilon_ab t_b d_Rs
   +bar Q_r^a (Y_u)_r s t_a* u_Rs+h.c.                  (Y1)
```

Products retain their order. In particular `t` is even, but its adjoint
still reverses the two constituent fermions. The last three terms are
four-fermion derivative contacts after the literal definition of `t` is
inserted. They are part of the action-coordinate change and cannot be
discarded on the ground that the original loop contained no quark line.
They are not additional loop integrations or independently determined
physical interaction strengths.

The first term in (Y1) has a further exact simplification. Orthogonality
`phi^dagger H=0` removes both contributions to `ell` proportional to `phi`.
It does not remove the differentiated scalar contribution:

```text
overline(ell_alpha)(Y_e)_alpha p H e_Rp
 =−(i/3)(F^dagger Y_e)_i p
      [(D_mu phi)^dagger H]
      overline(A_N N_Li)gamma^mu e_Rp.                 (Y2)
```

The adjoint of the explicit `i` fixes this minus sign. Generically
`(D_mu phi)^dagger H` is nonzero, even though
`(D_mu phi)^dagger H+phi^dagger D_mu H=0`. A constant neutral scalar
projection with vanishing gauge background would miss these interactions.
A charged weak connection can make this contraction nonzero even for
constant scalar components.

The same terms can be obtained from the literal variational sources. The
spectator additions are

```text
(S_L)_alpha^a=(Y_e)_alpha p H^a e_Rp,
(S_phi)^b=sum_alpha p a bar L_alpha^a(Y_e)_alpha p epsilon_ab e_Rp
           +sum_r s a bar Q_r^a(Y_d)_r s epsilon_ab d_Rs
           −sum_r s bar u_Rs(Y_u)_r s* Q_r^b.
```

Thus the change in the full Hermitian spectator action is

```text
delta_S L_spec=overline(ell)S_L−sum_b t_b(S_phi)^b+h.c. (Y3)
```

The up-type part of (Y3) is the conjugate of its displayed term in (Y1).
Both descriptions include it once. These source additions must also be
included wherever the source-dependent `x,y` field shifts or the sourced
representative are evaluated. Equations (Y1)–(Y3) do not replace that
complete source evaluation or the kinetic/gauge/scalar action variations.

## Local-normalizing coordinates and the bridge

The separately defined local-normalizing map has

```text
delta_C L=−phi F(−partial^2+m_H^2)N_L/2,
delta_C phi=delta_C H=0,
delta_C Q=delta_C e_R=delta_C u_R=delta_C d_R=0.
```

Its charged-lepton Yukawa variation is identically zero by
`phi^dagger H=0`; the quark Yukawa variations are zero because none of
their fields changes. This pointwise statement holds for the displayed
local map without a scalar equation, integration by parts, or a
constant-background assumption. It does not assert that the whole
electroweak action is unchanged.

Writing `Xi_old=Xi_S+kappa v_S=Xi_C+kappa v_C` gives
`Xi_S=Xi_C+kappa(v_C−v_S)`. The scalar part of this bridge is
`delta_bridge phi=+t*`, `delta_bridge H=−epsilon t`, and its lepton part
is `delta_bridge L=ell−phi F(−partial^2+m_H^2)N_L/2`.
Applying the bridge to `L_spec` produces exactly the negative of (Y1).
This is the explicit spectator cancellation needed to relate the two
coordinate conventions. Applying both maps with their original signs
would count the source-map variation twice. The bridge acts on `S_0` at
this order; its variation of an already order-`kappa F` packet is beyond
the linear-seed truncation. Other EFT insertions need their own power
counting.

## Reproduction and limits

Run `python -I -B -O verify_spectator_yukawa.py` in this directory. The
standalone SymPy checker uses exact Gaussian rational coefficients and
exterior Grassmann polynomials. Three generic weak/scalar-derivative and
off-shell momentum choices retain two active, two sterile and two quark
flavors, two Weyl components, and noncommuting complex Yukawa matrices.
The actual `ell`, `t`, and local-normalizing lepton row are constructed;
the Lagrangian is then substituted in each slot. A separate variational
source contraction checks (Y3), including the up-type adjoint ordering.
Nonzero contact witnesses reject missing spectators, incorrect scalar
shift signs, an extra half and an unreversed adjoint. Coefficient
polynomials test these action identities; they are not physical amplitudes
or an exhaustive invariant-basis calculation. Color is a spectator copy.

The contact coefficients here describe only the displayed field
substitutions. The source-dependent sterile shifts, other terms in the
matched packet, and further reductions can contribute to the same final
operator classes. Four-dimensional algebra does not define finite
evanescent terms. Complete matching, running and observable predictions
remain unresolved.
