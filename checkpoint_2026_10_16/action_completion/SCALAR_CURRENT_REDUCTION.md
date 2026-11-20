# Scalar gauge-current descendants and their sourced transport

16 October 2026. Input: candidate
`5f9889b76661faa155c681a5ecfac5d036398bc8`, especially the selected
[`ELECTROWEAK_REDUCTION.md`](../../checkpoint_2026_10_15/interacting_reduction/ELECTROWEAK_REDUCTION.md)
representative and its
[`CURRENT_DESCENDANTS.md`](../../checkpoint_2026_10_15/interacting_reduction/CURRENT_DESCENDANTS.md).
This calculation reduces the **scalar part of that gauge-current
subpacket**, retaining its fermion sources and giving an explicit local
field map for the additional EOM step. It does not integrate a new loop.

## Fields, sources and the three algebraic scalar tensors

Keep `D=partial+i Acal`, `phi=tilde H`, `T^I=tau^I/2`,
`y_L=y_phi=-1/2`, `rho=phi^dagger phi`, `A_N=i slash partial`,
`K bar L=-i(D_mu bar L)gamma^mu`, and
`F=C_seed Y_nu*` in that flavor order. `N=N^c` has the same real positive
diagonal retained mass `M` as the input action. A constant real symmetric
mass also works with the displayed ordered products.

Define the even weak row and singlet current

```text
b_a^mu=sum_(alpha i)F_(alpha i)bar L_alpha^a gamma^mu N_Li,
j^mu=b_a^mu phi^a,     j_I^mu=b_a^mu T^I_ab phi^b.
```

For one derivative on the scalar fields, a useful algebraic basis is

```text
O_1=rho b_mu D^mu phi,
O_2=j_mu(phi^dagger D^mu phi),
O_3=j_mu((D^mu phi)^dagger phi),
O_norm=j_mu partial^mu rho=O_2+O_3.                       (1)
```

The brackets in each scalar product fix its weak-index ordering. The
last equality is a pointwise covariant-product identity because `rho` is
neutral. It is not an integration-by-parts or EOM substitution.

There are exactly **two** weak invariant contractions with the derivative
on `phi` (`O_1,O_2`), and **one** with it on `phi^dagger` (`O_3`), when
`b_mu` is treated as one even anti-doublet vector. With the derivative on
`phi`, the independent scalar slots transform as `2 tensor 2 tensor 2*`,
which contains two copies of the output doublet. With it on
`phi^dagger`, the two identical undifferentiated `phi` fields furnish
`Sym^2(2)=3`, and `3 tensor 2*` contains one output doublet. The three
contractions are therefore algebraically independent before any derivative
is moved onto a fermion.

An exact coefficient-space calculation supplements that argument: the
most general derivative-on-`phi` output vector has 16 coefficients; its
48 weak-covariance equations have rank 14. The conjugate-derivative vector
has 12 coefficients; its 36 equations have rank 11. The explicit tensors
in (1) span those complete nullspaces. This is completeness of the
declared **scalar-derivative weak contractions**, not of every fermion,
Lorentz or derivative operator in the retained-sterile EFT.

## Literal scalar-current coefficients

The scalar currents of the same gauge-invariant reference action are

```text
J_Y^mu(phi)=i y_phi[phi^dagger D^mu phi-(D^mu phi)^dagger phi],
J_I^mu(phi)=i phi^dagger T^I D^mu phi
            -i(D^mu phi)^dagger T^I phi.
```

Both include their connection terms. The Pauli completeness tensor gives
the ordered, pointwise identity

```text
sum_I j_I,mu J_I^mu(phi)=i O_1/2-i(O_2+O_3)/4.
```

Consequently the selected scalar part of the prior sourced representative is

```text
P_scalar=-(g'^2 y_L j_mu J_Y^mu(phi)+g^2 j_I,mu J_I^mu(phi))/6
        =i[-2g^2 O_1+(g^2-g'^2)O_2+(g^2+g'^2)O_3]/24
        =i[(g^2+g'^2)O_norm-2g^2 O_1-2g'^2 O_2]/24.     (2)
```

Here `y_L y_phi=1/4`; the negative hypercharges do not give an additional
minus in their product. Equation (2) is exact tensor algebra. Its
normalization is fixed by the already selected gauge-EOM representative;
it is not an additional matching calculation.

## Integration by parts with the full fermion derivatives

Let `[bar U V Z]_F=sum_(alpha i)F_(alpha i)bar U_alpha V Z_i`.
Covariant product differentiation, without imposing an EOM, gives

```text
partial_mu j^mu=b_mu D^mu phi
  +i[(K bar L)phi N_L]_F-i[bar L phi A_N N_L]_F.         (3)
```

The row and column connections cancel in exactly that weak matrix order.
In particular the row source has the positive `i K bar L` sign, whereas
`gamma partial N_L=-i A_N N_L`. Applying (3) to `partial_mu(rho j^mu)`
yields

```text
O_norm=partial_mu(rho j^mu)-O_1
       -i rho([(K bar L)phi N_L]_F-[bar L phi A_N N_L]_F).
```

It would be incorrect to replace this by `O_norm=-O_1` while retaining
an interacting fermion action. Define

```text
d=(g^2+g'^2)/24,
e_barL=K bar L-bar S_L,
e_R=A_N N_L-M N_R-S_R.
```

The complete scalar subpacket has the following chosen representative:

```text
P_scalar,red= -i(3g^2+g'^2)O_1/24-i g'^2 O_2/12
       +d rho([bar S_L phi N_L]_F
                    -[bar L phi(M N_R+S_R)]_F).         (4)
```

Its exact off-shell certificate is

```text
P_scalar-P_scalar,red
 =i d partial_mu(rho j^mu)
  +d rho([e_barL phi N_L]_F-[bar L phi e_R]_F).          (5)
```

Thus two scalar-derivative contractions remain after this specified
IBP/fermion-EOM step, together with correlated non-derivative source
operators. This is an EOM equivalence of the action, not a pointwise
equality after dropping residuals. No scalar equation, free sterile
equation, stationary background or division by `rho` or a mass is used.
Sources are the total variational sources of the same retained reference
action; optional renormalizable spectators add their actual contributions
to (4)-(5).

## Explicit Yukawa and Majorana scalar descendants

Retain the neutrino and charged-lepton reference Yukawa interactions
before simplifying weak contractions:

```text
bar S_L,alpha^a=bar N_Rj Y_nu*_(alpha j)(phi^a)*
               +bar e_Rr Y_e*_(alpha r)(H^a)*
               +bar S_L,other,alpha^a,
S_Ri=sum_beta(Y_nu)_(beta i)X_beta^c+S_R,other,i,
X_alpha=phi^dagger L_alpha,     bar X_alpha=bar L_alpha phi,
X_beta^c=sum_a phi^a(L_beta^a)^c.
```

The explicit non-derivative part of (4) is

```text
+d rho^2 (Y_nu^dagger F)_(ji) bar N_Rj N_Li
-d rho (F M)_(alpha i) bar X_alpha N_Ri
-d rho (F Y_nu^T)_(alpha beta) bar X_alpha X_beta^c
+d rho (Y_e^dagger F)_(ri)bar e_Rr(H^dagger phi)N_Li
+d rho([bar S_L,other phi N_L]_F
                         -[bar L phi S_R,other]_F).      (6)
```

The charged-lepton monomial is retained until its displayed contraction:
`H^dagger phi=H^dagger epsilon H*=0` for the single commuting scalar
doublet. It therefore vanishes algebraically in this scalar subpacket.
This does not delete the charged-lepton currents, sources or operators
elsewhere in the theory. Additional scalar species or different source
orientations require their own contraction check.

Both scalar spin contractions in the first and third lines are **flavor
symmetric** after the Grassmann exchange:

```text
bar N_Rj N_Li=bar N_Ri N_Lj,
bar X_alpha X_beta^c=bar X_beta X_alpha^c.
```

Their spin contraction matrix is antisymmetric, in contrast to the
symmetric `sigma C` matrix of the prior strength operators. With
`H_N=Y_nu^dagger F`, `K=F Y_nu^T`, the ordered scalar coefficients are
therefore `+d(H_N+H_N^T)/2` and `-d(K+K^T)/2`. If sums use independent
off-diagonal flavor pairs, these coefficients double; diagonal terms
retain their displayed coefficients. For the selected symmetric seed,
`Y_nu^dagger C_seed Y_nu*` is already symmetric. No mixed active/sterile
`1/2` is appended to the mass term, whose flavor order remains `F M`.
These conventions specify literal ordered monomials; another convention
that includes a scalar-operator `1/2` must convert its coefficients.

After factoring out `F`, the derivative monomials in (1), the
`rho^2 bar N_R N_L` term and the active `rho bar X X^c` term have canonical
dimension seven. The `rho bar X N_R` monomial has dimension six, with
coefficient `F M` of dimension `-2` but fixed inverse-scale order
`Lambda^-3` under retained-mass counting. Yukawa and gauge couplings are
dimensionless. The common counterterm factor remains
`kappa=1/(16 pi^2 epsilon)`, and the action contains the Hermitian
conjugate once.

## An explicit additional local field map

The residual in (5) can be removed against the gauge-invariant
renormalizable action `S_0` by a further first-order substitution. All
fields on the right are new fields. Define

```text
ell_alpha^a=d rho sum_i F_(alpha i)phi^a N_Li,
x_i=-d rho sum_(alpha a)F*_(alpha i)(phi^a)*L_alpha^a,
L_old,alpha^a=L_alpha^a-kappa ell_alpha^a,
N_old=N-kappa(x+x^c).                                   (7)
```

Equivalently the sterile shift in (7) is positive `kappa d rho` times the
Majorana-real completion of `sum F* phi^dagger L`. No scalar or gauge
field is shifted in this additional map. The lepton adjoint is its actual
Hermitian conjugate. `ell` is a charged left-chiral doublet; `x` is a
neutral left-chiral polynomial. Their dimensions are the field dimensions.
The integrated action variation cancels the second line of (5) and its
conjugate. Majorana completion uses `x+x^c` without an extra half, because
variation of the two Majorana quadratic slots cancels their action factor
`1/2`. Thus

```text
S_0[old]+kappa integral(P_scalar[old]+h.c.)
 =S_0[new]+kappa integral(P_scalar,red[new]+h.c.)
  +O((kappa F)^2),                                      (8)
```

modulo the boundary term. Equation (7) is a specified additional
coordinate map. It is not silently composed with the earlier source or
normalization maps. Varying another retained EFT operator additionally
produces its appropriate `kappa F C_d` contribution; its inverse-scale
order must be counted separately rather than called quadratic in kappa.

## Relation to the existing coefficient sectors and reproduction

The third line of (6) uses the same **raw flavor-symmetric source matrix
`K_S`** that appeared in the flavor-symmetric weak strength channel, but
now multiplies a scalar dimension-seven monomial. The first line uses the
sterile matrix `Y_nu^dagger F`; the second uses the existing ordered mass
matrix `F M`. These are different operator coordinates, with the specific
correlation set by (2)-(6). This step introduces no field strength and
does not by itself change the previously reported active/retained-mass
tensor coefficients or four-fermion-current coordinates. Further
derivative reductions and other portions of the full selected action can
contribute to the same final scalar coordinates and must be combined in a
common declared basis. The coefficients in (6) are only this scalar
gauge-current contribution.

Run `python3 -I -S -B -O verify_scalar_current_reduction.py`. The
standalone standard-library checker proves the full weak polynomial
identities and covariance ranks, tests (3)-(5) on eight generic
noncommuting off-shell connection/matter/source jets, expands all displayed
Yukawa/Majorana scalar sources in an exact exterior algebra with two active,
three sterile and two charged-lepton flavors, and checks the real Majorana
completion of (7). Its completed optimized run passes **75 positive exact
checks** and rejects **12 deliberate mutations**. The saved JSON records
the full categories, rank systems and source/IBP controls.

This is a concrete source-retaining reduction of the selected scalar
gauge-current subpacket. It does not prove complete derivative-EFT basis
closure, supply the missing loop sectors or counterterms, canonically
normalize a full interacting theory, define a subtraction/evanescent
prescription, or yield a physical Wilson prediction, finite amplitude,
running coefficient or observable.
