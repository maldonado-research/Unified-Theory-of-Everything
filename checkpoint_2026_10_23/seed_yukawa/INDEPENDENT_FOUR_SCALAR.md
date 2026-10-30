# Independent four-scalar derivative-seed UV audit

23 October 2026, Pacific time. Input commit
`f1f9a892e15545c70e4fa05827ea7fefc5a438eb`.

This audit independently integrates the primitive one-loop four-scalar
pole with one original derivative seed and two conjugate neutrino Yukawa
vertices. Its large-momentum/angular route and literal Wick calculation
were fixed before comparing with the producer's Feynman-parameter result.
They agree on the complete displayed local polynomial and action sign.
No producer code, saved coefficients or repository modules are imported.

## Action, field orientation and result

Use the reference [literal action](../../checkpoint_2026_10_15/operator_map/ACTION_SIGN.md),
with `eta=(+---)`, `d=4-2epsilon`, `kappa=1/(16 pi^2 epsilon)`,
`phi=tilde H`, `P0=-Box` on the **whole** neutral composite, and

```text
L_seed=(1/2) C_ab (bar L_a phi) P0(phi^T L_b^c)+h.c.,
C=C^T,
L_Y=-bar L Y phi N_R+h.c.,
L_mass=-(bar N_L M N_R+h.c.)/2.
```

The first integration may be done with real positive diagonal `M`.
The final flavor expressions extend to an arbitrary constant complex
symmetric retained mass. Define

```text
rho=phi^dagger phi,
A=rho (partial_mu phi)^dagger partial^mu phi,
B=(phi^dagger partial_mu phi)(phi^dagger partial^mu phi),
K1=Tr(C Y* M Y^dagger),
K3=Tr(C Y* M M^dagger M Y^dagger).
```

The primitive **effective-action** pole from this seed orientation is

```text
Gamma_div,S = integral d^4x kappa [K3 rho^2+K1(B-A)],
Gamma_div = Gamma_div,S+its Hermitian conjugate.             (1)
```

Its cancelling counterterm action is the negative of (1). This is a
local off-shell pole for the specified diagrams and original field
coordinates. It is not a complete Higgs potential counterterm,
physical Wilson coefficient, beta function or finite threshold.

The two terms have different dimensions. `K1` has dimension minus two
and multiplies dimension-six scalar derivative operators. `K3` is
dimensionless and multiplies the dimension-four quartic monomial.
Keeping the sterile mass is essential for both. The sterile field is
retained in the action; the local UV pole does not integrate it out.

## One-loop route and the actual numerator

Let the two incoming `phi` momenta be `p,p'=P-p` and the two incoming
`phi^dagger` momenta be `k,k'=-P-k`. The two active propagators have
momenta `q` and `P-q`; choose the sterile denominator to be
`(q+k)^2-m^2`. All scalar momenta are generic and off shell.

At the seed, two internal active fields contract with the light fields
of the two conjugate Yukawa vertices. The remaining internal line is
the anomalous sterile mass propagator. It has no kinetic slash
numerator between its same-chirality endpoints. Contracting the two
active sigma matrices gives the two-gamma scalar trace. For one of the
two identical-scalar seed pairings the stripped integrand is

```text
N(q)/[D1 D2 D3],
N(q)=2[q^2-q.P](q-p)^2,
D1=q^2,
D2=(q-P)^2,
D3=(q+k)^2-m^2.                                         (2)
```

The factor `(q-p)^2` is the momentum of the differentiated neutral
composite. Replacing it by `q^2` drops the scalar and mixed derivatives
and changes the answer. The other seed pairing replaces `p` by `p'`.

With the normal loop measure, the UV part of (2) is `i kappa R`, where

```text
R(p,k)=2m^2-P^2-P.k-P.p+2p^2+2p.k.                       (3)
```

Its mass term is present even when the numerator mass is already
factored out. It comes from the massive propagator's denominator
expansion and ultimately supplies the cubic ordered mass in `K3`.

## Independent auxiliary-mass/angular integration

Introduce a positive auxiliary squared mass `a` and `Q=q^2-a`.
The denominator shifts in (2) are

```text
D1=Q+a,
D2=Q+[a-2q.P+P^2],
D3=Q+[a+2q.k+k^2-m^2].
```

Expand each inverse geometrically. Terms of total geometric order
three or greater are UV integrable: the largest numerator degree is
four, each shift adds at most one power of `q`, and the denominator
then gives superficial degree at most `2-3=-1`. Retaining orders zero,
one and two therefore exhausts the simple-pole coefficient.

The radial UV residue, stripped of `i/(16 pi^2 epsilon)`, is computed
from the Gamma pole rather than assigned in advance:

```text
Res int (q^2)^r/(q^2-a)^n
 =(-1)^(r-n) a^(2+r-n) Gamma(r+2)
    Res Gamma(n-r-2+epsilon)/[Gamma(2)Gamma(n)].
```

The isotropic contractions are performed in `d` dimensions before
extracting the simple pole:

```text
q_mu q_nu -> eta_mu nu q^2/d,
q_mu q_nu q_rho q_sigma
 ->(eta_mu nu eta_rho sigma+eta_mu rho eta_nu sigma
                         +eta_mu sigma eta_nu rho)
       (q^2)^2/[d(d+2)].
```

After summing every retained geometric term, `a` cancels exactly and
the coefficient of the external `k^2` also cancels. Evaluation of the
remaining simple-pole prefactors at `d=4` gives (3). The adjacent JSON
records all seven nonzero expansion contributions before setting `d=4`.
The auxiliary mass is an IR device, distinct from the physical `m`.
No scaleless integral is set to zero and then used to infer absence
of a UV divergence. No finite or evanescent-scheme contribution is
extracted from this simple-pole calculation.

A second independent algebraic reduction checks the integrated result.
The exact numerator identity

```text
2q.(q-P)=q^2+(q-P)^2-P^2
```

reduces (2) to two quadratic-numerator bubbles minus `P^2` times the
logarithmic triangle moment. The bubble identity

```text
Res int (l+b)^2/[l^2((l+Qext)^2-m^2)]
 =m^2-Qext.b+b^2
```

follows from the tadpole residue, vector moment `-Qext/2`, and scalar
bubble residue. Substitution of `(Qext,b)=(P+k,P-p)` and `(k,-p)` gives
exactly (3). This check uses no Feynman parameters.

## Literal Wick signs, Majorana factor and seed half

The overall sign in (1) is fixed from the same action. In a Weyl basis
write `epsilon=i sigma^2`, `n=N_L`. The four-component convention
`C_spin=i gamma^2 gamma^0` gives

```text
N=(n,-epsilon n^dagger,T),
L_seed=-(C/2) phi phi bar l epsilon bar l^T P0,
B_Y=-Y* phi^dagger l^T epsilon n,
<n n^T>=-i M epsilon/(q^2-m^2).
```

The ordered active contractions are
`<bar l l>=-i sigma(q)/q^2`. In the sequence

```text
bar l_a bar l_b l_c n_e l_f n_g
```

the pairings `(02)(14)(35)` and `(04)(12)(35)` have crossing signs
`+` and `-`. The antisymmetric seed spin tensor makes their complete
spin contractions equal: each gives `2q.(P-q)`. The checker verifies
all sixteen Lorentz coefficient pairs using explicit two-component
matrices.

Including the seed coefficient `-1/2`, both negative Yukawa action
coefficients, the third-order expansion `i^3/2!`, and the three ordered
propagator phases `(-i)^3`, the phase factor before the sum of the two
spin contractions is `-1/4`. Thus the unlabelled functional integrand
is `q.(q-P)(q-p)^2`, half of the stripped (2). Differentiating the
identical scalar fields restores the labelled kernel below. The
quadratic Majorana action's half is already accounted for in its
canonical propagator. Appending another half is incorrect.

The closed-loop sign is already contained in this complete Wick
calculation. Appending an additional closed-loop minus would reverse
the answer. That alternative fails the independent constant-background
check as well as the literal phase calculation.

For that check take real one-flavor fields and couplings. The Weyl
mass matrix in a constant scalar background, to first seed order, is

```text
Mcal(q^2)=[[-C phi^2 q^2, y phi],
          [y phi,            M]].
```

At linear `C`, the determinant of `q^2-Mcal^2` has term
`2C q^2 y^2 M phi^4`. To order `C y^2`, its contribution to `log det`
is `2C y^2 M phi^4/(q^2-M^2)`. The Majorana effective action
`-i integral log det` therefore gives
`+2 kappa C y^2 M^3 phi^4`, which is exactly (1) plus its conjugate
in this real limit. This check fixes the action sign and mass-cube
normalization; it does not replace the momentum-dependent integration.

## Local-action reconstruction and conjugation

Let scalar weak indices be `a,b` and conjugate-scalar indices `c,d`.
Differentiating (1) gives the displayed-seed four-point kernel

```text
i V_S,UV=+i kappa {
 delta_ac delta_bd [2K3+K1 R0(p,k)]
+delta_ad delta_bc [2K3+K1 R0(p',k)] },
R0=R-2m^2.                                             (4)
```

The local monomials are independently differentiated by all assignments
of the two identical scalar and two conjugate-scalar fields. For the
direct weak contraction, `A` gives
`-[p.k+p'.k']`, while `B` gives `-2p.p'`. Their difference is exactly
`R0(p,k)`. The crossed weak contraction gives `R0(p',k)`. The potential
`rho^2` gives `2(delta_ac delta_bd+delta_ad delta_bc)`, fixing its
relative coefficient. Five generic rational momentum choices and all
sixteen weak-index choices check the complete labelled polynomial,
identical-scalar exchanges and Hermitian-conjugate derivative placement.
The exact symbolic reconstruction, rather than those examples alone,
fixes the local map.

In particular `R(p,k)=R(P-p,-P-k)`, as required by exchanging the two
conjugate scalars and reversing the loop. The two weak channels must
both be retained.

The conjugate of `B` has the derivatives on `phi^dagger`. Consequently
one must retain

```text
K1 B+K1* B^dagger-(K1+K1*)A,
```

rather than replace `K1(B-A)+h.c.` by `2Re(K1)(B-A)`.
For a general complex invariant the latter loses a nonzero real action
contribution proportional to the imaginary part and the antisymmetric
derivative placement. The checker supplies such a complex-field
example, verifies reality, and rejects the premature real-part
replacement. It also checks global unitary weak rotations. No scalar
EOM, vacuum condition or current-source deletion is used here.

For general complex symmetric mass, the anomalous propagator's large-
momentum mass series has ordered coefficients `M` and `M M^dagger M`.
Their Yukawa/flavor contractions are exactly `K1,K3` above. Under the
passive convention

```text
Y'=V^dagger Y U*, M'=U^dagger M U*, C'=V^dagger C V*,
```

both traces are invariant. The checker verifies this with rectangular
complex Yukawas, symmetric complex masses, singular masses and zero
mass, without inverses. Replacing the ordered palindrome by `M^3`,
conjugating the wrong mass, or omitting the Yukawa conjugation fails.

## Reproduction and scope

Run `python -I -B -O independent_four_scalar.py` with pinned SymPy.
The [receipt](INDEPENDENT_FOUR_SCALAR_RESULTS.json) records **397 positive
controls and 19 rejected mutations**. Many controls are entries of
Lorentz or weak-component comparisons; they are not 397 separate loop
integrations. Explicit exceptions keep failures active under optimization.

This result concerns the original seed, two neutrino Yukawa vertices,
four external scalars and zero external gauge field. The two-fermion/
two-scalar matter sector and other interaction orders have separate
calculations. Equation (1) is not a complete scalar/fermion EFT map,
physical parameter running, finite threshold, observable prediction,
new fundamental mathematics or external novelty determination. These
independent implementations are AI-assisted internal verification,
not independent human peer review.
