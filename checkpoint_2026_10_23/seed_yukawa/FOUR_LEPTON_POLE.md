# A four-lepton pole from the derivative seed and two Yukawa vertices

23 October 2026, Pacific. Input: `f1f9a892e15545c70e4fa05827ea7fefc5a438eb`.

This calculation integrates a previously uncomputed one-loop four-lepton
triangle with one original derivative seed and two conjugate neutrino-Yukawa
vertices. A retained Majorana mass makes the same-orientation sterile
propagator possible. The result is a mass-dependent dimension-six contact
pole, linear in the dimension-seven seed. No sterile field is removed from
the EFT, and this isolated graph contribution is not a total physical Wilson
coefficient, complete beta function or finite scattering prediction.

## Literal action and ordered result

Use the [literal reference action](../../checkpoint_2026_10_15/operator_map/ACTION_SIGN.md),
`eta=(+---)`, `d=4-2 epsilon`, `kappa=1/(16 pi^2 epsilon)`, and
`S=T exp(i integral L_int)`. The terms needed here are

```text
L_seed = (1/2) C_alpha beta (bar L_alpha^a phi^a)
                       (-Box)(phi^b L_beta^{c b}) + h.c.,
L_Y = -bar L_alpha^a Y_alpha i phi^a N_Ri + h.c.,
phi=tilde H, C=C^T, M=M^T,
L_M = -(bar N_L M N_R+h.c.)/2.
```

The derivative acts on the entire second neutral composite. Define the
symmetric active-flavor mass tensor

```text
R_gamma delta = (Y* M Y^dagger)_gamma delta,
R=R^T.
```

For an unambiguous ordered-spinor expression use a chiral gamma basis,
`epsilon_spin=[[0,1],[-1,0]]`, `L=(l,0)` and `bar L=(0,bar l)`, and set

```text
P_alpha beta gamma delta
  = sum_(a,b) (bar l_alpha^a epsilon_spin bar l_beta^{b T})
              (l_gamma^{a T} epsilon_spin l_delta^b).
```

All displayed spinor components anticommute. The local UV part of the
one-loop action and its cancelling counterterm are

```text
L_loop^UV = +(kappa/2) C_alpha beta R_gamma delta
                         P_alpha beta gamma delta + h.c.,
L_ct      = -(kappa/2) C_alpha beta R_gamma delta
                         P_alpha beta gamma delta + h.c.                 (1)
```

The amplitude from this action carries an additional overall `i`. The
Hermitian conjugate in (1) is the conjugate-seed graph, rather than an
extra copy of the displayed complex coefficient. The graph does not have
a closed fermion loop: the sterile line joins two open light-fermion
sources. The Wick ordering sign along that line is included below.

For four-dimensional physical external spinors define
`J_alpha gamma^mu=bar L_alpha^a gamma^mu L_gamma^a`. Explicit Pauli
completeness with Grassmann reordering gives

```text
J_alpha gamma^mu J_beta delta,mu = -2 P_alpha beta gamma delta.
```

Thus the same physical-spinor projection has

```text
L_ct = +(kappa/4) C_alpha beta (Y* M Y^dagger)_gamma delta
                    (bar L_alpha gamma_mu L_gamma)
                    (bar L_beta gamma^mu L_delta) + h.c.                 (2)
```

Equation (1) is the primary ordered result; (2) is its four-dimensional
Fierz translation. This calculation does not specify a complete
`d`-dimensional evanescent-operator basis or the finite terms associated
with changing that basis.

## Charge conjugation and the complex mass orientation

In the chiral basis the reference `C_spin=i gamma^2 gamma^0` is
`diag(epsilon_spin,-epsilon_spin)`. Consequently

```text
L^c=(0,-epsilon_spin bar l^T),
bar L_alpha L_beta^c=-bar l_alpha epsilon_spin bar l_beta^T.
```

Let `n=N_L`. The conjugate Yukawa term can be ordered as

```text
B = -Y*_(gamma i) phi^dagger l_gamma^T epsilon_spin n_i
  = -Y*_(gamma i) phi^dagger n_i^T epsilon_spin l_gamma.
```

The equality uses Grassmann exchange as well as the antisymmetric spin
matrix. The left-Weyl holomorphic mass is `M*` in this reference action.
For a real diagonal mass the same-orientation propagator is

```text
<n_i n_j^T>(q) = -i M_i delta_ij epsilon_spin/(q^2-M_i^2+i0).
```

For general constant complex symmetric `M`, its leading large-momentum
term is `-i M_ij epsilon_spin/q^2`. This ordering follows directly by
expanding the inverse four-component quadratic kernel. Its first
mass insertion is

```text
P_L slash q (M P_R+M^dagger P_L) slash q P_L/q^2
   = M P_L,
<N_L N_L^T> = P_L S_N C_spin^T P_L^T.
```

The lower expression supplies the minus spin epsilon. No full matrix
propagator denominator needs to be commuted through `M`. Later terms in
its high-momentum expansion are two powers more convergent and cannot
alter the logarithmic four-lepton residue. Contraction of the two Yukawa
flavor slots therefore gives `Y* M Y^dagger`, not its adjoint or a mass
inverse. The result is valid for singular and zero masses without using
an inverse; at `M=0` this anomalous triangle is absent.

## Actual triangle integral and UV/IR separation

Two scalar propagators join the two seed scalar slots to the two `B`
vertices. The remaining sterile endpoints join through the anomalous
propagator. The other four fields are external leptons. There is no
free weak-doublet trace: each scalar propagator identifies one seed weak
index with the light lepton at its Yukawa endpoint.

The two internal scalars belong to opposite neutral composites of the
seed. If their momenta scale as `q` and `-q`, the derivative on the second
whole composite contributes `(q+p_external)^2`. Its leading term is
`q^2` with coefficient one. Replacing the seed by derivatives only on its
external lepton would remove this term and erase the actual pole.

For a real diagonal sterile mass, let `D0,D1,D2` denote the two scalar and
one sterile quadratic denominators, with generic shifted off-shell
external momenta. The leading integral is

```text
I = integral d^d q/(2pi)^d q^2/(D0 D1 D2),
1/(D0 D1 D2)=2 integral_simplex dx dy
                         1/[(l^2-Delta+i0)^3],
q=l-(linear combination of external momenta).
```

The parameter simplex has `x>=0`, `y>=0`, `x+y<=1` and area `1/2`.
Translation of the quadratic numerator leaves `l^2` with coefficient
one; the remaining term linear in `l` integrates to zero, and its
constant terms multiply a UV-finite three-denominator integral. The
radial masters follow from

```text
integral d^d l/(2pi)^d 1/(l^2-Delta+i0)^n
 = i(-1)^n Gamma(n-2+epsilon) Delta^(2-epsilon-n)
        /[(4pi)^(2-epsilon) Gamma(n)].
```

In particular

```text
[l^2/(l^2-Delta+i0)^3]_UV
 = [1/(l^2-Delta+i0)^2 + Delta/(l^2-Delta+i0)^3]_UV
 = i kappa.
```

The complete Feynman-parameter residue is therefore

```text
I_UV = 2 (simplex area 1/2) i kappa = i kappa.            (3)
```

The checker integrates exact simplex monomials and Gamma-pole
recurrences, including generic shifted numerator coefficients. No
anomalous dimension or preassigned four-lepton coefficient is used as
an integration input.

Retain generic nonexceptional off-shell momenta and positive internal
squared masses while extracting the UV pole. Equivalently a controlled
positive auxiliary mass can protect the IR without changing this
logarithmic residue; it is distinct from the retained physical matrix
`M`. The UV result analytically continues to the signed scalar mass
parameter of the reference action. A scaleless integral is not set to
zero and then used to infer a UV cancellation. Mass and momentum
independence of (3) follows from the degree-zero bound, not from dropping
physical `M` in the anomalous numerator.

## Literal expansion factors, Wick sign and scalar assignments

Keep both the seed half and the expansion coefficient `1/2!` of two
identical `B` vertices. In the chiral basis the seed action has the
additional minus from `bar L L^c=-bar l epsilon_spin bar l^T`.
For one ordered scalar pairing the literal factors are

| Factor | Value |
| --- | --- |
| Three vertex factors | `i^3` |
| Seed spin-conjugation sign | `-1` |
| Two Yukawa action signs | `(-1)^2` |
| Anomalous sterile propagator, spin epsilon stripped | `-i M` |
| Two scalar propagator phases | `i^2` |
| Wick crossing between the two sterile endpoints | `-1` |
| Seed and two-Yukawa expansion factors | `1/(2 times 2!)=1/4` |

The Wick sign is explicit. With the `B` bilinears ordered as
`n_r epsilon_rs l_s`, the relevant fermion sequence is

```text
[bar l_alpha, bar l_beta, n_at_B1, l_gamma, n_at_B2, l_delta].
```

The sterile contraction crosses one external `l_gamma`. Its spin matrix
contraction satisfies

```text
sum_(r,t) epsilon_(r s) epsilon_(t u) epsilon_(r t)
 = epsilon_(s u).
```

Multiplying all factors in the table gives **`+1/4`**, before the `i
kappa` from (3). There is no additional closed-loop minus or Majorana
half to append. The two scalar Wick pairings give

```text
(1/4) C_alpha beta R_gamma delta sum_(a,b)
 (bar l_alpha^a epsilon bar l_beta^{b T})
 [(l_gamma^{a T} epsilon l_delta^b)
  +(l_gamma^{b T} epsilon l_delta^a)].                    (4)
```

Using `C=C^T`, `R=R^T` and the symmetry of an odd-spinor epsilon bilinear
under joint exchange of its two field labels, the two terms in (4) agree
after flavor summation. This gives exactly (1). The factor two is two
identified scalar assignments; it is not a free weak trace. In an
external-amplitude calculation, differentiating the four Grassmann
fields supplies their fermion permutations from (1), without introducing
another action-level half or duplicating its Hermitian conjugate.

For real one-flavor `C,Y=y,M`, (1) including its conjugate gives
`L_loop^UV=kappa C y^2 M P`. At one fixed weak component,
`P=4 bar l_1 bar l_2 l_1 l_2`, which is nonzero. Thus identical-flavor
Grassmann statistics do not erase this four-lepton pole.

## Executable controls and interpretation

Run `python -I -B -O verify_four_lepton.py` in the root pinned environment.
The standalone program uses exact SymPy matrices and an exterior algebra
with independent barred/unbarred lepton components. Its adjacent
[receipt](FOUR_LEPTON_RESULTS.json) records the check totals. The program
uses explicit exceptions, imports no repository code and needs no
optimization-sensitive assertions.

The controls include actual simplex/radial residues, complete shifted
quadratic numerator coefficients, complex-mass inverse-propagator
ordering, the spin epsilon contractions, both scalar pairings, the
literal phases and action cancellation, exact component Fierz and
Hermitian-conjugation identities, a nonzero one-flavor component, and
active/sterile flavor covariance. Mutations reverse the Wick/CT sign,
omit the crossed scalar assignment, append an extra Majorana half,
delete the seed half, drop scalar derivatives, use the mass adjoint, or
apply a commuting-spinor Fierz sign.

With `L_old=V L_new`, `N_L,old=U N_L,new`, the prior convention gives
`Y'=V^dagger Y U*`, `M'=U^dagger M U*` and
`C'=V^dagger C V*`; the new source satisfies `R'=V^T R V`.
The different transformations of the barred and unbarred flavor slots
are essential for the Hermitian contact.

This supplies a specific primitive one-loop local pole and its
cancelling action. Earlier equation-of-motion reductions and field
changes can generate additional terms in the same four-lepton class.
Those contributions must be collected consistently before quoting a
total coefficient in a chosen physical basis. No such total, finite
threshold matching, observable, external novelty or empirical
confirmation is claimed here.
