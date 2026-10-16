# Pure-neutrino-Yukawa one-loop two-point poles

15 October 2026, Pacific client date; the directory retains the 16 October
UTC workspace date. Input candidate:
`0770c304c7e05f213dd796d3fd5a194107f38858`.

This note actually integrates the renormalizable pure-neutrino-Yukawa
one-loop two-point functions for `L`, retained Majorana `N`, and
`phi=tilde H`. It supplies missing kinetic and mass counterterm inputs.
It contains no seed insertion and does not calculate dimension-seven
running. Gauge couplings, other Yukawas, scalar quartic interactions and
all higher-dimensional insertions are excluded from this coupling sector.

## Literal action, amplitudes and result

Use the reference conventions of
[`ACTION_SIGN.md`](../../checkpoint_2026_10_15/operator_map/ACTION_SIGN.md):

```text
eta=(+---), d=4-2epsilon, kappa=1/(16 pi^2 epsilon),
sigma^{mu nu}=i[gamma^mu,gamma^nu]/2,
C=i gamma^2 gamma^0, N=N^c,
phi=tilde H, D=partial+i Acal,
L_Y=-(bar L_alpha^a Y_(alpha i)phi^a N_Ri+h.c.),
L_kin=bar L i slash D L+(Dphi)^dagger Dphi
      -m_H^2 phi^dagger phi+bar N(i slash partial-M)N/2.
```

The primary retained mass `M=diag(M_i)` is real positive. `m_H^2` is the
signed parameter in the literal action, and `Y` is an arbitrary complex
active-by-sterile matrix. Set gauge backgrounds to zero for these loops.
Define

```text
H_L=Y Y^dagger,     H=Y^dagger Y,
T_Y=tr(H),         T_M=tr(H M^2)=sum_i M_i^2 H_ii.
```

`i Gamma_loop` below means the **amputated 1PI correction to the inverse
kernel**, whose tree kernels are `i slash p` and `i(p^2-m_H^2)`. It is
not a self-energy symbol with an extra conventional minus. Columns have
the incoming wave `exp(-ip.x)`; their adjoints have `exp(+ip.x)`.
The poles are

```text
i Gamma_L^UV=+i kappa (H_L)_(alpha beta) slash p P_L/2
                         times delta_ab,
i Gamma_N^UV=+i kappa slash p[H P_R+H^T P_L],
i Gamma_phi^UV=+i kappa[T_Y p^2-2T_M] delta_ab.           (1)
```

The full Majorana expression in (1) has both conjugate chiral channels.
There is **no proper chirality-flipping fermion mass pole** in this
restricted sector. The scalar has both a kinetic pole and a nonzero
retained-mass pole. No weak multiplicity is appended to its fixed external
component.

## Phases, contractions and multiplicities

The propagators are `i(slash q+M_i)/(q^2-M_i^2+i0)`,
`i slash q/(q^2+i0)` for the massless active fermion, and
`i/(q^2-m_H^2+i0)` for `phi`. The two Yukawa vertices are
`-i Y P_R` and `-i Y^dagger P_L` in their open-chain orientations.
An open two-point chain therefore has literal phase

```text
(-i)^2 i^2=+1.
```

Its loop integral contributes the positive `i/(16 pi^2 epsilon)` master
pole. The scalar graph has the same vertex/propagator phase and one
additional **closed fermion minus**, giving the opposite sign before
the spin trace is integrated. In the ordered product
`(bar L N)(bar N L)`, the nested Wick pairings have one reversed
contraction `<bar L L>=-<L bar L>`; interchanging the two vertices gives
the same sign. Their two expansion assignments cancel `1/2!`.

There is no extra closed-loop `1/2` for this mixed active/sterile loop.
The `1/2` in the Majorana quadratic action is already removed by
differentiating its two fermion slots, so it gives the canonical Majorana
propagator. The two vertices in the scalar loop have distinct Yukawa and
conjugate-Yukawa orientations fixed by the external `phi,phi^dagger`.
Appending another Majorana half would halve a graph that has already
been counted once.

The weak contractions are also different in the three graphs. A fixed
external active component receives one internal scalar component; a
fixed external scalar component receives one active component. In the
sterile graph the internal doublet index is summed, giving **two**.
This is exactly the factor multiplying the vector-bubble `1/2` below.

For the external Majorana two-point kernel, the ordinary oriented chain
gives `H slash p P_R`. Conjugating the full open chain, including
`p -> -p` and the spin/flavor transpose, gives `H^T slash p P_L`:

```text
C[slash(-p)P_R]^T C^(-1)=slash p P_L.
```

These are the two required chiral channels, rather than two copies of
either one. In the counterterm action their common quadratic factor
`1/2` cancels when the Majorana slots are varied; it does not divide (1).

## Actual Feynman-parameter integrations and UV/IR separation

Take `p^2<0` and first take positive internal squared masses where present.
The off-shell bubbles are IR safe; the scalar bubble has one massless
active line, with an integrable parameter endpoint. The local UV
coefficients extend analytically to the signed scalar mass parameter.
No scaleless integral is set to zero and then interpreted as having no
UV pole. The computation extracts the UV part directly from the complete
off-shell bubble.

For `Delta>0`, dimensional regularization gives

```text
I_n(Delta)=mu^(2epsilon) integral d^d l/(2pi)^d
                              1/(l^2-Delta+i0)^n
 =i(-1)^n Gamma(n-2+epsilon)Delta^(2-epsilon-n)
                              /[(4pi)^(2-epsilon)Gamma(n)].
```

Thus the simple-pole masters are

```text
I_1^UV=i kappa Delta,
I_2^UV=i kappa,
integral l^2/(l^2-Delta+i0)^2 |_UV=2i kappa Delta.        (2)
```

The last result follows from `l^2=(l^2-Delta)+Delta`; it must not be
replaced by only one of its two terms. Gamma recursion gives
`Res Gamma(epsilon)=1` and `Res Gamma(epsilon-1)=-1`, fixing the signs
in (2).

For the active graph, assign `q` to the sterile line and `q-p` to the
scalar. After combining denominators with weight `x` on the sterile one,

```text
Delta_L=x M_i^2+(1-x)m_H^2-x(1-x)p^2,
q=l+(1-x)p.
```

The projector chain is
`P_R(slash q+M_i)P_L=slash q P_L`: its mass numerator is exactly zero.
The odd `l` integral vanishes by translation-invariant dimensional
regularization, and

```text
integral_0^1 dx (1-x)=1/2.
```

This proves the first pole in (1); the denominator masses do not enter
its logarithmic kinetic coefficient. For the sterile graph, use the
massless active numerator and
`Delta_N=(1-x)m_H^2-x(1-x)p^2`. The same parameter moment gives `1/2`,
the weak doublet gives two, and the Majorana-conjugate chain supplies
the other chirality. The massless internal active propagator has no
chirality-flipping numerator. Same-oriented scalar contractions are
absent in the unbroken complex-scalar vacuum, so they do not generate an
additional sterile mass graph or an anomalous `phi phi` two-point pole.

For the scalar graph at fixed weak component, the numerator is actually
traced:

```text
tr_spin[P_R(slash q+M_i)P_L(slash q-slash p)]
 =2 q.(q-p).
```

There is no linear mass numerator. With
`Delta_phi=x M_i^2-x(1-x)p^2` and `q=l+(1-x)p`, this becomes
`2[l^2-x(1-x)p^2]` plus an odd term. Equations (2) give its pole

```text
i kappa integral_0^1 dx [4Delta_phi-2x(1-x)p^2]
 =i kappa integral_0^1 dx [4x M_i^2-6x(1-x)p^2]
 =i kappa(2M_i^2-p^2).                                  (3)
```

Multiplying the closed fermion minus and summing `|Y_(alpha i)|^2`
gives exactly the third pole in (1). In particular the scalar mass term
comes from the massive denominator in the radial integral even though
the traced mass numerator itself vanishes. Four-dimensional open-chain
and two-gamma projector identities suffice for these simple poles; no
anomalous closed `gamma5` trace occurs here.

## Cancelling action and bare-field/parameter map

The literal zero-gauge counterterm action that cancels (1) is

```text
L_ct,2pt= -kappa bar L H_L i slash partial L/2
          -kappa bar N i slash partial[H P_R+H^T P_L]N/2
          -kappa T_Y (partial phi)^dagger partial phi
          +2kappa T_M phi^dagger phi.                    (4)
```

The Majorana expression is written once; a duplicate Hermitian conjugate
is not appended. In the gauge-invariant reference action the field
counterterms multiply the existing covariant kinetic terms. Only their
zero-gauge two-point coefficients are calculated here; this is not an
evaluation of additional mixed gauge/Yukawa vertices.

Define the fields and dimensionful parameters by

```text
L_0=(1+delta Z_L/2)L,
phi_0=(1+delta Z_phi/2)phi,
N_R0=(1+delta Z_R/2)N_R,
N_L0=(1+delta Z_R*/2)N_L,
m_H,0^2=m_H^2+delta m_H,param^2,
M_0=M+delta M_param.
```

The field poles are

```text
delta Z_L=-kappa H_L/2,
delta Z_R=-kappa H,
delta Z_phi=-kappa T_Y.                                 (5)
```

The coefficient of `-phi^dagger phi` in the actual counterterm action
is **not** the parameter pole by itself. Expansion of the bare mass term
gives

```text
delta m_H,coeff^2=delta m_H,param^2+m_H^2 delta Z_phi
                =-2kappa T_M,
delta m_H,param^2=kappa[m_H^2 T_Y-2T_M].                 (6)
```

For the Majorana mass orientation
`L_mass=-(bar N_L M N_R+h.c.)/2`, the bare-field expansion produces
the total displayed chiral mass counterterm

```text
delta M_coeff=delta M_param
              +delta Z_R^T M/2+M delta Z_R/2.
```

Its proper loop pole is zero, so

```text
delta M_param=kappa(H^T M+M H)/2,
delta M_coeff=0.                                       (7)
```

This distinction prevents confusing absence of a chirality-flipping
proper pole with absence of Majorana mass-parameter running.

## General complex symmetric mass and flavor covariance

For a constant complex symmetric `M` in that same chiral mass action,
the right-handed propagation is controlled by `M^dagger M`, as follows
by applying the paired chiral equations twice to `N_R`. The scalar
mass trace in (1)-(7) becomes

```text
T_M=tr(H M^dagger M)=tr(Y M^dagger M Y^dagger)
   =sum_(alpha i)|(Y M^dagger)_(alpha i)|^2 >=0.          (8)
```

It is `M^dagger M`, rather than `M M^dagger`, for this orientation.
The kinetic matrices and ordered symmetric mass pole in (5),(7) retain
their forms. With field coordinates `L'=V L`, `N_R'=U N_R`,

```text
Y'=V Y U^dagger, H'=U H U^dagger,
M'=U* M U^dagger,
delta M'=U* delta M U^dagger.
```

Both traces and the action counterterms transform accordingly. The real
diagonal input basis is enough for (1); running in fixed flavor
coordinates generally produces off-diagonal mass entries. Maintaining a
diagonal Takagi mass basis during running requires an additional scale-
dependent flavor-coordinate choice and the correlated Yukawa rotation.

## Minimal-subtraction scale relation

For clarity, the following restricted renormalizable scale inputs are
derived rather than inferred from the coefficient of a selected off-shell
pole. Use `Y_0=mu^epsilon(Y+delta Y)`, `d=4-2epsilon`, minimal subtraction,
and the convention `gamma_f=Z_f^(-1)mu dZ_f/dmu/2`, so
`mu df/dmu=-gamma_f f` at fixed bare field. At this order the tree term
`mu dY/dmu=-epsilon Y+...` differentiates every quadratic Yukawa residue
with factor `-2epsilon`. Therefore (5) gives

```text
16pi^2 gamma_L=H_L/2,
16pi^2 gamma_NR=H,
16pi^2 gamma_NL=H^T,
16pi^2 gamma_phi=T_Y.                                  (9)
```

For the dimensionful parameters, holding the bare quantities fixed and
using the same `-2epsilon` homogeneity on (6)-(7) gives

```text
16pi^2 beta_(m_H^2)=2m_H^2 T_Y-4T_M,
16pi^2 beta_M=H^T M+M H.                               (10)
```

These are only the pure-neutrino-Yukawa renormalizable contributions in
the declared fixed field coordinates and subtraction convention. The
proper Yukawa vertex and the field factors must be combined before
reporting the corresponding `beta_Y`; that is a separate calculation.
Neither (9) nor (10) is a dimension-seven physical beta function.

## Reproduction and limits

Run `python3 -I -S -B -O verify_yukawa_two_point.py`. The independent
standard-library implementation checks exact Gamma-pole recurrences,
symbolic Feynman-parameter polynomials, every two-gamma trace coefficient,
open-chain projectors, literal phases, closed Wick sign, weak multiplicity,
Majorana half normalization, complete kinetic flavor conjugation, and
complex-mass ordering/covariance. The optimized run passes **94 positive
exact checks** and rejects **12 deliberate mutations**. These counts
describe implementation controls of the displayed loop calculation, not
94 independent integrations. The adjacent JSON records all categories.

The UV integrals retain generic off-shell momentum until after extraction.
No finite threshold, on-shell pole mass, scalar vacuum/stability claim or
phenomenological bound is calculated. Gauge, other Yukawa, scalar
quartic and higher-operator contributions remain outside this restricted
sector. These new renormalizable inputs do not close the complete
dimension-seven counterterm vector or establish an observable prediction.
