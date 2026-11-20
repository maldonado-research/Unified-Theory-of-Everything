# Independent active field-content projection audit

16 October 2026. Input: candidate.9, `4ed8560e`, and the source-retaining
selected packet. This derivation uses the declared **neutrino-Yukawa reference
action**, optionally including the scalar quartic, with constant real symmetric
retained Majorana mass. It was derived independently of the producer's new
projection. It performs no loop integration or sterile-field elimination.

The complete literal zero-retained-N projection of that representative is
nonzero. It contains a neutral-composite two-derivative scalar bilinear,
retained-mass and scalar-potential terms, a lepton-number-conserving composite
derivative bilinear, and the already specified matrix-curvature source.
The original unreduced selected packet has zero literal projection. These
facts are compatible because its field redefinition mixes active and sterile
coordinates: the two zero-N subspaces are different subspaces of field space.

## Definition of the projection and exhaustive source substitution

Let `Pi_0` set `N_L,N_R` and all their derivatives to zero. It is an
algebraic field-content substitution. It neither solves an equation of motion
nor performs a path integral, matching step, or Schur elimination.

Use `phi=tilde H`, `X_alpha=phi^dagger L_alpha`, `rho=phi^dagger phi`,
`A=i slash partial`, `P0=A^2=-Box`, and `F=C_seed Y*`. The full weak
connection is `D=partial+i Acal`. Its curvature is `Fcal`, distinct from
the flavor coefficient `F`. With the stated neutrino-Yukawa sources,

    Pi_0(S_L)=Pi_0(bar S_L)=Pi_0(S_phi)=0,
    Pi_0(S_L^N)=Y^dagger X,
    Pi_0(S_R)=Y^T X^c.                                      (P1)

The sterile sources in the second line must not be discarded. Thus the
actual polynomials in the selected representative become

    Pi_0(R)=Y^T X^c,
    Pi_0(T)=M Y^dagger X+Y^T A X^c,
    Pi_0(U)=M^2 Y^T X^c+M Y^dagger A X+Y^T P0 X^c.           (P2)

All Yukawa and mass matrices are constant. Their flavor order is displayed;
no commuting-flavor assumption is made.
The separately reviewed [complex-mass extension](COMPLEX_MASS_COVARIANCE.md)
replaces `M^2` by `M M^dagger` in these projected mass terms and retains
`F M Y^dagger` in the LNC term. The new standalone audit here deliberately
keeps its declared real symmetric mass convention.

Every line of the representative is accounted for:

| Original family | Projection |
| --- | --- |
| `5/6[bar L phi U]_F` | Three terms from the last line of (P2) |
| `(2m_H^2/3+lambda rho/3)[bar L phi R]_F` | Active scalar source |
| `1/3[(K bar S_L)phi R]_F` | Zero, since `bar S_L` and its derivatives vanish |
| `-1/6[bar L S_phi R]_F` | Zero, since `S_phi` vanishes |
| `-1/6[bar S_L phi T]_F` | Zero |
| `m_H^2/2[bar S_L phi N_L]_F` | Zero |
| Higgs-derivative curvature pair | Zero, from its explicit `N_L` slot |
| Gauge-current products | Zero, from the explicit `N_L` in `j_mu,j_mu^I` |
| `-1/12[bar L Fcal phi sigma R]_F` | Active curvature source |

This exhaustion is complete for (P1). It is not an exhaustion of an
unspecified spectator theory. For example, `bar S_L` can contain
`bar e_R Y_e^dagger H^dagger`, which does not vanish merely by setting
retained N to zero; derivatives can also prevent an `H^dagger phi=0`
simplification. Any such retained spectator sources require their own
additional projected terms.

## Complete projected representative and normalizations

Define the active flavor matrices

    K=F Y^T,       K2=F M^2 Y^T,       V=F M Y^dagger,
    K_S=(K+K^T)/2,                    (K2)_S=(K2+K2^T)/2.

Before choosing independent flavor pairs, the projected non-Hermitian
orientation is

    Pi_0(P_red)=5/6 bar X K2 X^c
               +5/6 bar X V A X
               +5/6 bar X K P0 X^c
               +(2m_H^2/3+lambda rho/3)bar X K X^c
               -1/12 sum_alpha beta K_alpha beta
                    bar L_alpha Fcal_mu nu phi
                                  sigma^{mu nu}X_beta^c.       (P3)

Append the Hermitian conjugate once in the action, together with the same
common `kappa`. The derivative is on the **whole neutral composite** `X^c`,
including the scalar. Replacing it with a derivative on the lepton alone
changes (P3).

The scalar charge-conjugation spin contraction is antisymmetric. Grassmann
exchange therefore gives a flavor-symmetric scalar bilinear. The first
term in (P3) uses `(K2)_S` pointwise. The fourth uses `K_S` pointwise.
For the third, exchanging the two fermion slots and integrating the even
derivative twice shows

    integral bar X K P0 X^c=integral bar X K_S P0 X^c.       (P4)

Its ordered-flavor coefficient is **`+5 K_S/6` for `P0=-Box`**. If the
operator is defined with literal `+Box`, the corresponding coefficient is
`-5 K_S/6`. This is a coordinate coefficient of the projected representative,
not a physical Wilson prediction. In an independent-pair convention the
off-diagonal scalar/Box coefficients double, whereas the diagonal entries
retain their ordered coefficients. Hermitian conjugation does not introduce
an additional factor into this fixed chiral orientation.

The lepton-number-conserving derivative term has a different projection:

    integral [bar X V A X+h.c.]
       =integral bar X(V+V^dagger)A X.                     (P5)

It involves a Hermitian coefficient after including its conjugate, not a
flavor-symmetric transpose projection. Its derivative acts on `X`, and all
scalar derivatives implicit in that composite remain.

For the curvature term, write

    Fcal_mu nu=g'y_L B_mu nu 1_2+g T^I W^I_mu nu,
    T^I=tau^I/2,       y_L=-1/2.

The B channel selects `K_A=(K-K^T)/2`. The W channel requires simultaneous
weak/flavor exchange. With `P_ab^I=(T^I phi)_a phi_b`, its relevant tensor is

    K_A tensor P_S^I+K_S tensor P_A^I.                    (P6)

The symmetric flavor channel in (P6) cannot be omitted. In particular it
has a nonzero one-flavor charged component on `phi=(u,0)`. The source
coefficient remains `-1/12` with the original full nonlinear curvature.
No derivative EOM has been used to rearrange other terms in (P3) into this
explicit curvature coordinate.

## Why a nonzero active projection is consistent with the original action

The original selected packet is linear in retained N or its derivatives,
so `Pi_0(P_orig)=0` identically. The source-reducing map does not preserve
that slice. In its integrated residual certificate, `ell,t,j,j^I` vanish
on the new zero-N slice, but the sterile coefficients do not.

Let `W_i=sum_alpha F_alpha i bar X_alpha`, and
`K0 W=-i(partial_mu W)gamma^mu`. On this slice the actual row covectors are

    bar x_i=5/6[(W M^2)_i+K0^2 W_i]
       +(2m_H^2/3+lambda rho/3)W_i
       -1/12 sum_alpha F_alpha i bar L_alpha
                                 Fcal_mu nu phi sigma^{mu nu},
    bar y_i=5/6(K0 W)_j M_ji.                             (P7)

The sterile residuals are `e_R=-S_R`, `e_L=-S_L^N`. Substitution in the
off-shell certificate gives, modulo its already declared boundary term,

    Pi_0(P_red)=bar x_i S_Ri+bar y_i S_Li^N.               (P8)

Equation (P8) reproduces every term of (P3), including the derivative on
the complete source and the ordered mass factors. The full Majorana
variation with `N_old=N-kappa(w+w^c)` supplies (P8) and its conjugate.
There is no extra Majorana half factor. Thus the discrepancy between the
two literal projections is exactly a variational source effect of the
coordinate change, not a new active interaction inferred by integrating
out a sterile particle.

## Projection does not commute with the additional scalar-current EOM step

The original scalar-current subpacket contains `j_mu` or `j_mu^I`, each
linear in N, and therefore projects to zero. The additional reduction in
`SCALAR_CURRENT_REDUCTION.md` has the exact residual

    P_scalar-P_scalar,red=i d partial_mu(rho j^mu)
        +d rho([e_barL phi N_L]_F-[bar L phi e_R]_F),
    d=(g^2+g'^2)/24.

Taking the projection **after** this reduction yields

    Pi_0(P_scalar,red)=-d rho bar X K X^c
                     =-d rho bar X K_S X^c.               (P9)

Taking it before the reduction yields zero instead. The projected residual
is `+d rho bar X K X^c`, precisely the negative of (P9), so nothing has been
lost from the action identity. In that further representative the explicit
coefficient of `rho bar X K_S X^c` is `lambda/3-d` rather than `lambda/3`.
The two coefficients belong to different field-coordinate choices.

There is also a direct source-action witness. The additional scalar-current
map has, on the zero-N slice,

    delta N_L=+kappa d rho F^dagger X,
    delta N_R=(delta N_L)^c,
    delta L=0.

Inserting it into the neutrino Yukawa monomial and its conjugate gives
exactly (P9) and its conjugate; the symmetric Grassmann scalar identity
changes the apparent `Y F^T` ordering to the equivalent `F Y^T` orientation.
The kinetic and retained-Majorana variations have an explicit old N factor
and vanish on the slice. This is an actual variational certificate for the
noncommuting operations, not a free-field EOM argument.

## Independent exact reproduction and boundary of the result

Run with SymPy 1.14.0:

```sh
python -I -B -O independent_active_audit.py > INDEPENDENT_ACTIVE_RESULTS.json
```

The core imports no repository modules and uses no Python assertions. It
constructs explicit chiral Dirac/Pauli matrices and exact sparse Fourier
polynomials with exterior Grassmann coefficients. The active fields have
two opposite generic four-momenta, both weak components and both spin
components. Each scalar component has three Fourier modes. There are 32
independent Grassmann generators. Scalar derivatives in `P0 X^c` are thus
tested explicitly, rather than hidden by a constant background.

The test uses actual `F=C_seed Y*`, complex noncommuting flavor matrices,
and real positive symmetric `M=[[2,1],[1,3]]`. Generic Hermitian weak
curvature jets test matrix placement. It constructs (P2), evaluates every
line in the exhaustion table, and compares the full resulting density with
(P3). Separate controls test scalar/Box exchange, flavor-pair factors,
the Hermitian LNC term, B/W combined-slot exchange, the actual sterile
covector certificate (P8), and the literal Yukawa variation giving (P9).

The saved optimized run records **56 exact checks** and **15 rejected
negative controls**. Before integration the active density has 2,016
Fourier/Grassmann monomials; its zero-total-momentum component has 176.
Negative controls remove derivative, retained-mass, quartic, LNC or weak
source terms, commute flavor matrices, lose scalar derivatives, discard
diagonal scalar bilinears, equate different zero-N slices, or force
projection and scalar-current EOM reduction to commute.

This is complete field-content projection for the stated source theory,
not a complete derivative/EOM quotient, a sterile matching result, or a
physical beta function. A source-dressed Schur calculation could give an
additional implementation check using all mixed and sterile blocks and
transported probes, but the general cancellation is already established by
the prior congruence/Schur identities. The active block alone cannot provide
that calculation. No new physical cancellation or observable is claimed.
