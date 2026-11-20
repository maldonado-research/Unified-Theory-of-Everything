# Active field-content projection of the selected sourced packet

16 October 2026. Input: candidate `4ed8560e` and the literal sourced
representative in
[`ELECTROWEAK_REDUCTION.md`](../../checkpoint_2026_10_15/interacting_reduction/ELECTROWEAK_REDUCTION.md).
This note projects that representative onto monomials containing **zero
retained sterile-field factors**, with literal neutrino-Yukawa sources and
the optional scalar quartic. It also states and checks the additional
remainder from explicitly retained Standard Model spectator Yukawas.

This is a field-content projection in specified action coordinates.
It neither integrates out the sterile field nor applies its equation of
motion. In particular, the original mixed counterterm has zero such
projection, whereas its source-reduced representative generally does not.
The difference comes from the field map, not from a newly computed loop
or a physical anomalous dimension.

## 1. Definition and exhaustive source substitution for the declared reference

Keep `phi=tilde H`, `F=C_seed Y_nu*`, `C_seed^T=C_seed`,
`X_alpha=phi^dagger L_alpha`, `bar X_alpha=bar L_alpha phi`,
`A=i slash partial`, `P0=−partial^2`, `rho=phi^dagger phi`, and
`c=2m_H^2/3+lambda rho/3`. The mass `M` here is the retained real positive
diagonal matrix of the reference action. No mass inverse is introduced.
The separate [complex-mass extension](COMPLEX_MASS_COVARIANCE.md) describes
the ordered products needed outside that mass convention.

Define `Pi_0N` by setting every component of `N_L,N_R` and all their
derivatives to zero in the polynomial action. Conjugate sterile components
are included in the same removal. The active composites remain arbitrary
spacetime-dependent fields; no neutral background or scalar equation is
imposed.

For the declared neutrino Yukawa action the sources are

    S_L=Y_nu phi N_R,
    bar S_L=bar N_R Y_nu^dagger phi^dagger,
    S_phi=−bar N_R Y_nu^dagger L,
    S_L^N=Y_nu^dagger X,
    S_R=Y_nu^T X^c.                                      (A1)

Thus the first three vanish under `Pi_0N`, while the last two do not.
The composites in the existing representative become

    R -> Y_nu^T X^c,
    T -> M Y_nu^dagger X+A Y_nu^T X^c,
    U -> (M^2+P0)Y_nu^T X^c+M A Y_nu^dagger X.          (A2)

The ordinary derivatives act on the whole neutral `X` or `X^c`. Neither
`A^2=P0` nor the substitution in (A2) removes derivatives on `phi`.

For completeness, every remaining class in the sourced packet projects
as follows.

| Original class | Reason its zero-`N` projection vanishes |
|---|---|
| `(K bar S_L)phi R` | `bar S_L` and its covariant derivatives contain `N` |
| `bar L S_phi R` | `S_phi` contains `N`, even when the term has four fermions |
| `bar S_L phi T` | Its row source contains `N` |
| `m_H^2 bar S_L phi N_L` | Both the source and explicit sterile slot contain `N` |
| `Q_Vh,Q_Eh` | Each has an explicit `N_L` |
| `j_mu J_Y^mu,j_I,mu J_I^mu` | The first current contains `N_L`, including when the other current is fermionic |

This table is exhaustive for the literal **neutrino-Yukawa-plus-quartic
reference**. Section 5 explains why its source conclusions must be extended
when spectator Yukawas are retained. Counting fermions alone would be
insufficient: the explicit generator test also counts how many belong to
the retained sterile field.

## 2. Complete projected packet in these coordinates

Introduce constant active-flavor matrices

    K=F Y_nu^T,
    K2=F M^2 Y_nu^T,
    G=F M Y_nu^dagger.                                  (A3)

The complete projection of this declared sourced packet is

    Pi_0N P_red =
        (5/6)bar X K P0 X^c
       +(5/6)bar X K2 X^c
       +(5/6)bar X G A X
       +(2m_H^2/3+lambda rho/3)bar X K X^c
       −(1/12)sum_(alpha beta a b)
          K_(alpha beta)bar L_alpha^a(Fcal_mu nu)_ab phi^b
                                         sigma^{mu nu}X_beta^c.   (A4)

The full counterterm action is
`kappa[(Pi_0N P_red)+(Pi_0N P_red)^dagger]`, with the original
`kappa=1/(16 pi^2 epsilon)`. The field-strength term is exactly the prior
active source descendant; its weak placement and full raw matrix `K`
must be retained. Its already established independent strength coordinates
are described in
[`ACTIVE_TENSOR_BASIS.md`](../../checkpoint_2026_10_15/operator_basis/ACTIVE_TENSOR_BASIS.md).
The other four structures in (A4) do not exhaust a general active EFT.

The canonical operator dimensions differ. `bar X X^c` has dimension five;
`bar X A X` has dimension six; `bar X P0 X^c`, `rho bar X X^c`, and the
strength term have dimension seven. Their coefficients carry the
corresponding retained mass powers, while all originate at the selected
inverse-scale order `Lambda^-3` with fixed retained-mass counting.

## 3. Statistics, integration by parts, and literal action normalization

For anticommuting fields the scalar contraction is symmetric:

    bar X_alpha X_beta^c=bar X_beta X_alpha^c.

The `rho` factor is commuting and leaves this property unchanged. For the
neutral Box operator its even integrated transpose gives

    integral bar X_alpha P0 X_beta^c
      =integral bar X_beta P0 X_alpha^c.                 (A5)

Consequently (A4)'s scalar and Box terms depend on

    K_S=(K+K^T)/2,       K2_S=(K2+K2^T)/2.              (A6)

Equation (A5) is an action identity modulo a boundary term, not a
pointwise exchange of arbitrary derivatives. The strength term has
different spin statistics and weak slot structure; replacing its entire
`K` by `K_S` would discard a real channel.

For the lepton-number-conserving neutral derivative term,

    integral[(bar X G A X)+h.c.]
       =integral bar X(G+G^dagger)A X,
    G_H=(G+G^dagger)/2.                                 (A7)

The displayed Hermitian action on the right is included **once**.
Appending another Hermitian conjugate would double it. Derivatives act on
the scalar inside `X` as well as on the lepton.

These flavor quotients have the following coordinate counts within each
specified operator shape: a complex symmetric scalar or integrated Box
coefficient has `n_L(n_L+1)/2` complex entries, and a Hermitian LNC
coefficient has `n_L^2` real entries. Separate exact component ranks at
`n_L=1,2,3` check these statements. They are not an assertion that every
matrix is generated by this selected seed or that these shapes remain
independent after a further full EFT EOM quotient.

The literal input action uses `+(1/2)C7 bar X P0 X^c+h.c.`. To state
the normalization unambiguously, define analogous half-normalized scalar
coordinates and one Hermitian LNC coordinate as follows.

| Chosen action term | Raw coefficient projection from (A4), divided by `kappa` |
|---|---|
| `+(1/2)delta C7_box bar X P0 X^c+h.c.` | `(5/3)K_S` |
| `+(1/2)delta C5 bar X X^c+h.c.` | `(5/3)K2_S+(4m_H^2/3)K_S` |
| `+(1/2)delta C7_rho rho bar X X^c+h.c.` | `(2lambda/3)K_S` |
| `+bar X delta C6_X A X`, Hermitian action included once | `(5/6)(G+G^dagger)=(5/3)G_H` |

This table is a **coefficient projection in the chosen sourced
coordinates**, not a physical renormalization constant or beta function.
In particular the original mixed packet projects to zero before its field
map. The table must not be interpreted as a newly calculated closed
active-sector anomalous dimension. Other action terms and field changes
can populate these same coordinates.

## 4. An independent check from the actual source-map field substitution

At `N_new=0`, the source-reducing map has no lepton, scalar, or gauge
shift, but its sterile shift remains nonzero. Set

    U_i=(F^dagger X)_i,
    x_i=(5/6)(P0+M_i^2)U_i+c U_i
          −(1/12)sum_(alpha a b)F*_(alpha i)(phi^b)*
                     (Fcal_mu nu)_ba sigma^{mu nu}L_alpha^a,
    y_i=(5/6)M_i A U_i,
    delta N_Ri=−(y_i+x_i^c).                             (A8)

Then directly substituting the old sterile field into the reference
Yukawa action gives

    Pi_0N delta_vS S0
       =−bar X Y_nu delta N_R+h.c.
       =(Pi_0N P_red)+h.c.,                             (A9)

where equality includes the scalar statistics and integrated derivative
identities above. The sterile kinetic and mass actions have no
first-order term at `N_new=0`; both their slots start at the shift and
would be second order. No sterile mass inverse appears. The checker
substitutes (A8) explicitly, including charge conjugation and the curvature
term, and compares it with (A4).

Equation (A9) explains the coordinate dependence constructively. The
surface `N_S=0` maps to nonzero old sterile composites. It is not the same
field-content surface as `N_old=0`. The physical equivalence of full
actions does not require their coordinate-specific projections to agree.

## 5. Explicit spectator-Yukawa remainder

Setting `N=0` does not delete Standard Model matter sources. With the
spectator Yukawas of
[`SPECTATOR_YUKAWA_TRANSPORT.md`](../action_completion/SPECTATOR_YUKAWA_TRANSPORT.md),
the additional sources are

    (S_Le)_alpha^a=(Y_e)_alpha p H^a e_Rp,
    (S_phi,spec)^b=
        sum_(alpha p a)bar L_alpha^a(Y_e)_alpha p epsilon_ab e_Rp
       +sum_(r s a)bar Q_r^a(Y_d)_r s epsilon_ab d_Rs
       −sum_(r s)bar u_Rs(Y_u)_r s* Q_r^b.               (A10)

Here `H=−epsilon phi*` is the same single scalar doublet. The complete
additional projected source remainder is

    Delta P_spec=(1/3)[(Krow bar S_Le)phi S_R]_F
                   −(1/6)[bar L S_phi,spec S_R]_F,
    Krow bar psi=−i(D_mu bar psi)gamma^mu.              (A11)

Indeed `bar S_Le phi=bar e_R Y_e^dagger H^dagger phi=0`, which removes
the undifferentiated row-source contraction. It does not remove
`(Krow bar S_Le)phi`: a derivative acting on a zero product must include
the derivative of the other scalar factor. The first term in (A11)
therefore survives in general. The second includes active four-fermion
terms from all three displayed spectator Yukawas.

The checker constructs explicit `e_R,d_R,u_R,Q_L` exterior fields and the
actual constrained `H`, substitutes these sources in the literal
representative, and compares its difference from (A4) with (A11). It also
adds their source terms to (A8) before checking the reference-Yukawa
variation. This verifies the displayed spectator remainder without
claiming exhaustion for arbitrary additional interactions or reducing it
into a complete active derivative/four-fermion basis.

## 6. Projection does not commute with the additional scalar-current map

Let `d=(g^2+g'^2)/24`. The previously specified further scalar-current
map is

    delta L=−d rho Phi F N_L,
    delta N_L=+d rho F^dagger X,
    delta N_R=(delta N_L)^c,
    delta phi=delta B=delta W=0.                         (A12)

Before this map the entire scalar-current subpacket contains `N` and
has zero `Pi_0N` projection. After substituting (A12) in the **reference
Yukawa action**, however, its first-order zero-`N` term is

    Pi_0N delta S_Y
      =−d rho bar X Y_nu F^T X^c+h.c.
      =−d rho bar X K_S X^c+h.c.                        (A13)

The second equality uses scalar statistics: `Y_nu F^T=K^T` has the same
contraction as `K`. This agrees exactly with the active source term in
the already reduced scalar-current subpacket. It is nonzero in general,
even with one active flavor. The extra change therefore shifts the raw
half-normalized `C7_rho` coordinate in the table by
`−2 kappa d K_S`, while leaving the fact of action equivalence intact.

Equation (A13) accounts for the failure to commute by an explicit source
field variation. It supplies no threshold matching, sterile elimination,
or physical running. Applying `Pi_0N` before composing a field map would
lose that variation.

## 7. Exact verification and remaining work

Run `python -I -S -B -O verify_active_source_projection.py` here. The
standalone checker uses only the standard library, exact Gaussian rational
Fourier coefficients, explicit Weyl gamma matrices, constrained Majorana
completion, and exterior products through four fermions. It imports no
earlier calculation modules or saved coefficient outputs.

The original side substitutes literal sources in **every class** of
`P_red` before counting retained-`N` generators. The compact side is built
independently from (A4). Nonconstant scalars contain opposite non-null
Fourier modes, so omitting scalar derivatives is detectable under the
integral. Noncommuting weak connections generate the curvature explicitly.
Additional checks establish integrated scalar symmetry, LNC Hermiticity,
literal action normalization, source-map Yukawa transport, the spectator
remainder, and the noncommuting scalar-current projection.

Negative controls remove retained-mass or source pieces, insert an extra
seed half, delete derivatives on the scalar, replace the Box coefficient
by its antisymmetric part, discard the weak flavor channel, omit Hermitian
conjugation, assume that all spectator sources vanish at `N=0`, or reverse
or omit the induced scalar-map source. The adjacent JSON records actual
counts, retained-`N` and total fermion degrees, and operator ranks.

The completed isolated optimized replay passes **48 positive controls**
and rejects **19 mutations**, for **67 checks total**. Its recorded runtime
is approximately 142 seconds in the execution environment. The scalar and
integrated Box ranks are `1,3,6`, and the Hermitian LNC ranks are `1,4,9`,
at one, two, and three active flavors respectively. These counts describe
algebraic checks within this calculation, not additional loop evaluations.

The resulting projection is a useful exact action-coordinate inventory.
Complete derivative/EOM basis matching, contributions from other retained
coefficients and their field maps, finite evanescent conventions, and
physical running remain separate tasks. No new loop or observable is
calculated here.
