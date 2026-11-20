# Neutral projection and canonical-coordinate checks

15 October 2026. Restricted first-order quadratic calculation for the selected
packet. This is an application of standard field redefinition and matrix
identities, not new fundamental mathematics, complete EFT matching, physical
running, or a dipole prediction.

## 1. The actual selected packet has a simple neutral quadratic projection

Use the literal reference action and all six placements of `O_J` in
`checkpoint_2026_10_15/operator_map/ACTION_SIGN.md`, with the cancelling action

    delta L = kappa F_(alpha i) (O_J+O_D)_(alpha i) + h.c.,
    F=C_seed Y_nu*,     C_seed^T=C_seed.

Here `kappa` is a real first-order bookkeeping parameter; for this selected
counterterm its coefficient is `1/(16 pi^2 epsilon)`. For this quadratic
projection only, set the gauge backgrounds to zero and put
`phi=tilde H=(u,0)^T`, with constant real `u`. This is a background choice,
not an assertion that `u` solves the scalar equation. Write the active neutral
left field as `nu_L`; the retained Majorana left component is `N_L`.

`O_D` contains the gauge field strength or its derivative and gives no
quadratic term on this background. `-D^2 phi=0`, and integration by parts
places each remaining row derivative on the sterile field. On a common
momentum `p`, the six `O_J` placements give

    (5/6) slash p p^2,
    (1/3) p^2 slash p,
    0,
    -(1/6) slash p p^2,
    (m_H^2/2) slash p,
    (m_H^2/2) slash p.

Thus the original, unreduced packet projects to

    delta L^(2) = kappa u F_(alpha i)
                  bar nu_Lalpha A(-Box+m_H^2) N_Li + h.c.       (N1)

The coefficient is one: `5/6+1/3-1/6=1` and `1/2+1/2=1`.
There is no additional factor `1/2` on this mixed active-sterile bilinear.
The `1/2` in the four-component Majorana kinetic term becomes one ordinary
left-Weyl kinetic term; it does not halve the mixed block. This projection
was also checked independently by the separate action audit.

For `q=(nu_L,N_L)^T`, with all entries left-Weyl fields, the kinetic matrix is

    K(p^2)=I+kappa Z(p^2),
    Z(p^2) = [ 0                       u(p^2+m_H^2)F ]
             [ u(p^2+m_H^2)F^dagger     0             ].        (N2)

`Z=Z^dagger` for real `p^2`, and `tr Z=0`. In position space its real
polynomial multiplier is `u(-Box+m_H^2)`. This is a local derivative
operator, not an independent constant wavefunction coefficient. The
one-derivative part alone is `Z_as(0)=u m_H^2 F`; using only that part while
discarding the induced higher-derivative mass/source terms would change the
selected quadratic problem.

The signed scalar parameter remains exactly `m_H^2` from the reference
action. With the optional potential `-m_H^2 rho-lambda rho^2`, a nonzero
stationary background obeys `m_H^2+2 lambda u^2=0`. No such condition is
imposed in (N1)-(N2). In the literal action with `lambda=0`, a generic
nonzero `u` is not a vacuum.

## 2. Declare the mass orientation before transforming flavor matrices

Use the holomorphic left-Weyl convention

    L_mass=-1/2 [q^T epsilon Mcal q+h.c.],  Mcal=Mcal^T,
    epsilon=i sigma_2.

With the charge-conjugation convention of `ACTION_SIGN.md`, this is
`-1/2[bar q_L^c Mcal q_L+h.c.]`. The reference Dirac monomial is instead
`-bar nu_L (uY_nu) N_R+h.c.`. Therefore

    Mcal_ref = [ 0          uY_nu* ]
               [ uY_nu^dagger M*  ].                         (N3)

Here the retained `M` is the real positive diagonal matrix of the actual
checkpoint, so `M*=M`. More generally a symmetric complex four-component
mass in `-1/2[bar N_L M N_R+h.c.]` becomes the lower holomorphic block `M*`.
The antiholomorphic four-component mass block in
`-1/2[bar q_L M_four q_L^c+h.c.]` is
`M_four=[[0,uY_nu],[uY_nu^T,M]]`; `Mcal=M_four*`. Conjugating `Y_nu` in
(N3) is essential. Subsequent formulas use the declared `Mcal` and should
not be transferred to the other convention by silently keeping the same
matrix names.

The seed itself also gives an active higher-derivative Majorana monomial
on this background. If the reference quadratic problem is retained exactly,
it must appear in its momentum-dependent `Mcal(p^2)` and in the generalized
kernel below. For a calculation linear in the seed, (N2) is already linear
in `C_seed`, so its product with that seed monomial is second order in
`C_seed`. The induced shifts below may consequently use the renormalizable
matrix (N3) at that declared linear-seed order. This is order counting, not
an assertion that the tree seed vanishes.

## 3. Canonical kinetic normalization induces correlated mass/source shifts

For arbitrary Hermitian `Z` and symmetric `Mcal`, let

    q=T q_c,      T=I-kappa Z/2.

At first order,

    T^dagger K T=I+O(kappa^2),
    Mcal_c=T^T Mcal T
           =Mcal-kappa/2 (Z^T Mcal+Mcal Z)+O(kappa^2).        (N4)

The mass uses a transpose, not an adjoint, and remains symmetric. This
holds for the actual polynomial `Z(p^2)` point by point; in position space
the transpose includes integration by parts. Its derivatives have even
order, so the displayed momentum polynomial is unchanged. The kinetic
term can thus be made canonical by a local perturbative derivative
redefinition, but (N4) then contains momentum-dependent Majorana mass
terms. Those terms have not disappeared and must be retained in a pole or
matching calculation.

More explicitly, write `Z=[0 H;H^dagger 0]` and
`Mcal=[A B;B^T D_s]`. The coefficients of `kappa` in its block shifts are

    delta A   =-(H* B^T+B H^dagger)/2,
    delta B   =-(H* D_s+A H)/2,
    delta D_s =-(H^T B+B^T H)/2.                           (N4a)

For the actual renormalizable reference, `A=0`, `B=uY_nu*`,
`D_s=M*`, and `H=u(p^2+m_H^2)F`. Therefore

    delta A   =-u^2(p^2+m_H^2)(F*Y_nu^dagger+Y_nu*F^dagger)/2,
    delta B   =-u(p^2+m_H^2)F* M*/2,
    delta D_s =-u^2(p^2+m_H^2)(F^T Y_nu*+Y_nu^dagger F)/2. (N4b)

These correlated entries are coordinate coefficients. In particular a
nonzero `delta A` cannot by itself be read as a new physical active mass:
the full congruence preserves the mass matrix rank and the original
generalized pole equation. A pre-existing zero mode of a pure kinetic
perturbation remains a zero mode for invertible kinetic normalization.

If a row source vertex is `J^dagger G q+h.c.`, then

    G_c=G T=G-kappa G Z/2.                                 (N5)

For a holomorphic spinor source `J^T epsilon G q+h.c.`, the same row
flavor rule holds. A column coefficient multiplying `q^dagger` instead
transforms with `T^dagger`; one must declare which vertex is being used.
For derivative `T`, transport the derivatives as well. A mass shift
without the correlated vertex/source shift is not an equivalent quadratic
theory.

Canonical normalization leaves flavor-coordinate freedom. The general
first-order solution is

    T=I+kappa(-Z/2+Omega),      Omega^dagger=-Omega.

It adds `kappa(Omega^T Mcal+Mcal Omega)` to (N4) and `kappa G Omega` to
(N5), while preserving the kinetic term. Coordinates related this way
give the same poles and source-dressed response. An individual mass or
Wilson-matrix entry need not be invariant.

## 4. Pole and external-response certificates require no mass inverse

At rest, use `E` as momentum and a doubled Weyl/charge-conjugate flavor
kernel (the two-component epsilon factors have been absorbed in the
charge-conjugate component):

    D(E) = [ E K       -Mcal^dagger ]
           [ -Mcal      E K^T       ].                     (N6)

Define `T_N=diag(T,T*)`. Then, for constant or real-polynomial coefficients,

    D_c=T_N^dagger D T_N,
    det D_c=|det T|^4 det D
           =(1-2 kappa tr Z)det D+O(kappa^2).               (N7)

The factor is a perturbatively invertible coordinate Jacobian. It does
not move physical low-energy poles. For the actual off-diagonal `Z` in
(N2), `tr Z=0`, so the determinants agree through first order without any
rescaling. A momentum-dependent change of variables is understood as a
formal derivative expansion in its regime of validity. Solving the
truncated polynomial nonperturbatively and interpreting high-energy extra
roots as new particles is not justified.

For `K` perturbatively invertible, eliminating the conjugate equation also
gives the smaller generalized squared-pole kernel

    Q(s)=s K-Mcal^dagger K^(-T) Mcal,       s=E^2,
    Q_c=T^dagger Q T,
    det Q_c=(1-kappa tr Z)det Q+O(kappa^2).                 (N8)

Only a kinetic inverse is used: `K^(-T)=I-kappa Z^T+...`. No sterile or
active mass inverse is needed, so massless fields and degenerate masses
are allowed. For momentum-independent `Mcal`, the canonical squared-mass
matrix has first-order perturbation

    delta(Mcal_c^dagger Mcal_c)
       =-Mcal^dagger Z^T Mcal
        -1/2 {Z,Mcal^dagger Mcal}.                         (N9)

This formula also evaluates a momentum-dependent problem at each `s`,
but its roots then follow the full `s`-dependent equation, rather than a
single constant squared-mass matrix. In a degenerate subspace, diagonalize
the restricted perturbation or use (N6)/(N8); do not divide by zero mass
splittings. Exact invertible mass congruence preserves rank: a null vector
`v` transports as `T^(-1)v`. A finite numerical matrix retaining only the
first-order mass shift can have spurious higher-order rank changes;
rank claims must respect the truncation.

For an arbitrary doubled source row matrix `W`, let `W_c=W T_N`. Away
from a pole,

    W_c D_c^(-1) W_c^dagger = W D^(-1) W^dagger.            (N10)

This is a useful stronger check than comparing spectra: it detects omitted
external/source normalization even when masses alone happen to agree.
The inverse here is the off-shell propagator kernel, which can be
invertible at generic `E` even when `Mcal` is rank deficient.

## 5. Zero momentum: Schur reduction has an explicit domain

Write a symmetric holomorphic mass matrix and row source as

    Mcal=[A B; B^T D_s],       G=(G_a,G_s).

Only if `D_s` is invertible, define

    S=A-B D_s^(-1) B^T,
    G_eff=G_a-G_s D_s^(-1) B^T,
    C_J=G_s D_s^(-1) G_s^T.                                (N11)

Completing the square with
`q_s=y_s-D_s^(-1)B^T q_a` block-diagonalizes the mass by congruence.
The reduced source has the vertex `G_eff` and the associated local
source-source contact `C_J/2`. Thus

    det Mcal=det D_s det S,
    rank Mcal=rank D_s+rank S.

These statements need no inverse of `S`. For the checker's actual
two-active/one-sterile renormalizable matrix, `S` is singular and one
massless mode remains; no fictitious inverse of that block is taken.
If `S` is also invertible, the zero-momentum response satisfies

    G Mcal^(-1) G^T = C_J+G_eff S^(-1) G_eff^T.             (N12)

The checker verifies (N12) only for an explicitly declared invertible
active-mass spectator benchmark. It is not a coefficient computed from
the selected packet. Omitting either `C_J` or the shifted source vertex
fails this identity. If `D_s` is rank deficient, (N11) is unavailable:
retain the full kernel or declare a separately justified restricted
subspace; do not substitute an undeclared pseudoinverse. With kinetic
mixing, the split into an ``active'' mass block and a ``sterile'' mass
block is coordinate dependent. The source-dressed response identifies
what an external probe measures.

## 6. Relation to the full source-retaining field redefinition

Equations (N1)-(N12) canonically normalize the original packet's quadratic
neutral projection. They do not certify that the Hessian of its complete
source-reduced representative is obtained solely by a linear fermion
Jacobian on an arbitrary background. In the full residual construction,
the scalar shift starts at fermion degree two. Its contraction with a
nonzero scalar tadpole contributes at fermion degree two as well.
Consequently that term must be retained when comparing quadratic actions
away from a stationary background. The separate field-redefinition derivation
independently identified this caveat. At a stationary scalar background,
the tadpole vanishes, but that condition must follow from a declared
scalar potential and background equation.

The selected counterterm is one restricted insertion. A calculation of
physical running additionally requires the full retained coefficient set,
all relevant diagrams and field/parameter counterterms, an independent
operator basis and its lower-order closure, and a declared treatment of
dimension-dependent/evanescent and finite terms. Neither canonical
normalization here nor the previously obtained representative coefficient
`-1/12` supplies those missing ingredients.

## 7. Reproduction and discriminating controls

Run

    python -I -B -O verify_canonical.py

The script imports only SymPy and the Python standard library, and no
repository calculation modules. It uses exact complex rational matrices,
the literal six-placement Fourier kernel, the actual `F=C_seed Y_nu*`
orientation, a rank-deficient two-active/one-sterile reference mass, three
off-shell energies, and a separate nonzero-trace Hermitian control.
Determinant derivatives are computed by column multilinearity, so they
require no mass inverse. The doubled kernel, the independently eliminated
Weyl equation, and source-dressed propagators supply distinct checks.

Negative controls discard row derivatives, induced mass terms, derivative
mass terms, external/source transport, the Majorana transpose, and the
Schur source contact or vertex. These checks establish algebraic
consistency and discriminate concrete errors; they are not independent
loop computations or proofs of EFT completeness. Results are saved in
`CANONICAL_RESULTS.json`.
