# Restricted gauge-leg closure and the neutral photon projection

15 October 2026. This note collects the power-counting and projection
consequences of the evaluated two- and three-gauge residues. It does not
add a further loop integration or supply a finite scattering amplitude.

## Exactly which hierarchy is closed

Fix one literal derivative-seed insertion, one conjugate sterile Yukawa,
the `bar L phi N_L` matter orientation, one loop, and the lowest total
gauge-coupling order `n` for `n` labelled external gauge fields. The
connections are those of the canonical fermion/scalar kinetic actions.
There are no further operator insertions, spectator interactions, internal
gauge exchanges, external matter radiation or reducible gauge-tree
attachments in this stated hierarchy.

Let `f` count ordinary fermion insertions, `s` ordinary scalar insertions,
and `a` scalar seagulls. Then

    n=f+s+2a,
    V=2+f+s+a,       I=2+f+s+a,       L=I-V+1=1.

The crossed seed supplies two powers of loop momentum, fermion propagator
numerators supply `f+1`, and ordinary scalar vertices supply `s`.
Consequently the superficial UV degree is

    omega=4+(f+s+3)-2(f+s+a+2)=3-n.

The direct seed has no loop momentum and instead has degree `1-n`.
It is individually UV finite for `n>=2`. At `n=2`, the crossed pole has
external-momentum degree one; at `n=3`, degree zero. A scalar squared-mass
insertion lowers the UV degree by two, so neither new pole has a scalar
mass term. With `n>=4`, each member of this primitive one-loop matter-line
family is UV finite. There are no proper loop subdivergences in these
single-cycle diagrams. Original infrared singularities are not removed
by this UV statement.

For labelled gauge fields, the number of diagrams is

    N_n = sum_(f+s+2a=n) n! (s+a)! / [s! 2^a a!].

To obtain this expression, choose the fermion labels, scalar single
labels and unordered seagull pairs, then order the fermion vertices and
the scalar-line objects separately. The anticommutator in a seagull
already includes its two gauge-field assignments. This gives `N_2=7`
and `N_3=33`: 24 three-single-vertex graphs plus nine graphs containing
one seagull. The two independent integrations retain every diagram.

The zero-/one-gauge results, the newly evaluated two-gauge weak/mixed
kernel, and the independently integrated three-gauge kernel therefore
cover every potentially divergent multiplicity in this declared
hierarchy. This is a closure of one external-gauge hierarchy, not closure
of the EFT coefficient basis, all external matter sectors, higher gauge
orders or physical running.

## Why the single-photon cancellation does not remove every photon contact

Use `T^I=tau^I/2`, `y_L=-1/2`, and the declared neutral rotation
`e=g s_W=g' c_W`. Define

    Q=T^3+y_L I=diag(0,-1),
    W^+=(W^1-iW^2)/sqrt(2),   W^-=(W^1+iW^2)/sqrt(2),
    T^+=T^1+iT^2,             T^-=T^1-iT^2,
    X_A=e Q,                  X_+=g T^+/sqrt(2),
    X_-=g T^-/sqrt(2).

For constant `phi=(u,0)^T` and an external neutral lepton, the weak
matrix entry is `11`. A pure photon generator annihilates that column.
The previously calculated one-photon local pole vanishes; the newly
checked two-/three-gauge pure-photon kernels also have vanishing neutral
entries. In the cubic formula, commuting photon generators additionally
make both commutator sums vanish identically.

The presence of charged weak fields changes the matrix products. From
the actual three-gauge loop kernel, with `E_A=slash epsilon_A` and
similarly for `E_+,E_-`, the neutral entry is

    (R_AW+W-)_11 = e g^2 {
        [E_A,E_-] E_+/8
        +[(epsilon_A.epsilon_-)E_+
          +(epsilon_A.epsilon_+)E_-
          -2(epsilon_+.epsilon_-)E_A]/6 }.

The common `u F` and pole normalization remain as in the source action.
This is a labelled local vertex kernel; its permutation factors are
already included.

An independent matrix reduction explains the result. In the sigma term
only `X_+[X_A,X_-]` has a nonzero `11` entry, equal to `-e g^2/2`.
The nested commutators with remaining `X_+` and `X_-` each have `11`
entry `-e g^2/2`, while `[X_A,[X_+,X_-]]=0`. These give exactly the
displayed coefficients. Both the action audit and the direct 33-loop
sum check the projection. For an explicit nonzero algebraic example,
take `epsilon_A=e_1`, `epsilon_+=epsilon_-=e_2`, with spatial orthogonal
basis vectors; the kernel is `e g^2 slash(e_1)/12`.

This new contact does not determine a finite process, an on-shell
amplitude, a transition moment or an experimental signal. Connected
amplitudes also require the appropriate propagator-linked graphs,
external-state normalization, finite terms and input conventions.
The distinction is useful precisely because the one-photon UV result
alone cannot support a statement about every photon-containing vertex.

## Remaining scientific boundary

The loop matching here uses four-dimensional tensor/Clifford identities
for simple poles. A finite or higher-loop use needs a declared
dimensional/evanescent prescription. The earlier source-retaining field
change was derived for a selected hypercharge representative; it must
not silently be promoted to the full interacting non-Abelian physical
coefficient map. Other matter amplitudes, lower-order insertions,
correlated field/parameter/contact counterterms and physical matching
remain to be assembled.
