# One hypercharge leg of the routed derivative seed

15 October 2026. Local research candidate; not a physical beta function.

This calculation extends the public routed seed by attaching one abelian
gauge boson to each of its two charged internal lines. The two loop
integrands and their relative charges are declared explicitly below. Their
ultraviolet residues determine a transverse local polynomial that was
invisible to the previous zero-gauge three-point restriction. Common
action-to-amplitude phases and flavor normalization must be retained when
embedding this result in the full counterterm action.

## Conventions and restricted graph sector

Use Minkowski metric (+---), dimension d=4-2 epsilon, and strip the common
i/(16 pi^2 epsilon) from loop integrals. Set p=r+h+k: p is the incoming
neutral sterile momentum, and r,h,k are the outgoing lepton, Higgs and
hypercharge-boson momenta. All momenta are initially off shell. Let m2 be
the signed scalar mass-squared parameter in q^2-m2; the internal lepton
is massless. The sterile does not circulate in these graphs, so its mass
does not belong in either loop denominator.

The seed is the symmetric-flavor singlet-composite derivative operator
with reduced four-point polynomial
epsilon(ac)epsilon(bd)(r_a+h_c)^2 +
epsilon(ad)epsilon(bc)(r_a+h_d)^2. Its contraction with the Yukawa vertex
has the already audited relative weights 2:1. Write its coefficient as
C_(7,seed), defined by this positive momentum-squared polynomial. The
literal action +J^T W1 Box J/2 instead has C_(7,seed)=-W1 in the public
Fourier convention. The relevant flavor factor is C_(7,seed) Y_nu^* in
the public four-component convention, with its conjugate channel also
required.

Only one W1 seed, one sterile Yukawa vertex, and one ordinary hypercharge
vertex occur. The charged internal lepton and Higgs have Q_L=-Q_H. No
direct seed gauge vertex is introduced: J=epsilon LH is a singlet, and
the derivative acts on the complete J. These are the two 1PI internal
attachments to this seed bubble. External-leg radiation is reducible and
is not included in this 1PI vertex. No extra scalar quartic, Yukawa pair,
dimension-six insertion, or gauge loop is included.

Define t=p-k, v=h+k, u=2p-k, and slash a=gamma_mu a^mu. The open gamma
chain uses the public ordering inherited from the zero-gauge seed. Its
rightmost chiral projector is suppressed in the algebra below and can
be restored unchanged. The dimension-five source chain was freshly checked against [Zhang, Appendix A, Eq. (35)](https://arxiv.org/html/2405.18017v3#A1.E35) on 15 October 2026. Its `-3/2` reference pole and both sourced descendants in Eqs. (37)–(38) agree with the earlier calibration. This check does not fix the full four-point action-level sign. The source dimension-five open chain reads, from
the outgoing lepton toward the sterile, seed P_R, lepton propagator,
Yukawa P_L. The propagator adjacent to the seed carries q-k after photon
emission; the propagator adjacent to the Yukawa vertex carries q. This
fixes (slash q-slash k) gamma^nu slash q, rather than its reverse. A
reversed chain can still pass the longitudinal contraction, so that test
alone cannot fix the order. There is no closed chiral trace. Only the simple
UV residue is extracted; no finite or evanescent matching follows from
the four-dimensional residue identities.

## Two routed integrals and their exact Ward contraction

The fermion attachment is

    I_L^nu = integral [2t^2+(q-h-k)^2]
               (slash q-slash k) gamma^nu slash q
             / [q^2 (q-k)^2 ((q-p)^2-m2)].

The scalar attachment is

    I_H^nu = integral [2t^2+(q-h)^2]
               slash q (2p-2q-k)^nu
             / [q^2 ((q-p)^2-m2) ((q-t)^2-m2)].

For the unweighted zero-gauge bubble define

    T(p,h) = integral slash q [2p^2+(q-h)^2]
                         / [q^2 ((q-p)^2-m2)],
    B1(p)  = integral slash q / [q^2 ((q-p)^2-m2)].

Let delta=2(t^2-p^2). Numerator cancellation and a loop shift give the
exact regulated-integral identities

    k_nu I_L^nu = T(t,h)-T(p,h+k)-delta B1(p),
    k_nu I_H^nu = T(t,h)-T(p,h)-delta B1(p).

The fermion identity follows from
(slash q-slash k) slash k slash q =
q^2(slash q-slash k)-(q-k)^2 slash q.
The scalar identity follows from
k.(2p-2q-k)=((q-p)^2-m2)-((q-t)^2-m2).
The common T(t,h) and delta B1 pieces cancel at Q_L=-Q_H. Thus the
charge-weighted sum is Q_H times the difference I_H-I_L, whose Ward
contraction is T(p,h+k)-T(p,h). This cancellation is an integrand-level
control of the declared relative signs, not an independently assumed
physical Ward identity for every EFT counterterm.

## Residues from denominator cancellation and rank-two/three triangles

Write R_L,H for the stripped residues. The zero-gauge residue is

    T_UV(p,h) = [p^2+m2-(2/3)p.h+h^2/2] slash p
                           +[p^2/6-m2/2] slash h.

The rank-two triangle residue is g^(ab)/4. For triangle shifts a,b its
rank-three residue is [g^(ab)(a+b)^c + permutations]/12. To avoid using
rank-four trace coefficients in a dimensionally ambiguous way, cancel
the q^2 term of each seed numerator before reduction. The remaining
triangle terms have rank at most three. The needed bubble tensors are

    B1(t;0,m2) = t/2,
    B2^(ab)(t;0,m2) = t^a t^b/3 + (m2/4-t^2/12)g^(ab).

For the fermion q^2 term, shift q=l+k; its residue is

    slash t gamma^nu slash t/3
      +(t^2/6-m2/2)gamma^nu+slash t gamma^nu slash k/2.

For the scalar q^2 term, shift q=l+p in the equal-mass bubble; its
residue is -k^nu slash k/6+(k^2/6-m2)gamma^nu. Keeping the remaining
triangle terms gives

    R_L^nu = slash t gamma^nu slash t/3
      +slash t gamma^nu slash k/2
      -[slash v gamma^nu slash(p+k)
             +slash(p+k) gamma^nu slash v]/6
      +slash k gamma^nu slash v/2
      +[-5t^2/6-v^2/2+v.(p+k)/3-m2/2]gamma^nu,

    R_H^nu = -k^nu slash k/6-u^nu slash h/6+h^nu slash u/3
      +[k^2/6-m2+h.u/3-h^2/2-t^2]gamma^nu.

An independent derivation shifts the full uncontracted triangle and
integrates the Feynman simplex moments. It agrees with both expressions,
including their mass terms. The two methods do not infer rank-four
finite terms from a pole tensor. A separate general tensor-integral
engine supplies a further computational check.

## Transverse polynomial beyond the old representative

Normalize the combined reduced residue as R^nu=R_H^nu-R_L^nu. Define

    Gamma(a,nu,b) = slash a gamma^nu slash b
                   -a^nu slash b-b^nu slash a+(a.b)gamma^nu.

This is the completely antisymmetric three-gamma contraction in the
displayed order. With r=p-h-k, direct Clifford reduction gives

    R^nu = G_J^nu + D^nu,

    G_J^nu = [-(2r+k)^nu/3+(2h+k)^nu/6] slash p
                           +(p^2/6-m2/2)gamma^nu,

    D^nu = [r^nu slash k-(r.k)gamma^nu+Gamma(k,nu,r)]/2
                          +[k^nu slash k-k^2 gamma^nu]/3.

Equivalently, D^nu=[slash k,gamma^nu]slash r/4
+[k^nu slash k-k^2 gamma^nu]/3. The commutator form is checked separately.

G_J is precisely the one-B polynomial of the public six-term O_J
representative with the corresponding explicit charge/sign convention.
It satisfies k.G_J=T_UV(p,h+k)-T_UV(p,h), while k.D=0. The two source
charge normalizations must not be identified silently: G_J was stripped
of g' y_L, whereas the charged integral combination was divided by
Q_H=-y_L. The common loop/counterterm sign is fixed only after the same
gauge-vertex and seed-action conventions are imposed on both sides.
Here their *polynomials* are compared in the common normalization
defined by their equal Ward contractions.

All m2 dependence lies in G_J. There is no extra mass-dependent
transverse polynomial in this particular pair of integrals. The
transverse remainder has two derivatives, as required for a
dimension-seven four-field vertex, and need not vanish.

For an on-shell outgoing massless lepton, ubar(r) slash r=0. For a real
transverse gauge boson, k^2=0 and k.epsilon=0. Then

    ubar(r) D.epsilon
       = ubar(r)[(r.epsilon)slash k-(r.k)slash epsilon].

For example, p=(2,0,0,0), r=(1,0,0,1),
h=k=(1/2,0,0,-1/2), epsilon=(0,1,0,0) obey conservation and all light
on-shell conditions. This matrix element is nonzero, so the new
transverse residue cannot be discarded solely by a longitudinal Ward
test or the previous zero-gauge matching.

The optimized standard-library producer checks the independently
contracted primitive moments, both individual Ward identities including
their auxiliary bubbles, the combined Ward identity, the full local
decomposition, the compact commutator, the scalar-mass cancellation,
and a spin sum of 4 in the example above. Deliberate omission of the
entire D polynomial still passes Ward but fails full matching. Dropping
either attachment or assigning equal signed charges fails Ward;
dropping either transverse substructure fails full off-shell matching.

## What this establishes and what remains

The new project result is the explicitly matched transverse polynomial
of the declared two 1PI attachments, including a check that its retained
scalar-mass dependence cancels. The zero-gauge bubble, general primitive
tensor identities, field-shift Schur invariance, and flavor-kernel lemma
were already known and are used here only as controls.

This is not a complete independent Green basis, a renormalized physical
amplitude, an anomalous dimension, or a new unification result. The
common full Feynman-rule phase, literal action coefficient, complex
flavor order and Majorana convention must be carried consistently into
the full action. Operators related to D by light-field EOM generate
additional source contacts. The derivative-of-field-strength part also
requires the sourced gauge EOM if reduced. Correlated field and
parameter shifts, direct four-field contact diagrams, other external
field sectors at their specified coupling orders, and a subtraction
scheme remain necessary before a physical running coefficient is
claimed.
