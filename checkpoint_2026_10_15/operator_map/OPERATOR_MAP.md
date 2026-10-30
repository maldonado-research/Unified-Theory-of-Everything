# Local transverse operators and their sourced descendants

15 October 2026. Restricted four-dimensional simple-pole algebra; this is
not a complete dimension-seven EFT basis or a physical running coefficient.

## Conventions and the exact restricted space

Write phi=tilde H, eta=diag(+---), D=partial+i g' y B, y_L=y_phi=-1/2,
D bar L=partial bar L-i g'y_L bar L B, and g_y=g'y_L. Outgoing bar L,
phi and B have waves exp(+ir.x), exp(+ih.x), exp(+ik.x); incoming N has
exp(-ip.x), p=r+h+k. A=i slash partial, sigma^{mu nu}=i[gamma^mu,gamma^nu]/2.
The action vertex is i times the Fourier kernel of the monomial. Flavor
indices on each operator are bar L_alpha ... N_beta; the Hermitian
conjugate channel is also required.

At fixed chirality, consider local Lorentz-covariant one-B kernels that
are homogeneous of degree two in the independent off-shell momenta r,h,k,
have one free vector index, and use an odd Clifford chain followed by P_L.
There are no other external vectors. Define

    V_a^nu = a^nu slash k - (a.k) gamma^nu,    a=r,h,k,
    E_a^nu = slash k gamma^nu slash a - k^nu slash a
             -a^nu slash k + (k.a) gamma^nu, a=r,h.

The transverse space is exactly five-dimensional, with independent basis
(V_r,V_h,V_k,E_r,E_h). This is an off-shell tensor-space statement in this
declared sector, not a count of the full EFT operator basis.

For completeness: grade one starts with 6 gamma^nu(a.b) and 9 a^nu slash b
structures, a,b in {r,h,k}. The Ward map has rank 12, leaving the three
V_a. Grade three has E(r,nu,h), E(k,nu,r), E(k,nu,h); its first tensor has
a nonzero generic Ward contraction, leaving E_r,E_h. At fixed chirality
grade three is a dual-vector structure, and these five remain independent.
The attached exact check proves the rank-12 statement by polynomial
coefficient elimination, checks the five projected tensors have rank five
on exact off-shell data, and checks their Ward identities. The grade-one
and grade-three Ward images cannot cancel: their vector coefficients lie,
respectively, in the span of r,h,k and its dual orthogonal direction.

A term proportional to a scalar mass squared can only be m^2 gamma^nu
at this mass dimension. It is not transverse for generic k, so no local
mass-squared transverse direction exists in this restricted space.

## Five explicitly placed gauge-invariant operators

Let G^{mu nu rho}=gamma^{[mu}gamma^nu gamma^{rho]} with unit-weight
antisymmetrization. Define the following monomials, including their
normalization by g_y:

    Q_Vr = -g_y (D^nu bar L) phi gamma^mu N_L B_{mu nu},
    Q_Vh = -g_y bar L (D^nu phi) gamma^mu N_L B_{mu nu},
    Q_Vk = -g_y bar L phi gamma^mu N_L partial^nu B_{mu nu},
    Q_Er = -g_y/2 (D_rho bar L) phi G^{mu nu rho} N_L B_{mu nu},
    Q_Eh = -g_y/2 bar L (D_rho phi) G^{mu nu rho} N_L B_{mu nu}.

The electroweak doublet indices contract between bar L and phi. All five
have mass dimension seven and vanish at B=0. Their one-B Fourier kernels
are exactly g_y(V_r,V_h,V_k,E_r,E_h)P_L. Thus an action term c Q_Vr has
vertex i c g_y V_r P_L, and similarly for the other four. The derivative
connection terms supply the correlated two-gauge vertices; they must be
kept in the gauge-invariant operators even though this calculation matches
only one B. Operators containing nonabelian covariant derivatives similarly
carry their required mixed gauge vertices.

With gamma5=i gamma0 gamma1 gamma2 gamma3, epsilon^{0123}=+1 and
tilde B^{mu nu}=epsilon^{mu nu rho sigma}B_{rho sigma}/2, one has
G^{mu nu rho}=i epsilon^{mu nu rho sigma}gamma_sigma gamma5. On N_L,
Q_Eh=i Q_tildeVh, where Q_tildeVh is Q_Vh with B replaced by tilde B.
The two Higgs-derivative ordinary/dual operators are independent; their
sum is the particular combination selected below.

Two useful derivative-on-N definitions are

    Q_Vp = +g_y bar L phi gamma^mu partial^nu N_L B_{mu nu},
    Q_Ep = +g_y/2 bar L phi G^{mu nu rho} partial_rho N_L B_{mu nu}.

Integration by parts (IBP) gives

    Q_Vp = Q_Vr+Q_Vh+Q_Vk,       Q_Ep = Q_Er+Q_Eh,       modulo a total derivative.

The second relation uses the abelian Bianchi identity
partial_[rho B_mu nu]=0. These relations exhibit why six arbitrary
field-derivative ordinary/dual monomials do not make six independent
directions once divergence and IBP are treated consistently.

## Coefficient map of the computed reduced tensor

The supplied reduced residue is

    D^nu = [slash k,gamma^nu]slash r/4 + V_k^nu/3
         = (V_r^nu+E_r^nu)/2 + V_k^nu/3.

Consequently the unique coefficient vector in the five-dimensional
declared tensor basis is (1/2,0,1/3,1/2,0), represented by

    O_D = (Q_Vr+Q_Er)/2 + Q_Vk/3
        = i g_y/4 (D_rho bar L) phi sigma^{mu nu}gamma^rho N_L B_{mu nu}
          -g_y/3 bar L phi gamma^mu N_L partial^nu B_{mu nu}.

Its one-B action vertex is i g_y D^nu P_L. This is a reduced-normalization
map, independent of a loop sign convention. If the correlated action
audit gives iM_loop=-i g_y F D/(16 pi^2 epsilon), the canceling action is
+F O_D/(16 pi^2 epsilon)+h.c. If its overall pole sign is reversed, this
counterterm sign reverses. The map itself does not infer that sign.

## Complete light-lepton and hypercharge source reduction

Specify sources from the same reference action, with canonical kinetic
terms -B^2/4, bar f i slash D f and (D H)^dagger(D H):

    i slash D L_alpha = (S_L)_alpha,
    -i(D_mu bar L_alpha)gamma^mu = (bar S_L)_alpha,
    partial_mu B^{mu nu} = g' J_Y^nu.

At renormalizable order,

    (S_L)_alpha = (Y_e)_{alpha a} H e_Ra
                 +(Y_nu)_{alpha j} phi N_Rj,
    J_Y^mu = (i/2) H^dagger <->D^mu H
             + sum_generations [ -1/2 bar L gamma^mu L
               -bar e_R gamma^mu e_R +1/6 bar Q gamma^mu Q
               +2/3 bar u_R gamma^mu u_R -1/3 bar d_R gamma^mu d_R ].

Color and weak components are summed. A sterile neutrino has zero
hypercharge. The covariant scalar current includes its gauge connections.
Additional light fields contribute with their own hypercharges if present.
In a larger reference EFT these symbols denote TOTAL sources, including
functional derivatives of every retained lower operator. In particular
J_Y receives -(1/g') delta S_lower/delta B_mu, with the variational
derivative including derivatives on B. One must not import a current
formula written with D=partial-i g'YB without translating its sign and
the order of the two field-strength indices.

The Clifford identity

    [gamma^mu,gamma^nu]gamma^rho
      = gamma^rho[gamma^mu,gamma^nu]
        +4(eta^{nu rho}gamma^mu-eta^{mu rho}gamma^nu)

gives the exact sourced relation

    Q_Er ~= Q_Vr -g_y/2 (bar S_L) phi sigma^{mu nu}N_L B_{mu nu},

and hence

    O_D ~= Q_Vr -g_y/4 (bar S_L) phi sigma^{mu nu}N_L B_{mu nu}
            +g'^2 y_L/3 (bar L phi gamma^mu N_L) J_{Y mu}.

Here ~= means EOM equivalence at first order in the inserted coefficient,
not equality of off-shell Green functions. Crucially,
partial^nu B_{mu nu}=-g'J_{Y mu}; this antisymmetry sign is responsible
for the positive current coefficient in this derivative-on-L form.
The scalar current and every fermion current survive in this equation.

For an arbitrary complex flavor matrix F multiplying O_D,

    F_{alpha beta}(bar S_L)_alpha phi sigma N_Lbeta B
       = (Y_nu^dagger F)_{j beta}(H^dagger H)
           bar N_Rj sigma N_Lbeta B + lower-source terms.

The charged-lepton Yukawa term vanishes here by H^dagger tilde H=0;
it was retained until this explicit algebraic contraction. No matrices
have been commuted. For Majorana N, the tensor bilinear is antisymmetric
in j,beta, so only the antisymmetric part of Y_nu^dagger F contributes.
For the particular F=C_seed Y_nu^* with C_seed^T=C_seed, that matrix is
symmetric and this particular renormalizable sterile-tensor contact
vanishes. It does not vanish for arbitrary F or automatically for lower
EFT sources. This simplification does not remove the gauge-current terms.

## Equivalent Higgs-derivative form and full sterile source packet

It is useful to expose a second exact IBP route, rather than confuse a
choice of source packet with absence of the other one. Define

    Q_AN = g_y bar L phi sigma^{mu nu} A N_L B_{mu nu}.

Direct IBP on the compact O_D, followed by the Bianchi identity, gives

    O_D = -(Q_Vh+Q_Eh)/2-Q_Vk/6-Q_AN/4, modulo a total derivative.

Equivalently, use Q_Vp+Q_Ep=-Q_AN/2 together with the preceding lepton
source identity and IBP relations. The explicit bar S_L terms then cancel
between those intermediate relations. Thus it would be inconsistent to
append that source again to this Higgs-derivative form.

With the sterile source equation A N_L=M N_R+S_R^N, this becomes

    O_D ~= -(Q_Vh+Q_Eh)/2
           -g'^2 y_L/6 (bar L phi gamma^mu N_L)J_{Y mu}
           -g_y/4 bar L phi sigma^{mu nu}(M N_R+S_R^N)B_{mu nu}.

The antisymmetric-three-gamma form is primary; its Clifford/IBP identities
do not require replacement by a four-dimensional dual tensor. Only for
the specified four-dimensional residue and gamma5 conventions may one
substitute Q_Eh=i Q_tildeVh. Since Q_tildeVh=-g_y Q_Hdual when
Q_Hdual=bar L(D^nu phi)gamma^mu N_L tilde B_{mu nu} is unnormalized,
the unnormalized Higgs pair is +(Q_H+i Q_Hdual)/2 in O_D/g_y.
No finite evanescent relation is inferred from that substitution.

For a general mass matrix the mass coefficient is F M in that order.
In the positive diagonal mass convention this is still a matrix product,
not a license to pull M through F. For the renormalizable Majorana
neutrino Yukawa action the paired equation has
S_L^N=Y_nu^dagger phi^dagger L and S_R^N=(S_L^N)^c; the source equations
fix phases. The last displayed term contains both the mass times the
dimension-six dipole and the two-lepton/two-scalar/field-strength contact.
Dropping S_R^N while retaining the interacting Yukawa theory is not an
EOM reduction. Lower EFT terms add their variational sources at whatever
inverse-scale order is being retained.

There is a useful further flavor check. Write X_alpha=phi^dagger L_alpha,
so S_R^N=Y_nu^T X^c in these source conventions. The last source operator
has coefficient F Y_nu^T multiplying bar X_alpha sigma X_gamma^c B.
Only its antisymmetric part survives the Grassmann/charge-conjugation
tensor identity. In the selected symmetric-seed sector, set
H_Y=Y_nu Y_nu^dagger (this H_Y is a flavor matrix, not the Higgs field).
Then

    [F Y_nu^T]_antisym = [C_seed H_Y^T-H_Y C_seed]/2.

It is generally nonzero even though the lepton-route coefficient
Y_nu^dagger C_seed Y_nu^* is symmetric. The two source packets therefore
cannot be equated by dropping both contacts. If H_Y=U diag(lambda_i)U^dagger
and C'=U^dagger C_seed U^*, the active-source coefficient transforms to
(lambda_j-lambda_i)C'_{ij}/2. It vanishes exactly when C' has support only
within equal-eigenvalue blocks. This conditional linear-algebra statement
is distinct from a full operator-basis or anomalous-dimension result.

The apparent current coefficients +1/3 and -1/6 belong to different
derivative representatives: the conversion of Q_Vr carries the rest.
They are not competing predictions. For example, in unnormalized raw
operators Q_H=bar L(D^nu phi)gamma^mu N_L B_{mu nu},
Q_N=bar L phi gamma^mu partial^nu N_L B_{mu nu},
Q_div=bar L phi gamma^mu N_L partial^nu B_{nu mu},

    O_D/g_y ~= Q_H+Q_N-2 Q_div/3-(bar S_L)phi sigma N_L B/4,

which displays the same IBP/current bookkeeping before the sterile EOM.

## Power counting and what a real B wave misses

For k^2=k.epsilon=0, V_k.epsilon=0. Thus the Q_Vk direction is invisible
to a real external transverse B wave, but its sourced reduction contains
(bar L phi gamma N) times the scalar current and every hypercharge
fermion current. These are genuine higher-multiplicity local operators,
not zero operators. The E_r-to-V_r relation additionally uses the lepton
EOM, so its discarded row-slash term represents the spinor source above.
The nonzero real-wave tensor alone does not determine any of these
off-shell/source coefficients.

If F has dimension -3 and this is a first-order 1/Lambda^3 calculation,
using renormalizable total sources is consistent at that inverse-scale
order. A source from an independent dimension-d operator gives an extra
factor 1/Lambda^(d-4), hence order 1/Lambda^(d-1) when multiplied by F;
it must be included if that mixed order is retained. This statement
assumes ordinary canonical power counting and no assigned compensating
enhancement. In sequential matching, retained lower coefficients and
the actual target order must be specified before truncating the sources.

The one-B vertex has coefficient g'F, whereas the gauge-current source
packet has coefficient g'^2 F. A projection selecting one power of g'
in four-field 1PI graphs does not prove that those higher-multiplicity
contacts vanish. Similarly, the Yukawa-source packets add another Yukawa
factor and the mass dipole adds M. These are correlated changes of
representative, not newly calculated loop amplitudes in those sectors.

## Reproduction and limits

Run `python verify_operator_map.py` or `python -O verify_operator_map.py`.
The standalone script uses only the standard library, exact Gaussian
rationals and polynomial coefficient elimination. Its checks test the
five-dimensional restricted tensor space, direct Fourier placement,
IBP/Bianchi identities, both sourced algebraic routes, current weights
and flavor order. Deliberate omissions/sign mutations must be rejected.
These checks do not calculate another loop, establish a full EFT basis,
fix a subtraction scheme, or identify a physical beta function.
