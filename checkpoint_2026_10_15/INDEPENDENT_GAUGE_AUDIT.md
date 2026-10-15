# Independent audit of the one-hypercharge-leg reduced UV tensor

15 October 2026. This is a separately authored AI-assisted derivation, not external human peer review. The calculation below was performed by homogeneous large-loop-momentum expansion independently of the producer's parameter/denominator reduction. It accepts the declared two ordered integrals as its starting point and establishes their reduced UV residue. It does not establish a physical anomalous dimension, a complete operator basis, or an unqualified action-level sign.

## Verdict and conventions

The proposed reduced sum `R = R_H - R_L` equals `G_J + D`, with

    G_J^nu = (-2 p^nu/3 + h^nu + k^nu/2) slash(p)
              + (p²/6 - m_H²/2) gamma^nu,
    D^nu = [slash(k), gamma^nu] slash(r)/4
              + (k^nu slash(k) - k² gamma^nu)/3,
    r = p - h - k.

The full transverse remainder is nonzero and independent of the signed scalar mass parameter. Both terms in D separately satisfy `k_nu D^nu = 0`. The coefficients are meaningful for the specified integrals and representative. They must not be identified with a physical beta coefficient or a convention-independent Wilson-coefficient sign.

Use metric (+---), d=4-2epsilon and C_UV=i/(16 pi² epsilon). All displayed integrals have this common UV pole stripped. Generic nonexceptional off-shell momenta define the integrals before any on-shell projection. Angular UV extraction below is understood with an auxiliary IR separation; an unregulated scaleless expansion must not be assigned zero and mistaken for a vanishing UV residue.

## Routing, ordered numerators, and Ward identity

Write `t=p-k`, and

    D_p = (q-p)² - m_H²,    D_t = (q-t)² - m_H².

The fermion and scalar integrands are respectively

    F^nu = [2t²+(q-h-k)²] (slash(q)-slash(k)) gamma^nu slash(q)
             / [q² (q-k)² D_p],
    S^nu = [2t²+(q-h)²] slash(q) (2p-2q-k)^nu / [q² D_p D_t].

Following the open chain from the external spinor row at the seed toward the incoming sterile spinor puts the seed-side momentum q-k on the left, the gauge matrix next, and the Yukawa-side momentum q on the right. The crossed composite momentum and the direct composite momentum are q-h-k and t for F, and q-h and t for S. The relative SU(2) weights are the original 1 and 2 because hypercharge is proportional to the identity in doublet space. No extra identical-line symmetry factor appears.

Let T(p,h) be the previously computed bubble polynomial and I(p)=slash(p)/2 its rank-one residue. Direct contraction and a legitimate dimensional-regularization shift give

    k.F = T(t,h) - T(p,h+k) + 2(p²-t²) I(p),
    k.S = T(t,h) - T(p,h)   + 2(p²-t²) I(p).

For example, the first identity follows from

    (slash(q)-slash(k)) slash(k) slash(q)
       = q² (slash(q)-slash(k)) - (q-k)² slash(q).

The scalar identity uses `k.(2p-2q-k)=D_p-D_t`. Therefore `k.(S-F)=T(p,h+k)-T(p,h)`. With y_H=-y_L, the declared charge-weighted algebraic combination y_L F+y_H S is y_H(S-F).

There is a critical limit to the Ward test: reversing the fermion-chain order to `slash(q) gamma^nu (slash(q)-slash(k))` preserves this contracted identity but changes the transverse residue. Thus Ward agreement alone cannot establish the ordering. The checker includes this exact negative control.

## Independent UV extraction

For either triangle, expand

    1/[q² ((q-a)²-ma²) ((q-b)²-mb²)]
      = (q²)^(-3) {1 + 2q.(a+b)/q²
          - (a²+b²-ma²-mb²)/q²
          + 4[(q.a)²+(q.b)²+(q.a)(q.b)]/(q²)² + ...}.

If N_j is the numerator's homogeneous piece of degree j in q, the logarithmic terms are N_2 times the leading denominator, N_3 times its linear correction, and N_4 times its quadratic correction. Lower numerator degrees and higher denominator corrections do not contribute to this pole. The d=4 angular contractions at the simple-pole coefficient are

    <q_mu q_nu>/q² = g_mu_nu/4,
    <q_mu q_nu q_rho q_sigma>/(q²)²
       = (g_mu_nu g_rho_sigma + g_mu_rho g_nu_sigma
          + g_mu_sigma g_nu_rho)/24.

Cancelling explicit q² factors reduces the analytic calculation to these two identities. The independent script instead retains all polynomial factors and also uses the rank-six angular coefficient 1/(4*6*8), checking its traces exactly. Setting smooth rational functions of d to d=4 is valid for this simple UV residue; their epsilon corrections can affect finite terms. No dual-epsilon or gamma5 replacement of the open triple-gamma chain is used in the analytic tensor.

Define s=h+k and v=p+k. The fermion numerator pieces are

    N4 = q² slash(q) gamma^nu slash(q),
    N3 = -q² slash(k) gamma^nu slash(q)
           -2(q.s) slash(q) gamma^nu slash(q),
    N2 = 2(q.s) slash(k) gamma^nu slash(q)
           +(s²+2t²) slash(q) gamma^nu slash(q).

For the scalar numerator, put b=2p-k:

    N4 = -2q² q^nu slash(q),
    N3 = q² b^nu slash(q) + 4(q.h) q^nu slash(q),
    N2 = -2(h²+2t²) q^nu slash(q) -2(q.h)b^nu slash(q).

The resulting separate residues are

    R_L^nu = slash(k) gamma^nu (slash(h)-slash(p))/2 + C_L gamma^nu
      + [2p^nu slash(p)-h^nu(slash(p)+slash(k))
          -(p^nu+k^nu)slash(h)]/3,

    C_L = -h²/2 -(h.k)/3 +2(p.h)/3 -7p²/6 -k² +2(p.k) -m_H²/2,

    R_H^nu = C_H gamma^nu
      +h^nu(2slash(p)-slash(k))/3
      -(2p^nu-k^nu)slash(h)/6 -k^nu slash(k)/6,

    C_H = -h²/2 +2(p.h)/3 -(h.k)/3 -p² +2(p.k) -5k²/6 -m_H².

Subtracting produces the compact independent expression

    R_H^nu-R_L^nu = [(p²+k²)/6-m_H²/2] gamma^nu
       +h^nu slash(p)-2p^nu slash(p)/3
       +k^nu slash(h)/2-k^nu slash(k)/6
       +slash(k) gamma^nu (slash(p)-slash(h))/2.

Subtracting the explicitly gauged G_J gives exactly D above. This agrees with the producer's separately obtained expression.

## Representative, action conventions, and scope

The public representative uses D_L=partial+i g' y_L B and the row derivative D barL=partial barL-i g' y_L barL B, with y_tildeH=y_L. Expanding its six fixed derivative placements gives `g'y_L G_J`. Its contraction is the same bubble difference as R. This comparison deliberately normalizes each reduced tensor to that Ward difference: the triangle charge factor y_H and the representative factor y_L have opposite signs. Full propagator, vertex, Fourier and Majorana conventions must be included before comparing unstripped amplitudes or assigning a counterterm sign. The reviewed draft states this limitation explicitly.

For a literal action +(1/2) J^T W1 Box J with exp(-ipx), the matching convention using positive C_(7,seed) q² requires C_(7,seed)=-W1. The singlet J=epsilon LH is neutral; partial derivatives on the whole J introduce no independent seed gauge vertex. With exactly the specified seed, Yukawa and single hypercharge vertex, the two internal attachments are the 1PI triangles in this external-field sector. External-leg emissions are reducible and matter for a complete amplitude, not for this selected 1PI residue. Broader coupling sectors, sourced EOM contacts and physical matching remain open.

The chirality projector can stay at the right of the odd open gamma chain. There is no closed anomalous trace. Explicit four-dimensional matrices verify the pole's retained tensors and the witness below, but do not determine evanescent finite terms or a complete subtraction scheme.

## Nonzero on-shell witness and reproducibility

For a massless external lepton and real transverse gauge field, the `slash(r)` in the commutator expression cannot simply be dropped: it is to the right of the gamma matrices, not adjacent to the external spinor row. Correct anticommutation gives

    ubar(r) D(epsilon)
       = ubar(r)[(r.epsilon) slash(k)-(r.k) slash(epsilon)].

Take p=(2,0,0,0), r=(1,0,0,1), h=k=(1/2,0,0,-1/2) and epsilon=(0,1,0,0). Momentum is conserved; r,h,k are null and k.epsilon=0. With sterile mass 2, standard spin sums and the right P_L projector give the squared stripped D matrix element 4, or 2 after averaging over the two initial sterile spins. This witnesses a nonzero tensor matrix element, not a separately observable contribution from the chosen subtraction.

Run `python3 -I -S -B -O verify_gauge_asymptotic.py`. The standard-library script imports no producer code, reconstructs the literal numerators with exact Gaussian-rational matrices, and extracts their UV poles by homogeneous angular expansion. It checks 24 generic off-shell rational kinematic choices (1,536 matrix entries), the Ward identity, exact transversality, Clifford and angular-trace controls, omitted-attachment and omitted-transverse mutations, a Ward-blind chain-reversal control, and the nonzero spin-summed witness. Explicit failure exceptions remain active under optimized Python. This is reproducible internal algebra verification, not proof-assistant formalization or independent human validation.
