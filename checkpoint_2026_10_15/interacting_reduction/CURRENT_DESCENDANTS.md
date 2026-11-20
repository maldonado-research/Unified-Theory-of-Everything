# Weak-current descendants and neutral/charged field coordinates

15 October 2026. These are exact weak-index and action-coordinate consequences
of the selected source reduction, not additional loop integrations. Inputs are
candidate `d822e9f4d70ad5b455436f40abc3c2b45909bbef` and the
[matrix source reduction](ELECTROWEAK_REDUCTION.md). Write
`F=C_seed Y_nu*`, `phi=tilde H`, `T^I=tau^I/2`, and `y_L=-1/2`.
The Lorentz index is kept on every current; the following products contract it.

## Ordered current contractions

Define the Grassmann-even row
`b_a^mu=sum_(alpha,i) F_(alpha i) bar L_alpha^a gamma^mu N_Li`,
`j^mu=b^mu phi`, and `j_I^mu=b^mu T^I phi`. With the declared Dplus
connection the reference-action currents are defined by

    e_B^mu = partial_nu B^(nu mu) - g' J_Y^mu,
    e_W^(I mu) = (D_ad)_nu W^(I nu mu) - g J_I^mu.

Thus the gauge-current part of this representative is

    P_current = -(g'^2 y_L j_mu J_Y^mu + g^2 j_I,mu J_I^mu)/6.

Every retained charged field contributes to these currents. For a retained
fundamental fermion `f`, `J_I^mu(f)=bar f_c gamma^mu T^I_cd f_d`. The
fundamental completeness identity gives the explicit ordered descendant

    sum_I j_I,mu J_I^mu(f)
      = (1/2) sum_(a,b) b_a,mu phi_b (bar f_b gamma^mu f_a)
        -(1/4) (b_mu phi) sum_b bar f_b gamma^mu f_b.          (C1)

No fermion is moved past another fermion here. This is a weak-index
contraction, not a spinor Fierz transformation. It does not apply unchanged
to higher weak representations. Weak singlets contribute to `J_Y` but not
to `J_I`.

The scalar's weak current, expressed using `phi`, is

    J_I^mu(phi)=i phi^dagger T^I D^mu phi
                    -i (D^mu phi)^dagger T^I phi.

For `rho=phi^dagger phi`, the same completeness tensor gives

    sum_I j_I,mu J_I^mu(phi)
       = (i/2) rho b_mu D^mu phi -(i/4) j_mu partial^mu rho.  (C2)

The norm derivative must be retained for a varying scalar background.
It is an ordinary derivative because `rho` is neutral. The hypercharge
scalar current is `i y_L[phi^dagger D phi-(D phi)^dagger phi]`, and supplies
its own term in `J_Y`. Equation (C2) keeps the scalar connection already
inside `D phi`; deleting it changes the action.

## Constant neutral scalar projection

Choose `phi=(u,0)^T`, with real constant `u`, solely to display weak indices.
No scalar equation or stationary-background condition is imposed. Define

    j_0^mu=u F_(alpha i) bar nu_Lalpha gamma^mu N_Li,
    j_c^mu=u F_(alpha i) bar e_Lalpha gamma^mu N_Li.

Then `j=j_0`, `j_1=j_c/2`, `j_2=i j_c/2`, `j_3=j_0/2`.
The source-reducing real field substitutions, with all right-hand fields
new and `kappa` first order, give

    delta B   = -kappa g' (j_0+j_0^dagger)/12,
    delta W^3 =  kappa g  (j_0+j_0^dagger)/12,
    delta W^+ =  kappa g j_c/(6 sqrt(2)),
    delta W^- =  kappa g j_c^dagger/(6 sqrt(2)).             (C3)

Here `W^+=(W^1-i W^2)/sqrt(2)` and `W^-=(W^1+i W^2)/sqrt(2)`.
Let `g_Z=sqrt(g^2+g'^2)`, `s_w=g'/g_Z`, `c_w=g/g_Z`,
`A=c_w B+s_w W^3`, and `Z=-s_w B+c_w W^3`. It follows that

    delta A=0,
    delta Z=kappa g_Z (j_0+j_0^dagger)/12.                  (C4)

These equations refer to this field-coordinate direction on this scalar
background. They do not set the photon field strength, higher-point photon
contacts, or physical photon amplitudes to zero.

For the full retained currents define the coupling-weighted neutral source
`mathcal J_Z=c_w g J_3-s_w g' J_Y`. The same projection of `P_current` is

    P_current = -g_Z j_0,mu mathcal J_Z^mu/12
                -g^2 j_c,mu (J_1^mu+i J_2^mu)/12.           (C5)

The complete counterterm includes its Hermitian conjugate. The charged
source in (C5) is lost by a free weak-field EOM substitution. The currents
include their scalar/gauge-background dependence; (C5) does not supply a
four-fermion-only theory or an observable coupling shift.

## Exact controls

Run `python -I -B -O verify_current_projection.py` in this directory with
the repository's pinned SymPy environment. The standalone checker proves
all 16 components of the Pauli completeness tensor, checks (C1) and (C2)
at eight generic complex weak-index choices, verifies finite non-diagonal
weak rotations, and checks the exact symbolic neutral/charged rotations.
It passes **61 positive controls** and rejects **20 mutations**. No project
code is imported, and Python assertions are not required. This tests the
stated algebra; completeness of the physical coefficient basis and running
remains a separate task.
