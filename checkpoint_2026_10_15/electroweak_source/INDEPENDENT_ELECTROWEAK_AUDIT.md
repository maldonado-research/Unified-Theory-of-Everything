# Independent weak-index, UV-pole and source-packet audit

15 October 2026. This audit starts from the literal action in the previous
`operator_map/ACTION_SIGN.md`. It imports no repository calculation code.
It concerns the selected coefficient linear in `C_seed Y_nu*` and one
external electroweak gauge coupling, at the simple UV pole. It does not
claim a physical anomalous dimension or a complete dimension-seven basis.

## Result

Let `F=C_seed Y_nu*`, `phi=tilde H`, `T^I=tau^I/2`, and
`D=partial+i g T^I W^I+i g' y_L B` on `L` and `phi`, with `y_L=-1/2`.
The full one-W pole in the literal action and source-chain convention is

    i M_W^UV = -i g/(16 pi^2 epsilon)
               F_(alpha i) (T^I)_(ac)
               bar u_(L alpha)^a R^nu P_L u_(N i) epsilon_nu^I*.

Here the external scalar carries the `phi^c` index, `p=r+h+k`, and

    R^nu = G_J^nu + D^nu,
    G_J^nu = [-(2r+k)^nu/3+(2h+k)^nu/6] slash p
             +(p^2/6-m_H^2/2) gamma^nu,
    D^nu = [slash k,gamma^nu] slash r/4
           +(k^nu slash k-k^2 gamma^nu)/3.

Thus the summed simple UV pole has the same momentum polynomial as the
previous hypercharge result, with `g' y_L delta_ac` replaced by `g T^I_ac`.
That replacement is false diagram by diagram: its validity for the sum
requires the trace and UV cancellation demonstrated below. The cancelling
local action has the opposite amplitude sign. The established factors
`(i seed)(i fermion)(-i Yukawa)(i scalar)=-1` and the additional gauge
insertion `i(-i)=+1` are unchanged; no new Majorana factor is inserted.

## Literal index route and conjugation

Writing the seed in `phi` components produces direct and crossed tensors
`delta_ac delta_bd` and `delta_ad delta_bc`, respectively. The conjugate
Yukawa contracts the internal charge-conjugate fermion and scalar indices.
The internal fermion is `chi=L^c` with generator `-T^{I T}`; the internal
scalar is `phi` with generator `T^I`. Explicit index summation gives

| Pairing | Scalar insertion | Fermion insertion |
|---|---|---|
| Direct | `tr(T^I) delta_ac` | `-tr(T^I) delta_ac` |
| Crossed | `T^I_ac` | `-T^I_ac` |

The seed momentum factor before the weak-index sum is `t^2`, not `2t^2`,
for the direct pairing, where `t=p-k`. The familiar direct numerator
`2t^2` already includes the zero-gauge doublet trace `tr(1)=2`. For
hypercharge, `tr(y_L 1)=2y_L` restores the old direct weight. For each
SU(2) generator the direct factor is zero. Assigning `2T^I` to that pairing
would therefore double-count a zero-gauge multiplicity as a gauge matrix.

The original epsilon-index route is checked separately. With
`E=[[0,1],[-1,0]]`, pseudoreality is

    T^I E + E T^{I T}=0.

Both `L^c` and `H*` carry the anti-doublet generator `-T^{I T}`. Contracting
the direct/crossed seed tensors `E_ac E_bd`, `E_ad E_bc` with the Yukawa
epsilon gives direct zero, crossed scalar `+(T^I E)_ac`, and crossed
fermion `-(T^I E)_ac`. Converting the external `H*` to `phi=E H*` agrees
with the table. All three generators are checked, including the imaginary
second generator, for which silently omitting the transpose is unsafe.

The Lorentz charge-conjugation identity remains
`C gamma^{nu T} C^{-1}=-gamma^nu`; transposition reverses the ordered
propagator product. The chain in this route is consequently
`(slash q-slash k) gamma^nu slash q`, with the seed-side propagator on
the left. The representation conjugation supplies `-T^T` once; adding
another charge or momentum reversal would count that conjugation twice.

## Independent integration and Ward test

Only the crossed numerator remains for one W. The two integrands are

    I_Fx^nu = integral (q-h-k)^2 (slash q-slash k) gamma^nu slash q
                        /[q^2 (q-k)^2 ((q-p)^2-m_H^2)],
    I_Sx^nu = integral (q-h)^2 slash q (2p-2q-k)^nu
                        /[q^2 ((q-p)^2-m_H^2)((q-t)^2-m_H^2)].

The physical group-weighted combination is `T^I(I_Sx-I_Fx)`. Independently
expanding their full numerators at large loop momentum gives the displayed
`R=G_J+D`. An analytic shortcut explains why it agrees with hypercharge.
For the removed direct numerator `2t^2`, the rank-two triangle residue
`q^mu q^nu -> eta^{mu nu}/4` gives

    Res I_F,direct^nu = 2t^2 gamma_mu gamma^nu gamma^mu/4
                     = -t^2 gamma^nu,
    Res I_S,direct^nu = 2t^2 (-2) gamma_mu eta^{mu nu}/4
                     = -t^2 gamma^nu.

Both individual poles are nonzero but cancel in their difference. Therefore
`Res I_Fx=R_L+t^2 gamma` and `Res I_Sx=R_H+t^2 gamma`, using the old
hypercharge notation. This cancellation is a simple-pole result; it does
not equate the two finite direct diagrams.

Define the crossed bubble `T_x(p,h)=T(p,h)-2p^2 B1(p)`. Numerator
cancellation in the unintegrated crossed triangles gives

    k.I_Fx = T_x(t,h)-T_x(p,h+k),
    k.I_Sx = T_x(t,h)-T_x(p,h).

Consequently `k.(I_Sx-I_Fx)=T_x(p,h+k)-T_x(p,h)` equals the old full
bubble difference: its direct piece is independent of `h`. The SU(2)
group route is therefore consistent with the same covariant `G_J` Ward
response, while `k.D=0`. A longitudinal Ward test by itself still cannot
determine `D` or the ordered fermion chain.

The standalone checker uses a rational Clifford algebra and a separately
implemented asymptotic angular UV functional. For three denominators
shifted by `a_i`, writing `c_i=a_i^2-m_i^2`, the terms required for a
degree-four numerator are

    1 + 2 q.sum(a_i)/q^2 - sum(c_i)/q^2
      +4 sum_(i<=j) (q.a_i)(q.a_j)/(q^2)^2.

It extracts the logarithmic coefficient with angular moments
`eta/4`, `sym(eta eta)/24`, and `sym(eta eta eta)/192`. Its starting
numerators, rather than the earlier reduced `R_L,R_H` formulas, are the
integration inputs. Four-dimensional angular and Clifford coefficients
are sufficient for this simple pole; their dimension-dependent
corrections affect finite terms and are not claimed here.

The degree-two local momentum polynomial is exhausted using the 12
diagonal and 66 pairwise-sum substitutions in `(p,h,k)` for each of the
four polarization basis vectors, plus four independent scalar-mass
tests and a generic rational off-shell case. These are evaluations of
the UV polynomial defined at generic momenta. They do not instruct one
to set exceptional original massless integrals to zero as scaleless, or
identify an IR pole with a UV pole. Signed scalar masses are retained.

## Neutral photon and Z projections

Use `A=s_W W^3+c_W B`, `Z=c_W W^3-s_W B`, and
`g s_W=g' c_W=e`. The common pole polynomial then multiplies

    photon: e (T^3+y_L 1)=e diag(0,-1),
    Z:      g c_W T^3-g' s_W y_L 1.

For the constant neutral scalar background `phi=(v/sqrt(2),0)^T` and
the matching upper lepton component, the one-photon local UV pole
vanishes. The neutral Z coefficient is instead
`(g c_W+g' s_W)/2=e/(2s_W c_W)`, which is nonzero. The lower component
has photon coefficient `-e`, providing a nonzero charged control.

This statement concerns only this matched pole and the specified constant
neutral-Higgs projection. It is not a statement that a complete broken
electroweak finite amplitude or every neutrino transition moment vanishes.
The external VEV projection uses the corresponding zero-Higgs-momentum
local polynomial; it does not supply a complete treatment of electroweak
symmetry breaking, other insertions, thresholds or source descendants.

## Review of the selected hypercharge source packet

An independent inspection of `SELECTED_PACKET_REDUCTION.md` finds the
combined explicit mass/source dipole coefficient `-g_y/12` correct in its
declared representative. The `O_L/3` lepton-curvature contribution is
`+g_y/6`, while the Higgs-derivative IBP route for the transverse packet
contributes `-g_y/4`; their sum is `-g_y/12`. Applying
`A N_L=M N_R+S_R` retains both terms with this common coefficient.

The displayed residual certificate has the consistent signs for the
scalar residual, row residual, sterile source tower, and
`-g_y/6 (bar L phi gamma_mu N_L) e_B^mu`. The scalar-mass coefficient
is `m_H^2/2+m_H^2/6=2m_H^2/3`, and the ordered mass matrix multiplying
the explicit dipole is `F M`.

This coefficient is not a physical-basis result. The same packet in its
alternative lepton-source representative displays `+g_y F M/6` as the
explicit mass dipole, with compensating derivative and current/source
terms. Derivatives of sources remain in both presentations. Neither
isolated coefficient can be read as a physical anomalous dimension or
transition moment. The source note explicitly restricts its calculation
to the U(1) connection; it must not be presented as the full sourced
SU(2) reduction supplied by this one-W calculation.

## Reproduction and remaining scope

Run `python3 -I -S -B -O verify_su2_independent.py`. The frozen result is
recorded in `SU2_INDEPENDENT_RESULTS.json`. All checks use explicit
exceptions and remain active with Python assertions disabled. There are
317 direct integrations/comparisons of the crossed pole, 317 checks of
each direct attachment and their cancellation, 169 contracted Ward
checks, and exact weak-index/projection checks with mutation controls.

The negative controls include the false `2T` direct group weight, the
wrong charge-conjugate crossed sign, removal of the transverse remainder,
equal crossed charges, and the B-only neutral-photon inference. The
direct-weight mutation is tested before summation because its spurious
direct poles would themselves cancel and could evade a summed-pole test.

Neither this audit nor the one-W result independently matches two-W or
mixed B-W diagrams, their seagulls and nonabelian contacts. A covariant
lift fixes a set of correlated vertices, but that prediction still needs
the applicable matching test. Full sourced weak-current reduction,
independent operator-basis reduction, correlated lower-dimensional
counterterms, parameter/field shifts, and a subtraction prescription
remain outside this checkpoint. No physical beta function follows.
