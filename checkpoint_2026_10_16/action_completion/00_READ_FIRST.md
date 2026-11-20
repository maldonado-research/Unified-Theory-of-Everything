# Selected action-map composition and scalar-current completion

16 October 2026. Provisional continuation of candidate
`5f9889b76661faa155c681a5ecfac5d036398bc8`, prepared with AI assistance.
The inputs are the previously matched electroweak local packet, its sourced
field changes, and its separately specified local-normalizing coordinates.
This checkpoint performs further classical action algebra. It integrates
no additional loop and does not complete physical matching or running.

## Two coordinate maps, explicitly related

The [field-map composition](FIELD_MAP_COMPOSITION.md) relates the original
source-reducing and local-normalizing substitutions for every declared
field. If `Xi_old=Xi_S+kappa v_S=Xi_C+kappa v_C`, the bridge is
`Xi_S=Xi_C+kappa(v_C−v_S)` through first order. Applying the bridge to the
full reference action transports its sources along with its fields.
The induced coefficient tables of the two original maps cannot be added.

The explicit neutral quadratic calculation compares the actual six-placement
packet, its sourced representative, and direct substitutions in the
kinetic, Yukawa, Majorana and scalar-potential action. The nonlinear scalar
shift contributes even when only a fermion Hessian is being computed:

```text
phi=(u,0), t_1=F_alpha i bar nu_Lalpha A_N N_Li/6,
delta_S L_potential^(2)
  =u(m_H^2+2lambda u^2)(t_1+t_1^dagger).
```

This term is generally nonzero off shell and is canceled by the scalar
part of the coordinate bridge. A fermion-Jacobian-only comparison misses
it. The [independent composition audit](INDEPENDENT_COMPOSITION_AUDIT.md)
checks the actual actions with exterior Grassmann Fourier modes, including
noncommuting flavor matrices and non-diagonal real retained masses. At a
stationary scalar background the tadpole term vanishes. Neither that
specialization nor the classical chain rule determines a transformed
quantum mean-field effective action.

## Scalar currents with their sources retained

The [scalar-current reduction](SCALAR_CURRENT_REDUCTION.md) resolves the
remaining explicit scalar gauge-current part of the selected representative.
Write `rho=phi^dagger phi`,
`b_mu^a=sum F_alpha i bar L_alpha^a gamma_mu N_Li`, and
`j_mu=b_mu phi`. Two derivative structures can be retained after an IBP
and full fermion-EOM reduction. In the declared convention their
derivative contributions are

```text
−i(3g^2+g'^2)rho b_mu D^mu phi/24
−i g'^2 j_mu(phi^dagger D^mu phi)/12.
```

The associated source term is essential:

```text
+(g^2+g'^2)rho F_alpha i [
   bar S_Lalpha phi N_Li−bar L_alpha phi(M N_R+S_R)_i]/24.
```

For the literal neutrino Yukawa action it induces active scalar,
retained-mass and additional-sterile scalar monomials. An exact off-shell
residual identifies the discarded EOM terms; interaction sources are not
set to zero. The [independent scalar audit](INDEPENDENT_SCALAR_AUDIT.md)
checks the weak contractions, complete-current divergence and source signs.
These are contributions from this explicit subpacket. Other derivative,
source and coordinate terms can enter the same final operators.

## Charged-lepton and quark spectators

The [spectator Yukawa transport](SPECTATOR_YUKAWA_TRANSPORT.md) substitutes
the actual lepton and scalar shifts in the retained charged-lepton and
quark interactions. The source-reducing map induces derivative contacts,
including a nonzero charged-lepton term proportional to
`(F^dagger Y_e)[(D_mu phi)^dagger H]`. The specified local-normalizing map
gives zero spectator-Yukawa variation by `phi^dagger H=0` and its unchanged
scalar/quark fields. The bridge cancels the source-map spectator terms.
This calculation prevents silent deletion of the sources when changing
coordinates; it does not generate separate loop predictions for quarks.

## Reproduction, literature and remaining scope

Run the five standalone scripts in the root pinned Python environment:

```sh
python -I -B -O verify_field_map_composition.py
python -I -B -O independent_composition_audit.py
python -I -B -O verify_scalar_current_reduction.py
python -I -B -O independent_scalar_audit.py
python -I -B -O verify_spectator_yukawa.py
```

Every script uses explicit failures, imports no project calculations, and
keeps its checks active under optimization. Adjacent receipts state the
actual controls and their scope. The root runner executes them on a
disposable copy, preserving historical scientific files and receipts.
These are internal AI-assisted calculations and implementation audits,
not independent human peer review or empirical validation.

The [literature comparison](LITERATURE_ACTION_COMPLETION.md) records four
focused searches, 32 distinct paper records, and primary-body comparisons.
A 2025 background-source benchmark supports retaining gradient-dependent
quadratic terms. Another distinguishes a classical field substitution
from the nonlinear transformation of a quantum effective action. Paper
coefficients are not imported into the selected packet.

The selected all-field bridge is now explicit, and its full neutral
quadratic action is checked off shell. This does not canonically normalize
all interacting derivative terms. The remaining scalar/derivative and
additional-sterile sectors, total coefficients in an independent
derivative/EOM basis, lower-order insertions, other loop/counterterm sectors,
dimensional/evanescent definitions, gauge fixing/ghost and composite-source
quantum transport, finite matching and one physical observable remain
unresolved. No completed unification, new fundamental mathematics,
physical dipole prediction, empirical validation or external novelty is
established. The checkpoint remains a GitHub review candidate; Zenodo
continues to identify its unchanged earlier baseline.
