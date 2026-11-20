# One weak field and the combined sourced packet

15 October 2026. Provisional calculation prepared with AI assistance.
Continue from the [literal action](../operator_map/ACTION_SIGN.md) and
[two-B checkpoint](../two_gauge_leg/00_READ_FIRST.md).

There are two new calculations, with distinct scopes.

The one-W loop calculation derives the weak-index contractions from the
seed action, with `T^I=tau^I/2` and `T_chi^I=-T^{I T}`. The direct seed
pairing vanishes by `tr T^I=0`. Its hypercharge counterpart has equal
fermion and scalar UV residues that cancel in their opposite-charge sum.
The actual crossed-triangle integrations therefore give

\[
i\mathcal M_{W^I}^{\rm UV}
=-\frac{i g F}{16\pi^2\epsilon}\,T^I_{ac}
\bar u_L(G_J^\nu+D^\nu)P_Lu_N\,\varepsilon_\nu^*,
\qquad F=C_{7,\mathrm{seed}}Y_\nu^*.
\]

The same Lorentz
polynomial was obtained for one B, with `g T^I` replaced by `g' y_L`.
Combining B and W3 in the constant neutral-Higgs projection gives
`e Q_phi`, with `Q_phi=y_L 1+T3=diag(0,-1)`. The selected neutral-photon
**simple UV pole vanishes**. Z and charged weak coefficients can remain
nonzero. This does not set the finite photon amplitude to zero, compute
a broken-phase decay rate, or establish a physical dipole beta function.

The sourced calculation instead combines the already matched hypercharge
operators `P_F=F(O_J+O_D)`. It supplies a fully expanded source-retaining
representative and an exact off-shell equation-residual certificate.
In that representative the explicit mass/source tensor term is

\[
-\frac{g'y_L}{12}\,[\bar L\widetilde H\sigma^{\mu\nu}
(M N_R+S_R^N)]_F B_{\mu\nu},
\]

before the common counterterm loop factor. The coefficient is `1/6-1/4`:
the curvature term in `O_J` must be combined with the isolated transverse
packet. All derivative source terms and current contacts are retained.
Another equivalent representative displays a different explicit dipole
coefficient, so `-1/12` is not a basis-invariant observable. This reduction
does not perform another loop integration and is restricted to U(1).

Read the evidence:

- [One-W derivation](ONE_W_DERIVATION.md) and
  [independent angular audit](INDEPENDENT_ELECTROWEAK_AUDIT.md).
- [Selected source packet and residual certificate](SELECTED_PACKET_REDUCTION.md).
- [Primary-source and recent-literature audit](ELECTROWEAK_SOURCE_CHECK.md):
  generator normalization, sourced weak currents, a newer field-redefinition
  benchmark and bounded 2025–2026 phenomenology leads.

From this directory, the weak checkers require only the standard library:

```sh
python3 -I -S -B -O verify_one_w.py
python3 -I -S -B -O verify_su2_independent.py
```

The source checker uses the pinned root SymPy environment:

```sh
python3 -I -B -O verify_selected_packet.py
```

The checkers import no repository modules and use explicit failures that
remain active under optimization. They report 695 weak checks, 1461
independent angular/index checks, and 142 source identities with 20 rejected
mutations. The independent angular route exhausts the declared quadratic
momentum coefficients and scalar-mass coefficient; these are evaluations
of a UV polynomial defined at generic off-shell kinematics, not scaleless
original-integral prescriptions. JSON files preserve the focused results.

The remaining task is an action-level field redefinition and complete
physical basis, canonical normalization, correlated counterterms and a
subtraction prescription. No two-W or mixed BW loop matching is supplied
here. Finite matching, a physical RGE and an observable remain open. Four
new bounded literature searches do not establish novelty or exhaustive
coverage. The baseline DOI does not archive this provisional continuation.
