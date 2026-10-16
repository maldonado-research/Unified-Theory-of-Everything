# Derivative placement, complete source projection and complex retained masses

16 October 2026. Provisional continuation of candidate
`4ed8560e019dfe74ed6c5178b0f990fe3ccb3008`, prepared with AI assistance.
The input is the previously matched electroweak packet and its
source-preserving representative `P_red`. This checkpoint performs further
classical action algebra, with explicit coefficient transport and independent
component audits. It evaluates no additional loop or physical observable.

## Curvature derivatives and their correlated sources

The [curvature derivation](CURVATURE_DERIVATIVE_BASIS.md) fixes the complete
four-dimensional Lorentz/weak invariant space for the declared placement
`bar L gamma N_L (D phi) B_mu nu` or its weak-strength counterpart.
There are two complex Lorentz directions per strength, with one weak
intertwiner each. This is algebraic completeness at this derivative
placement, before the full EFT IBP/EOM quotient.

In the declared gamma/epsilon convention the actual selected combination
obeys

```text
−(Q_Vh+Q_Eh)/2=H(Fcal_+),
Fcal_+=(Fcal+i tilde Fcal)/2,   tilde Fcal_+=−i Fcal_+,
H(Fcal)=bar L Fcal_mu nu(D^nu phi)gamma^mu N_L.
```

Thus the packet occupies one specific chiral strength projection in this
placement. The Hermitian-conjugate action has the conjugate projection.
The B/W coefficients remain tied to the same source matrix and declared couplings. This is not a general selection rule for physical running.

An exact identity moves the scalar derivative to an explicitly defined
row-derivative representative while retaining its boundary term, both
fermion interaction sources and the complete gauge currents. If that
identity is applied to the displayed part of `P_red`, the explicit
curvature-source coefficient changes from `−1/12` to `+1/6`, and its
current coefficient from `−1/6` to `+1/3`. A new lepton tensor-source term
belongs to the same transport. These related coefficients describe
equivalent sourced action representatives. They are not different
predictions for a dipole or a measured coupling.

The actual symmetric seed makes the induced pure-neutrino sterile tensor
source vanish by its flavor exchange. A charged-lepton weak source can
survive: `H^dagger phi=0` does not imply `H^dagger T^I phi=0`.
The [independent curvature audit](INDEPENDENT_CURVATURE_AUDIT.md) derives
the Clifford/dual signs and the source-retaining identity separately.

## Collecting all terms with no retained sterile field

The [active source projection](ACTIVE_SOURCE_PROJECTION.md) collects the
complete no-retained-sterile-field part of the original `P_red` for the
declared neutrino-Yukawa action plus optional scalar quartic. It exhausts
the literal source terms instead of projecting each previously reported
coefficient table in isolation. With `X=phi^dagger L`, `P0=−partial^2`,
real symmetric retained `M`, and

```text
K=F Y_nu^T, K_M=F M^2 Y_nu^T, G=F M Y_nu^dagger,
c=2m_H^2/3+lambda rho/3,
```

the result is

```text
P_active=(5/6)bar X(K_M X^c+K P0 X^c+G A X)
         +c bar X K X^c
         −K_alpha beta bar L_alpha Fcal_mu nu phi
                       sigma^{mu nu}X_beta^c/12,
A=i slash partial.
```

The full action includes its Hermitian conjugate. Scalar/Box monomials
use symmetric flavor matrices; the LNC derivative term uses a Hermitian
matrix after adding its conjugate; the weak tensor retains the combined
weak/flavor exchange. Literal factors for a seed defined with `1/2` are
given in the note. They are selected coefficient projections, not
physical counterterms or beta functions inferred from this block alone.

Setting the retained sterile field to zero is a field-content operation.
It is not an EOM, a functional integration, a low-energy matching step or
an on-shell condition. The original mixed packet has zero active-only
projection, while its source-reduced coordinates have the nonzero block
above. Field transformations mix the retained sterile and light fields.
For the additional scalar-current map, applying that map before the
projection produces the active term `−d rho bar X K X^c`,
`d=(g^2+g'^2)/24`; projecting the original scalar-current subpacket first
would give zero. The explicit sterile Yukawa variation accounts for the
difference. The [independent active audit](INDEPENDENT_ACTIVE_AUDIT.md)
checks the complete projection and its Grassmann/derivative symmetries.

The main displayed exhaustion assumes the stated neutrino sources.
Retained Standard Model spectators have their additional source remainder,
described separately in the projection note. In particular, the
differentiated charged-lepton source cannot be dropped merely because
`H^dagger phi=0`.

These projection formulas start from the original `P_red`. They are not
the active coefficient table of the alternative row-derivative
representative introduced above. Changing derivative placement or
performing the additional scalar-current map requires the same full
action/source transport before coefficient tables can be combined.

## General complex Majorana masses

The [mass-covariance extension](COMPLEX_MASS_COVARIANCE.md) removes the
real-mass restriction from the selected chiral source recurrences and
integrated field-map rows. For constant `M=M^T`, the opposite equation
contains `M^dagger`. The ordered powers are `M M^dagger` and
`M M^dagger M`; the active matrix becomes `K_M=F M M^dagger Y_nu^T`.
Exact checks retain complex unitary active/sterile changes, singular and
zero masses and independent nonzero sources. No mass is inverted.

## Reproduction and scope

Run the five standalone scripts in the root pinned environment:

```sh
python -I -B -O verify_curvature_derivative.py
python -I -B -O independent_curvature_audit.py
python -I -B -O verify_active_source_projection.py
python -I -B -O independent_active_audit.py
python -I -B -O verify_complex_mass_covariance.py
```

Each imports no project calculation and uses explicit failures that remain
active under optimization. Adjacent receipts distinguish exact matrix,
component and exterior-polynomial controls from rejected mistakes. Root
replay operates on a disposable copy. These are internal AI-assisted
derivation/implementation audits, not independent human peer review.

The [literature comparison](LITERATURE_DERIVATIVE_PROJECTION.md) records a
bounded search and primary-source checks. A recent helicity/RG benchmark
has applicability assumptions not established by this selected off-shell
packet. Its theorem is not used to infer a beta-function zero here.

The full retained-sterile derivative/EOM basis, total physical coefficient
map, other loop/counterterm sectors, lower-order insertions, dimensional
and evanescent definitions, quantum probe/gauge-fixing transport, finite
matching and an observable remain unfinished. The new identities do not
establish unification, new fundamental mathematics, external novelty or
an empirical result. The candidate is for GitHub review; Zenodo continues
to identify the unchanged earlier baseline.
