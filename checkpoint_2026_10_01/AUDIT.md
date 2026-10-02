# Independent skeptical review of the local sourced-EOM checkpoint

1 October 2026 (America/Los_Angeles). This is an independent AI-assisted algebra review, not external
human peer review, experimental validation, or proof-assistant certification.

## Verdict

No correctness blocker was found in the bounded claims of
`GREEN_REPRESENTATIVE.md` and `SOURCED_STERILE_REDUCTION.md`. The six-term
operator image reproduces the declared public three-point UV polynomial.
Its sourced reductions demonstrate interactions omitted by a free-field
replacement, and a correlated leading Schur identity is valid in the declared
positive diagonal sterile-mass basis. None of these results determines the
full physical dimension-seven mixing vector or a unification result.

## Independent checks

The reviewer first reconstructed the relevant identities analytically, then
read the producer scripts and ran their optimized controls. The following
summarizes analytic and computational controls. The included
`independent_clifford_checks.py` separately checks the Clifford algebra,
noncommuting-connection curvature, gauge-field trace and complex flavor
orientation with exact SymPy arithmetic.

- **Fourier signs and chirality.** The incoming `i partial` gives `p`, the
  outgoing row `-i D` gives `r`, and the negative Laplacians give the positive
  squared momenta. Since `bar L=bar ell P_R`, an odd gamma chain ending in
  `N_L` is allowed. Every monomial has dimension seven.
- **Light scalar source.** The equation from the stated potential gives
  `-D² tilde H = m_H² tilde H + 2 lambda (H†H) tilde H - S_tildeH`.
  The scalar-mass multiple of `slash p` is therefore `1/2+1/6=2/3`.
- **Gauge curvature.** With `D=partial+i Gamma` and `[D,D]=i F`,
  `(i slash D)²=-D²-sigma F/2`. Dirac adjunction yields the positive
  `bar L F sigma/2` in the row-Laplacian reduction. An explicit constant,
  noncommuting SU(2)-matrix connection verifies that row identity; reversing
  its curvature sign fails.
- **Sourced sterile identities.** Iterating the two conjugate EOM gives all
  three source terms in the cubic derivative. Defining nonzero residuals
  reproduces the exact off-shell identity in the independent audit.
  Merely substituting `A³N_L=M³N_R` is invalid for nonzero sources.
- **Complex flavor orientation.** A nonreal rectangular example verifies
  `delta C5 + delta(Y M^-1 Y^T)=0`. Replacing the required transpose by a
  Hermitian conjugate produces a nonzero residue. The signs follow from the
  negative Yukawa and positive half-normalized composite contact in the
  declared Weyl action. This check assumes invertible positive diagonal `M`;
  arbitrary-basis formulas need the appropriate conjugation rules.
- **Gauge ambiguity is physical at the vertex level.** The displayed
  `bar L tilde H gamma^mu i partial^nu N_L B_mu_nu` has a vanishing
  gauge-free three-point kernel. At the stated four-leg kinematics its
  transverse physical-polarization kernel has spin-summed square `4`.
  Independently computed gamma matrices and an exact trace reproduce that
  value. Its Ward contraction vanishes. Thus adding its coefficient can
  change a gauge-boson amplitude while preserving the matched three-point
  polynomial. This does not establish a contribution to the same scalar
  beta coefficient at the same coupling order.

## Test assessment

`python -I -S -B -O verify_green_representative.py` passes 605 exact controls;
`python -I -S -B -O verify_sourced_sterile_identity.py` passes 640.
Both use explicit exceptions, so optimization does not remove their checks.
The wrong-row-sign, dropped scalar descendant, wrong-curvature,
dropped-source, wrong-complex-transpose and missing finite-quadratic-contact
mutations distinguish the particular failures they advertise.

The Fourier tests compare two polynomial representations, not a new loop
integration. The formal source-tower test is strengthened by the independent
script's residual identity and explicit solutions of both sourced equations.
The finite matrix tests verify a field-redefinition invariant, not running.
No test establishes a complete gauge basis, reconstructs all source-contact
amplitudes, or substitutes for full counterterm matching.

## Review correction and remaining limits

The review requested an explicit conjugate-channel coefficient map. The
independent audit now states that its undotted Weyl coefficient is
`Y = Y_nu*` relative to `-bar L Y_nu tilde H N_R+h.c.`, with the epsilon
convention translated consistently. The corresponding `F` is likewise
conjugated across channels. This corrected a potential ambiguity in comparing
the documents; it did not change the sourced algebra or any computed result.

The six monomials are a **representative**, with IBP/EOM relations not yet
quotiented. Gauge covariance does not imply basis completeness. The example
with an external gauge field establishes nonuniqueness of completion, not a
computed Wilson coefficient. The single-seed Yukawa-shaped pole matrix is
not an anomalous dimension. The cancellation of the leading Schur variation
holds for any coefficient of the removable cubic derivative and therefore
cannot select an unknown physical loop coefficient.

The explicit spinor, chirality and physical-polarization controls use
four-dimensional matrices. They do not classify evanescent structures in
`d=4-2 epsilon`, choose a complete gamma5 prescription, or verify finite
terms that can arise from dimensional continuation. The formal Clifford
identity is compatible with continuation, but that fact does not extend the
four-dimensional test package into a complete renormalized operator basis.
A dimensional-regularization calculation at a claimed physical order must
state and consistently use its physical/evanescent basis and subtraction
prescription.

The optional independent symbolic implementation is included for reproduction:

```sh
python -m pip install -r checkpoint_2026_10_01/requirements-review.txt
python -B checkpoint_2026_10_01/independent_clifford_checks.py
```

Run these commands from the repository root in its virtual environment. The
review pins SymPy 1.14.0 and mpmath 1.3.0. Its four exact control groups are
separate from the 1,245 standard-library producer checks and from the root
runner's numerical replay. The shipped JSON records the symbolic result.

The next physics calculation must retain one fixed full action and power
counting, match the additional external-leg Green functions and source
contacts, and include field/parameter renormalization and the subtraction
scheme. Finite matching and an observable require further work.
