# Independent exact cubature audit of the three-gauge UV pole

15 October 2026. This audit independently integrates the 33 specified matter-loop
1PI topologies, with the literal single-seed and one-source orientation. It is
restricted to this Green-function sector and is not a complete EFT calculation,
an anomalous dimension, a finite amplitude, or a physical dipole prediction.

## Independent construction and result

The auditor first derived the kernels and completed its coefficient table without
reading the producer implementation. It uses exact Gaussian-integer angular
cubature rather than the producer's monomial Wick-pairing tensor functional.
No producer or repository module is imported by the standalone checker.

Let `F=(f1,...,fnf)` and `S=(s1,...,sns)` list ordinary attachments from Yukawa to
seed. The crossed seed has the leading polynomial `q^2`. Cancelling it against
the unshifted fermion denominator leaves the ordinary leading numerator

`(-2)^ns product(q.e_s) /q E_fnf /q ... E_f1 /q`.

It has angular degree four. There are six permutations for each `nf=0,1,2,3`,
hence 24 ordinary diagrams. The scalar seagull gives three diagrams with the
remaining ordinary vertex on the fermion and six with it on the scalar, for
nine more diagrams. Their remaining numerators have angular degree two, with
one additional factor `e_i.e_j`.

The conjugate fermion representation is `-T^T`. Contracting the literal crossed
seed and source tensors yields

`scalar_chain * transpose(fermion_chain) = (-1)^nf reverse(S) F`.

A seagull replaces a scalar pair by `{T_i,T_j}`. Its phase relative to an
ordinary insertion is `i*i=-1`, whereas an ordinary propagator-vertex pair has
`i*(-i)=+1`. The baseline seed/fermion/Yukawa/scalar phase is `-1`. These literal
phases also reproduce the already fixed commuting two-gauge weights: `-1` for
one ordinary attachment on each line and `-2` for the scalar seagull. No
three-gauge target is used to normalize the integration.

After integration, all 64 polarization basis triples agree with the separately
specified cubic covariant-action kernel. With `E_i=/epsilon_i`, free color
letters `X_i`, and `k` the unused label for each `i<j`, the common residue is

```
R3 = -1/4 sum_(i<j) X_k [X_i,X_j] [E_i,E_j] E_k
     -1/3 sum_(i<j) [X_k,[X_i,X_j]]
          ((epsilon_k.epsilon_i) E_j - (epsilon_k.epsilon_j) E_i).
```

The loop pole has the already calibrated overall negative amplitude sign;
the local counterterm has the opposite sign. Couplings can be placed in each
`X_i`, or factored as `g^3` when all gauge couplings are equal. In particular,
for `X=(-I/2, sigma1/2, sigma2/2)` and polarizations `(e0,e0,e1)`, the result is
`i gamma^1 sigma3/8`, with the mixed coupling `gprime*g^2` stripped. This witness
retains the imaginary second Pauli generator. The commuting three-gauge residue
vanishes; the mixed BWW and noncommuting weak structures do not.

## Exact angular integration and UV/IR separation

Use the 24 vectors with two nonzero entries `+1` or `-1`, and replace spatial
components by `i` times those entries. Before normalization each vector has
Minkowski norm two. Dividing by `sqrt(2)` gives the complex Wick-rotated 24-cell
sphere points. Uniform averaging reproduces the full second and fourth sphere
moments exactly:

`<q_mu q_nu> = g_mu_nu/4`,

`<q_mu q_nu q_rho q_sigma> = (g_mu_nu g_rho_sigma + g_mu_rho g_nu_sigma + g_mu_sigma g_nu_rho)/24`.

The checker exhausts all coordinate monomials through degree four. It evaluates
at the unnormalized Gaussian-integer points, then divides a degree-four average
by four and a degree-two average by two. Every operation uses integer pairs or
exact rational numbers; there is no floating-point quadrature.

The original integrals are defined at generic off-shell momenta. A common
positive auxiliary mass in the asymptotic denominators separates the radial UV
logarithm from IR behavior. The radial simple-pole normalization is
`i/(16*pi^2*epsilon)`. Mass corrections and external momentum shifts lower the
UV degree and are finite here. Replacing the exact cancellation
`q^2/(q^2-Maux^2)` by its leading one changes only UV-finite terms. No original
massless integral is set to zero as scaleless.

For `nf` ordinary fermion, `ns` ordinary scalar, and `ng` scalar-seagull
vertices, the crossed seed superficial UV degree is

`4 + (nf+ns+3) - 2*(nf+ns+ng+2) = 3-n`,

where `n=nf+ns+2ng` is the external gauge count. Consequently the three-gauge
pole has no external momentum or mass dependence. The direct seed has degree
`1-n` and is UV finite here. The same restricted one-matter-loop family with
`n>=4` is UV finite by power counting. This statement excludes additional
internal gauge lines, interactions, loops, or operator insertions.

## Coverage and reproducibility

Run:

```
python3 -I -S -B -O verify_three_gauge_cubature.py --coefficients THREE_GAUGE_COEFFICIENTS.json
```

`THREE_GAUGE_CUBATURE_RESULTS.json` records **10,407 explicit checks**. This number includes low-level
arithmetic and moment assertions; it is not a count of physical discoveries or
independent experimental tests. The substantive coverage is:

- 33 literal topology color contractions, with exact dense 3x3 generators whose
  six cubic word matrices have rank six.
- 64 polarization basis triples, all six free noncommutative color words, and
  the full Clifford algebra. Multilinearity and zero external momentum degree
  make this an exhaustive coefficient test of the specified local tensor.
- 128 free-word Bose relabeling checks and 64 commuting-color zero checks.
- Deliberate omission/reversal controls: seagull omission and phase reversal
  each fail at 40 basis triples; reversed fermion color order fails at 24;
  omitted one three-fermion ordering fails at 36; omitted sigma-curvature and
  adjoint-curvature terms fail at 60 and 36 respectively.

After the independently generated table was complete, a separate cross-check
compared every raw diagram with the producer before comparing totals.
`THREE_GAUGE_COMPARISON_RESULTS.json` records **2,112 individual diagram agreements**
(`33*64`) and **64 full six-free-word Clifford tensor agreements**. Reproduce
this receipt by putting both integration cores and the comparison harness in
one directory, then running:

```
python3 -I -S -B -O compare_three_gauge_routes.py
```

The comparison harness explicitly imports both adjacent integration cores and
writes its JSON result to standard output. It writes no files. Each integration
core remains standalone and does not import the other route. The harness was
verified in an isolated scratch directory containing copies of the three files,
with the same layout intended for publication.

The result is an independently checked non-Abelian completion of the selected
local pole through its final UV-divergent external-gauge multiplicity. It does
not establish complete counterterm closure across other external-field sectors
or yield physical running.
