# Action sign, local operators and sourced contacts

15 October 2026. A continuation of the one-hypercharge-leg candidate,
prepared with AI assistance. This packet uses the previously evaluated loop
residues and advances their action interpretation and sourced reduction.
It does not calculate a complete anomalous dimension or an observable.

## Results

1. **The action sign is fixed in a literal convention.** The displayed
   `bar L ... N_L` channel propagates `L^c` and `tilde H`, with hypercharges
   `-y_L` and `y_L`. Their weighting is `y_L(R_H-R_L)`. This corrects the
   earlier physical charge labels without changing any integral tensor.
   With `F=C_(7,seed)Y_nu*`, the declared zero/one-B poles are canceled by

       L_ct = F[O_J + g' y_L Q_T]/(16 pi² epsilon) + h.c.,
       Q_T  = Q_sigma/4 + Q_div/3,
       Q_sigma = i(D_mu bar L) sigma^{ab} gamma^mu tilde H N_L B_ab,
       Q_div = bar L gamma^mu tilde H N_L partial^nu B_{nu mu}.

   Read [ACTION_SIGN.md](ACTION_SIGN.md) for the literal action, Fourier
   rules, slot normalization and charge-conjugated open line. Its two-B
   contacts are predictions of the chosen covariant lift, not evaluated
   two-B loop amplitudes in this action audit. The subsequent
   [seven-diagram continuation](../two_gauge_leg/00_READ_FIRST.md) now
   matches those contacts by actual integration with independent checks.
2. **The declared off-shell transverse space has five directions.** It
   assumes independent generic momenta, locality, degree two, the specified
   odd Clifford chain and constant coefficients. The loop occupies a
   particular vector in that space. An independent mass-squared transverse
   direction is excluded under these premises. This is not the full EFT
   operator basis or an EOM quotient.
3. **Both source-reduction routes are explicit.** The local operators give
   Higgs-derivative ordinary/dual field-strength structures, a retained-mass
   dipole, a sterile-Yukawa source and the full hypercharge current. The
   current includes the Higgs and all charged SM fermions. Dropping a term
   because it vanishes for a real gauge wave loses its sourced contacts.
   See [OPERATOR_MAP.md](OPERATOR_MAP.md).
4. **The two flavor source tensors behave differently.** For symmetric C,
   the sterile tensor `Y† C Y*` is symmetric and kills that Majorana dipole
   source. The active tensor instead has antisymmetric part

       K = (C H^T-H C)/2,    H=Y Y†.

   It can be nonzero. In an eigenbasis of H it vanishes precisely when C
   joins only equal-eigenvalue blocks. [FLAVOR_SOURCE.md](FLAVOR_SOURCE.md)
   proves the exact kernel/rank and norm statements. This is a conditional
   source coefficient, not a net physical transition moment.

## Evidence and limits

The producer, independent critic and exact Ward-kernel classifier test
Fourier placement, actual action signs, charge conjugation, off-shell
independence, IBP/Bianchi relations, source/current signs, complex flavor
order and deliberate errors. See [INDEPENDENT_OPERATOR_AUDIT.md](INDEPENDENT_OPERATOR_AUDIT.md)
for independent review and its precise coverage. These are internal
AI-assisted checks, not independent human peer review.

[SOURCE_BASIS.md](SOURCE_BASIS.md) checks the established retained-sterile
operator family and sourced equations against three versioned primary
texts. Only request metadata is stored in [SOURCE_REQUEST_RECEIPTS.json](SOURCE_REQUEST_RECEIPTS.json);
full papers and private source files are excluded. No new fundamental
operator, mathematics or physical discovery is claimed.

The gamma-three representative and its polynomial identities can be kept
without dualization. The dual-field-strength map is explicitly
four-dimensional, with fixed epsilon/gamma5 conventions; it is not an
evanescent subtraction prescription. EOM reductions use the gauge-invariant
reference action, before gauge fixing. Coefficients of independent lower
EFT sources require their own declared inverse-scale counting.

The gauge-current descendant has an additional g', and Yukawa sources
have an additional Y. They belong to the correlated reduced packet even
though a four-point, single-g' or single-Y graph selection does not compute
those higher-multiplicity sectors separately. A complete coefficient needs
the remaining diagrams, wave-function/parameter counterterms, accepted
basis/scheme and finite threshold matching. The old withdrawn coefficient
and the baseline DOI remain unchanged.

## Reproduce

The root runner executes three new standard-library stages on disposable
copies: `verify_action_sign.py`, `verify_operator_map.py` and
`verify_transverse_space.py`. Use `python -B verify.py --full` from the
repository root for the complete replay.

The additional independent symbolic audit uses the existing pinned
SymPy review dependencies:

```sh
python -B checkpoint_2026_10_15/operator_map/verify_operator_critic.py
```

No service or new credentials are needed for the mathematical checks.
