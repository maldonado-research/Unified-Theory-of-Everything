# Two external hypercharge fields: seven-diagram continuation

15 October 2026. Provisional research candidate; continue from the literal
action in [ACTION_SIGN.md](../operator_map/ACTION_SIGN.md).

The seven actual one-loop 1PI diagrams at one derivative-seed insertion,
one conjugate neutrino Yukawa and two hypercharge couplings give

\[
R_{BB}=(\varepsilon_1\cdot\varepsilon_2)\not p
-\frac14\bigl([\not k_1,\not\varepsilon_1]\not\varepsilon_2
+[\not k_2,\not\varepsilon_2]\not\varepsilon_1\bigr).
\]

For `F=C_(7,seed)Y_nu*`, `y_L=-1/2`, `p=r+h+k1+k2`, metric `(+---)` and
the declared D-plus/Fourier convention, the pole is
`-i F (g' y_L)^2 R_BB P_L/(16 pi^2 epsilon)`.
It is cancelled by the two-B content of the previously chosen covariant
counterterm action. Both fermion orderings, both scalar orderings, both
mixed assignments and the canonical scalar seagull are included. The
residue is scalar-mass independent; the seagull and field-strength
connection are essential. The direct seed's `2t^2` term is UV finite in
this vertex, so this calculation does not independently retest its 2:1
pairing weight.

Read the evidence in this order:

1. [Derivation](TWO_B_DERIVATION.md): vertices, phases, seven routed
   integrands, simplex moments and individual residues.
2. [Independent audit](INDEPENDENT_TWO_B_AUDIT.md): large-momentum angular
   extraction, all individual residues, both Ward contractions and
   deliberate-error controls.
3. [Vertex audit](VERTEX_WARD_AUDIT.md): exact unintegrated Ward identity
   with its loop-momentum shift term and charge-conjugate chain.
4. [Source and method audit](SOURCE_METHOD.md): primary action/Ward
   references, precise local-polynomial uniqueness assumptions and the
   conditional three-B consequence.

From this directory, these standard-library commands use explicit failures
that remain active under optimization:

```sh
python3 -I -S -B -O verify_two_b.py
python3 -I -S -B -O verify_two_b_angular.py
python3 -I -S -B -O verify_two_b_vertices.py
```

The first two check all 256 coefficients of the degree-one external-momentum
polynomial, bilinear in polarizations, alongside generic rational cases.
These coefficient evaluations act on the UV polynomial defined at generic
off-shell momenta; they do not set exceptional scaleless integrals to zero.
JSON files record each checker's output; retrieval receipts record primary
source requests without republishing downloaded papers.

The local degree-one two-B completion is unique **given** the complete
lower-point kernels, the same Ward hierarchy and generic independent
off-shell momenta. No full EFT basis or physical anomalous dimension is
established. The three-B zero-pole consequence is conditional on locality,
power counting and that hierarchy; no three-B loop was integrated. Finite
terms, EOM reduction with all sourced contacts, other coupling orders and
observable predictions remain open. The baseline DOI does not archive
this provisional continuation.
