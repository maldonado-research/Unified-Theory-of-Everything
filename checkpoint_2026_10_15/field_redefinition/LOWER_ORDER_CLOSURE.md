# A lower-order insertion required at inverse-scale order three

15 October 2026. A declared algebraic benchmark, prepared with AI assistance.
This is a concrete closure test, not another evaluated loop or a new operator.

Let `t=1/Lambda`, `rho=H^dagger H`, and let `n` be a column of left-handed
sterile Weyl fields. Spinor contractions in `n^T M n` are implicit. Declare

    L_kin = (i/2)[n^dagger K barsigma^mu partial_mu n
                     -(partial_mu n^dagger) K barsigma^mu n],
    K = 1 + t^2 rho Z,                 Z=Z^dagger,
    L_mass = -(1/2)[n^T (M+t rho C+t^3 rho^2 D) n + h.c.],
    L_source = -[n^T Y J + h.c.],       M=M^T, C=C^T, D=D^T.

`J` is a declared gauge-singlet, Grassmann-odd left-Weyl light-field source,
independent of `n` and without a hidden `t` expansion in this benchmark,
with canonical dimension `5/2` (for example a lepton-Higgs
composite in a consistent charge-conjugation convention). Here `M` has
dimension one, `rho` dimension two, and `Z,C,D,Y` dimension zero. The
dimension-five sterile-scalar and dimension-six derivative operators are
optional additional inputs; their coefficients have NOT been obtained from
the selected seed calculation. In particular this benchmark is not a claim
that a nonzero `Z` or `C` is generated in that calculation.

The first canonical field change is

    n = (1 - t^2 rho Z/2) n_c + O(t^4).

This is a formal local expansion about `t=0`. A finite square-root
normalization would additionally require a positive invertible `K` in its
domain; no such global assertion for arbitrary `rho,Z` is used here.

It gives `S^dagger K S=1+O(t^4)`. For nonconstant `rho`, the symmetrized
kinetic action above also contains the connection
`(i/2)n_c^dagger[S^dagger K (partial_mu S)
 -(partial_mu S^dagger)K S]barsigma^mu n_c`. Its order-`t^2` term cancels;
indeed `S`, `K` and their derivatives here are polynomials in the same
constant Hermitian `Z`. This cancellation is special to this declared
normalization problem. It is not permission to discard derivative contacts
from a general field-dependent or derivative field redefinition.

Direct substitution into the entire mass and source action gives

    M_c(rho) = M + t rho C
          - (t^2 rho/2)(Z^T M + M Z)
          + t^3 rho^2 [D - (Z^T C + C Z)/2] + O(t^4),
    Y_c(rho) = Y - (t^2 rho/2) Z^T Y + O(t^4).

The transpose in these holomorphic bilinears is essential. Hermitian
conjugation instead belongs in the kinetic congruence. Both the mass-shaped
order-two term and the dimension-seven coefficient shift are symmetric.
No inverse of a mass or Yukawa matrix occurs; zero or singular matrices are
admissible.

Thus a dimension-six normalization acting on a dimension-five interaction
changes the order-`Lambda^-3` action by

    Delta D = -(Z^T C+C Z)/2.

Applying only the renormalizable equations to the dimension-six redundant
term would miss this contact. The checker exhibits explicit complex,
noncommuting examples in which that omission changes the order-three
action, as well as a rank-deficient mass example. This makes one previously
listed closure gap executable. It does not compute the missing coefficient
vector of the retained-sterile EFT or settle any physical running.

More generally, with `S=S4+t S5+t^2 S6+t^3 S7` and a field change
`Phi=Phi_c+t^2 f2+t^3 f3`, Taylor expansion of the action yields

    S7,c = S7 + delta_f2 S5 + delta_f3 S4.

The formula assumes no order-one field shift; if there is one, its mixed
and higher variations must also be included. The selected packet's own
order-three field change needs renormalizable sources through this order,
but that fact does not eliminate `delta_f2 S5` from the complete EFT.

Under a unitary sterile coordinate change `n=U n'`,
`Z'=U^dagger Z U`, `C'=U^T C U`, and
`Delta D'=U^T Delta D U`. These transformation laws are part of the
checker. The note concerns conventional field-coordinate algebra, with no
claim of mathematical novelty, a physical beta function or an observable.

Run the standalone checker with a Python environment containing SymPy:

```sh
python -I -B -O verify_lower_order_closure.py
```

The tests substitute the full polynomial field change before extracting
orders, compare it with the formulas above, check local kinetic connection
terms and complex basis covariance, and require omissions/adjoint errors to
fail. They use explicit exceptions that remain active under optimization.
