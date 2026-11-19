# Which sterile mass combinations can the on-shell cut distinguish?

The [off-shell test](../four_lepton/MASS_RESOLUTION.md) separated masses by their different thresholds. In the on-shell process the pair crossing the cut consists of two scalars of the same mass. All exchanged sterile masses therefore share a threshold. In one projection there are exact cancellations between three or more distinct masses over the entire cut. Tagging the weak components and retaining the angular dependence supplies a second independent test.

## Amplitudes and the two angular terms

Use massless external leptons in the unbroken theory. The complex scalar doublet has common squared mass a>0. For a concrete flavor-changing amplitude put the derivative coupling only on the outgoing e pair and the Yukawa couplings only on the incoming μ flavor. These labels specify doublet flavors, with the weak component designated separately. This selects a single seed orientation. In each group with positive Takagi mass m define a signed or complex weight `w_b=C_ee m Σ(Y* U)_{μi}²`, summed over the group with `b=m²`. Massless modes carry zero weight. The proof works componentwise for any fixed effective weights multiplying the same kernel.

Take physical `s>4a` and set

```text
β=√(1−4a/s),    N=s/2−a,    h=sβ/2,
J_b=log[(N+b+h)/(N+b−h)]/(2h).
E(s)= N Σ_b w_b J_b,       O(s)= Σ_b w_b [1−(N+b)J_b].
```

Here E and O denote the angle-even coefficient and the coefficient multiplying the cosine of the scattering angle. With the common external spin tensor and phase-space normalization stripped, the weak direct and exchanged routings give `E±cosθ O`. External legs of the same weak component select their sum `2E`. For an ordered mixed pair with initial components (1,2), the final ordered components (1,2) select the direct routing, while (2,1) select the exchanged routing. Equivalently use the direct routing at opposite angle cosines and strip its spin factor consistently. The massless left-handed spinor contraction has magnitude s. The weak assignments can therefore distinguish these kinematic coefficients at amplitude level. [Loop normalization](ON_SHELL_TRIANGLE.md) · [Cut and external states](INDEPENDENT_CUT.md)

## Exact degeneracy in the equal-component channel

Let `v=(N+h)/a=(1+β)/(1−β)>1`. The endpoints of the angular integral obey `(N+h)(N−h)=a²`. For a dimensionless squared mass `r=b/a`, put

```text
L_r(v)=2h J_b=log[(r+v)/(r+1/v)].
L_r(v)+L_(1/r)(v)=2 log v=2 L_1(v).
```

Thus a group at b and one at a²/b with equal weights produce exactly the same even cut as twice that weight at b=a. The identity is exact for every s>4a.

For a finite spectrum there is also a converse. The derivative with respect to v is the rational function

```text
dL_r/dv=1/v + 1/(v+r) − 1/(v+1/r).
```

Write W_r for the weight at b=ar, taking missing weights to be zero. If the even cut vanishes on an open energy interval, this rational derivative of the weighted sum vanishes identically. Its simple-pole coefficients at v=0 and v=−r give respectively `Σ W_r=0` and `W_r=W_(1/r)` for every r≠1. These conditions suffice because the derivative then vanishes and all L_r equal zero at v=1. They classify all equal-component null combinations for finitely many positive masses. They include a pair of reciprocal masses compensated by a self-reciprocal mass and combinations involving several reciprocal pairs. For a pair of distinct masses alone they force both weights to vanish.

If the coupled spectrum lies entirely at or above the scalar mass, every group with b>a lacks a reciprocal partner. The even cut then has no nontrivial null combination, including at b=a. The same argument applies to spectra entirely on the other side of a. A nonzero reciprocal cancellation must straddle the scalar mass. In the model used here a coupled sterile lighter than the scalar permits a scalar decay to a sterile and a massless lepton. The cut with two free scalar lines is a well-defined fixed-order result; interpreting a finite-width line shape or a measured rate requires the decay channels. [Perturbative assumptions](CHANNEL_SCOPE.md)

## How tagging the weak components removes the degeneracy

Suppose the even coefficient and the odd coefficient both vanish on the same open interval. From the even result `Σw_b=Σw_b J_b=0`. The odd coefficient reduces to `−Σ b w_b J_b`, so the same residue argument applies with weights `b w_b`. For a non-self-reciprocal pair we obtain

```text
w_b = w_(a²/b),    b w_b = (a²/b) w_(a²/b).
```

The two equations force both weights to zero. Only b=a could remain and the zero sum removes its weight as well. Conversely zero weights give zero in both projections. Thus, given the two amplitude functions over an interval and known a, there is no exact ambiguity between distinct finite positive-mass group weights in this model. Degenerate states within a group enter only through their summed coupling. The proof uses access to complex amplitude coefficients or signed real interference terms. Rates, finite-precision samples and detector-level extraction require additional information.

## A realizable example and an angular check

In units with a=1, choose the three Takagi masses `(1/2,1,2)`, `C_ee=1` in the chosen units and the conjugated incoming Yukawa row `(√2,i√2,1/√2)`. The mass-weighted contributions are `(1,−2,1)`. A common rescaling of all Yukawa couplings preserves the cancellation. Taking squared magnitudes before summing would incorrectly give a positive ultraviolet weight.

The equal-component cut is zero for all s>4. With final components (1,2) and `cosθ=1/2`, direct four-dimensional angular integration instead gives **β(E+cosθ O)=−0.04010098542 at s=9**. Exchanging the weak routing with the same angular convention gives +0.04010098542. The normalization strips the common spin tensor and the cut factor `1/(16π)`. These numbers are signed amplitude contributions.

There is an analytic check near threshold. Since `(N+b)J_b=1+h²/[3(N+b)²]+O(h⁴)`, the constant cancels in each odd contribution. For these weights `Σ w_b/(1+b)²=9/50` and `h=2β/(1−β²)` in the chosen units. Hence

```text
β(E+cosθ O)=−(6/25) cosθ β³+O(β⁵).
```

This suppression matters if one tries to use the angular signal close to threshold. A null angle-integrated cut alone does not establish the absence of the sterile exchange. The [checker](verify_mass_identifiability.py) tests the exact identities, the null space of a sample residue matrix, direct angular integrals, complex coupling weights and the threshold limit. [Numerical receipt](MASS_IDENTIFIABILITY_RESULTS.json)
