# When a zero four-lepton pole conceals a finite momentum response

At the kinematics of the [finite triangle calculation](FINITE_TRIANGLE.md), the logarithmic pole weights each retained sterile mass by the same loop residue. The finite difference weights each mass by a different function of momentum. This gives a precise cancellation test and a counterexample to inferring a finite zero from the vanishing pole.

## Cancellation by mass group

Use a positive common scalar squared mass `a`, spacelike `s=Q²`, and a finite collection of Takagi masses. Write `M=U diag(m_i) Uᵀ` with `m_i≥0`. Let `y_i` be column i of `Y* U`. In each group of equal positive mass m define a complex symmetric flavor tensor

```text
A_m = m Σ_{i: m_i=m} y_i y_iᵀ.
R_UV = Σ_m A_m,
R(s) = Σ_m A_m F(a,m²,s).
```

The transpose follows from the anomalous Majorana propagator. Relative complex phases can therefore cancel a group. A real orthogonal change of Takagi basis within a positive degenerate group preserves its tensor. Massless modes contribute zero to this numerator. The common spin tensor, seed C, and loop normalization can be restored from the [ordered residue](../../checkpoint_2026_10_23/seed_yukawa/FOUR_LEPTON_POLE.md).

**For the ordered graph with soft seed leptons, R(s) vanishes at every spacelike momentum if and only if each positive-mass group tensor vanishes.** It is enough to have a zero on an open spacelike interval. A componentwise statement follows by applying the same argument to any fixed linear projection of the tensor. At a single momentum additional accidental cancellations are possible.

To prove the converse, continue the parameter integral along the real axis toward its first normal threshold. Put `A=√a`, `T=(A+m)²`, `x*=m/(A+m)`, `d=ax+m²(1−x)`, and `c=x(1−x)`. For `s>−T`, the denominator `d+cs` is positive on the compact parameter interval, and the kernel is real analytic. The identity

```text
d − cT = (A+m)² (x−x*)²
```

locates its first zero in the interior. For `s=−T+δ` with `δ↓0`, the log term of F remains integrable. In the rational term, write the integrand as `ax/d − ax/(d+cs)`. Rescaling `x−x*=√δ v` gives

```text
lim_{δ↓0} √δ F(a,m²,−T+δ)
      = −π a x* / [(A+m)√(x*(1−x*))].
```

The coefficient is nonzero for a>0, m>0. Away from x* the rescaled contribution goes to zero. Near x*, c is bounded below by a positive constant, so a multiple of `1/(1+v²)` bounds the rescaled integrand. Every heavier group's threshold lies strictly farther away and its kernel stays finite here. If the sum is zero on an open spacelike interval, real analyticity extends that zero to the first threshold. Multiplication by √δ then forces the lightest group tensor to vanish. Removing that group extends the argument to the next positive mass. Repeating over the finite set proves necessity. Sufficiency follows directly from the sum.

Any fixed finite polynomial in s stays regular at these thresholds. Adding such a local term cannot remove the first nonzero group's singularity. Local operators can change the Taylor coefficients about s=0, including the slope. This argument tests the mass dependence of the specified graph in fixed coordinates. A threshold of this off-shell function is not by itself a predicted production rate.

## An exact pole zero and a nonzero finite value

Choose units in which the scalar mass is one. Take one active flavor, `C=1` in the chosen units, two sterile masses 1 and 4, and the conjugated Yukawa row `(1,i/2)` in their Takagi basis. Their pole weights are `1` and `4(i/2)²=−1`. The summed four-lepton pole of this graph vanishes exactly. The ordered finite difference, stripping the common spin tensor and loop prefactor, is `F(1,1,s)−F(1,16,s)`.

For positive b, differentiation of the convergent kernel gives

```text
∂ₛF(a,b,0) = −∫₀¹ b x(1−x)² / [ax+b(1−x)]² dx.
```

Its magnitude is strictly smaller than `1/(2b)` because the integrand is smaller than `x/b` in the interior. For b=a it equals `1/(12a)` in magnitude. In the example the slope difference is thus strictly below `−1/12+1/32=−5/96`. This proves a finite response near the subtraction point without having to resolve two nearly equal numbers numerically. The numerical slope is −0.06247262904, and the ordered difference at s=1 is **−0.06184365542**. These are reduced values in the stated units. The mirrored graph, with the seed on the two legs carrying momentum, has its own routing.

As a control, choose both masses equal to one and take `(1,i)`. Then the only group tensor vanishes, so the whole ordered difference vanishes for every s. Replacing the squared complex couplings by their squared magnitudes destroys both pole cancellations. The [standalone checker](verify_mass_resolution.py) tests the threshold identity, three numerical approaches to the branch point, both flavor examples and the bound. [Numerical receipt](MASS_RESOLUTION_RESULTS.json)

## Both orientations at the same four external momenta

In the one-flavor example the weights and the seed are real. Reversing the four odd fields gives one minus sign for each antisymmetric spin bilinear, so the conjugate has the same sign in the common ordered pair tensor. Keep the momenta of all labelled external legs fixed. In the mirror routing the propagators carry the same loop momentum and the whole-composite numerator contains the external momentum. Its subtracted contribution is `s B(a,b)`, where `B(a,b)=∫₀¹ x/[ax+b(1−x)] dx`. Thus this mirror adds a linear polynomial to each mass group's soft-seed response. The nonanalytic thresholds of the soft-seed orientation persist.

The first derivative at s=0 of the sum of the two kinematic functions for mass squared b admits the radial representation

```text
D(a,b) = ∫₀∞ [t/(t+a)²] [1/(t+b) − b t/(t+b)³] dt.
```

The first term is the mirror derivative and the second comes from the soft-seed angular integral. For b>0 the bracket is positive and its derivative in b is `−(2t²+b²)/(t+b)⁴`. The integral converges at both ends, and differentiation is locally dominated for a,b>0. **The paired slope strictly decreases with b.** Hence in this one-flavor projection any two unequal positive masses with equal and opposite real pole weights leave a nonzero paired slope, provided their common weight and real seed are nonzero. This slope is defined for the isolated graphs in the declared coordinates. Independent derivative contacts could shift the slope of the EFT vertex.

For the masses 1 and 4, the mirror slope difference is `1/2−(16 log(16)−15)/225`. Including both orientations yields a **paired slope of +0.30703217293** and a **paired response at s=1 of +0.30766114655**. Omitting the mirror would even give the opposite sign in this benchmark. The common spin tensor and factor `1/(32π²)` have been stripped throughout. The subtraction is at s=0. Neither an absolute four-lepton coefficient nor a scattering rate follows from these numbers.
