# Independent routed mixed-bubble UV calculation

6 September 2026 Pacific. **The declared derivative-seed bubble does not have a UV pole proportional to the whole dimension-five answer times one external p².** The crossed Wick pairing produces additional independent off-shell momentum structures. This completes a concrete kinematic calculation requested by N00AK-r1; it does not complete dimension-seven Green-basis matching or fix a physical beta-function coefficient.

This derivation was obtained independently of the parent task's exact-arithmetic implementation. The original massless result and the finite-Higgs-mass extension were sent to the parent before reading its result file. The final signed polynomial agrees with its result after the explicit loop orientation change below.

## 1. Reference normalization and declared seed

Zhang's Appendix A, Eq. (35), uses an internal numerator `slash k`, denominator `k²[(k+p)²−m_H²]`, and the two SU(2) pairings. Its dimension-five pole is proportional to `−(3/2) slash p`. Eq. (37) then uses a **sourced** sterile equation of motion. The published calculation is dimension five; the derivative-insertion calculation in this memo is our own extension, not a result quoted from Zhang. [Primary reference, Appendix A](https://arxiv.org/html/2405.18017v3#A1).

Let p be the incoming sterile momentum, r the outgoing lepton momentum, and `h=p−r` the outgoing Higgs momentum. At the LLHH vertex use all incoming assignments

`r_a=−r`, `h_c=−h`, `r_b=−k`, `h_d=k+p`.

Their sum vanishes. For the declared composite-derivative seed, the pairings carry

`(r_a+h_c)²=p²` and `(r_a+h_d)²=(k+h)²`.

With `ε₁₂=1`, contraction of the internal Yukawa epsilon gives

`Σ_bd ε^(ad) ε^(bc) ε^(bd)=ε^(ac)`,

`ε^(ac) Σ_bd ε^(bd) ε^(bd)=2ε^(ac)`.

Thus the common spinor/flavor factors multiply the integral

\[
T^\mu=\mu^{2\varepsilon}\int\frac{d^dk}{(2\pi)^d}
\frac{k^\mu[2p^2+(k+h)^2]}
{(k^2+i0)[(k+p)^2-m_H^2+i0]},\qquad d=4-2\varepsilon.
\]

The overall derivative-vertex Fourier sign and precise map from the normalized seed to `W₁` are external conventions; they multiply the entire polynomial uniformly. They cannot remove its independent momentum structures. Chiral projectors and `W₁Y_A*` also remain outside this purely tensor calculation.

## 2. UV tensor integrals

Remove the common pole `C_UV=i/(16π² ε)` and use `[I]_UV/C_UV` for the following coefficients. Feynman parametrization shifts `l=k+xp`, with

`Δ=x m_H²−x(1−x)p²`.

The elementary pole residues are `∫1/(l²−Δ)² → 1` and `∫l_μl_ν/(l²−Δ)² → Δg_μν/2`. The factor 1/d in the tensor reduction may be evaluated at d=4 for the simple-pole coefficient; its O(ε) part affects finite terms. Exact integration of the parameter powers gives

\[
I_0=1,\qquad I_1^\mu=-\frac12p^\mu,
\]
\[
I_2^{\mu\nu}=\frac13p^\mu p^\nu
+\left(\frac{m_H^2}{4}-\frac{p^2}{12}\right)g^{\mu\nu},
\]
\[
I_3^{\mu\nu\rho}=-\frac14p^\mu p^\nu p^\rho
+\left(\frac{p^2}{24}-\frac{m_H^2}{6}\right)
\left(p^\mu g^{\nu\rho}+p^\nu g^{\mu\rho}+p^\rho g^{\mu\nu}\right).
\]

Two useful contraction checks are

`g_μν I₂^(μν)=m_H²`,

`g_νρ I₃^(μνρ)=−m_H² p^μ`.

The second follows directly from cancellation of k² in the numerator, shifting `l=k+p`, and retaining the massive tadpole. **It is zero in the massless limit but must not be dropped before that limit when m_H² is retained.** These identities also fix the sign of the mass terms.

The UV coefficient can be calculated at nonexceptional Euclidean/spacelike external momentum, so no on-shell infrared pole is confused with the UV pole. The result is a local polynomial and can then be analytically continued. In the massless limit the original nonexceptional bubble is infrared safe. It is not legitimate to set every external invariant and every mass to zero first and infer the UV result from a scaleless integral.

## 3. The complete polynomial for this declared bubble

The crossed contribution alone is

\[
\frac{[J^\mu(h)]_{\rm UV}}{C_{\rm UV}}
=\left[-m_H^2+\frac23p\cdot h-\frac12h^2\right]p^\mu
+\left[\frac{m_H^2}{2}-\frac{p^2}{6}\right]h^\mu.
\]

Adding the weight-two direct contribution gives

\[
\boxed{\frac{[T^\mu]_{\rm UV}}{C_{\rm UV}}
=-\left[p^2+m_H^2-\frac23p\cdot h+\frac12h^2\right]p^\mu
+\left[\frac{m_H^2}{2}-\frac{p^2}{6}\right]h^\mu.}
\]

For `m_H²=0`, `J(0)=J(p)=0`; both are exact cancellation checks. Removing the derivative factors entirely gives `3I₁=−3p/2`, reproducing the stated dimension-five reference normalization. Multiplying that answer by p² would instead give `−3p²p/2`, which is generally different from the boxed result.

In external-lepton variables `r=p−h`, an equivalent spinor form is

\[
\frac{[\not T]_{\rm UV}}{C_{\rm UV}}
=-\left(\frac56p^2+\frac12m_H^2+\frac13r^2+\frac16h^2\right)\not p
+\left(\frac16p^2-\frac12m_H^2\right)\not r.
\]

In the massless-Higgs limit, this exhibits four independent cubic off-shell structures: `p² slash p`, `r² slash p`, `h² slash p`, and `p² slash r`. Keeping the Higgs mass also gives `m_H² slash p` and `m_H² slash r`. A single scalar coefficient multiplying `p² slash p` cannot represent them off shell.

For example, take Minkowski signature `(+---)`, `p=(3,0,0,0)`, `h=(1,1,0,0)`, `m_H²=0`. The original-k vector is `T/C_UV=(−45/2,−3/2,0,0)`, whereas the common-p² ansatz is `(−81/2,0,0,0)`. This is an algebraic off-shell counterexample, not a physical decay prediction.

## 4. Sign agreement with the parent calculation

The parent uses `ℓ=−k` and denominator `ℓ²[(ℓ−p)²−m_H²]`, and records the polynomial for **positive** numerator `ℓ^μ[2p²+(ℓ−h)²]`. Since the original numerator becomes `−ℓ^μ`, that recorded polynomial is `−T`:

\[
J_+^\mu/C_{\rm UV}
=\left[p^2+m_H^2-\frac23p\cdot h+\frac12h^2\right]p^\mu
+\left[\frac{p^2}{6}-\frac{m_H^2}{2}\right]h^\mu.
\]

There is no disagreement in the mass or mixed-momentum terms. The exact coefficient vectors in the basis

`[p²p^μ, (p·h)p^μ, h²p^μ, p²h^μ, m_H²p^μ, m_H²h^μ]`

are `[-1,2/3,−1/2,−1/6,−1,1/2]` in the original-k convention and its negative in the parent's positive-ℓ convention.

The independent rational script [routed_loop_independent.py](routed_loop_independent.py) derives the parameter moments, checks the SU(2) weights and tensor contractions, assembles the weighted polynomial, and verifies this signed agreement. It passed isolated optimized Python; the receipt is [ROUTED_LOOP_INDEPENDENT_CHECK.json](ROUTED_LOOP_INDEPENDENT_CHECK.json). The parent's separate [815-check receipt](routed_loop/ROUTED_POLE_CHECKS.json) also includes its independent large-loop angular expansion. Neither check count represents a proof of the missing EFT matching.

## 5. What this does and does not settle

This resolves the first concrete kinematic gap in N00AK-r1: the routed UV pole of the **declared single derivative seed** is now explicitly evaluated, and a uniform-p² replacement is false for that off-shell graph. The prior release's caution is strengthened from an unevaluated routing objection to a computed polynomial with independent checks.

If one merely sets `r²=h²=0` and sandwiches between a free massless external lepton and a sterile spinor, the surviving cubic coefficient is `−5p² slash p/6`, versus the naive `−3p² slash p/2`. This ratio `5/9` is a narrow free-on-shell projection. **It is not a replacement beta-function coefficient 5/3**, and it cannot be used to reinstate a physical Schur cancellation.

The remaining required work is substantive:

1. Fix the derivative seed's complete Lagrangian normalization and Majorana/Fourier conventions, and its relation to the full sequential tree packet rather than one pole numerator alone.
2. Match the full off-shell result to a gauge-covariant Green basis. Three-point zero-gauge data do not fix all field-strength and covariant-derivative completion terms.
3. Apply sourced sterile, lepton and Higgs equations of motion, including every induced contact/operator term. Setting external free EOM factors to zero prematurely discards precisely the information needed to distinguish these descendants.
4. Include the other operator insertions, field/parameter variations, graphs and threshold pieces required by the chosen physical question. A UV pole is not a finite matching coefficient or an observable.
5. Only then compare `dot C₅`, `dot Y_A`, the survivor-map variation and any candidate cancellation in one consistent scheme.

No claim is made here that the eventual physical coefficient must change, that the eventual cancellation fails, or that the massless-neutrino rank theorem is overturned. The result is bounded progress in the TOE's actual unresolved loop calculation.
