# Independent review of the routed dimension-seven seed bubble

**Date:** 2026-09-06. **Verdict:** the pairing-dependent momenta, relative Wick weights and stripped ultraviolet tensor are correct for the declared single derivative seed. A complete amplitude must explicitly state its common Fourier/operator sign and normalization. This does not establish a physical dimension-seven beta function or a replacement for the historical coefficient three.

Reviewed: `verify_routed_pole.py`, the N00AK-r1 technical note §5, and the primary dimension-five example in [Zhang, arXiv:2405.18017v3, Appendix A, Eq. (35)](https://arxiv.org/html/2405.18017v3#A1). The proof checks below were performed independently of the script’s numerical tests; the script was read rather than merely rerun. This is a separate AI-assisted derivation review, not external human peer review.

## 1. Seed and combinatorics

Take a symmetric flavour matrix `W1=W1^T` and the singlet Weyl composite

\[
J_\alpha=\epsilon_{ij}L^i_\alpha H^j.
\]

Spinor contraction is understood. Use a seed proportional to

\[
\frac12 J^T W_1\Box J+\text{h.c.},
\]

with its eventual common sign declared as in §3 below. In the four-point vertex let lepton legs carry indices `a,b` and momenta `r_a,r_b`, and Higgs legs `c,d` and momenta `h_c,h_d`; every momentum in this paragraph is outgoing.

The two inequivalent Higgs pairings give the polynomial

\[
\epsilon^{ac}\epsilon^{bd}(r_a+h_c)^2
+\epsilon^{ad}\epsilon^{bc}(r_a+h_d)^2,
\]

up to the common normalization and Fourier sign. Exchanging the two fermionic composite slots gives the same bilinear after the spinor antisymmetry and Grassmann exchange are included; the conventional `1/2` accounts for these identical-slot assignments. The symmetry of `W1` is essential to using this form.

Momentum conservation makes each complementary composite momentum the negative of the displayed one. Consequently, placing `Box` on the other composite reproduces the same squared momentum. There is no extra relative factor from that exchange. The cross term in `(r+h)^2` already includes the factor two from differentiating the product; it must not be added again.

Contraction with the Yukawa tensor `epsilon^(bd)` gives

\[
\sum_{b,d}\epsilon^{ac}\epsilon^{bd}\epsilon^{bd}
=2\epsilon^{ac},
\qquad
\sum_{b,d}\epsilon^{ad}\epsilon^{bc}\epsilon^{bd}
=\epsilon^{ac}.
\]

For example, at `a=1,c=2` the first sum has two nonzero equal contributions, while the second has only `b=1,d=2`. Antisymmetric external indices reproduce the same relative `2:1` weights with their common epsilon sign. These are the dimension-five reference’s two tensor contractions with the derivative momenta retained. No new relative Wick factor or separate diagram symmetry factor was found: the internal lepton and Higgs are distinct fields.

## 2. Routing check at both vertices

Let the sterile momentum `p` enter the bubble, and external lepton and Higgs momenta `r,h` leave it, so `p=r+h`. At the derivative seed take

\[
r_a=r,\quad r_b=-\ell,\quad h_c=h,\quad h_d=\ell-p.
\]

Their sum is zero. The composite momenta are

\[
r_a+h_c=r+h=p,
\qquad
r_a+h_d=r+\ell-p=\ell-h.
\]

The complementary momenta are `-p` and `-(ell-h)`, respectively. At the Yukawa vertex the outgoing internal momenta are `ell` and `p-ell`, while the outgoing sterile momentum is `-p`; conservation holds there too.

The two propagator denominators are therefore

\[
D=\ell^2[(\ell-p)^2-m_H^2],
\]

with the usual `i0` prescriptions. With a stripped internal numerator `ell_mu`, the integrand is

\[
\boxed{
\frac{\ell_\mu[2p^2+(\ell-h)^2]}{
\ell^2[(\ell-p)^2-m_H^2]}.
}
\]

The direct term has weight two and fixed `p^2`; the crossed term has weight one and the loop-dependent `(ell-h)^2`. Replacing their sum by `3p^2` is not a permissible routing identity.

## 3. Required common sign and normalization

Zhang’s Eq. (35) uses an internal variable `k`, a numerator proportional to `k_slash`, and the scalar denominator `(k+p)^2-m^2`. Its dimension-five pole is proportional to `-3 p_slash/2`. Taking `ell=-k` produces the present denominator but introduces an overall minus in that numerator. The positive stripped integral here therefore gives `+3 p_mu/2` at dimension five; the full diagram still needs that common routing/Feynman-rule sign.

For Fourier convention `exp(-i q.x)` and Minkowski signature `(+---)`, `Box` acts as `-q^2`. Hence:

- If the reduced dimension-seven vertex is **defined** by `C5 -> +W1 q^2` relative to the reference vertex, its full amplitude has the additional overall minus from `k=-ell`.
- For a literal **positive** seed `+(1/2) J^T W1 Box J` with the same dimension-five operator normalization, the Fourier minus and the routing minus cancel.
- Changing the seed’s overall sign or using a coefficient such as `Ctilde7=W1/2` changes the common map; it does not change the relative momentum polynomial or six tensor coefficients.

Thus `i/(16 pi^2 epsilon)` alone is not a complete specification of the physical amplitude. The report should state whether the displayed object is the stripped integral, a reduced vertex with `+q^2`, or the literal action convention. It should not silently identify all three. Nothing in this review fixes a convention-independent sign for a later beta function.

## 4. Independent ultraviolet reduction

Strip `i/(16 pi^2 epsilon)` in dimensional regularization with `d=4-2epsilon`. Work first at generic off-shell momentum and nonexceptional mass, so ultraviolet poles are not erased by setting a scaleless integral to zero.

Feynman parametrization and the shift `ell=t+x p` give

\[
\Delta(x)=x m_H^2-x(1-x)p^2.
\]

The primitive ultraviolet residues are

\[
\left[\int\frac{\ell_\mu}{D}\right]_{\rm UV}
=\frac12p_\mu,
\]

\[
\left[\int\frac{\ell_\mu\ell_\nu}{D}\right]_{\rm UV}
=\left(\frac{m_H^2}{4}-\frac{p^2}{12}\right)g_{\mu\nu}
+\frac13p_\mu p_\nu.
\]

To see the metric coefficient independently, the shifted numerator contributes `g_mu_nu/d` times the pole of `t^2/(t^2-Delta)^2`, namely `2 Delta`; at the pole this is `Delta/2`. Integrating `Delta/2` over `x` gives `m_H^2/4-p^2/12`. The longitudinal coefficient is `integral_0^1 x^2 dx=1/3`.

The cubic trace is best reduced by cancelling its massless denominator before integration:

\[
\left[\int\frac{\ell_\mu\ell^2}{D}\right]_{\rm UV}
=\left[\int\frac{\ell_\mu}{(\ell-p)^2-m_H^2}\right]_{\rm UV}
=m_H^2p_\mu.
\]

Expanding the seed polynomial yields

\[
T_\mu=
\left[\int\frac{
\ell_\mu\ell^2-2\ell_\mu(\ell\cdot h)
+(h^2+2p^2)\ell_\mu}{D}\right]_{\rm UV}.
\]

Substitution gives

\[
\boxed{
T_\mu=p_\mu\left[p^2+m_H^2-\frac23p\cdot h+\frac12h^2\right]
+h_\mu\left[\frac16p^2-\frac12m_H^2\right].
}
\]

This independently confirms the script’s coefficient vector

\[
(1,-2/3,1/2,1/6,1,-1/2)
\]

in its stated basis. All terms have momentum degree three, as required for this seed. With `m_H=0`, `p=(3,0,0,0)` and `h=(1,1,0,0)`, the result is `(45/2,3/2,0,0)`, whereas the naive common-factor answer is `(81/2,0,0,0)`.

The source script’s simultaneous reversal of `p,h` checks odd parity. That test should not be described as a general loop-momentum-shift test. The analytic reduction above does use legitimate shifts in dimensional regularization; a separately encoded shift comparison would be additional verification.

## 5. Free light-leg projection and its limit

If `r^2=h^2=m_H^2=0`, `p=r+h` implies `p.h=p^2/2`. For an external free massless lepton, `ubar(r) r_slash=0`, so inside this matrix element `ubar(r) h_slash=ubar(r) p_slash`. Thus

\[
\overline u(r)\gamma^\mu T_\mu
=\frac56p^2\overline u(r)\not p.
\]

The coefficient `5/6` and its ratio `5/9` to the naive `3/2` follow. This is a projection of this one seed diagram only. It is not an operator identity and is not permission to replace the historical physical coefficient three by `5/3`, `5/6`, or another inferred value.

## 6. Scope and attempted falsification

The following potential failures were checked and did not invalidate the stripped result: momentum conservation at both vertices; complementary-composite derivative placement; identical-fermion and Higgs pairings; relative SU(2) factors; denominator mass assignment; polynomial degree; the shifted tadpole; and the independent `h_mu` structure that rules out a common `p^2` factor.

`J=epsilon LH` is a gauge singlet with zero hypercharge. Therefore a plain derivative on the whole composite is appropriate; an omitted gauge seagull is not automatically required for this particular `W1 Y* g^0` calculation. This fact does not complete the gauge-dependent and gauge-coupled sectors of a full anomalous-dimension calculation.

Still required for physical claims are: an explicit action-to-vertex normalization; the full relevant Green-operator basis; all operators and diagrams at the claimed order; sourced lepton, Higgs and sterile equations of motion; integration-by-parts and field-redefinition bookkeeping; correlated contact and parameter variations; and conversion from pole counterterms to physical running in a declared scheme. Free on-shell equations cannot substitute for those sourced equations.

**Conclusion:** accept the stated tensor as the UV residue of the explicitly routed derivative seed. Do not promote it to complete dimension-seven matching, a full EOM cancellation, finite threshold correction, or observable prediction.
