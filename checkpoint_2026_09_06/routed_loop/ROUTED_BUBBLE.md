# TOE: a routed dimension-seven seed bubble

6 September 2026 Pacific. This checkpoint computes the missing **off-shell UV momentum polynomial for one declared derivative-seed mixed bubble**. It advances the N00AK-r1 scope correction from a routing objection to an evaluated integral. It does not complete the dimension-seven operator packet, its sourced equation-of-motion reduction, or a physical beta function.

## 1. Fixed conventions and diagram scope

Use Minkowski signature (+---), Fourier convention exp(-ipx), d=4-2 epsilon, and dimensional regularization. Strip the common logarithmic pole factor C_UV=i/(16 pi² epsilon). Keep generic external momenta during integration, with nonexceptional p², so a scaleless on-shell integral is not mistaken for an absent UV pole. The internal lepton is massless. The scalar propagator has denominator (ell-p)²-m_H²+i0; m_H² denotes its signed mass-squared parameter, not a fitted cosmological scale.

Let p enter from the sterile-neutrino leg, and let r and h leave on the lepton and Higgs legs, respectively, so p=r+h. The derivative vertex has all-outgoing momenta

| Slot | Field | Momentum |
|---|---|---|
| a | external lepton | r |
| b | internal lepton | -ell |
| c | external Higgs | h |
| d | internal Higgs | ell-p |

They sum to zero. The paired composite momenta are r+h=p and r+(ell-p)=ell-h. Their complementary composites have the opposite momenta and hence identical squares.

Define the **momentum-space seed convention** by replacing the dimension-five vertex tensor with

\[
iC_{7,\mathrm{seed}}^{\alpha\gamma}
\{\epsilon^{ad}\epsilon^{bc}(\ell-h)^2
+\epsilon^{ac}\epsilon^{bd}p^2\}P_R.
\]

The flavor tensor is symmetric. C_(7,seed) denotes the coefficient of the positive momentum-squared polynomial in this definition. If it is identified with N00AK's W1, that identification uses the declared momentum-space kernel convention K(z)=W0+zW1+.... A literal position-space Box gives -q² with our Fourier convention. Therefore an action written with a different relative Box sign must be translated first; it must not silently inherit the amplitude sign here. This checkpoint fixes a vertex convention, not every historical action sign.

Contracting the Yukawa epsilon^(bd) gives weights **one for the crossed term and two for the direct term**. The symmetry factor1/2 of a symmetric composite bilinear is compensated by exchanging the fermion slots, as for the dimension-five vertex; it does not introduce an extra relative factor between the two pairings.

The reference dimension-five diagram is Zhang, Appendix A, Eq.(35), whose final pole is proportional to -3 slash-p/2. That reference fixes the comparison of momentum flow and the standard tensor integral; it does not provide the dimension-seven answer derived below. [Primary source](https://arxiv.org/html/2405.18017v3#A1).

## 2. The evaluated tensor integral

The reduced object is

\[
J^\mu=\mathrm{Res}_{C_{UV}}
\int\frac{d^d\ell}{(2\pi)^d}
\frac{\ell^\mu[2p^2+(\ell-h)^2]}
{(\ell^2+i0)[(\ell-p)^2-m_H^2+i0]}.
\]

Feynman-parametrize and shift ell=q+xp. The denominator becomes (q²-Delta+i0)² with Delta=x m_H²-x(1-x)p². The normalized UV residues are

\[
I^\mu=\frac{p^\mu}{2},\qquad
I^{\mu\nu}=\left(\frac{m_H^2}{4}-\frac{p^2}{12}\right)g^{\mu\nu}
+\frac{p^\mu p^\nu}{3}.
\]

These follow from integrating x, x² and Delta over[0,1], and from the pole of the rotational tensor q^mu q^nu. Cancelling ell² in the cubic numerator before integrating gives a shifted massive tadpole:

\[
\mathrm{Res}_{C_{UV}}\int\frac{d^d\ell}{(2\pi)^d}
\frac{\ell^\mu\ell^2}{\ell^2[(\ell-p)^2-m_H^2]}
=m_H^2p^\mu.
\]

Therefore

\[
\boxed{J^\mu=p^\mu\left[p^2+m_H^2-\frac23p\cdot h+\frac12h^2\right]
+h^\mu\left[\frac16p^2-\frac12m_H^2\right].}
\]

Equivalently, in terms of the outgoing lepton momentum,

\[
\boxed{\not J=
\left(\frac56p^2+\frac13r^2+\frac16h^2+\frac12m_H^2\right)\not p
+\left(\frac12m_H^2-\frac16p^2\right)\not r.}
\]

The four cubic off-shell structures are p² slash-p, r² slash-p, h² slash-p and p² slash-r; there are also two Higgs-mass descendants. With our vertex convention, Zhang's loop variable k=-ell supplies an overall minus, so the amplitude pole is

\[
i\mathcal M_{7,\mathrm{seed}}^{UV}
=-\frac{i}{16\pi^2\epsilon}\epsilon^{ac}
(C_{7,\mathrm{seed}}Y_\nu^*)_{\alpha\beta}
\bar u_\ell(r)\not J P_Lu_N(p).
\]

Replacing the complete derivative vertex with the dimension-five vertex yields -C_UV(3/2)slash-p and reproduces the reference normalization. No gamma5 trace or anomalous chiral trace occurs in this open fermion chain.

## 3. An independent derivation

The second route expands the integrand at large loop momentum and keeps only terms homogeneous of degree-4. The needed angular averages at d=4 are g/4 and the three symmetric metric products/24. The cubic numerator contributes m_H² p; the quadratic term contributes -(2/3)(p.h)p+(p²/6-m_H²/2)h; and the linear term contributes (p²+h²/2)p. Their sum reproduces J.

This operation extracts the UV residue. It is not an instruction to integrate an unregulated scaleless expansion and set every term to zero. Equivalently, an auxiliary IR regulator may be used to separate the UV pole before removing it. The original nonexceptional integral and the massive-tadpole calculation keep UV and IR roles distinct.

The separate reviewer used Zhang's original k-plus routing and independently derived the rank-one, rank-two and rank-three tensor poles. Its result equals -J after k=-ell, including the Higgs-mass terms. The seed-vertex reviewer separately checked all momentum slots and SU(2) weights. See the review files alongside this note.

## 4. A concrete falsification of the old shortcut

The shortcut would give J_naive=(3/2)p²p. Set p=(3,0,0,0), h=(1,1,0,0), m_H²=0, in arbitrary consistent momentum units. The actual result is

\[
J=(45/2,3/2,0,0),\qquad J_{naive}=(81/2,0,0,0).
\]

The actual vector is not even parallel to p. No change of a single common coefficient can reproduce its full off-shell dependence.

As a separate diagnostic, after the integration only, take free massless light legs r²=h²=0 and bar-u(r)slash-r=0. Then the surviving coefficient is(5/6)p² slash-p, versus(3/2)p² slash-p in the shortcut. Their ratio 5/9 is a **three-point kinematic projection**, not a replacement beta-function coefficient. Substituting5/3 for the old proposed3 in the N00AK code would skip the unresolved operator analysis and is not justified by this result.

## 5. What remains for a physical coefficient

This diagram is the no-external-gauge-field component of one derivative seed. Before extracting the running of the survivor Yukawa or Weinberg coefficient, one must:

1. Match the polynomial to a declared dimension-seven Green basis, with full action and Majorana normalization fixed.
2. Supply its gauge completion. Three-point momentum data alone do not fix every independent operator containing gauge fields.
3. Reduce derivatives on N, L and H with their **sourced** equations of motion. Squared light-field derivatives generate interactions and field-strength terms; using only their free on-shell equations loses descendants.
4. Include all correlated contact, field and mass terms relevant at the same order, and compare the complete survivor matching combination.
5. Compute the finite threshold and remaining operator classes before claiming an observable.

The correct next target is now an explicit finite polynomial and its operator image, rather than an unspecified routing question. The old N00AK-r1 correction remains valid. Exact tree-kernel identities are unchanged; no completed dimension-seven physical cancellation, new neutrino mass prediction, or unified theory is asserted.

## Verification and status

`python3 -I -S -B -O verify_routed_pole.py` passes 815 exact rational checks: SU(2) contractions, Feynman-parameter residues,400 comparisons of tensor reduction with the separate UV angular expansion, momentum-reversal controls, and the nonproportional counterexample. These controls validate the displayed algebra under the declared integral identities. They are not independent expert peer review or proof-assistant formalization. The checks use explicit failures and remain active under optimized Python.
