# Exact finite-pole tree kernels and a conditional Schur-cancellation audit

**Ricardo Maldonado**  
**TOE-N00AK-r1 · 4 September 2026**  
**Research note and correction to the N00AK interpretation**

## Abstract

This release preserves the exact finite-pole identities for the lepton-number-violating tree kernel in a sequential Type-I seesaw calculation. It also corrects the interpretation of the N00AK dimension-seven one-loop claim. The proposed counterterm matrices satisfy a Schur-cancellation identity, but their numerical and algebraic tests do not establish that they are the complete counterterms of the stated effective field theory. In particular, extending the published dimension-five mixed bubble by a common external-momentum-squared factor has not accounted for the momentum routing of crossed field pairings. The dimension-seven coefficient, its all-odd extension, and its rationally summed loop interpretation remain conditional pending a routed off-shell calculation. The tree-kernel results and their truncation errors do not depend on that unresolved step. The public package presents this corrected note and self-contained verification code; historical N00AK notes and receipts are not included as current release evidence.

## 1. Continuity and correction

The preceding published version, verified on Zenodo on 4 September 2026, is *UNIFIED ε-Lattice Locks and Modular Closure v2026.05.16: Determinant-Valuation Geometry, CKM Ultrametric Defect, Minimal Anomaly No-Go, and Spectator Completion Candidate*, version `v2026.05.16-math-audit`, DOI [10.5281/zenodo.20241403](https://doi.org/10.5281/zenodo.20241403). N00AK belongs to a later sequential Type-I effective-field-theory research branch within the same program. This note supplies no derivation connecting the earlier ε-lattice construction to a complete ultraviolet theory.

The original N00AK report described a restricted mass-enhanced dimension-seven bubble/EOM sector as closed. That physical interpretation is withdrawn here pending the missing calculation. What survives is an exact tree-level result and an exact algebraic implication **conditional on the proposed counterterms**. Statements in original N00AK source material describing the dimension-seven loop sector, its all-odd family, or its rational resummation as established are superseded by this correction. No nonzero physical coefficient has been disproved by this audit, and no replacement coefficient is asserted.

## 2. Exact tree kernel

Use a positive Takagi mass basis for the removed sterile modes, indexed by \(D\). Their masses are \(M_d>0\), and their dimensionless lepton-flavour Yukawa columns are \(y_d\). Strip the convention-dependent overall propagator sign and retain the lepton-number-violating mass-numerator part of tree exchange. Define

\[
K_D(z)=\sum_{d\in D}\frac{y_dy_d^T/M_d}{1-z/M_d^2},
\qquad
W_n=\sum_{d\in D}\frac{y_dy_d^T}{M_d^{2n+1}}.
\]

Here \(z\) has mass dimension two, \(K_D\) has mass dimension minus one, and \(W_n\) has mass dimension \(-(2n+1)\). The transpose is not a Hermitian conjugate. This is the tree propagator expansion used in Elgaard-Clausen and Trott; their Eq. (3.2) supplies its physical setting. Their dimension-seven seed uses a documented normalization equivalent to \(\widetilde C^{(7)}=W_1/2\). The seed alone does not represent the complete sequential operator packet. [On expansions in neutrino effective field theory](https://arxiv.org/abs/1703.04415)

For \(|z|<\min_d M_d^2\), the geometric series gives

\[
K_D(z)=\sum_{n=0}^{\infty}z^nW_n,
\qquad
R_N(z)=K_D(z)-\sum_{n=0}^{N}z^nW_n
=\sum_d\frac{y_dy_d^T}{M_d}
\frac{(z/M_d^2)^{N+1}}{1-z/M_d^2}.
\]

The finite-remainder identity also holds away from the poles as an algebraic identity; convergence of the infinite local series requires the stated disk. With \(M_{\min}=\min_dM_d\), \(Y_D\) the matrix of removed columns and \(\rho=|z|/M_{\min}^2<1\), submultiplicativity yields

\[
\|R_N(z)\|_2\le
\frac{\|Y_D\|_2^2|z|^{N+1}}
{M_{\min}^{2N+3}(1-\rho)}.
\]

This is an absolute bound. Matrix cancellations can make a relative error large when the exact kernel is small.

Combine residues at identical masses, deleting a node only when its combined matrix residue is exactly zero. Let the number of remaining distinct nodes be \(c\), assumed known. If

\[
P(z)=\prod_{a=1}^{c}(1-z/M_a^2)=\sum_{k=0}^{c}p_kz^k,
\qquad p_0=1,
\]

then \(Q(z)=P(z)K_D(z)\) is a matrix polynomial of degree at most \(c-1\). Consequently,

\[
\sum_{k=0}^{c}p_kW_{n-k}=0\quad(n\ge c),
\qquad
Q(z)=\sum_{n=0}^{c-1}z^n\sum_{k=0}^{n}p_kW_{n-k}.
\]

Thus the first \(c\) moments reconstruct the rational tree kernel when its nodes are known. This is an application of classical rational generating-function algebra. Unknown nodes require additional information; nearly coincident nodes can make inverse recovery ill-conditioned. These identities do not determine a loop mixing coefficient.

## 3. Frozen benchmark and what its errors measure

The N00AK benchmark uses \((M_1,M_2,M_3)=(2,4,8)\times10^{12}\) GeV and the supplied full-precision Yukawa matrix. Define the relative error as \(\|K_{\rm approx}-K_D\|_F/\|K_D\|_F\), evaluated at \(z=M_{\rm probe}^2\). For a single nonzero pole with \(r=M_{\rm probe}/M_D<1\), the errors are exactly \(r^2\) when retaining \(W_0\), and \(r^4\) when retaining \(W_0+zW_1\).

| Removed modes | Probe | \(W_0\) error | \(W_0+zW_1\) error |
|---|---|---:|---:|
| \(N_3\) | \(N_2\) | 25% | 6.25% |
| \(N_3\) | \(N_1\) | 6.25% | 0.390625% |
| \(N_2\) | \(N_1\) | 25% | 6.25% |
| \(N_3,N_2\) | \(N_1\) | 22.1643102559% | 5.62145822425% |

The first three rows follow directly from the mass ratios. The last row is the complex-matrix benchmark recorded in the original N00AK certificate; its finite-pole reconstruction from \(W_0,W_1\) had a recorded relative residual of approximately \(7.72\times10^{-17}\). That residual is quoted as historical provenance; fresh public verification records its own numerical results. These are errors in the selected tree kernel. They are not errors in a complete neutrino-mass prediction or a sequential one-loop observable.

## 4. Conditional dimension-seven cancellation

Let \(A\) index surviving modes, with positive diagonal \(M_A\). Define

\[
S_A=Y_AM_A^{-1}Y_A^T,\qquad
P_{A,2}=Y_AM_A^2Y_A^\dagger,\qquad
\mathcal A_2=P_{A,2}W_1+W_1P_{A,2}^T,
\]

and \(\dot X=16\pi^2dX/d\ln\mu\). The N00AK proposal is

\[
\dot Y_A=-3W_1Y_A^*M_A^3,
\qquad
\dot C_{5,\mathrm{EOM}}=+3\mathcal A_2.
\]

**Assuming these expressions**, holding \(M_A\) fixed in this channel and using the symmetry of \(W_1\), direct variation gives

\[
\dot S_A=\dot Y_AM_A^{-1}Y_A^T
+Y_AM_A^{-1}\dot Y_A^T=-3\mathcal A_2,
\qquad
\dot C_{5,\mathrm{EOM}}+\dot S_A=0.
\]

This implication is exact. It also holds if every displayed factor three is replaced by the same arbitrary complex scalar \(\kappa\): the variation uses a transpose, not a conjugate transpose. Therefore, cancellation alone cannot establish that the field-theory calculation yields three, or establish the counterterm ansatz. Tests of complex-matrix orientation, signs, covariance and precision remain useful checks of the implementation within that ansatz.

For \(N_3\) removed and \(N_1,N_2\) active, the original benchmark recorded \(\|W_1\|_F=1.0912763574195238\times10^{-41}\) GeV\(^{-3}\) and \(\|\mathcal A_2\|_F=5.564319499800995\times10^{-19}\) GeV\(^{-1}\). With the proposed factor three, each cancelling component had norm approximately \(1.6693\times10^{-18}\) GeV\(^{-1}\); the original receipt recorded a relative residual of \(1.54\times10^{-16}\). These are historical values of the ansatz, separately checked by the public verification code. Nonzero components and a small residual demonstrate numerical consistency, not a completed loop derivation.

For the higher moments, an assumed family

\[
\dot y_a^{[n]}=-\kappa M_a^{2n+1}W_ny_a^*,
\qquad
\dot C_5^{[n]}=\kappa\big(P_{A,2n}W_n+W_nP_{A,2n}^T\big)
\]

likewise cancels under variation of \(S_A\), where \(P_{A,2n}=\sum_aM_a^{2n}y_ay_a^\dagger\). Below all removed poles this assumed family sums algebraically to \(\dot y_a=-\kappa M_aK_D(M_a^2)y_a^*\), with the corresponding opposite contact term. The exact rational tree identity is established; identifying this summed expression with the actual loop response is conditional. Summing an ansatz does not establish the missing momentum-dependent counterterms.

## 5. Why the loop interpretation remains open

Zhang's Appendix A gives a dimension-five redundant mixed lepton–Higgs–sterile bubble. Its vertex contains direct and crossed SU(2) field pairings; see Eq. (35) in the HTML rendering of arXiv:2405.18017v3. N00AK extended the resulting dimension-five pole by multiplying the whole answer by an external \(p^2\), based on differentiation of the composite \(J_\alpha=\epsilon_{ij}L_\alpha^iH^j\). [Threshold Effects on the Massless Neutrino in the Canonical Seesaw Mechanism, Appendix A](https://arxiv.org/html/2405.18017v3#A1)

The obstruction follows from retaining the derivative in each Wick pairing. With all-incoming lepton momenta \(r_a,r_b\) and Higgs momenta \(h_c,h_d\), the derivative vertex has pairing-dependent structure, up to a common Fourier sign and normalization,

\[
\epsilon^{ad}\epsilon^{bc}(r_a+h_d)^2
+\epsilon^{ac}\epsilon^{bd}(r_a+h_c)^2.
\]

If \(a,c\) are external while \(b,d\) are internal, the first composite momentum contains loop momentum and the second is fixed externally. Contracting with the Yukawa SU(2) tensor gives respective weights one and two, acting on different momentum polynomials. Thus extracting one common \(3p^2\) is unsupported. This is an inference from the derivative operator and its Wick pairings, not a quoted dimension-seven result from Zhang. The momentum polynomials must be retained before the tensor integrals are combined.

This identifies a missing derivation, not a computed replacement answer. The coefficient might survive a complete calculation or change after the required terms are assembled. Neither outcome is established here. Graph-connectivity counts and sourced-equation-of-motion algebra do not fill that gap. No complete dimension-seven anomalous-dimension matrix, finite hard threshold, gauge-independent sequential observable, or experimentally tested prediction is supplied by this note.

## 6. Next calculation: a routed off-shell bubble

The immediate research target is a reproducible off-shell calculation of the one-insertion \(W_1Y_A^*\) mixed bubble and its correlated local descendants:

1. Fix the Lagrangian, field and Majorana conventions, regulator, subtraction scheme, Fourier signs and Green-operator basis. Record the precise map to the dimension-five reference normalization.
2. Assign independent external lepton and Higgs momenta \(p_L,p_H\), with total sterile momentum \(p_N=-(p_L+p_H)\) when all are incoming. Enumerate every direct and crossed SU(2) pairing and write its differentiated-composite momentum polynomial before integration.
3. Retain off-shell kinematics and reduce the loop tensor integrals for each routing. Separate ultraviolet poles from possible infrared singularities with a documented prescription. Reproduce the published dimension-five bubble first as a normalization check.
4. Decompose the dimension-seven ultraviolet pole into all independent local momentum structures. Match a complete Green-basis counterterm list for this sector before applying equations of motion or integrations by parts.
5. Apply the sourced survivor equations of motion consistently, including contact terms and any required mass or field variations. Compare the resulting \(\dot C_5\), \(\dot Y_A\) and survivor-map variation; test the proposed coefficient only after these steps.
6. Reproduce the routed result with a separate implementation and check allowed loop-momentum shifts and equivalent basis reductions. Archive the integrands, reductions and conventions, including any noncancelling terms.

This calculation has not been completed in N00AK-r1. A later full sequential prediction also requires the remaining operator classes, correlated lepton-number-conserving and time-ordered terms, threshold matching and running to an observable. The broader unification program therefore remains a research objective.

## Reproducibility and AI assistance

During preparation, a 4 September 2026 rerun of the original archived code passed its core, separate math, separate EOM-tower, exact-rational finite-pole, mutation and archive-contract checks. Normal and optimized runs produced identical outputs within that runtime. Several regenerated numerical receipts differ byte-for-byte from the historical receipts; no bitwise reproduction of every historical output is claimed. Nested-parent replay and the full off-shell loop calculation were not rerun or completed by these checks. These preparation checks are distinct from the new public code's verification, which records its own environment and results.

The public package includes the source matrix and self-contained code with explicit conditional labels. It omits the historical notes and receipts, whose selected values above serve as provenance. Historical labels such as “independent” or “clean-room” describe project test implementations; they do not document external expert review. This correction separates checks of supplied algebra from derivation of physical counterterms. Checksum agreement establishes file identity, not scientific validity.

This research note and release review were prepared with AI assistance using OpenAI Codex/GPT-6 Astra. AI-assisted drafting, code review and parallel checks are disclosed as such and are not represented as independent human peer review. Ricardo Maldonado is the named author of this research release.
