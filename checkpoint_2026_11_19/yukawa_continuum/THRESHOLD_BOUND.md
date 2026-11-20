# The full exclusive cut at the lightest three-particle threshold

In the benchmark with $M/m=10$, the leading onset of the complete $\phi^-\mu^-N_1$ intermediate-state cut is fixed by scalar decay on both sides. Write $C_{DD}$ for the signed decay–decay contribution and $C_{\rm full}$ for sewing the complete trees in this exclusive channel. With $\Delta=s/m^2-4$, each external scattering angle obeys
$$
 |C_{\rm full}-C_{DD}| < \tfrac25\Delta\, |C_{DD}|,
 \qquad 0<\Delta\le 1/100 .
$$
The bound also applies to the angular mean. Since $C_{DD}$ starts as $\Delta^{3/2}$ with a nonzero negative reduced coefficient, the additional exchange terms are $O(\Delta^{5/2})$ or smaller. This argument uses norms throughout and is independent of the relative exchange sign. [Threshold coefficient and three-body phase space](THREE_BODY.json) · [Exact arithmetic certificate](THRESHOLD_BOUND.json)

## Which amplitudes can interfere

Use $m=1$, fix the spin phase $B_{34}^*B_{12}=s$, and strip the common coupling phase and magnitudes. The three-body state contains a charged scalar of momentum $k$, a Majorana fermion of momentum $n$ and a massless muon of momentum $p$. Set $k^2=n^2=1$ and $p^0=\omega$. There are three external lepton slots at the three Yukawa vertices. Charge fixes the outgoing scalar to one of the two incoming-lepton slots. In the **decay** tree, the external $N_1$ meets the outgoing muon, and both initial leptons exchange a singlet with an anomalous mass numerator. On the electron side this exhausts the amplitude. On the muon side the external $N_1$ can instead meet the other incoming lepton. The remaining two vertices then exchange a singlet with a chirality-preserving momentum numerator. There are two such **exchange** terms corresponding to the two initial leptons. Assigning $k$ and $n$ to the same vertex would disconnect the diagram. The direct muon tree has coupling $g\sum_j M_j(Y_{\mu j}^*)^2$ with the propagator for each summand understood. The exchange tree has coupling $g^*\sum_j |Y_{\mu j}|^2$ with its own propagators. Both carry the common factor $g^*|g|^2$ in this texture. [Flavor texture and charged-tree phases](../uv_completion/MATCHING.md)

Let $E=\sqrt{s}/2$ be either incoming lepton energy and set $Q_j=-(p_j-k)^2$. The dimensionless anomalous kernels are $F_\mu(Q)=15/[(1+Q)(16+Q)]$ and $F_e(Q)=-30Q/[(100+Q)(400+Q)]$. Each decay tree is a scalar factor times the same outgoing two-fermion spin vector $V$. For the vector norm sum over both $N_1$ spins with the outgoing muon helicity fixed. In the original Weyl normalization $\|V\|^2=2n\cdot p=(n+p)^2-1$. Including the scalar decay propagator and $|B_{12}|=\sqrt{s}$ gives the exact stripped muon decay norm
$$
 \|D_\mu\|=\frac{\sqrt{s}[F_\mu(Q_1)+F_\mu(Q_2)]}{\sqrt{2n\cdot p}} .
$$
The electron decay tree uses the two electron-axis transfers in $F_e$. Wherever the four $Q_j$ are positive, $F_\mu>0$ and $F_e<0$. Hence the decay–decay product after phase removal has constant negative sign and magnitude $\|D_\mu\|\,\|D_e\|$ pointwise. [Charged-scalar spin sum](SPECTRUM.md)

## Uniform bound on the exchange terms

The invariant mass of the system recoiling against either massive particle is at least one. Thus $E_k,E_n\le E<1.002$, $|\boldsymbol k|,|\boldsymbol n|\le\sqrt{\Delta}/2\le0.05$. The invariant mass recoiling against the massless muon is at least two, giving $\omega\le\Delta/(2\sqrt{s})\le\Delta/4$. These inequalities hold at every decay angle. They imply
$$
 .89<Q_j<1.12,\quad 2n\cdot p\le2(1.052)\omega,\quad
 |(p_i-n)^2-1|=2p_i\cdot n\ge1.9 .
$$
For the exchanged singlet let $q_j=p_j-k$. Both Weyl matrices $q_j\cdot\sigma$ and $q_j\cdot\bar\sigma$ have operator norm bounded by $|q_j^0|+|\boldsymbol q_j|\le(1.002-1)+1.002+0.05=1.054$. The normal-propagator sum after stripping $|g|^2$ has magnitude
$$
 H_\mu(Q)=\frac1{1+Q}+\frac1{4(16+Q)}.
$$
In either exchange, the vertex with the external Majorana spin has norm at most $\sqrt{2E(E_n+|\boldsymbol n|)}$. This follows from $\|u_i\|^2=2E$ for a massless Weyl spinor and the largest eigenvalue $E_n+|\boldsymbol n|$ of the massive Weyl completeness matrix. The other fermion chain has two massless external spinors and one momentum numerator, and is bounded by $2\sqrt{E\omega}\,\|q_j\cdot\sigma\|$. The antisymmetric spin matrix has unit operator norm. Multiplying these bounds with both propagator factors gives
$$
 \|X_i\|\le
 \sqrt{2E(E_n+|\boldsymbol n|)}\,2\sqrt{E\omega}\,
 \frac{1.054\ H_\mu(.89)}{1.9}<.89\sqrt{\omega},\qquad
 \|D_\mu\|\ge\frac{2[2F_\mu(1.12)]}{\sqrt{2(1.052)\omega}}>
 \frac{1.13}{\sqrt{\omega}} .
$$
The two exchange terms together therefore satisfy $\|X_1+X_2\|/\|D_\mu\| < 2(.89)\omega/1.13 < (2/5)\Delta$ at interior phase-space points. Cauchy–Schwarz bounds their overlap with $D_e$. The pointwise decay–decay product has a fixed phase and sign, so integration against the positive phase-space measure proves the stated absolute bound relative to the integrated cut. The same sign also permits angular averaging. The certificate checks the positive rational inequalities after squaring, including the beam energy and recoil bounds. The interval applies only close to $s=4m^2$. [Reproduction script](check_threshold.py)
