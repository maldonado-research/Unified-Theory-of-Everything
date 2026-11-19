# Complete charged three-particle cut in the four-singlet benchmark

At $4<s/m^2<25$ the only three-particle intermediate state in the order $\gamma |g|^2$ Hermitian cut of $\mu^-_L\mu^-_L\to e^-_Le^-_L$ is $\phi^-\mu^-N_1$. The two additional exchange routings in the muon-initiated tree increase the magnitude of its $J=0$ contribution by **9.03%, 15.33% and 20.40%** at $s/m^2=6,8,10$ relative to the scalar-line insertions. Their sign follows analytically at every angle. Here $J=0$ denotes the external-cosine mean after choosing $B_{34}^*B_{12}=s$. The two-scalar cut of the Yukawa correction and the real part of the complete two-loop amplitude are separate inputs to a four-lepton cross section. [Benchmark couplings](../uv_completion/MATCHING.md) · [Three-body insertion and tree counting](THREE_BODY.md)

## An effective kernel after the decay-angle integral

Use $m=1$ throughout. Write $k,p,n$ for the outgoing pole scalar, massless muon and $N_1$ momenta. Let $Q=p+n$, $U=Q^2$, and take the beam momenta $p_1,p_2$ along $\pm z$. Denote a positive transfer by $q_j=-(p_j-k)^2$. The positive muon mass-exchange kernel and the normal-propagator sum are
$$
 F_\mu(q)=\frac{15}{(1+q)(16+q)},\qquad
 H(q)=\frac{1}{1+q}+\frac{1}{4(16+q)}.
$$
The weights in $H$ are squared magnitudes $|Y_{\mu j}|^2/|g|^2=(1,1/4)$. The anomalous propagator uses $M_j(Y_{\mu j}^*)^2/(g^*)^2=(1,-1)$. The two structures carry the same total weak phase after sewing. The full state contribution follows by replacing each $F_\mu(q)$ in the scalar-insertion cut by
$$
 \boxed{F_{\rm full}(q,U)=F_\mu(q)+D(U)\frac{q}{U+q}H(q),\qquad
 D(U)=\frac{U\ln U}{U-1}-1.}
$$
In particular $D(U)>0$ for $U>1$, and $q>0$ throughout the physical sphere at finite $s$. The per-transfer fractional increase of the positive muon kernel is $D(U)q(13+q)/[12(U+q)]$. The electron tree is negative for both its transfers. Hence the **signed** contribution of the extra diagrams to this off-diagonal cut is strictly negative at every external angle for $s>4$.

To see completeness, at this coupling order the connected $2\to3$ tree on either side contains three Yukawa vertices. Counting their sterile slots requires an odd number of external steriles. Charge $-2$ then fixes one neutral sterile and two charged particles. Counting the scalar slots requires an odd number of external scalars, hence one charged scalar and one charged lepton. Below $s=25$ the only accessible sterile is $N_1$. On the electron side its Yukawa vertex must carry the outgoing muon. This fixes a scalar of momentum $Q$ decaying into $p,n$. On the muon side either incident vertex can emit $N_1$ as well. The other incident vertex then emits $k$. Emitting both $N_1$ and $k$ at the same vertex disconnects the tree. A disconnected splitting of an on-shell external massless lepton into a massive pair has no phase space. These observations exhaust the three-particle state and its tree routings in the specified gaugeless theory.

## Spin sum and sign of the exchange term

Take $i\ne j$, attach $n$ to incoming leg $i$, and attach $k$ to incoming leg $j$. The internal scalar denominator is $d_i=(p_i-n)^2-1=-2p_i\!\cdot n$. The sterile momentum is $r_j=p_j-k$. Let $\chi$ be a left-handed spinor with $\chi_p\chi_p^\dagger=p^0-\boldsymbol\sigma\cdot\boldsymbol p$ and $B_{ab}=\chi_a^T\epsilon\chi_b$, $\epsilon=i\sigma_2$. After stripping common external phases, the anomalous muon numerator contains $B_{12}\bar u(p)P_R v(n)$. Its spin sum against the decaying electron-side scalar includes $\sum|\bar u(p)P_R v(n)|^2=U-1$. For the exchange orientation the two chains are $(\bar u(n)P_L u_i)(\bar u(p)P_R\not r_j P_Lu_j)$. The Majorana relation $v=C\bar u^T$, with $C=\mathrm{diag}(\epsilon,-\epsilon)$ in the chiral basis, gives
$$
 \sum_{\lambda_{N_1}}
 (\bar u(n)P_Lu_i)(\bar u(p)P_Rv(n))^*=M_1 B_{ip},\qquad M_1=1.
$$
This follows from the mixed-chirality massive completeness relation $\sum u_L(n)u_R(n)^\dagger=M_1\mathbf1$. The normal numerator retains $r_j^0+\boldsymbol\sigma\cdot\boldsymbol r_j$. For the fixed external ordering the two exchange diagrams enter with opposite signs, which are absorbed by $B_{ij}$ below. One can check the relative sign directly by ordering the sterile fields as $\bar N_1\bar N_2N_3$, with subscripts here denoting the vertices. The anomalous left-left block is $(C^{-1}(\not r+M))_{LL}=-M\epsilon$. Contracting the reversed normal pair $\bar N_2 N_3$ also contributes a minus. Contracting $\bar N_1 N_3$ instead adds one Wick interchange. Thus the first exchange and anomalous spin products have a common sign before swapping the incoming labels.

Relative to the decay sewing the ratio of the missing scalar denominator $U-1$ to the decay spin trace $U-1$ is one. In the positive $F_\mu$ convention the unaveraged extra kernel for this route is $H(q_j)K_{ij}(p)/d_i$, where
$$
 K_{ij}(p)=\frac{B_{ip}}{B_{ij}}\,\chi_p^\dagger
 (r_j^0+\boldsymbol\sigma\cdot\boldsymbol r_j)\chi_j.
$$
The normalization uses $B_{ij}$ with the same ordering as the two incident legs of each route. Put $\eta_i=+1$ for the $+z$ leg and $-1$ for the $-z$ leg. The massless equation on $\chi_j$ removes the $p_j$ term in $r_j$, leaving
$$
 K_{ij}(p)=-p\cdot k-\eta_i\big[p^0k_z-k^0 p_z-i(\boldsymbol p\times\boldsymbol k)_z\big].
$$
The result is linear in $p$. In particular $K_{ij}(p_i)=0$ and $K_{ij}(Q)=-q_j$ since $Q=P-k$.

## The decay average

In the $Q$ rest frame, define $z_i=p_i\cdot Q/\sqrt U$, $\omega=(U-1)/(2\sqrt U)$, $\nu=(U+1)/(2\sqrt U)$ and $t=\hat p_i^*\cdot\hat p^*$. Then $d_i=-2z_i(\nu+\omega t)$. Averaging the spatial components perpendicular to the $p_i$ axis gives zero. The two remaining integrals are
$$
 J=\int_{-1}^1\frac{dt/2}{\nu+\omega t}
   =\frac{\ln U}{2\omega},\qquad
 J_1=\int_{-1}^1\frac{t\,dt/2}{\nu+\omega t}
   =\frac{1-\nu J}{\omega}.
$$
Restoring covariance expresses the averaged four-vector as
$$
 \left\langle\frac{p^\mu}{d_i}\right\rangle_{p^*}
 =-\frac{\omega}{2z_i}\left[(J-J_1)\frac{Q^\mu}{\sqrt U}
                               +J_1\frac{p_i^\mu}{z_i}\right].
$$
Here $\omega(J-J_1)=\sqrt U\,J-1=D(U)$ and $2p_i\cdot Q=U+q_j$. Applying the linear functional $K_{ij}$ gives $\langle K_{ij}(p)/d_i\rangle=D(U)q_j/(U+q_j)$. This establishes the boxed replacement, including its sign, before projecting the external angle.

## Angular projection, threshold and checks

Let $h=(s-1-U)/2$, $v=\sqrt{h^2-U}$, $q_\pm=h\pm vx$, and $S_{\rm full}=F_{\rm full}(q_+,U)+F_{\rm full}(q_-,U)$. With $S_e$ and the sphere axes as in the insertion calculation, the full state contribution is
$$
 A_{3}(s,c)=\int_{1}^{(\sqrt{s}-1)^2}\frac{2v\,dU}{256\pi^3 U}
 \int\frac{d\Omega_k}{4\pi}\, S_{\rm full}(x,U)\, S_e(y).
$$
Both sides of the cut carry the upper-rim external-state convention of the leading box. The formula gives the coefficient of $\gamma |g|^2$ after taking the Hermitian discontinuity. Rotational invariance still factorizes the $J=0$ sphere average into $4\bar F_{\rm full}\bar F_e$ at each $U$. For $J_b=(2v)^{-1}\ln[(b+h+v)/(b+h-v)]$, the exchange average uses
$$
 \left\langle\frac{q}{(U+q)(b+q)}\right\rangle_x
 =\frac{U J_U-b J_b}{U-b}
 \ \longrightarrow\ J_b-\frac{b}{(b+h)^2-v^2}\quad(U=b).
$$

| $s/m^2$ | Full $A_{3}(s,0)$ | Full $\langle A_3\rangle_c$ | Exchange contribution to $\langle A_3\rangle_c$ |
|---:|---:|---:|---:|
| 6 | $-4.357691833\times10^{-7}$ | $-4.356892150\times10^{-7}$ | $-3.606903594\times10^{-8}$ |
| 8 | $-1.458437323\times10^{-6}$ | $-1.457577494\times10^{-6}$ | $-1.937814101\times10^{-7}$ |
| 10 | $-2.964545555\times10^{-6}$ | $-2.961200295\times10^{-6}$ | $-5.017304897\times10^{-7}$ |

For $\Delta=s-4\downarrow0$ the exchange term has a simple coefficient. At $U=1+t$ and $q=1$, $D=t/2+O(t^2)$ and $H(1)/F_\mu(1)=7/6$. Thus the leading per-transfer fractional correction is $7t/24$. The two transfers cancel the odd terms in $v$ before the angular integral, and the threshold weight $\sqrt{\Delta-2t}\,dt$ has mean $\langle t\rangle=\Delta/5$. Consequently at each fixed $c$, and after the $J=0$ projection,
$$
 \frac{A_3}{A_{\rm scalar\ insertion}}-1
   =\frac{7}{120}\Delta+O(\Delta^2).
$$
The leading full threshold coefficient is $-1.097860094\times10^{-7}\Delta^{3/2}$.

The [replay program](check_full_threebody.py) compares the invariant-mass integral against direct recoil-momentum integration with both unintegrated spin chains reduced to the Pauli trace above. The direct route boosts the decay momenta, samples their full sphere, forms the two sterile and scalar transfers, and samples the production sphere and external angles. At $s=8$ and the five momenta in the receipt, the ratio of full to scalar-only spin overlaps is $1.257438753877467$ before the decay-angle integration. A separate MadGraph5_aMC tree calculation agrees at that point. The offline checker also verifies the absolute complex products at two asymmetric points, including their imaginary parts. Using the generated helicity amplitudes at $s=8$ gives $-1.457577493715894\times10^{-6}$ for the integrated full $J=0$ cut. [Numerical receipt](FULL_THREE_BODY.json) · [Generator, conventions and reproduction](EXACT_MG.md) · [MadGraph figures](MADGRAPH.json)
