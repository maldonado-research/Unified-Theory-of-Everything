# The subthreshold cusp of the polarized four-lepton amplitude

The complete leading and pure-Yukawa next-order cuts fix the first two nonanalytic terms of the $J=0$ amplitude on the real axis below $s=4m^2$. Set $m=1$, $\eta=\sqrt{4/s-1}>0$, and strip $\gamma=y_e^2(g^*)^2$ with $B_{34}^*B_{12}=s$. At $\mu=m$, in each coupling sector the real reduced amplitude has the local form
$$
 \mathcal M_0(s)/\gamma=R(-\eta^2)-\kappa_1\eta+\kappa_3\eta^3+O(\eta^5),
$$
where $R$ is regular in $1-4/s$ near zero. The coefficients of $R$, including the finite two-loop threshold value, require dispersive information beyond these local cuts. The $J=0$ convention is the uniform external-cosine mean. [Complete cut](ABSORPTIVE.json) · [Cusp receipt](CUSP.json)

| Coupling factor in $\mathcal M_0/\gamma$ | $\kappa_1$ | $\kappa_3$ |
|---|---:|---:|
| $1$ at leading order | $-5.20101351102\times10^{-5}$ | $-1.15377808597\times10^{-4}$ |
| $E=|y_e|^2$ | $-1.10563778623\times10^{-5}$ | $-2.45292152906\times10^{-5}$ |
| $G=|g|^2$ | $+1.31099848240\times10^{-6}$ | $+2.29468997904\times10^{-6}$ |

## Expanding the physical cut

Above threshold let $\beta=\sqrt{1-4/s}$. The entries in the table are the coefficients of $\kappa_1\beta+\kappa_3\beta^3+O(\beta^5)$ in the Hermitian discontinuity divided by $\gamma$. The common weak phase and the real relative residues identify this with the imaginary part of the reduced upper-rim amplitude.

Write $q=s/2-1-(s\beta/2)x$ for the positive spacelike transfer on the two-scalar cut. Every production kernel $K$ is analytic in $q$ near 1. Since $q=1-2x\beta+2\beta^2+O(\beta^3)$, its uniform mean has the expansion
$$
 \langle K(q)\rangle=K(1)+\beta^2\mathcal D K+O(\beta^4),\qquad
 \mathcal D K=2K'(1)+\tfrac23 K''(1).
$$
For any two single-transfer kernels $U,V$, sewing with $s\beta/(8\pi)$ gives
$$
 \kappa_1[U,V]=\frac{U(1)V(1)}{2\pi},\qquad
 \kappa_3[U,V]=\kappa_1[U,V]+\frac{(\mathcal D U)V(1)+U(1)\mathcal D V}{2\pi}.
$$
At leading order $U=F_e$, $V=F_\mu$. At the next order replace the tree in the selected flavor block by its sterile and proper-vertex corrections and add the combined external residue times the leading pair. All factors here adopt the scalar pole mass, four positive $\overline{\rm MS}$ sterile masses $(10,20,1,4)$ and full-theory $\overline{\rm MS}$ Yukawas at $\mu=m$. [Sterile correction](STERILE.md) · [Vertex correction](VERTEX.md) · [External residues](EXTERNAL.md)

The vertex derivatives can be taken before numerical integration. For the crossed-transfer density $\rho_b(u)$ and squared sterile mass $b$, the coefficient of $(q-1)^n$ in the spacelike triangle is $(-1)^n\int_1^\infty\rho_b(u)\,du/(u+1)^{n+1}$. All three moments converge absolutely. For $b=1$ the undifferentiated triangle from the independent Feynman parameters reduces to $\int_0^1(1-x)\operatorname{atanh}(x)\,dx/x=\pi^2/8-\ln2$. [Crossed-transfer spectral density](INDEPENDENT.md)

The coincident three-particle channel contains a scalar, a massless muon and $N_1$. Its leading cut equals $F_e(1)F_\mu(1)(s-4)^{3/2}/(96\pi^3)$ per $G$. The exchange correction starts at $(s-4)^{5/2}$. Using $s-4=4\beta^2/(1-\beta^2)$ gives
$$
 \kappa_{3,G}^{(3)}=\frac{F_e(1)F_\mu(1)}{12\pi^3}
 =-8.78288075127\times10^{-7}.
$$
It adds to the two-scalar value $+3.17297805417\times10^{-6}$ and leaves $\kappa_{1,G}$ fixed. This coincidence uses the benchmark's renormalized light sterile mass $M_1=m$ in the fixed-order cut. A pole-mass conversion or resummation that resolves a splitting of the thresholds calls for a new threshold expansion. [Three-particle threshold coefficient](../../checkpoint_2026_11_19/yukawa_continuum/THREE_BODY.md) · [Complete exclusive cut](../../checkpoint_2026_11_19/yukawa_continuum/FULL_THREE_BODY.md)

For $J=2$, each production moment $\langle P_2(x)K(q)\rangle$ starts at $(4/15)K''(1)\beta^2$ because $\langle P_2 x^2\rangle=2/15$. Sewing two such moments costs four powers of $\beta$. For example the leading $J=2$ moment divided by $\beta^5$ tends to $5.13300654963\times10^{-8}$. The three-particle measure starts at $\beta^3$, and its recoil momentum is uniformly $O(\beta)$, giving $O(\beta^7)$ for its $J=2$ moment. Higher even moments require higher powers of the recoil momentum, and odd moments vanish by scalar exchange for the two-body cut and by the external identical-lepton assignments for both cuts. Thus the first two cusp coefficients apply at each fixed interior scattering cosine in the same spin phase.

## Physical sheet and an independent real-box check

Put $w=1-4/s$. On the upper rim the locally nonanalytic amplitude is $i\kappa_1\sqrt w+i\kappa_3 w^{3/2}+\cdots$. Continuing from above through $\mathrm{Im}\,s>0$ gives $\sqrt w=+i\eta$ below threshold and the signs displayed in the first equation. Only the reduced amplitude has these real coefficients when $\gamma$ is complex.

One can recover the continuation from dispersion without choosing a square-root ansatz. Integrate the local cut up to $w'=b^2$ using $x=w'$. At $w=-\eta^2$, the two integrals with measure $dx/[\pi(x+\eta^2)]$ and numerators $x^{1/2},x^{3/2}$ give respectively
$$
 \frac2\pi[b-\eta\arctan(b/\eta)],\qquad
 \frac2\pi[b^3/3-b\eta^2+\eta^3\arctan(b/\eta)].
$$
Since $\arctan(b/\eta)=\pi/2-\arctan(\eta/b)$, their nonanalytic parts are $-\eta$ and $+\eta^3$. Transforming the original $ds'/(s'-s)$ kernel to $w'$ adds $dw'/(1-w')$, independent of $s$. Local subtraction terms, the part of the right-hand cut away from this threshold, and crossed-channel terms analytic near this positive $s$ enter the regular background. The three-particle endpoint integral has a smooth expansion after rescaling $U-1$ by its interval length and factoring its square-root phase space. Together with the spacelike analyticity of the two-particle production functions, this gives the half-integer series through the orders used here.

Run `python -I -B -O check_cusp.py --output CUSP.json`. The checker compares the differentiated crossed-channel integrals against the unexpanded complete cut at $\beta=0.02,0.04,0.06,0.08$. A cubic fit in $\beta^2$ to the cut divided by $\beta$ recovers $\kappa_1$ within $10^{-9}$ relative and $\kappa_3$ within $10^{-7}$ in each sector. Numerical subthreshold dispersion checks both signs with two local cutoffs. For a separate test using the full leading box, eight evaluations on the physical real axis at $\eta=0,0.0005,\ldots,0.0035$ resolve its linear and cubic terms with the smooth background fitted freely. They agree with the predicted coefficients within $3\times10^{-13}$ and $2\times10^{-10}$ relative, respectively. [Cusp checker](check_cusp.py) · [Leading full box](../../checkpoint_2026_11_19/uv_completion/DISPERSION.md)
