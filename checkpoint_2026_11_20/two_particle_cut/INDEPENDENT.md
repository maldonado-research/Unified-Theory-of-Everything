# Reconstructing the proper vertex from crossed-channel scattering

Sewing two physical tree amplitudes in the crossed lepton–scalar channel fixes the sign and the finite part of the proper Yukawa vertex. With the scalar pole mass as unit and all four singlets active in $\overline{\rm MS}$ at $\mu=m$, the resulting complete scalar-pair correction to the $J=0$ discontinuity at $s=8m^2$ is $\gamma[-2.65731240\,10^{-5}|y|^2+3.67404321\,10^{-6}|g|^2]$. Here $\gamma=y^2(g^*)^2$ and $B_{34}^*B_{12}=s$. [Independent computation](check_independent.py) · [Numerical record](INDEPENDENT.json)

## Tree sewing in the sterile transfer

Consider one lepton-flavor block. Write $a=m^2$, $B_j=M_j^2$, $w_j=M_j c_j^2$, where $c_j$ is its Yukawa or its conjugate for the chosen external orientation. Strip the ordered external spin determinant and its common minus sign. The anomalous and normal tree kernels per transfer are
$$
 F(z)=\sum_j\frac{w_j}{B_j-z},\qquad
 G(z)=\sum_i\frac{|c_i|^2}{B_i-z}.
$$
Crossing to $\ell^-\phi^+\to \ell^+\phi^-$ gives two anomalous assignments with the same sign multiplying the same ordered spin determinant. The normal elastic $\ell^-\phi^+\to \ell^-\phi^+$ tree carries $+G(z)$ multiplying the slash of the total momentum. In the cut product, a direct anomalous assignment dresses the sterile propagator. Pairing the elastic tree with a crossed anomalous exchange instead gives the proper vertex on one endpoint. The direct contraction sums both weak components. The crossed assignment ties the internal component to the charged external one. The conjugate intermediate state gives the other endpoint, with the same kinematic and weak factor in a single-flavor block. [Ordered-tree conventions](../../checkpoint_2026_11_19/uv_completion/MATCHING.md) · [Sterile insertion](STERILE.md)

Let $u=q^2>a$ be the crossed energy squared, $q=p+k$, $p^2=0$, $k^2=a$, and $l$ the active intermediate momentum. In the rest frame of $q$ put $r=(1-\hat{\boldsymbol p}\cdot\hat{\boldsymbol l})/2$. The external and intermediate lepton energies are both $E=(u-a)/(2\sqrt u)$. With the outgoing lepton along $+z$, its left-handed spinor is proportional to $(0,1)^T$. For incoming spinor $v$ and active spinor $v_l$, the unitarity product divided by the ordered spin determinant is
$$
 \sqrt u\,\frac{v_p^T\epsilon v_l v_l^\dagger v}{v_p^T\epsilon v}.
$$
The massless projector is $v_l v_l^\dagger=E(1-\hat{\boldsymbol l}\cdot\boldsymbol\sigma)$. Its azimuthal mean retains the upper diagonal component in this ratio, giving $2\sqrt u E r=(u-a)r$. The virtual anomalous transfer is
$$
 t_r=(l-k)^2=2a-u+\frac{(u-a)^2}{u}\,r,\qquad
 B_j-t_r\ge B_j-a^2/u>0 .
$$
The angular mean $\int_0^1dr$ and two-body velocity $1-a/u$ then give the following kernel in a convention where $16\pi^2$ is extracted:
$$
 \rho_j(u)=\frac{\mathrm{Im}\,I_j(u+i0)}{\pi}
   =\left(1-\frac au\right)^2\int_0^1\frac{r\,dr}{B_j-t_r},\qquad
 I_j(z)=\int_a^\infty\frac{\rho_j(u)\,du}{u-z}.
$$
Indeed the ordered phase space integrated over angles is $(1-a/u)/(8\pi)$ and the discontinuity divided by $2i$ supplies another half. Dividing the numerator by $u$ leaves the second velocity. Both cut particles are distinct. Near $u=a$ the density is $O((u-a)^2)$ for $B_j>a$ and $O(u-a)$ for $B_j=a$. At infinity it is $O(\ln u/u)$, giving a convergent integral. For each individual mass the density is positive. For $Q_2>Q_1\ge0$ this gives $(a+Q_1)/(a+Q_2)<I_j(-Q_2)/I_j(-Q_1)<1$ by averaging the ratio $(u+Q_1)/(u+Q_2)$ over the positive measure.

Restoring the numerator $z$ and the two possible endpoints yields the signed result
$$
 \delta F_{\rm proper}(z)=\frac{2z}{16\pi^2}G(z)\sum_j w_j I_j(z).
$$
Thus for a real single-column Yukawa the spacelike correction is negative relative to its positive tree kernel. In the massless-scalar limit the pole endpoint sign also agrees with the factor $1-M b K^T M$ of Buchmüller and Plümacher, with $b_j=-I_j/(16\pi^2)$. [Buchmüller and Plümacher, Eqs. (15), (17), (18), (57)](https://arxiv.org/abs/hep-ph/9710460)

## Integral checks and complete scalar-pair projection

Combining the active denominator $l^2$, scalar denominator $(q-l)^2-a$ and anomalous denominator $(l-k)^2-B_j$ with parameters $x,b,c$ satisfying $x+b+c=1$ shifts $l$ by $bq+ck$. The open chain reduces to $z(p\cdot l)/(p\cdot q)$ after the vector loop integral, giving numerator $z(b+c)$. The simplex denominator is $b(a-xz)+c(B_j-ax)$. Its only singularities in $z$ start at $a$ for $B_j\ge a$. At spacelike infinity the integrand tends to zero under an integrable bound from any fixed spacelike $z$, which fixes the unsubtracted dispersion relation. Integrating the relative scalar/sterile parameter gives
$$
 I_j(z)=\int_0^1(1-x)\,
 \frac{\ln[(B_j-ax)/(a-xz)]}{B_j-ax-(a-xz)}\,dx.
$$
At real $u>a$ this gives $\mathrm{Im}\,I_j/\pi=\int_{a/u}^1 (1-x)dx/[B_j-a+x(u-a)]$. Both the timelike density and the spacelike or subthreshold value agree numerically with the spin-summed reconstruction for all four benchmark masses. A separate calculation disperses in the scalar invariant while holding the sterile transfer fixed. [Scalar-invariant calculation](VERTEX.md)

For the complete scalar-pair check the standalone program integrates the finite sterile bubble and both external pole residues. In units $a=\mu^2=1$ their parameter integrals per unit squared Yukawa are
$$
 \begin{aligned}
 b_N(z)&=-\int_0^1(1-x)\ln[(1-x)(1-xz)]\,dx,\\
 h_L(B)&=-\int_0^1(1-x)\ln[1+x(B-1)]\,dx,\\
 h_\phi(B)&=-\int_0^1\ln[x(B-1+x)]\,dx
       -(B-1)\int_0^1\frac{1-x}{B-1+x}\,dx.
 \end{aligned}
$$
The last product has continuous limit zero at $B=1$. The sterile addition to a single $F$ is $4zG(z)F(z)b_N(z)/(16\pi^2)$. With $\beta=\sqrt{1-4/s}$, average each single-transfer kernel uniformly over $z=1-s/2+(s\beta/2)\cos\theta$. The leading projection is $A_0/\gamma=s\beta\bar F_e\bar F_\mu/(8\pi)$ for unit first Yukawas in both rows. Correcting either side uses its mean corrected kernel. For an extra Yukawa pair from block $f$, external residues multiply the leading cut by $-\sum_{i\in f}|c_i|^2[h_L(B_i)+2h_\phi(B_i)]/(16\pi^2)$. Each scalar tree carries two scalar legs and each flavor supplies two external lepton legs to the sewn cut. [Sterile self-energy](STERILE.md) · [External residues](EXTERNAL.json)

| $s/m^2$ | Proper vertex, electron | Proper vertex, muon | Two-scalar total, electron | Two-scalar total, muon |
|---:|---:|---:|---:|---:|
| 6 | $+1.10989314\,10^{-6}$ | $+6.46226760\,10^{-7}$ | $-1.38316851\,10^{-5}$ | $+1.78723358\,10^{-6}$ |
| 8 | $+2.13277672\,10^{-6}$ | $+1.58008833\,10^{-6}$ | $-2.65731240\,10^{-5}$ | $+3.67404321\,10^{-6}$ |
| 10 | $+3.20052146\,10^{-6}$ | $+2.81830430\,10^{-6}$ | $-3.98654872\,10^{-5}$ | $+5.83741910\,10^{-6}$ |

All entries divide out $\gamma$ and the extra squared Yukawa modulus of the indicated flavor. The independent program uses SciPy adaptive integration over the full crossed cut and 32- and 48-point Gauss–Legendre rules for the physical angular mean. Run `python -I -B -O check_independent.py --output INDEPENDENT.json`. [Checker](check_independent.py) · [Numerical record](INDEPENDENT.json)
