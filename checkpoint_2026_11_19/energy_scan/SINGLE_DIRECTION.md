# A four-energy lower bound with a muon beam alone

For the leading amplitude of $\mu^-_L\mu^-_L\to e^-_Le^-_L$ in the specified two-mass model, the **equal-exposure mean of four direct cross sections** has a sharp positive floor. The four center-of-mass energies obey $s/a=5,6,9,11$, where $a>0$ is the common scalar squared mass. After minimizing over one complex Weinberg coefficient and two independent complex local four-lepton coefficients, the floor is
$$
 a\,\frac14\sum_{i=1}^4\sigma_{\mu\mu\to ee}(a\xi_i)
 \;\geq\;\kappa |g_7|^2,\qquad
 \kappa\simeq 1.33764501629727\times10^{-9}.
$$
The determinant formula below defines the exact $\kappa$. Rational interval arithmetic certifies $1.33764501629\times10^{-9}<\kappa<1.33764501631\times10^{-9}$. One common choice of the three nuisance coefficients attains the minimum of the four-energy mean. The coefficient uses polarized incident leptons, integration over the full solid angle, and a factor of $1/2!$ for identical final electrons. [Numerical receipt](SINGLE_DIRECTION_RESULTS.json) · [Numerical certificate](CERTIFICATE.md)

## Physical inputs and the fitted amplitudes

Work in the gaugeless, unbroken theory with massless external leptons, the charged weak component, only the muon neutrino-Yukawa row occupied, and electron-flavor seed entries $C_{7,ee}$ and $C_{5,ee}$. In the positive Takagi basis choose $m_j/\sqrt a=(1,4)$ and $w_j=m_j(Y^*_{\mu j})^2=(w_0,-w_0)$ with $w_0\ne0$. The opposite weights remove the primitive dimension-seven UV pole. The scalar is stable against decay to a massless lepton and either sterile state in this sector, with zero two-body phase space for the equal-mass state. Write
$$
 \xi=s/a,\quad n(\xi)=\xi/2-1,\quad
 g_7=a C_{7,ee}w_0,\quad g_5=C_{5,ee}w_0,\quad
 f(\xi)=a[T(a\xi;a,a)-T(a\xi;a,16a)] .
$$
Here $T(s;a,b)$ is the positive-Euclidean triangle analytically continued to the upper lip. The positive-action convention gives $K_7=B_0+(s/2-a)T$ and $K_5=-T$. The finite bubble difference $g_7\Delta B_0$ can be absorbed into $d_0$. After removing the external ordered spin tensor $W$, the full $J=0$ reduced amplitude in this family is
$$
 r(\xi)=(g_7n(\xi)-g_5)f(\xi)+d_0+d_1\xi,\qquad
 \mathcal M_0=\frac{W}{8\pi^2a}r(\xi),\qquad |W|=s .
$$
The off-diagonal entries $d_0,d_1$ of two independent Hermitian local interactions may each be complex. Hermiticity supplies their conjugates in the inverse direction. Every coefficient above is held fixed throughout the scan. [Loop and cut normalization](../on_shell/ON_SHELL_TRIANGLE.md) · [Weinberg and derivative seed conventions](LOWER_OPERATOR.md)

For the angle-independent form, $\sigma=|\mathcal M_0|^2/(32\pi s)$, so $a\sigma=\xi |r|^2/(2048\pi^5)$. More generally, remove the external helicity phase by dividing the amplitude pointwise by $W/(8\pi^2a)$. For a reduced angle-dependent amplitude $R$, the exact polarized total rate is proportional to $\langle |R|^2\rangle_\Omega$, with $\langle X\rangle_\Omega=\int d\Omega\,X/(4\pi)$. Both pairs have helicity difference zero, and the reduced $J=0$ coefficient in this normalization is $r=\langle R\rangle_\Omega$. Orthogonality gives $\langle |R|^2\rangle_\Omega=|r|^2+\langle |R-r|^2\rangle_\Omega$. Thus the bound also applies if additional partial waves are present and the local $J=0$ term has the form $d_0+d_1\xi$. Setting the higher waves to zero attains the rate of the reduced model. [Flux and identical-particle measure](../on_shell/INVERSE_RATE_BOUND.md)

## Complex projection and the sharp coefficient

Let $h_i=\sqrt{\xi_i/4}$. Define a four-by-three matrix $X$ and a four-vector $u$ by
$$
 X_{i,:}=h_i(f_i,1,\xi_i),\quad u_i=h_i n_i f_i,\quad
 \eta=(-g_5,d_0,d_1)^T .
$$
The $J=0$ contribution to the mean rate, multiplied by $2048\pi^5a$, is $\|g_7 u+X\eta\|^2$. The norm is the complex Euclidean norm. Put $D_i=(-1)^{i+4}\det X_{\widehat i}$ for row labels $i=1,\ldots,4$, where the hat deletes that row. For $D\ne0$,
$$
 q_i=\frac{\overline{D_i}}{(\sum_j|D_j|^2)^{1/2}},\qquad
 q^\dagger X=0,\qquad q^\dagger u=\frac{\det[X\,u]}{(\sum_j|D_j|^2)^{1/2}}.
$$
The normalized conjugate cofactors give the Hermitian left null vector. Orthogonal projection onto it gives the minimum, since all three columns of $X$ are independent. Equivalently, by Cauchy–Binet,
$$
 \kappa=
 \frac{1}{2048\pi^5}\frac{|\det[X\,u]|^2}{\sum_i |D_i|^2}
 =\frac{1}{2048\pi^5}
 \frac{\det([X\,u]^\dagger[X\,u])}{\det(X^\dagger X)} .
$$
The minimizing nuisance vector is $\eta_*=-g_7(X^\dagger X)^{-1}X^\dagger u$. For example, the coefficients after absorbing the bubble are approximately
$$
 \frac{g_5}{g_7}=0.35043247-1.39776343i,\quad
 \frac{d_0}{g_7}=-0.27456776-0.99478889i,\quad
 \frac{d_1}{g_7}=0.04469979-0.11312283i .
$$
Each three-row nuisance minor is nonzero for this design. For rows $i<j<k$, divide the minor by $h_i h_j h_k$ to obtain $(\xi_k-\xi_j) f_i+(\xi_i-\xi_k) f_j+(\xi_j-\xi_i) f_k$. Inserting the certified triangle rectangles gives its real part between 0.29 and 0.61 for each of the four triples. Consequently any three of these four direct amplitudes can be tuned to vanish exactly, while a common choice cannot cancel all four for $g_7\ne0$. The exposure mean implies that at least one of the four cross sections reaches the stated lower bound on the mean. [Signed minors and optimum](SINGLE_DIRECTION_RESULTS.json)

## Evaluation on the physical lip

For $b/a=1,16$, the difference of the two angular cut logarithms reduces to
$$
 L(v)=\log\frac{1+kv}{1-kv},\qquad k=15/17,\qquad v=\sqrt{1-4/\xi}.
$$
Indeed, with $y=(1+v)/(1-v)$ the ratio inside the difference of the two logs simplifies to $(16y+1)/(16+y)$. The triangle spectral gap is $L(v)/\xi$. Its unsubtracted dispersion representation in $t=4a/s'$ is
$$
 f(\xi+i0)=\int_0^1\frac{L(\sqrt{1-t})}{4-(\xi+i0)t}\,dt .
$$
Writing $t=1-u'^2$ and subtracting $L(v)$ at the pole gives
$$
 \xi f(\xi+i0)=\int_0^1\frac{2u'[L(u')-L(v)]}{u'^2-v^2}\,du'
 +L(v)\left[\log\frac{1-v^2}{v^2}+i\pi\right].
$$
The subtracted integrand has a removable singularity. Using the hyperbolic arctangent difference formula avoids cancellation there. At 64 decimal digits, this integral and two independently evaluated Feynman-simplex logarithmic integrals give $f_i$ values that agree to absolute error below $5\times10^{-63}$ between methods. This comparison measures reproducibility. The interval calculation supplies the certificate. [High-precision checker](verify_single_direction.py) · [Interval checker](verify_scan_certificate.py)

| $\xi_i$ | $\operatorname{Re} f_i$ | $\operatorname{Im} f_i$ |
|---:|---:|---:|
| 5 | 0.688648027318 | 0.524315519178 |
| 6 | 0.495341316572 | 0.588478163543 |
| 9 | 0.205980396526 | 0.550610436224 |
| 11 | 0.115374883848 | 0.499763220179 |

For the weighted columns in the order specified above, $\det[X\,u]\simeq-0.075074694116-0.024958390436i$, $\det(X^\dagger X)\simeq7.466100749121$, and $|q^\dagger u|^2\simeq0.000838340006383398$. The two-norm condition numbers of $X$ and $[X\,u]$ are approximately 262 and 1790 in this column normalization. Direct QR minimization, cofactor projection and the Gram ratio agree. The QR fit works directly on the weighted amplitudes. [Reproducible determinants and singular values](SINGLE_DIRECTION_RESULTS.json)

## Why a continuous scan also has a positive minimum

Take any compact above-threshold interval of nonzero length with a positive continuous exposure density. In its weighted $L^2$ norm, the three nuisance functions span a finite-dimensional closed subspace. A zero minimum at $g_7\ne0$ would imply a single triple of complex constants with
$$
 F(z)=(g_7 n(z)-g_5) f(z)+d_0+d_1 z
$$
zero almost everywhere on the upper lip across that interval. Continuity makes the upper-rim values zero on the whole interval. The physical-sheet $f$ is holomorphic on $\mathbb C\setminus[4,\infty)$, as is clear from its convergent dispersion integral. Schwarz reflection across a small open segment is legitimate for $F$ because **its boundary values are zero**, even with complex coefficients. The reflected extension vanishes on an interior line segment, so the identity theorem gives $F=0$ in the upper half-plane and hence throughout the connected physical cut plane. Taking both boundary values of this identity would give $(g_7 n(\xi)-g_5)\,2\pi i L(v)/\xi=0$ everywhere above threshold. Here $L(v)>0$ and $n$ is nonconstant, contradicting $g_7\ne0$. This proves strict positivity for the stated continuous scans and complements the separate quantitative four-point calculation.

The result constrains cross sections computed from the retained leading amplitude and the stated local $J=0$ family. Other loop shapes, higher independent local powers, a different sterile spectrum, or additional amplitudes in a broken-phase model change the fitting problem. If more scalar-lepton operators contribute at the same order, the derivative signal depends on the net coefficient of the projected tree slope. Converting the result to an experimental rate would also require a physical energy scale and Wilson coefficients. [Operator scope and power counting](CONTACT_SCOPE.md)
