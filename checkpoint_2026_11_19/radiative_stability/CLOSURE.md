# A four-rate test of arbitrary complex electron-coupling detuning

The three energies $s/m^2=6,8,10$ admit an exact consistency relation for the leading four-singlet cross sections. The inputs are the three direct $\mu^-_L\mu^-_L\to e^-_Le^-_L$ rates and one inverse rate at $s/m^2=6$. A common positive normalization and any fixed complex change $Y_{e2}^2=-2y_e^2(1+\delta)$ cancel in this relation. The quartic interaction gives a small, nonzero first derivative of the relation in the tuned model. The masses and relative muon Yukawas must take their benchmark values throughout. [Rate coefficients](DETUNING.md) · [Quartic calculation](RESCATTERING.md)

## Construct the test from rates

At energy $i=s/m^2$ denote the direct and inverse scaled rates by $R_i^+=m^2\sigma_{\mu\mu\to ee}/a_i$ and $R_i^-=m^2\sigma_{ee\to\mu\mu}/a_i$. Here $a_i$ is the **unit-coupling** tuned scaled cross section, so the unknown factor $N=|y_e y_\mu|^4>0$ remains in each $R$. A common further multiplicative normalization may be absorbed into $N$. Put $\delta=x+iy$ and use the angular Gram coefficients $a_i,b_i,c_i$ from the detuning calculation. Then
$$
d_i=2\operatorname{Re}b_i/a_i,\quad e_i=c_i/a_i,\quad
f_i=-2\operatorname{Im}b_i/a_i,\qquad
R_i^\pm=N[1+d_ix+e_i(x^2+y^2)\pm f_i y].
$$
In particular $f_6=-32.01989927$ is nonzero because the two heavy-mass contributions have different absorptive phases. The following four inputs determine four real quantities:
$$
 T=\frac{R^+_6-R^-_6}{2f_6},\qquad
 r_i=R_i^+-f_i T,\qquad
 \begin{pmatrix}U\\ V\\ W\end{pmatrix}=
 \begin{pmatrix}1&d_6&e_6\\1&d_8&e_8\\1&d_{10}&e_{10}\end{pmatrix}^{-1}
 \begin{pmatrix}r_6\\r_8\\r_{10}\end{pmatrix}.
$$
If both directions are measured at each energy, their arithmetic means can be used for $r_i$. Each measured signed half-difference divided by $f_i$ then provides a consistency check on $T$. With just the four specified inputs the reconstruction above uses the predicted energy dependence of the detuning.

The three-by-three matrix $D$ in this equation has determinant $-345.940168$. In the leading model $(U,V,W,T)=N(1,x,x^2+y^2,y)$ at **all** finite $x,y$. Thus the dimensionless statistic
$$
 \mathcal I=\frac{W}{U}-\left(\frac VU\right)^2-\left(\frac TU\right)^2
$$
vanishes exactly. Conversely, $U>0$, $W\ge0$ and $\mathcal I=0$ allow $N=U$ and $\delta=(V+iT)/U$ to reproduce the four rates. If the inputs comprise only the three reciprocal means, positivity $UW\ge V^2$ is the consistency condition and the sign of the imaginary detuning is undetermined. For a known real detuning, three direct rates suffice and $T=0$. The executable checks the exact identity using rational arithmetic and symbolic $x,y,N$. It then reconstructs three independently integrated sets of direct and inverse rates at real, imaginary and mixed detunings. [Checker](check_closure.py) · [Results](CLOSURE.json)

## Response to the quartic interaction and numerical conditioning

Keep $\delta=0$ and the scalar pole mass fixed. The tuned quartic correction has one common Yukawa phase and equal reciprocal rates at the order calculated. Consequently $T=0$ and $r_i=N(1+\lambda\ell_i)+O(\lambda^2)$ within this perturbative sector, with $\ell_i$ from the two-loop interference. At $\lambda=0$, $U=N$ and $V=W=0$. Differentiation gives
$$
 \left.\frac{d\mathcal I}{d\lambda}\right|_0
 =(D^{-1}\ell)_3=5.21219\times10^{-5}.
$$
The derivatives $(U',V',W')/N$ are $(0.01998055,-0.01925191,0.0000521219)$. A positive infinitesimal quartic can therefore imitate a positive inferred $|\operatorname{Im}\delta|^2$ if only the three means are used. The fourth directional measurement tests that interpretation. The derivative above is local to the tuned point. It uses the specified quartic interference at fixed masses and Yukawas. Other corrections to the rates require their own projections.

The slope is the ratio $\det[\mathbf1,d,\ell]/\det[\mathbf1,d,e]$. Exact interval arithmetic over independent relative perturbations of up to one part per million in *every entry* of $d,e,\ell$ encloses it between $5.1636\times10^{-5}$ and $5.2609\times10^{-5}$. This interval tests coefficient sensitivity, including the precision at which the inputs are printed. Errors shared by the physical model or by its loop normalization are outside this numerical check.

For three separately measured reciprocal means with independent fractional errors bounded by $\epsilon$, the first-order error in $\mathcal I$ at the tuned point can reach $\epsilon\sum_i|(D^{-1})_{3i}|=0.1332293\,\epsilon$. At this point a common normalization error cancels. Resolving just the linear response in a worst-case comparison would require those differential errors to be smaller than about $3.91\times10^{-4}|\lambda|$. Using four direct/inverse inputs also requires propagating the error in $T$ into the reconstructed means. These figures quantify the numerical difficulty of separating the two effects. Experimental sensitivity depends on rates, luminosities, backgrounds and correlations.
