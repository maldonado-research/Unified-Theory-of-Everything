# Quartic scalar rescattering in the four-lepton amplitude

The leading four-lepton amplitude in this model contains one loop. Adding a single scalar self-interaction produces a two-loop correction whose coefficient can be calculated in full. Fix the physical scalar pole mass $m$, write $a=m^2$ and use $V\supset+\lambda(H^\dagger H)^2$. The external leptons are massless and their helicities are fixed. Keep the [four-singlet masses and Yukawas](../uv_completion/MATCHING.md), the causal upper rim, and the external phase $B_{34}^*B_{12}=s$. Write $\gamma=y_e^2(y_\mu^*)^2$. We calculate the monomial $\lambda\gamma$ in the amplitude. The two-loop terms with an additional pair of Yukawa couplings have different coupling dependence.

## Complete diagrams at this order

The two electron vertices must connect to one another by the electron-block Majorana propagator. The two muon vertices likewise connect through their own block. All four active lepton slots are external. There are two scalar ends of like charge at the muon vertices, two of opposite charge at the electron vertices, and two ends of each charge at the scalar quartic. Count the contractions according to the number of lines at the quartic vertex that return directly to that same vertex.

- **Zero self-contractions.** All four scalar lines connect to the central vertex. There are $2!2!=4$ labeled contractions for the external charged component. The scalar vertex is $-4i\lambda$. Each side has its own unconstrained loop momentum. After the external spin factor is removed, the diagram is the product of two scalar triangles and depends only on $s$.
- **One self-contraction.** The quartic supplies a mass tadpole on either scalar propagator of the original box. Summing the charged and neutral components of the doublet gives a momentum-independent self-energy. The mass counterterm cancels the entire insertion when the scalar pole mass is fixed. Its derivative, which sets the field residue correction, vanishes.
- **Two self-contractions.** The quartic is a disconnected vacuum diagram. The normalized scattering amplitude excludes vacuum diagrams.

A one-loop counterterm proportional to $\lambda$ alone cannot renormalize the sterile or lepton field, the sterile mass or the Yukawa vertex in this unbroken theory. Fixing a running scalar mass instead of the pole mass would also require converting the scalar mass used in the leading one-loop box. [Mass and vertex conventions](QUARTIC.md) · [Leading box and phase space](../uv_completion/CUT.md)

## The two-loop coefficient and its cut

For sterile squared mass $b$, define the positive Euclidean triangle continued to the upper rim,
$$
T_b(s)=\int_{x,y\ge0,\,x+y\le1}\frac{dx\,dy}{a(x+y)+b(1-x-y)-sxy-i0}.
$$
With the first mass in each pair labeled by $m$ or $M$, the dimension-minus-one functions are
$$
 S_\mu=m(T_a-T_{16a}),\qquad S_e=M(T_{M^2}-4T_{4M^2}).
$$
In Minkowski loop normalization each triangle is $-iT_b/(16\pi^2)$ before propagator and vertex factors. Restoring those factors and summing the scalar assignments gives
$$
\boxed{\mathcal M(s,c)=\mathcal M_0(s,c)+\lambda\gamma L(s)+\cdots,\qquad
 L(s)=-\frac{4s}{(16\pi^2)^2}S_e(s)S_\mu(s).}
$$
The spin convention is the same as the one used in the [finite box](../uv_completion/DISPERSION.md). Each $S$ uses its upper-rim value. There is no complex conjugation between them in the amplitude. $\mathcal M_0$ denotes the full one-loop box including the factor $\gamma$.

A separate normalization check cuts both scalar lines on either side of the central vertex. Let $\beta=\sqrt{1-4a/s}$, $h=s\beta/2$ and $N=s/2-a$. Put $J_b=(2h)^{-1}\ln[(b+N+h)/(b+N-h)]$. Then $\mathrm{Im}\,T_b=\pi\beta J_b$. For each flavor, $\overline F_\ell$ is the sum of these $J_b$ with the mass residues in $S_\ell$. Including the two identical intermediate scalars and taking the amplitude mean over the two scalar assignments, the corrected electron vertex equals $-2\lambda S_e/(16\pi^2)$ after its original $-2B_{34}^* y_e^2$ is stripped. The muon vertex obeys the same relation in its own conventions. Sewing the two trees with one corrected vertex gives
$$
 \frac{\mathrm{Disc}_s L}{2i}
 =-\frac{s\beta}{4\pi(16\pi^2)}
 [\overline F_\mu\,\mathrm{Re}\,S_e+\overline F_e\,\mathrm{Re}\, S_\mu].
$$
The terms proportional to $\overline F_\mu\mathrm{Im}\,S_e-\overline F_e\mathrm{Im}\,S_\mu$ cancel. The result equals the imaginary part of the product in the boxed expression with its stated sign. At zero external momentum the electron factor also gives the separate finite quartic correction to the lepton–scalar interaction. Adding that threshold a second time to the full diagrams would double count it. [Two-scalar unitarity normalization](../uv_completion/CUT.md) · [Triangle spectral representation](../energy_scan/LOWER_OPERATOR.md)

## Rate response at fixed masses and base Yukawas

Write $u(s,c)=\mathcal M_0(s,c)/\gamma$, and $r=m^2\sigma/|\gamma|^2$. With $\langle f\rangle_c=\int_{-1}^1 f(c)dc/2$ and the identical-final-state factor, the terms at the specified coupling orders are
$$
 r_0=\frac{m^2}{32\pi s}\langle |u|^2\rangle_c,\qquad
 r_1=\frac{m^2}{16\pi s}\,\mathrm{Re}[\langle u\rangle_c^* L],\qquad
 r=r_0+\lambda r_1+\cdots .
$$
Since $L$ is independent of the scattering direction in these phases, only the angular mean of the one-loop amplitude enters the interference. A square of the two-loop correction would be of higher order and by itself would omit other terms at that order.

| $s/m^2$ | $r_0$ | $r_1$ | $r_1/r_0$ |
|---:|---:|---:|---:|
| 6 | $4.00925370\times10^{-11}$ | $-4.82330547\times10^{-11}$ | $-1.20304322$ |
| 8 | $5.61821543\times10^{-11}$ | $-5.30510123\times10^{-11}$ | $-0.94426803$ |
| 10 | $7.17400709\times10^{-11}$ | $-5.72375648\times10^{-11}$ | $-0.79784650$ |

The ratio between the last and first settings has intercept $1.7893622161$ and linear coefficient $+0.7250436854$ in $\lambda$. Its fractional slope is $+0.4051967114$. Both individual rates fall for a small positive quartic, and the lower-energy rate falls faster. No finite value of $\lambda$ is assumed. The omitted powers $\lambda^2$ and the distinct two-loop Yukawa terms require their own calculations. The numerical responses refer to the chosen mass and coupling texture in a gaugeless theory.

The [program](check_rescattering.py) evaluates the triangles by Feynman parameters and unsubtracted dispersion and checks their upper-rim values against sixteen evaluations from an external library. It recomputes the box mean and checks the endpoints against the previous direct box calculation. It also tests the scalar Wick contractions and the cut normalization. [Numerical record](RESCATTERING.json) · [QCDLoop inputs and reproduction](QCDLOOP.md)

## Can an energy scan identify the scalar interaction?

Compare the quartic response with the independent tree-input change $Y_{e2}^2=-2y_e^2(1+\delta)$ at fixed masses. Let $\bar\sigma_i$ be the arithmetic mean of the two reciprocal cross sections at $i=s/m^2$. In the CP-aligned reference texture the forward and inverse rates coincide. Across a complex perturbation, the reciprocal mean removes the term linear in $\mathrm{Im}\,\delta$. Define $\ell_i=(\partial_\lambda\ln\bar\sigma_i)_0$ and $d_i=(\partial_{\mathrm{Re}\,\delta}\ln\bar\sigma_i)_0$. The quartic and detuning enter as independent perturbations around $\lambda=\delta=0$. The three $d_i$ are approximately $67.28565,52.52387,44.24096$. The endpoint ratio alone constrains the combination
$$
 \Delta\ln(\bar\sigma_{10}/\bar\sigma_6)
  \simeq 0.405197\,\lambda - 23.04468\,\mathrm{Re}\,\delta .
$$
For instance, a shift $\mathrm{Re}\,\delta\simeq-0.0176\lambda$ at zero quartic mimics the fractional endpoint-ratio change due to the quartic alone. Mixed perturbations, higher orders and other parameter changes are outside this linear response. [Exact leading rates for complex detuning](DETUNING.md)

To remove a common normalization and both linear detuning directions, define $r_{i6}=(\bar\sigma_i/\bar\sigma_6)/(\bar\sigma_i/\bar\sigma_6)_0$. Use the dimensionless statistic
$$
 X=\ln r_{8,6}-w\ln r_{10,6},\qquad
 w=\frac{d_8-d_6}{d_{10}-d_6}\simeq 0.64057189 .
$$
The weight makes the linear real-detuning coefficient vanish. The remaining quartic derivative is $(\ell_8-\ell_6)-w(\ell_{10}-\ell_6)\simeq-0.00078244$. Thus the two responses are linearly distinguishable at these three settings. Their energy dependences are so similar that removing the detuning also suppresses almost all the quartic response. At $\lambda=0$, the exact leading detuned rates give $X\simeq+0.00302$ at $\delta=+0.01$ and $+0.00673$ at $\delta=-0.01$. Even small quadratic effects matter for this statistic. These calculations motivate testing a wider range of collision energies or an independent measurement of the couplings before using an energy scan to infer $\lambda$.
