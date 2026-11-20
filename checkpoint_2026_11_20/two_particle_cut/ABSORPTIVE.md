# The finite absorptive amplitude through two Yukawa loops

In the unbroken, gaugeless four-singlet benchmark, two fixed-order on-shell cuts determine the order-six Yukawa discontinuity of polarized $\mu^-_L\mu^-_L\to e^-_Le^-_L$ below the second sterile threshold. The charged scalar pair contributes through its production amplitude at one loop. The charged scalar, muon and lightest sterile contribute through their tree production amplitudes. At $s=8m^2$, these channels contribute $+3.67404321\,10^{-6}$ and $-1.45757749\,10^{-6}$ to the reduced muon-Yukawa coefficient in the angular mean. [Numerical results](ABSORPTIVE.json) · [Complete three-particle trees](../../checkpoint_2026_11_19/yukawa_continuum/FULL_THREE_BODY.md)

## Inputs and channel count

Set the scalar pole mass $m=1$ in the formulas below. The reference full-theory $\overline{\rm MS}$ scale is $\mu=m$. The four positive sterile mass parameters are $(10,20,1,4)m$, with electron Yukawa row $(y,i\sqrt2 y,0,0)$ and muon row $(0,0,g,-ig/2)$. External charged leptons are massless, and external states have unit pole normalization. Define $E=|y|^2$, $G=|g|^2$, $\gamma=y^2(g^*)^2$ and use the spin phase $B_{34}^*B_{12}=s$. Holding the scalar pole fixed permits a finite residue at the coincident lowest lepton–sterile threshold. At this perturbative order the three-body tree uses the reference sterile masses. The expansion in couplings at fixed energy is nonuniform on the scale of the induced pole splitting $m|g|^2/(16\pi^2)$. A formulation with degenerate physical pole masses would require a different running input and a corresponding change in the finite coefficients. [Action and propagator conventions](../../checkpoint_2026_11_19/uv_completion/MATCHING.md) · [External states](EXTERNAL.md) · [Pole thresholds and scalar lifetime](POLE_THRESHOLD.md)

For $4<s<25$ write the anti-Hermitian part of the transition amplitude in the form

$$
 \frac{\mathcal M_{fi}-\mathcal M_{if}^*}{2i\gamma}
 =D^{(4)}(s,c)+E D_e^{(6)}(s,c)+G D_\mu^{(6)}(s,c)+\cdots .
$$

The superscripts count Yukawa factors before division by $\gamma$ and the additional squared coupling. The common phase makes these functions real. The ellipsis includes higher orders and other interactions. At degree six the two-body cut pairs a degree-two tree and a degree-four one-loop production amplitude, in either order. Each of the two incoming charged leptons must attach to a Yukawa vertex. Charge and the particle slots at the two tree vertices select $\phi^-\phi^-$. On the three-body cut both sides are degree-three trees. Charge and sterile-line counting leave $\phi^-\mu^-N_1$ in this energy window. A connected Yukawa tree into four particles needs at least four vertices, so that cut first contributes at degree eight. On-shell splittings of an external massless charged lepton into a sterile and a scalar are kinematically closed. The calculation treats the pure-Yukawa coefficient. The separate term linear in the Higgs quartic was computed in the preceding work. [Two-body counting and state factors](../../checkpoint_2026_11_19/uv_completion/CUT.md) · [Single-quartic sector](../../checkpoint_2026_11_19/radiative_stability/MATCHING_AUDIT.md)

On a one-loop production side, the two lepton and two scalar fields attach to three holomorphic Yukawa vertices. The fourth vertex is conjugate. At least one holomorphic vertex carries an external lepton and scalar together. Cutting its sole internal line makes that diagram reducible. Enumerating the contractions gives corrections to external legs, a sterile self-energy insertion and proper three-point insertions at the two endpoints. The conjugate orientation obeys the same counting. After minimal subtraction the ultraviolet poles in this sector are accounted for by field and parameter renormalization. The finite proper vertex survives at spacelike momentum. [Contraction census](../../checkpoint_2026_11_19/radiative_stability/MATCHING_AUDIT.md) · [Vertex numerator and signs](VERTEX.md) · [Sterile matrix inversion](STERILE.md)

## Angular dependence

Let $q(x)=s/2-1-(s\beta/2)x$, $\beta=\sqrt{1-4/s}$, and $x=\hat{\boldsymbol p}_1\cdot\hat{\boldsymbol k}$. With $y'=\hat{\boldsymbol p}_3\cdot\hat{\boldsymbol k}$ and external cosine $c=\hat{\boldsymbol p}_1\cdot\hat{\boldsymbol p}_3$, define the single-assignment kernels $F_f$, their normal-propagator weights $G_f$ and their sum over identical scalar assignments $S_f(x)$ by

$$
 F_f(q)=\sum_i\frac{w_{fi}}{q+b_{fi}},\qquad
 G_f(q)=\sum_i\frac{h_{fi}}{q+b_{fi}},\qquad
 S_f(x)=F_f(q(x))+F_f(q(-x)).
$$

In electron and muon order, the squared masses $b$, anomalous residues $w$ and conserving weights $h$ are $((100,400),(10,-40),(1,2))$ and $((1,16),(1,-1),(1,1/4))$. The leading cut is $s\beta\langle S_\mu(x)S_e(y')\rangle_{\boldsymbol k}/(32\pi)$. At the next order replace $F_f$ in one tree at a time by

$$
 U_f(q)=-\frac{2qG_f(q)}{16\pi^2}
     \left[F_f(q)\,B(q)+\sum_i w_{fi}I_i(q)\right],
 \quad
 B(q)=2+\frac1q-\frac{(1+q)^2}{q^2}\ln(1+q).
$$

Here $I_i$ is the proper vertex integral with an internal sterile of squared mass $b_{fi}$. The first term accounts for both weak components in the sterile bubble. The second sums the two endpoints. Each replacement contributes to the coefficient of its own $E$ or $G$. In addition the four external scalar factors and the two active-lepton factors of each flavor multiply the whole leading cut, giving $D^{(4)}(0.229679017971E-0.0170128006694G)$. These field factors include the loop normalization. In a convention with a pole-normalized scalar propagator inside the leading box, converting its four Yukawa endpoints to the running input couplings recovers the same pole-residue contribution. The three-body tree cut already includes the scalar continuum once. [External residues](EXTERNAL.md) · [Massive triangle and scalar-leg dispersion](VERTEX.md) · [Independent crossed-channel dispersion](INDEPENDENT.md)

The full two-body contribution follows by integrating this product over the intermediate direction. Adding the earlier three-body formula gives $D^{(6)}_f(s,c)$ at any physical angle in the open interval. A compact numerical evaluation uses Legendre moments $f_\ell=\tfrac12\int_{-1}^1dx\,P_\ell(x)F(q(x))$. At each two-body energy or each three-body recoil mass the intermediate-angle integral of two assignment sums gives $4\sum_{\ell\ \mathrm{even}}(2\ell+1) f_\ell g_\ell P_\ell(c)$. This follows by rotating either axis in the spherical-harmonic addition theorem. In the three-body term the muon kernel includes the previously derived, decay-angle-integrated exchange diagram. The numerical checker evaluates moments through $\ell=16$ and repeats the angular integral directly at $c=0,0.37,1$ at each of the three benchmark energies. The leading cut and the separate three-body angular mean reproduce the archived results. [Combined checker](check_absorptive.py) · [Three-body kernel](../../checkpoint_2026_11_19/yukawa_continuum/FULL_THREE_BODY.md)

## Energy dependence and rate information

Write $D_{f,0}=\langle D_f\rangle_c$ with the uniform external-cosine measure. The following dimensionless coefficients include all loop and symmetry factors.

| $s/m^2$ | Leading $D^{(4)}_0$ | $D^{(6)}_{e,0}$ | Two-body part of $D^{(6)}_{\mu,0}$ | Three-body part | $D^{(6)}_{\mu,0}$ |
|---:|---:|---:|---:|---:|---:|
| 6 | $-6.50420877\,10^{-5}$ | $-1.38316851\,10^{-5}$ | $+1.78723358\,10^{-6}$ | $-4.35689215\,10^{-7}$ | $+1.35154436\,10^{-6}$ |
| 8 | $-1.24868828\,10^{-4}$ | $-2.65731240\,10^{-5}$ | $+3.67404321\,10^{-6}$ | $-1.45757749\,10^{-6}$ | $+2.21646572\,10^{-6}$ |
| 10 | $-1.87155112\,10^{-4}$ | $-3.98654872\,10^{-5}$ | $+5.83741910\,10^{-6}$ | $-2.96120029\,10^{-6}$ | $+2.87621881\,10^{-6}$ |

At $s/m^2=9.84983756$ the negative three-body term cancels half the positive two-body muon coefficient. As the energy increases further, competition between the two terms gives their sum a local maximum of $3.85746\,10^{-6}$ near $s/m^2=17.50997$. These are results for the selected flavor coefficient, input texture and scale. The energy variation also appears in the signed ratio of angular means,

$$
 \frac{D_0(10)}{D_0(6)} =
 2.87744627\left[1+0.0003502794 E+0.005411434 G+\cdots\right].
$$

Constant external residue factors cancel in this ratio to the order shown. Its finite response still depends on the parameter convention. The explicit scale derivative of the reconstructed cut cancels the implicit running of the leading cut. At $c=0$ the scale derivatives for both flavor coefficients agree with the imaginary parts obtained earlier by differentiating the full box. [Scale and coordinate checks](EXTERNAL.md) · [Earlier full-amplitude logarithm](../../checkpoint_2026_11_19/yukawa_continuum/RG_LOG.md)

For a common phase, the absorptive contribution to the scaled polarized rate is $m^2\langle D^2\rangle_c/(32\pi s)$ with dimensionful $s$ restored. The checker evaluates its term linear in each new Yukawa pair by angular orthogonality. At $s=8m^2$ this component is $1.93873621\,10^{-11}+8.25158295\,10^{-12} E-6.88264614\,10^{-13} G$. The physical rate also contains the interference of the dispersive parts of the one- and two-loop amplitudes. The local real threshold terms are determined separately; the full real amplitude over the scan, the complete corrected cross section and a connection between these benchmark inputs and a unified theory remain open. [Threshold continuation](CUSP.md) · [Scan and rate components](ABSORPTIVE.json)
