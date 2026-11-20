# A dimension-five contact in the charged lepton-pair channel

For the charged process $\mu^-_L\mu^-_L\to e^-_Le^-_L$, a constant electron-pair contact and the derivative contact enter the two-scalar cut through the single complex combination $(s/2-a)C_7-C_5$. Here $a=m_\phi^2>0$. At one energy an independent $C_5$ can cancel the derivative contribution exactly. At two distinct energies above threshold the same cancellation requires $C_7=C_5=0$ if the sterile exchange has a nonzero projected weight at both energies.

Use the unbroken, gaugeless, zero-background theory, massless external leptons and the equal charged weak components. In original action coordinates, $X=\phi^\dagger L$, $\phi=\widetilde H$, and the electron contact terms in $\mathcal L$ are
$$
 \tfrac12\bar X_e C_5 X_e^c
 +\tfrac12\bar X_e C_7(-\Box)X_e^c+\mathrm{h.c.}
$$
Both coefficients denote their $ee$ entries, with mass dimensions $-1$ and $-3$. Only the muon row of $Y$ is active. At amplitude level retain one of these contacts and two Yukawas, with no other interactions in the cuts. Both orientations use causal propagators. With the ordered external spinor $W=-B_{34}^*B_{12}$, where $B_{ij}=u_i^T\epsilon u_j$, $|W|=s$. A positive Takagi mass $m_i$ contributes $w_i=m_i(Y_{\mu i}^*)^2$ and $b_i=m_i^2$. Each flavor direction has one seed orientation. [Positive dimension-five action](../../checkpoint_2026_11_04/common_action/COMMON_ACTION_MAP.md) · [Spin, flavor and derivative normalization](../on_shell/ON_SHELL_TRIANGLE.md)

## Electron tree and signed cut

Strip the spinors from the trees only after summing the two identical-scalar assignments. For $t_j=(p_j-k_1)^2$, $u_j=(p_j-k_2)^2$ on the scalar cut,
$$
 M_{34,\phi\phi}=-B_{34}^*[2C_5+C_7(t_3+u_3)],\qquad
 M_{\phi\phi,12}=B_{12}\sum_i w_i[(t_1-b_i)^{-1}+(u_1-b_i)^{-1}].
$$
The cut electron tree has intercept $C_5+aC_7$ and slope $-C_7/2$ in $s$ after dividing by $-2B_{34}^*$, since $t_3+u_3=2a-s$. The inverse electron tree uses the conjugated coefficient and conjugated spinor. In particular, matching this reduced on-shell amplitude at one energy fixes a single complex linear combination. Local changes of basis must transport both tree inputs and the induced interactions to compare this slope and intercept.

Set $h=s/2-a$, $v=s\beta/2$, $\beta=\sqrt{1-4a/s}$, and $x=\hat{\boldsymbol k}_1\cdot\hat{\boldsymbol p}_1$. Then $t_1=-h+vx$, $u_1=-h-vx$. Define the real transfer density for $s>4a$ by
$$
 P_b(s)=\beta\int_{-1}^1\frac{dx}{2(b+h-vx)}
       =\frac1s\log\frac{b+h+v}{b+h-v}.
$$
It vanishes below and at threshold. The product of trees contains two equal integrated sterile routings and the electron factor $2(C_5-hC_7)$. Sewing with $\tfrac12\times\tfrac1{2!}\int d\Phi_2$, where $d\Phi_2=\beta\,d\Omega/(32\pi^2)$, supplies
$$
 A_{fi}=\frac{M_{fi}(s+i0)-M_{if}(s+i0)^*}{2i}
      =\frac{W}{8\pi}\,(h C_7-C_5)\sum_i w_i P_{b_i}(s).              \tag{1}
$$
The first half comes from unitarity and the factorial counts the intermediate scalars. The minus on $C_5$ follows from the spacelike sterile transfer denominator. This convention reproduces the positive derivative density $hP_b$ at each mass. The external identical-lepton phase-space factorial enters only on forming a rate. [Sewing and state normalization](../on_shell/INDEPENDENT_CUT.md)

## Finite lower-operator form factor

Use the three denominators $D_1=k^2-a$, $D_2=(P-k)^2-a$, $D_N=(k-p_1)^2-b$, with $P=p_1+p_2$ and $+i0$. The equal-component derivative numerator averages to $(D_1+D_2)/2-h$. The lower-operator vertex has the same Wick factors and numerator one. In units of $i/(16\pi^2)$ the scalar triangle equals $-T_b$,
$$
 T_b(s)=\int_{x,y\geq0,\ x+y\leq1}\frac{dx\,dy}{d-sxy-i0},
 \qquad d=a(x+y)+b(1-x-y).
$$
The sign follows from the three denominator signs in the Wick rotation. The positive dimension-five action thereby gives $K_{5,b}=-T_b$ and the derivative action gives $K_{7,b}=B_b+hT_b$, with $B_b=B_0(0;a,b)$. The direct loop amplitude before adding a matched local four-lepton term is
$$
 M_{fi}^{\rm loop}=\frac{W}{8\pi^2}\sum_i w_i
    [ C_7 B_{b_i}+(h C_7-C_5)T_{b_i}(s) ].                           \tag{2}
$$
The reciprocal upper rim sends $(W,w_i,C_5,C_7)$ to their conjugates and evaluates the same $T_{b_i}(s+i0)$. This gives (1) even for complex couplings. A complex off-diagonal amplitude can have an ordinary imaginary part different from $A_{fi}$.

One real integral computes the finite $C_5$ primitive. Write $d_0=b+(a-b)x$, $d_1=a-sx(1-x)$, and $q=b-a+sx$. Integrating $y$ gives
$$
 \operatorname{Re}T_b(s)=\int_0^1\!dx\;\frac{\log(d_0/|d_1|)}{q},
 \qquad \operatorname{Im}T_b(s)=\pi P_b(s).
$$
Here the logarithmic endpoint values are integrated by a limit, and the ratio at a removable $q=0$ is $(1-x)/d_0$. For $b>0$, the support of $\delta(d-sxy)$ begins at $4a$. Its $y$ integral contributes $1/q$ for $x\in[(1-\beta)/2,(1+\beta)/2]$. The resulting density agrees with the angular integral. Since $1/(d-sxy)= (xy)^{-1}/[d/(xy)-s]$ in the simplex interior, the triangle also obeys the **unsubtracted** representation
$$
 T_b(s)=\int_{4a}^{\infty}\frac{P_b(\sigma)\,d\sigma}{\sigma-s-i0}.
$$
The real upper rim uses a principal value within the cut. At large $\sigma$, $P_b(\sigma)=O(\log\sigma/\sigma)$, so the upper endpoint converges. For $b=a$ and $s\ge4a$ the exact continuation is $T_a(s)=2[\pi/2+i\,\mathrm{arcosh}\sqrt{s/(4a)}]^2/s$.

The momentum-independent bubble has the ultraviolet pole $+1/\bar\epsilon$ with
$$
 B_b=\frac1{\bar\epsilon}
      -\int_0^1 dx\,\log\frac{ax+b(1-x)}{\mu^2}.
$$
A nonzero $C_7\sum_i w_i$ therefore requires a local four-lepton counterterm. The $C_5$ triangle is superficially convergent and has no such pole. Independent finite Hermitian local four-lepton amplitudes contribute to the dispersive part even when the pole cancels. To predict an absolute rate they require matching. Changing $C_5$ changes an on-shell scalar-pair tree and hence the cut. If the input coefficients are assigned $C_5=O(\Lambda^{-1})$ and $C_7=O(\Lambda^{-3})$, the reciprocal rates below refer to the square of the declared amplitude. Interpreting all three squared and interference terms as a consistent prediction to a common order in $1/\Lambda$ requires control of every amplitude that contributes at that order.

## Stable two-mass benchmark

Take $m_i=(\sqrt a,4\sqrt a)$ and $w_i=(w_0,-w_0)$ with $w_0\ne0$. The specified Yukawas give no scalar decay phase space for these masses. Both the pole and $\mu$ dependence cancel in the primitive. Define $\Delta F=F_a-F_{16a}$. The residual local kernel is $\Delta B=(16/15)\log16-1$. The transfer-gap identity
$$
 \Delta P=\frac1s\log\left[1+
    \frac{2v(15a)}{(a+a^2/(h+v))(16a+h+v)}\right]
$$
uses $h^2-v^2=a^2$ and gives $\Delta P>0$ at every finite $s>4a$. Just above threshold $a\Delta P=(15/34)\beta+O(\beta^3)$. At $s=8a$ the two upper-rim kernels and the spectral gap are

| Quantity | Value |
| --- | ---: |
| $a\Delta K_5$ | $-0.273493435030441184-0.574435895925932904\,i$ |
| $\Delta K_7$ (absolute primitive) | $2.777908275480423539+1.723307687777798713\,i$ |
| $a\Delta P$ | $0.182848624652067526$ |
| $\Delta B$ | $1.957427970389099987$ |

At any selected $s_*>4a$, setting $C_5=(s_*/2-a)C_7$ cancels the two-scalar cut. The absolute loop kernel then retains $C_7 w_0\Delta B$ after stripping $W/(8\pi^2)$. A finite matched four-lepton contact can change this remaining amplitude. At other energies (1) is proportional to $(s-s_*)C_7\Delta P(s)/2$. For $g_7=a C_7w_0$, $g_5=C_5w_0$, and $\xi=s/a$, the established reciprocal identity gives
$$
 a\bar\sigma(s)\ \ge\ \frac{\xi}{2048\pi^3}(a\Delta P)^2
                \left|(\xi/2-1)g_7-g_5\right|^2.
$$
The squared rates describe this specified flavor sector. Additional input interactions that modify the scalar-pair subamplitudes can change the spectral function. [Reciprocal rate and phase space](../on_shell/INVERSE_RATE_BOUND.md)

Run `python -I -B -O verify_lower_operator.py` in this directory. It evaluates the absolute triangle in parameter space and by unsubtracted dispersion, sews the two on-shell trees, and checks the tuned crossing and the complex reciprocal amplitudes. [Verifier](verify_lower_operator.py) · [Numerical results](LOWER_OPERATOR_RESULTS.json)
