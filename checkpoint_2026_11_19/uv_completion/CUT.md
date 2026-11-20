# Physical scalar-pair cut of the renormalizable boxes

At leading order in the specified Yukawa interactions, the polarized process
$\mu^-_L\mu^-_L\to e^-_Le^-_L$ has a nonzero absorptive $J=0$ amplitude
at every finite $s>4m^2$. At $M=10m$ and $s=8m^2$ it implies
$m^2\sigma_{\mu\mu\to ee}\geq 1.938734024\times10^{-11}|\gamma|^2$.
Here $\gamma=(Y_{\mu 1}^*)^2(Y_{e1})^2$, with each index 1 labeling the first
sterile in its own flavor sector. The result follows from an exact cut with
both electron propagators retained. [Numerical receipt](CUT.json)

## Flavor and topology

Work at zero background in the gaugeless, unbroken theory, with massless
charged left-handed leptons and a charged scalar of squared mass $a=m^2>0$.
Use positive Majorana masses $(m,4m)$ for the muon block and $(M,2M)$ for
the electron block. Each sterile couples only to its designated flavor.
For incoming leptons write $y_\ell=(Y_{\ell 1}^*)^2$. Relative to $y_\ell$,
the muon squared Yukawas are $(1,-1/4)$ and the electron squared Yukawas are
$(1,-2)$. The corresponding chirality-flipping residues divided by $y_\ell$
are $w_{\mu}=(m,-m)$ and $w_e=(M,-4M)$. The outgoing electron tree conjugates
$y_e$, giving $\gamma=y_\mu y_e^*$. All residues left after this factorization
are real. The normalization follows the ordered spinor convention
$W=-B_{34}^*B_{12}$, $B_{ij}=u_i^T\epsilon u_j$, $|W|=s$.
[External states and sewing conventions](../on_shell/INDEPENDENT_CUT.md)

At degree two in each flavor's Yukawas every active slot is external.
The two muon vertices join by one anomalous muon-block propagator, and the
electron vertices join by one anomalous electron-block propagator. The
two scalar fields on the muon side have the same charge and must join the
conjugate fields on the electron side. Their two pairings give the direct
and crossed scalar boxes. Every connected graph has these four internal
lines. The mass numerators are independent of loop momentum, so each scalar
box integral falls as $d^4k/k^8$ in the ultraviolet. There is no lower-order
flavor-changing four-lepton tree for an external leg renormalization to dress.

The partition into the two incoming muons and two outgoing electrons crosses
precisely the scalar pair. A partition cutting the two sterile lines corresponds
to crossed external kinematics, with momentum transfer $t$ or $u\leq0$.
Partitions isolating one external lepton carry invariant zero. Equivalently,
at this order the two tree amplitudes with the same on-shell intermediate
particles on the two sides of unitarity exist only for that scalar pair.
In particular a mixed $\mu e$ sterile-exchange tree would require one
sterile coupling to both flavors. For $M=10m$ the scalar has no open
lepton–sterile two-body decay, including the zero-phase-space $m$ mode.
The cut therefore gives the complete physical $s$-channel discontinuity of
the specified leading amplitude. All exchanged tree transfers on the cut
are spacelike for finite $s$.

## Exact angle dependence

Let $\beta=\sqrt{1-4a/s}$, $h=s/2-a$, $v=s\beta/2$, and
$q(x)=h-vx$, with $x=\hat{\boldsymbol p}_1\cdot\hat{\boldsymbol k}_1$.
Then $t_1=-q(x)$, $u_1=-q(-x)$ and $q(x)\geq h-v=a^2/(h+v)>0$.
For $y=\hat{\boldsymbol p}_3\cdot\hat{\boldsymbol k}_1$ use the same function
$q(y)$ on the electron side. Define
$$
\begin{aligned}
 F_\mu(q)&=m\left(\frac1{a+q}-\frac1{16a+q}\right),\\
 F_e(q)&=M\left(\frac1{M^2+q}-\frac4{4M^2+q}\right)
       =\frac{-3Mq}{(M^2+q)(4M^2+q)},\\
 S_\ell(x)&=F_\ell(q(x))+F_\ell(q(-x)).
\end{aligned}
$$
The two incoming trees to the scalar pair are $-B_{12}y_\mu S_\mu(x)$
and $-B_{34}y_e S_e(y)$. The plus between transfers includes both boson
assignments after the antisymmetric spinor product has been factored out.
The two-body ordered measure is $d\Phi_2=\beta\,d\Omega/(32\pi^2)$.
Unitarity contributes $1/2$ and completeness contributes $1/2!$ for the
identical intermediate scalar pair. For $c=\hat{\boldsymbol p}_1\cdot
\hat{\boldsymbol p}_3$ the resulting cut is
$$
 A(s,c)=\frac{\mathcal M_{fi}-\mathcal M_{if}^*}{2i}
 =-\frac{W\gamma\beta}{32\pi}
   \int\frac{d\Omega_{\boldsymbol k}}{4\pi}\,S_\mu(x)S_e(y).                 \tag{1}
$$
The inverse upper-rim matrix element in this definition swaps the complete
external states. Each ordering is counted once. There is no extra Majorana
loop factor. Equation (1) also fixes the sign at finite angle in a chosen
external spinor phase convention.

For an alternative one-dimensional evaluation of the angular integral, let
$b,d$ be squared muon and electron sterile masses. Combine the two distinct
denominators with a parameter $u$. The sphere average of
$(D-v\boldsymbol n\cdot\boldsymbol N)^{-2}$ is
$[D^2-v^2|\boldsymbol N|^2]^{-1}$, which gives
$$
 I_{bd}(c)=\int_0^1\frac{du}
 {(h+ub+(1-u)d)^2-v^2[1-2u(1-u)(1-c)]}.
$$
The integral in (1) equals
$2\sum_{i,j} w_{\mu i}w_{ej}[I_{b_i d_j}(c)+I_{b_i d_j}(-c)]$.
The factor two here arises from reversing the intermediate scalar direction
within each pairing. At $W=-s$, $\gamma=1$, $M/m=10$, the full $c=0$ cut
is $-6.506712938\times10^{-5}$ for $s=6m^2$ and
$-1.875707071\times10^{-4}$ for $s=10m^2$.

## Angular momentum and the directional cross section

Choose the scattering-plane helicity phases so $W=-s$, or divide out $W$
before projecting and restore it afterward. Initial and final helicity
differences are zero, so the $J=0$ coefficient in the angular-average
normalization is the uniform average over $c$. In (1) first average the
outgoing axis at fixed scalar direction. This yields
$\int_{-1}^1 S_e(y)\,dy/2=2\overline F_e$, which is independent of
the scalar direction. Averaging the muon tree then gives $2\overline F_\mu$.
Explicitly,
$$
 \begin{aligned}
 J_b&=\frac1{2v}\log\frac{b+h+v}{b+h-v},&
 \overline F_\ell&=\sum_i w_{\ell i}J_{m_{\ell i}^2},\\
 A_0(s)&=-\frac{W\gamma\beta}{8\pi}\overline F_\mu\,\overline F_e .       \tag{2}
 \end{aligned}
$$
Both $F$ functions have mass dimension $-1$. Since $F_\mu(q)>0$ and
$F_e(q)<0$ for $q>0$, the two averages have strict opposite signs for
every finite $s>4a$. The finite electron masses produce higher even
partial waves in (1). The projection in (2) captures the part needed for
the following floor.

In this model all boxes have one common complex factor $W\gamma$.
The direct and inverse reduced upper-rim amplitudes consequently have the
same scalar loop functions with real relative coefficients. Writing
$\mathcal M_{fi}=W\gamma(H+iD)$ with real $H,D$ gives
$A=W\gamma D$ and $|\mathcal M_{fi}|^2=|W\gamma|^2(H^2+D^2)$ at each angle.
The identical final electrons contribute another $1/2!$, the massless
flux is $2s$, and the specified incident helicities require no spin average
or initial identical-particle factorial. Therefore
$$
 \sigma_{\mu\mu\to ee}
 =\frac{\langle|\mathcal M_{fi}|^2\rangle}{32\pi s}
 \ \geq\ \frac{\langle|A|^2\rangle}{32\pi s}
 \ \geq\ \frac{|A_0|^2}{32\pi s}
 =\frac{s\beta^2}{2048\pi^3}|\gamma|^2\overline F_\mu^2\overline F_e^2 .
                                                                    \tag{3}
$$
The second inequality is angular Cauchy–Schwarz. The same bound holds
in the inverse direction. With independent weak phases in extra interactions,
only the reciprocal-mean identity follows from the discontinuity alone.
This is a bound on rates computed from the stated leading one-loop amplitude.
[Reciprocal identity and flux](../on_shell/INVERSE_RATE_BOUND.md)

At threshold $q\to a$, $m\overline F_\mu\to15/34$ and, for $M=10m$,
$m\overline F_e\to-30/(101\cdot401)$. The $J=0$ floor in
$m^2\sigma/|\gamma|^2$ opens as $6.726917812\times10^{-12}\beta^2$.
The checker constructs the two outgoing scalar four-vectors, sews the
raw propagator denominators on the sphere, projects over the external
angle, and compares the result with (2). It also compares the right-angle
cut with the one-parameter integral above and verifies nonzero higher
waves. Run `python -I -B -O check_cut.py` to print its JSON or supply
`--output PATH` to write it. [Checker](check_cut.py) · [Receipt](CUT.json)
