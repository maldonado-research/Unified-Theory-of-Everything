# The finite box at fixed spacelike transfer

At leading order in the four Yukawa vertices this model fixes the real and imaginary parts of the charged four-lepton amplitude. Retaining the heavy propagators also determines the analytic contribution that was free in the earlier contact model. The dispersion calculation below holds the external transfer fixed while integrating over the scalar-pair invariant. The crossed box is added at the physical point only after that integral.

## Which integral the action gives

Use the masses, spinor convention and flavor blocks in [tree matching](MATCHING.md) and the [cut calculation](CUT.md). Put $a=m^2$, $\gamma=Y_{e1}^2(Y_{\mu1}^*)^2$ with the first column in each flavor block understood, and choose external spin phases with $B_{34}^*B_{12}=s$. Define $b_i=(m^2,16m^2)$, $d_j=(M^2,4M^2)$, $w_i=(m,-m)$ and $v_j=(M,-4M)$. For massless external states with scattering cosine $c$, the external transfers are $t=-s(1-c)/2$ and $u=-s(1+c)/2$. The two internal lepton-number-changing propagators each contribute their mass. The two scalar pairings give
$$
\mathcal M(s,c)=\frac{s\gamma}{16\pi^2}\sum_{i,j} w_i v_j
 [D(s,t;b_i,d_j)+D(s,u;b_i,d_j)] .
$$

Here $D$ has mass dimension $-4$. With $x_k\ge0$, $[dx]=\prod_{k=1}^4dx_k\,\delta(1-\sum x_k)$, its Feynman definition is
$$
D(s,t;b,d)=\int \frac{[dx]}
 {[a(x_1+x_2)+b x_3+d x_4-s x_1 x_2-t x_3x_4-i0]^2} .
$$
The scalar denominators are at loop momenta $k,P-k$; the two sterile denominators are at $k-p_1,k-p_3$. These differences of reference momenta directly give the displayed products $s x_1x_2$ and $t x_3x_4$. The factor $1/(16\pi^2)$ is the four-dimensional scalar loop normalization. Sewing the two trees fixes the overall sign in the declared phases. Earlier electroweak work gives the same mass numerators and $D_0(s,t)+D_0(s,u)$ structure for the charged-Goldstone box in the reverse process. Its full electroweak answer also includes gauge bosons. [Cannoni, Kolb and Panella, Eq. (4)](https://arxiv.org/pdf/hep-ph/0209120v2)

## Positive spectral representation for one mass assignment

Set $x_3+x_4=v$, $x_3=v\lambda$, $x_1=(1-v)y$, and $w=v/(1-v)$. Write
$$
 b_\lambda=\lambda b+(1-\lambda)d,\quad
 d_y=y(1-y),\quad A=a+b_\lambda,\quad B=b_\lambda-t\lambda(1-\lambda).
$$
At $t\le0$ both $A$ and $B$ are strictly positive. Begin at $s<0$. The transformed integral is $\int d\lambda\,dy\,w\,dw/Q^2$ with $Q=a+Aw+Bw^2-s d_y$, $0\le w<\infty$, and $0\le\lambda,y\le1$. Integration by parts gives
$$
 \int_0^\infty\frac{w\,dw}{Q^2}=
 \int_0^\infty \frac{A\,dw}{(A+2Bw)^2 Q}.
$$
The boundary term $w/[(A+2Bw)Q]$ vanishes at both limits. Change variables to $s'=(a+Aw+Bw^2)/d_y$. The weight per $dy\,d\lambda$ becomes $A/[A^2+4B(s'd_y-a)]^{3/2}$ when $s'd_y\ge a$. Substituting $y=(1+x\sqrt{1-4a/s'})/2$ and integrating from $x=-1$ to $1$ yields
$$
 \begin{aligned}
 D(s+i0,t;b,d)&=\int_{4a}^{\infty}\!\frac{\rho(s',t)\,ds'}{s'-s-i0},\\
 \rho(s',t)&=\sqrt{1-4a/s'}\int_0^1\frac{d\lambda}{A^2+B(s'-4a)} .
 \end{aligned}
$$
This identity follows from the convergent Euclidean integral and analytic continuation. It fixes the dispersive part, including its additive normalization. Its imaginary part $+\pi\rho(s,t)$ reproduces the one-parameter average of the two exchanged propagators on the cut. The signed flavor sum can have either sign despite positivity of each mass assignment. For fixed negative $t$ this derivation also applies to values of $s'$ for which $1+2t/s'$ lies outside the physical angular interval.

## Computing on the upper rim

Define $\beta=\sqrt{1-4a/s}$ for $s>4a$, $L=A^2$, $H=4aB-L$. Substituting $q=\sqrt{1-4a/s'}$ in the spectral integral and taking the real principal value gives
$$
 \begin{aligned}
 D(s+i0,t;b,d)&=\frac{8a}{s}\int_0^1 d\lambda\,
 \frac{-\beta\operatorname{artanh}\beta+i\pi\beta/2 + E(H/L)}
      {L+H\beta^2},\\
 E(x)&=\int_0^1\frac{dq}{1+xq^2}.
 \end{aligned}
$$
For $x>0$, $E(x)=\arctan(\sqrt{x})/\sqrt{x}$. For $-1<x<0$, $E(x)=\operatorname{artanh}(\sqrt{-x})/\sqrt{-x}$. At $x=0$ its value is one. This evaluates the timelike pole analytically, with one smooth mass-parameter integral left per mass assignment. The program uses a logarithmic mass coordinate and repeats the external-angle and mass integrals at two quadrature orders.

For fixed incident helicities and identical final electrons, $\sigma=\langle |\mathcal M|^2\rangle_c/(32\pi s)$, with $\langle\cdot\rangle_c=\int_{-1}^1dc/2$. The scaled results for $M=10m$ are given below. Each rate column is $m^2\sigma/|\gamma|^2$. The absorptive rate uses the squared imaginary part after removing the common phase. Its $J=0$ part uses the angular average before squaring.

| $s/m^2$ | $\mathop{\rm Re}[\mathcal M(s,0)/\gamma]$ | $\mathop{\rm Im}[\mathcal M(s,0)/\gamma]$ | Full | Absorptive | Absorptive $J=0$ |
|---:|---:|---:|---:|---:|---:|
| 6 | $-1.4143753842\times10^{-4}$ | $-6.5067129383\times10^{-5}$ | $4.0092537001\times10^{-11}$ | $7.0135499582\times10^{-12}$ | $7.0135491265\times10^{-12}$ |
| 10 | $-1.9338112607\times10^{-4}$ | $-1.8757070705\times10^{-4}$ | $7.1740070858\times10^{-11}$ | $3.4842174581\times10^{-11}$ | $3.4842037117\times10^{-11}$ |

The dispersive calculation reproduces the angularly projected tree cut. Direct integration of the Feynman simplex checks both amplitude components and the full angular integral for each energy. [Dispersion program](check_bounds.py) · [Dispersion results](BOUNDS.json) · [Independent direct integration](BOX_CROSSCHECK.md) · [External QCDLoop comparison](QCDLOOP.md)

## A uniform error bound on the absorptive amplitude

On the two-scalar cut let $z=-q_{\rm transfer}^2>0$, $x=z/M^2$. For $4<s/m^2\le11$ and $M=10m$, kinematics gives $z<(s-2m^2)\le9m^2$. The exact electron kernel divided by its dimension-seven approximation on each route is $R(x)=1/[(1+x)(1+x/4)]$. Its dimension-nine approximation to this ratio is $1-5x/4$. Their relative errors, measured against the **exact** kernel, obey
$$
 \frac1{R(x)}-1=\frac{5x+x^2}{4}\le 0.114525,\qquad
 1-\frac{1-5x/4}{R(x)}=\frac{21x^2+5x^3}{16}\le 0.0108590625 .
$$
Every transfer on the electron side has the same sign. Both transfers of the muon tree have the opposite, fixed sign. The sewn angular integral therefore inherits these error bounds at each external angle and in its $J=0$ average. The bound controls the physical discontinuity in this interval. The real part samples the spectral density out to infinite energy and has its own matching contribution. [Exact tree factor](MATCHING.md) · [Rational bounds in the numerical receipt](BOUNDS.json)
