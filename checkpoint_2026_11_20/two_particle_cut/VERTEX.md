# Proper Yukawa vertices on the charged two-scalar cut

The finite vertex diagrams give positive contributions to both Yukawa coefficients of the two-scalar Hermitian cut, in the external phase $B_{34}^*B_{12}=s$. At $s=8m^2$ they contribute $1.58008833\times10^{-6}\,\gamma |g|^2+2.13277672\times10^{-6}\,\gamma |y|^2$ to the angular mean. Here $g$ and $y$ are the first muon and electron Yukawas and $\gamma=(g^*)^2y^2$. The numbers include both endpoints on the corrected tree, both identical scalar assignments, and the intermediate scalar symmetry factor. [Benchmark action and tree normalization](../../checkpoint_2026_11_19/uv_completion/MATCHING.md) · [Numerical receipt](VERTEX.json)

## Proper three-point graph

Take the incident lepton momentum $p$, the outgoing scalar momentum $k$, $Q=p-k$, $t=Q^2=-q<0$, $p^2=0$ and $k^2=a=m^2$. In one flavor block let $b_i=M_i^2$, $r_i=(Y_i^*/Y_1^*)^2$, and $w_i=M_i r_i$. The two blocks have $r_\mu=(1,-1/4)$, $r_e=(1,-2)$. With the weak phase removed, the tree per assignment is $-B_{12} F(q)$, where $F=\sum_i w_i/(b_i+q)$. Define $G(q)=\sum_i |r_i|/(b_i+q)$ for a normal sterile propagator. The coefficient $V$ in $F+|Y_1|^2 V$ from the proper vertices is

$$
 V(q)=-\frac{2q}{16\pi^2}\,G(q)\sum_j w_j I_j(q),\qquad
 I_j(q)=\int_{x,y\ge0,\,x+y\le1}
 \frac{(1-z)\,dx\,dy}{xa+yb_j+xzq-yza},\qquad z=1-x-y .
$$

To obtain the numerator, route the scalar loop momentum $l$ from the vertex bearing $p$ to the conjugate Yukawa vertex. The internal sterile carries $p-l$, the internal active lepton carries $Q-l$, and the sterile bridge carries $Q$. The sterile between the two like Yukawa vertices supplies $M_j$. In a continuous fermion-flow orientation the chain acting on the incident spinor is $M_j\not Q(\not Q-\not l)u_L(p)$. Assign parameters $x,y,z$ respectively to $l^2-a$, $(p-l)^2-b_j$, and $(Q-l)^2$. Shifting $l$ by $yp+zQ$ and using the massless Dirac equation leaves $M_jt(1-z)u_L(p)$. The product of three triangle propagators, the two extra Yukawa vertex factors relative to tree level, and the Minkowski integration $-i/(16\pi^2)$ for the positive simplex give $+1/(16\pi^2)$ relative to the tree numerator. More explicitly, the ordered tree contains $(-i)^2i=-i$, and the new graph contains $(-i)^4i^4=1$ before its triangle integration. The external spin ordering is the same in this comparison. The high-momentum numerator grows at most linearly in $l$ for fixed $Q$.

The triangle sterile joins two like chiralities, giving $(Y_j^*)^2$ at its ends. The other sterile joins a vertex to its conjugate and gives $|Y_i|^2$. The reversed fermion flow at the opposite endpoint transposes the spin chain. Each slashed momentum changes sign under $C\gamma^{\mu T}C^{-1}=-\gamma^\mu$, leaving the scalar projection $t(1-z)$ unchanged. Exchanging the incident labels flips the ordered spin determinant and the external fermion order together. Thus the two possible triangle locations add with factor two for a given transfer. The scalar assignments $t,u$ are subsequently summed as at tree level. Both loop particles that carry weak isospin are fixed by the two external charged components, giving one charged scalar and one charged active lepton in this triangle.

For a direct graph count, on the $LLHH$ side the four-vertex Yukawa term contains three $LHN$ vertices and one conjugate vertex. The conjugate vertex has three internal legs. When these land at three distinct like vertices, the vertex joined by the normal sterile line carries both an external lepton and a scalar. The two remaining vertices exchange an anomalous sterile, making precisely the triangle and its bridge above. If the scalar and lepton arrive at the same vertex, their bubble dresses the sterile bridge. Coincidence with the sterile attachment instead dresses an external active or scalar leg. Coincidence of all three attachments disconnects the graph. This counts every external assignment, including the graphs which in the crossed channel are corrections to the opposite endpoint. [Prior explicit Wick-contraction enumeration](../../checkpoint_2026_11_19/radiative_stability/MATCHING_AUDIT.md)

## Threshold and independent checks

Set $x=(1-z)(1-h)$ and $y=(1-z)h$. The positive denominator remaining after the Jacobian cancellation is $(1-h)(a+zq)+h(b_j-za)$. Integrating over $h$ gives
$$
 I_j(q)=\int_0^1 (1-z)\,
     \frac{\ln[(b_j-za)/(a+zq)]}{b_j-za-(a+zq)}\,dz .
$$
Coincident arguments give the reciprocal of their common value. For $b_j=a$ the $z=1$ endpoint behaves as $(1-z)\ln[1/(1-z)]/(a+q)$ and the upper-rim limit is real and finite. On the physical scalar-pair cut $q$ stays strictly positive. The explicit factor $q$ makes the proper vertex vanish in the zero-transfer limit.

Two independent discontinuities check the triangle.

- **Cut the internal lepton and sterile.** Continue the scalar invariant to $k^2=U+i0$ at fixed $q$, keeping the internal scalar mass $a$. Then $b_j-za$ in the one-dimensional formula becomes $b_j-zU-i0$ for positive $z$. Above $U=b_j$ the spectral density is $\rho_j=\mathrm{Im}\,I_j/\pi=\int_{b_j/U}^{1}(1-z)dz/[a-b_j+z(q+U)]$. The limit at large $|U|$ supports an unsubtracted dispersion integral $I_j(q)=\int_{b_j}^\infty\rho_j(q,U)dU/(U-a)$. In units $a=b_1=1$, this gives $\rho_1=(U-1)D(U)/[U(U+q)]$ with $D=U\ln U/(U-1)-1$. The spectral weight of a scalar propagator carrying the single endpoint correction is $-\mathrm{Im}\,V_{\rm one}/[\pi(U-1)]$. Its lightest term is $q G(q)D(U)/[16\pi^2 U(U+q)]$, exactly the signed exchange contribution per scalar placement in the previously calculated three-particle cut. The checker integrates the dispersion relation for all four masses at all three scan energies. [Three-particle calculation and independent helicity-amplitude comparison](../../checkpoint_2026_11_19/yukawa_continuum/FULL_THREE_BODY.md)
- **Cut the internal active lepton and scalar.** With the internal and external scalar masses set to zero, continue the sterile transfer to $t>0$. On its upper rim, $\mathrm{Im}\,I_j=\pi\int_0^1(1-z)dz/(b_j+tz)=(\pi/t)[(1+b_j/t)\ln(1+t/b_j)-1]$. Buchmüller and Plümacher use $b_k=-I_k/(16\pi^2)$, in the convention of their scalar loop functions. This gives their absorptive proper form factor, Eqs. (17)–(18), including its sign. Their Eqs. (57)–(58) insert it as $-M b K M$ or its weak-conjugate version in an on-shell decay coupling. Our projected endpoint likewise has the positive relative coefficient $t I_j/(16\pi^2)$ before continuation to $t=-q$. The checker also verifies $I_j(q=b_j)=1/b_j$ in this massless-scalar limit. [Buchmüller and Plümacher, Eqs. (13)–(18) and (57)–(58)](https://arxiv.org/abs/hep-ph/9710460)

## Angular mean

Write $\beta=\sqrt{1-4a/s}$, $q(c)=s/2-a-s\beta c/2$ and $\bar F=\int_{-1}^1F(q(c))\,dc/2$, and define $\bar V$ the same way. The $J=0$ projection averages each production axis independently. Each tree has two scalar assignments. The ordered two-particle measure, intermediate $1/2!$ and unitarity $1/2$ then give $A_{2,\mu}^{V}=s\beta\bar V_\mu\bar F_e/(8\pi)$ and $A_{2,e}^{V}=s\beta\bar F_\mu\bar V_e/(8\pi)$. The subscripts specify which flavor has the extra squared Yukawa magnitude stripped off. [Unitarity and spin conventions](../../checkpoint_2026_11_19/uv_completion/CUT.md)

| $s/m^2$ | $A_{2,\mu}^{V}$ | $A_{2,e}^{V}$ |
|---:|---:|---:|
| 6 | $+6.46226760331\times10^{-7}$ | $+1.10989313869\times10^{-6}$ |
| 8 | $+1.58008832854\times10^{-6}$ | $+2.13277671877\times10^{-6}$ |
| 10 | $+2.81830430106\times10^{-6}$ | $+3.20052145882\times10^{-6}$ |

At $s\downarrow4a$, $q\to a$ and the projected vertex coefficients open as $3.49057758\times10^{-7}\beta$ and $8.87028602\times10^{-7}\beta$, respectively. Run `python -I -B -O check_vertex.py` to reproduce the table. [Checker](check_vertex.py) · [Receipt](VERTEX.json)


## External scalar-integral comparison

QCDLoop 2.1.0 evaluates the scalar functions after tensor reduction. With the momenta defined above let $d_0=l^2-a$, $d_1=(p-l)^2-b$, $d_2=(Q-l)^2$, and write $C_0$ for the integral of their inverse product in the library's convention. Contracting with the massless external momentum gives $2p\cdot l=d_0-d_1+a-b$. In the shifted integral the same contraction gives $z(t-a)$ inside the scalar integral. The production vertex integral is consequently

$$
 I(q,b)=\frac{B_0(a;b,0)-B_0(t;a,0)+(2a-b-t)C_0(0,a,t;a,b,0)}{t-a},\qquad t=-q.
$$

The two bubble poles cancel. In quadruple precision the library returns a finite real result at each of sixteen points with $a=1$, $b=1,16,100,400$ and $q=1/8,1,4,9$. The largest relative difference from the direct parameter integral with 64 nodes was $6.7\times10^{-15}$. At the lightest threshold point $I(1,1)=0.5405533695762245$. The combined verifier compares its direct integral to each archived value. [External program](vertex_qcdloop.cpp) · [Reference values](QCDLOOP_VERTEX.json) · [Combined checker](check_absorptive.py)

To regenerate the library values, use the QCDLoop source and generated configuration header described in the [pinned build instructions](../../checkpoint_2026_11_19/uv_completion/QCDLOOP.md). From the parent of the QCDLoop source directory, after copying the external program there, run

```bash
g++ -std=gnu++11 -fext-numeric-literals -O2 -Ibuild -Iqcdloop/src \
  vertex_qcdloop.cpp qcdloop/src/{triangle,bubble,cache,topology,tools,types}.cc \
  -lquadmath -o build/vertex_qcdloop
./build/vertex_qcdloop
```
