# One external weak gauge field and the neutral photon pole

15 October 2026. Restricted continuation of the literal action and source
orientation in `checkpoint_2026_10_15/operator_map/ACTION_SIGN.md` and the
triangle routing in `checkpoint_2026_10_15/GAUGE_LEG_DERIVATION.md`.
This calculation determines a selected simple UV pole. It makes no claim
about a complete anomalous dimension, finite amplitude, or decay rate.

## Result

Let `phi=tilde H`, `chi=L^c`, `T^I=tau^I/2`, and
`F=C_(7,seed)Y_nu*`. Restore the scalar wavefunction `phi_c` and the
outgoing lepton weak index `a`. The two internal one-W attachments give

\[
\boxed{
i\mathcal M_{W^I}^{\rm UV}
=-\frac{ig}{16\pi^2\epsilon}F_{\alpha i}
\bar u_{L\alpha}^{\,a}\,(T^I)_{ac}R^\nu P_Lu_{Ni}
\phi_c\varepsilon_\nu^*,\qquad R^\nu=G_J^\nu+D^\nu .}
\]

Here `R` is exactly the reduced one-B polynomial already found, with
`p=r+h+k` and the same gamma ordering:

\[
G_J^\nu=\left[-\frac{(2r+k)^\nu}{3}
 +\frac{(2h+k)^\nu}{6}\right]\not p
 +\left(\frac{p^2}{6}-\frac{m_H^2}{2}\right)\gamma^\nu,
\]
\[
D^\nu=\frac14[\not k,\gamma^\nu]\not r
 +\frac13(k^\nu\not k-k^2\gamma^\nu).
\]

This equality is derived by keeping the two seed pairings and integrating
the crossed triangles. It does not assume the replacement `g'y_L B ->
gT^I W^I`. Tracelessness removes the direct pairing separately on each
SU(2) attachment, but the direct pairing had already cancelled between
the two hypercharge UV residues. Thus no new transverse coefficient
appears at this order.

For a constant neutral Higgs background `phi=(v/sqrt(2),0)^T` and an
external neutrino, the selected photon UV pole is zero:

\[
g'c_Wy_L+gs_WT^3_{11}
=e\left(-\frac12+\frac12\right)=0.
\]

The conclusion applies to the entire local pole, including `D`, and not
just its longitudinal contraction. It does **not** state that the finite
neutral photon vertex or physical `N -> nu gamma` amplitude vanishes.

## Literal weak contractions

Keep the literal symmetric-flavor seed

\[
\frac12C_{\alpha\beta}\bar X_\alpha(-\Box)X_\beta^c
+\mathrm{h.c.},\qquad X=\phi^\dagger L.
\]

In the `phi` basis its two scalar contractions, after the established
fermion-slot exchange removes the factor `1/2`, are

\[
K^d_{ab;cd}=\delta_{ac}\delta_{bd},\qquad
K^x_{ab;cd}=\delta_{ad}\delta_{bc}.
\]

`a,c` are the external `bar L,phi` indices; `b,d` are the seed's internal
`chi,phi` indices. At the source Yukawa vertex the weak contraction is
`delta_ef`, where `e` labels `bar chi` and `f` labels `phi*`. The matrices
in the kinetic vertices are

\[
V_\chi^I=-ig(-T^{I\,T})\gamma^\nu P_R,
\qquad V_\phi^I=-igT^I(s+s')^\nu.
\]

The minus and transpose in the conjugate fermion representation follow
from charge conjugation of the kinetic term. They are additional to the
common vertex phase `-i`, and must not be counted twice. The two exact
weak contractions are

\[
\begin{aligned}
G_L(K)_{ac}
 &=\sum_{bdef}K_{ab;cd}(-T^{I\,T})_{be}
    \delta_{df}\delta_{ef}
 =-\sum_{bd}K_{ab;cd}T^I_{db},\\
G_H(K)_{ac}
 &=\sum_{bdef}K_{ab;cd}\delta_{be}T^I_{df}\delta_{ef}
 =+\sum_{bd}K_{ab;cd}T^I_{db}.
\end{aligned}
\]

They give the following table before the momentum integrals are summed:

| Seed pairing | Fermion attachment | Scalar attachment |
| --- | --- | --- |
| Direct | `-delta_ac tr(T^I)=0` | `+delta_ac tr(T^I)=0` |
| Crossed | `-T^I_ac` | `+T^I_ac` |

The generator normalization is `tr(T^IT^J)=delta_IJ/2`. In particular,
`T^2` is imaginary; replacing `-T^T` with `-T` would change its contraction.

For hypercharge, take the scalar generator `M=y_L 1_2` and the conjugate
fermion generator `-M^T`. The same contractions give direct factors
`-2y_L,+2y_L` and crossed factors `-y_L,+y_L`. This recovers the previously
audited direct-to-crossed weight `2:1` without changing flavor order.

## Original epsilon basis

The scalar change of basis is `phi=epsilon H*`,

\[
\epsilon=\begin{pmatrix}0&1\\-1&0\end{pmatrix},
\qquad T^I\epsilon=-\epsilon T^{I\,T}.
\]

The seed tensor becomes
`epsilon_ac epsilon_bd` or `epsilon_ad epsilon_bc`; the source Yukawa
tensor is `epsilon_ef`. Both `chi` and `H*` carry the SU(2) matrix
`-T^{I T}`. Contracting these tensors explicitly gives zero for either
direct pairing and `-(T^I epsilon)_ac` or `+(T^I epsilon)_ac` for the
fermion and scalar crossed pairings respectively. Thus the open external
tensor in the original `H*` variables is `T^I epsilon`, as required by
`bar L T^I tilde H`.

The transpose in `T^I epsilon=-epsilon T^{I T}` is essential. A uniform
replacement by `-epsilon T^I` fails for `I=2`. The checker evaluates all
internal indices in both scalar bases and tests this negative control.

## Separate momentum integrals and direct UV cancellation

Use `t=p-k`, `v=h+k`, `u=2p-k`, the metric `(+---)`, and strip
`i/(16 pi^2 epsilon)` from the loop residues. Define **unit-pairing**
integrals, keeping weak weights out of the numerators:

\[
\begin{aligned}
I_{L,d}^\nu&=\int_q
\frac{t^2(\not q-\not k)\gamma^\nu\not q}
 {q^2(q-k)^2((q-p)^2-m_H^2)},\\
I_{H,d}^\nu&=\int_q
\frac{t^2\not q(2p-2q-k)^\nu}
 {q^2((q-p)^2-m_H^2)((q-t)^2-m_H^2)},\\
I_{L,x}^\nu&=\int_q
\frac{(q-v)^2(\not q-\not k)\gamma^\nu\not q}
 {q^2(q-k)^2((q-p)^2-m_H^2)},\\
I_{H,x}^\nu&=\int_q
\frac{(q-h)^2\not q(2p-2q-k)^\nu}
 {q^2((q-p)^2-m_H^2)((q-t)^2-m_H^2)}.
\end{aligned}
\]

The usual `+i0` prescriptions are understood. All external momenta are
generic and off shell when the UV functional is defined.

The direct integrals only need the rank-two triangle pole
`Res(q_mu q_nu)=g_mu_nu/4`. Therefore

\[
R_{L,d}^\nu=\frac{t^2}{4}\gamma_\mu\gamma^\nu\gamma^\mu
=-\frac{t^2}{2}\gamma^\nu,
\qquad R_{H,d}^\nu=-\frac{2t^2}{4}\gamma^\nu
=-\frac{t^2}{2}\gamma^\nu.
\]

With the hypercharge trace weight two, these are both
`-t^2 gamma^nu`, exactly the parts removed when passing from the old
`R_L,R_H` to the crossed residues. Their difference is zero. Using
four-dimensional gamma contractions here fixes the simple pole only;
terms of order `epsilon` in these contractions can affect the finite
part. Equality of the two complete integrals is neither used nor claimed.

## Crossed residues from primitive tensor moments

For a triangle with denominator shifts `a,b`, the further moment is

\[
\operatorname{Res}(q^\mu q^\nu q^\rho)
=\frac{g^{\mu\nu}(a+b)^\rho+g^{\mu\rho}(a+b)^\nu
 +g^{\nu\rho}(a+b)^\mu}{12}.
\]

For a bubble `q^2-ma2`, `(q-P)^2-mb2`, its moments are

\[
B_0=1,\qquad B_1^\mu=\frac{P^\mu}{2},\qquad
B_2^{\mu\nu}=\frac{P^\mu P^\nu}{3}
 +\left(\frac{ma2+mb2}{4}-\frac{P^2}{12}\right)g^{\mu\nu}.
\]

Expand `(q-v)^2=q^2-2q.v+v^2` and cancel its `q^2` against the fermion
denominator. With `q=l+k`, the bubble shifts are `0,t`, with masses
`0,m_H^2`; its numerator is `/l gamma^nu (/l+/k)`. The remaining triangle
has shifts `k,p` and rank at most three. Applying these moments gives

\[
\begin{aligned}
R_{L,x}^\nu={}&\frac{\not t\gamma^\nu\not t}{3}
 +\frac{\not t\gamma^\nu\not k}{2}
 -\frac{\not v\gamma^\nu(\not p+\not k)
 +(\not p+\not k)\gamma^\nu\not v}{6}
 +\frac{\not k\gamma^\nu\not v}{2}\\
&+\left[\frac{t^2}{6}-\frac{v^2}{2}
 +\frac{v\cdot(p+k)}{3}-\frac{m_H^2}{2}\right]\gamma^\nu.
\end{aligned}
\]

For the scalar crossed numerator, cancel the `q^2` in `(q-h)^2`. Put
`q=l+p`; the bubble has equal masses and shift `-k`, while its numerator
is `(/l+/p)(-k-2l)^nu`. The remaining triangle has shift sum `u`. This
gives

\[
R_{H,x}^\nu=-\frac{k^\nu\not k}{6}
 -\frac{u^\nu\not h}{6}+\frac{h^\nu\not u}{3}
 +\left[\frac{k^2}{6}-m_H^2+\frac{h\cdot u}{3}
 -\frac{h^2}{2}\right]\gamma^\nu.
\]

Thus `R_L=2R_L,d+R_L,x` and `R_H=2R_H,d+R_H,x` reproduce the earlier
hypercharge formulas. Independently subtracting the newly integrated
crossed expressions gives

\[
R_{H,x}^\nu-R_{L,x}^\nu=G_J^\nu+D^\nu,
\qquad k_\nu D^\nu=0.
\]

The longitudinal part equals
`T_UV(p,h+k)-T_UV(p,h)`, with the same zero-gauge bubble as before.
Scalar-mass dependence lies entirely in `G_J`.

## Action sign, flavor, and a local one-W representative

The zero-gauge source chain has phase `(i)(i)(-i)(i)=-1`.
Each ordinary gauge insertion introduces one extra propagator and one
vertex, whose scalar phase is `i(-i)=+1`. Consequently the stripped
one-W graph is

\[
-gF\left[G_L(K^d)I_{L,d}+G_H(K^d)I_{H,d}
 +G_L(K^x)I_{L,x}+G_H(K^x)I_{H,x}\right]
=-gF T^I(I_{H,x}-I_{L,x}).
\]

The loop normalization supplies the remaining `i`, establishing the
negative pole in the result. The local counterterm action has the
positive coefficient. The order `F_{alpha i}=sum_beta C_{alpha beta}
Y^*_{beta i}` is fixed by the source contraction; weak matrices act in a
different space and do not commute the flavor factors. There is no
additional identical-field factor `1/2` in the mixed counterterm.

Let `O_J` now use the full electroweak covariant derivative. At one weak
gauge field a covariant representative is obtained by adding

\[
\frac{F}{16\pi^2\epsilon}\left[
\frac g4 i(D_\rho\bar L)\sigma^{\mu\nu}\gamma^\rho
 T^I\widetilde H N_L W^I_{\mu\nu}
+\frac g3\bar L\gamma_\mu T^I\widetilde H N_L
 (D_\nu W^{\nu\mu})^I\right]+\mathrm{h.c.}
\]

The respective linear-W kernels are `g T^I [slash k,gamma^nu]slash r/4`
and `g T^I(k^nu slash k-k^2 gamma^nu)/3`. At this order only the linear
parts of the field strengths are probed. Their nonlinear completion and
multi-W Green functions require their own comparison; writing this
covariant representative does not constitute that comparison.

## Neutral projection and its limits

Choose `A=s_W W^3+c_W B`, `Z=c_W W^3-s_W B`,
`e=gs_W=g'c_W`, and insert only the constant background
`tilde H=(v/sqrt(2),0)^T`. The remaining neutrino index is `a=1`.
The `B` and `W^3` residues carry exactly the same Lorentz matrix `R`,
including both transverse structures. Their photon coefficient therefore
vanishes, while the neutral-Z coefficient is

\[
-g's_Wy_L+gc_WT^3_{11}=\frac{e}{2s_Wc_W}.
\]

This is a projection of the computed symmetric-phase local pole; it
sets `h=0` for the constant scalar only after the UV polynomial exists.
It is not a broken-phase loop calculation with resummed masses or mixing.

The electromagnetic scalar generator is

\[
Q_\phi=T^3+y_L1_2=\operatorname{diag}(0,-1).
\]

Its neutral crossed entry is zero, but `tr Q_phi=-1` is nonzero. The
direct electromagnetic pairing therefore need not vanish as a complete
finite graph. Only its opposite-charge simple UV residues cancel. This
explicitly identifies why the pole cancellation cannot be promoted to a
claim about the finite vertex, a transition moment, or a decay rate.
The hypercharge-only result likewise cannot be read as a physical photon
prediction before the weak projection.

## Reproduction

Run `python3 -I -S -B -O verify_one_w.py` in this directory. The checker
uses standard-library exact Gaussian-rational matrices and imports no
repository code. Its **695 explicit checks** include all weak-index
contractions in both scalar bases; the conjugate representation and
generator normalization; actual direct and crossed tensor integrations
at 24 generic rational off-shell points; the local decomposition, Ward
identity, and scalar-mass dependence; action phases and complex flavor
order; and photon/Z projections. It rejects wrong conjugate charges,
missing transpose, commuted flavor order, omission of the transverse
pole, and treating the hypercharge pole as the photon pole.

These checks use the explicitly stated primitive bubble/triangle moments.
They do not derive those primitive moment identities again.
A separately authored checker at
`verify_su2_independent.py` instead expands both crossed
numerators at large loop momentum and uses angular moments. It reports
**1461 exact checks**, including exhaustive quadratic momentum
coefficient tests and scalar-mass coefficients. Its basis substitutions
apply to the UV polynomial defined at generic kinematics, not to an
unseparated scaleless integral.
