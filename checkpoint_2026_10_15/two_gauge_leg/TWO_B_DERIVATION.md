# Two external hypercharge fields: explicit one-loop UV match

15 October 2026. This evaluates the next 1PI Green function of the literal
seed sector fixed in `checkpoint_2026_10_15/operator_map/ACTION_SIGN.md`.
It is a restricted one-loop counterterm check, not an anomalous dimension
or a complete EFT calculation.

## Result and scope

Let `F=C_(7,seed) Y_nu*`, `y_L=-1/2`, `K=k_1+k_2`, and
`p=r+h+K`. With all external momenta initially generic and off shell, the
seven diagrams derived below give

\[
\boxed{R_2=(\varepsilon_1\cdot\varepsilon_2)\not p
-\frac14\left([\not k_1,\not\varepsilon_1]\not\varepsilon_2
 +[\not k_2,\not\varepsilon_2]\not\varepsilon_1\right).}
\]

Consequently, in the action, Fourier, and source-chain conventions already
fixed by the literal reference action,

\[
i\mathcal M_{BB}^{\rm UV}
=-\frac{iF(g'y_L)^2}{16\pi^2\epsilon}
\bar u_L R_2P_Lu_N.
\]

The external electroweak epsilon contraction is suppressed. The local
counterterm has the opposite amplitude, and its action kernel is
`+F(g'y_L)^2 R_2/(16 pi^2 epsilon)`. This agrees with the two-B Fourier
kernel of the previously chosen covariant action. The calculation below
integrates the actual loop diagrams; that action kernel is used only as
the comparison target after integration.

This coefficient is exactly linear in `C_(7,seed) Y_nu*` and quadratic in
`g'`. No internal gauge propagator, sterile propagator, extra Yukawa pair,
external-leg radiation, or additional operator insertion is included.
The two-gauge scalar seagull is included. No gauge vertex is assigned to
the neutral composite in the literal seed.

## Routing, vertices, and phases

Use `(+---)`, `d=4-2 epsilon`, and the UV integral normalization
`i/(16 pi^2 epsilon)`. The internal charged lines are
`chi=L^c` with charge `q_chi=-y_L` and `phi=tilde H` with charge
`q_phi=y_L`. Momentum `q` starts at the Yukawa vertex on the fermion line;
the scalar momentum there is `p-q`. All gauge momenta are outgoing.

Their emitted-B vertices are

\[
V_\chi(\varepsilon)=-ig'q_\chi\not\varepsilon,
\qquad
V_\phi(s,s';\varepsilon)=-ig'q_\phi(s+s')\cdot\varepsilon,
\qquad s-s'=k.
\]

The scalar kinetic term contains
`(g' q_phi)^2 B_mu B^mu phi^dagger phi`. Differentiating it in the two
gauge fields gives

\[
V_{\phi BB}=+2i(g'q_\phi)^2
(\varepsilon_1\cdot\varepsilon_2).
\]

The zero-B source chain supplies the phase
`(i seed)(i fermion)(-i Yukawa)(i scalar)=-1`. Each ordinary
gauge insertion adds `i(-i)=+1`, and therefore leaves that common sign
unchanged. For the seagull, the extra scalar propagator times the
two-gauge vertex gives `i(2i)=-2`; its full prefactor is hence
`+2F(g'y_L)^2(e_1.e_2)`. With the common `-F(g'y_L)^2` extracted,
the physical diagram weights are

\[
R_2=R_{FF}^{12}+R_{FF}^{21}+R_{SS}^{12}+R_{SS}^{21}
-R_{FS}^{12}-R_{FS}^{21}
-2(\varepsilon_1\cdot\varepsilon_2)R_{SG}.
\]

The minus for each mixed diagram is the product
`q_chi q_phi/y_L^2=-1`. The identical gauge legs give the two ordinary
orderings and the factor two in the seagull rule; there is no further
global factor of `1/2`. These signs follow from vertices and propagators
before any Ward test.

Set

\[
t=p-K,
\qquad D_a=(q-a)^2-m_H^2+i0,
\qquad A_f(q)=2t^2+(q-h-f)^2.
\]

Here `f` is the total gauge momentum emitted from the fermion line before
it reaches the seed. The seed polynomial follows from its direct and
crossed scalar pairings with the already audited relative weights `2:1`.
For each `(i,j)=(1,2),(2,1)`, define `E_i=slash epsilon_i` and
`Q=slash q`. Order `ij` on a single line means `i` is closer to Yukawa,
and `j` is closer to the seed. The two mixed diagrams have photon `i` on
the fermion and photon `j` on the scalar. The seven unweighted integrals
are

\[
I_{FF}^{ij}=\int_q
\frac{A_K(q)(\not q-\not K)E_j(\not q-\not k_i)E_i\not q}
{q^2(q-k_i)^2(q-K)^2D_p},
\]

\[
I_{SS}^{ij}=\int_q
\frac{A_0(q)\not q\,
 [(2p-2q-k_i)\cdot\varepsilon_i]
 [(2p-2q-2k_i-k_j)\cdot\varepsilon_j]}
{q^2D_pD_{p-k_i}D_t},
\]

\[
I_{FS}^{ij}=\int_q
\frac{A_{k_i}(q)(\not q-\not k_i)E_i\not q\,
 [(2p-2q-k_j)\cdot\varepsilon_j]}
{q^2(q-k_i)^2D_pD_{p-k_j}},
\]

\[
I_{SG}=\int_q\frac{A_0(q)\not q}{q^2D_pD_t}.
\]

Every open chain is written from the outgoing lepton at the seed to the
incoming sterile fermion at Yukawa. The suppressed rightmost `P_L` is
restored in the amplitude. These orders agree with the actual
charge-conjugate propagator identity in the prior action audit.

## UV reduction without an assumed action match

For the boxes the numerator has degree at most five in `q`. Only its
degree-four part and the degree-five part with one denominator shift
contribute to the logarithmic UV pole. For the triangle seagull, only
degrees two and three contribute. In particular,

\[
A_f(q)=q^2-2q\cdot(h+f)+\underbrace{(h+f)^2+2t^2}_{\text{UV-finite here}}
\]

can be replaced by its first two terms for this pole. Scalar masses
first enter the denominator expansion two orders below the leading
term, so cannot enter this degree-one residue. This argument establishes
mass independence for arbitrary signed `m_H^2`, without fixing it to a
special value.

The terms proportional to the direct seed weight `2t^2` are individually
UV finite in these diagrams. Therefore this new Green function does
**not** independently retest the seed's `2:1` pairing weight; that weight
is already tested in the earlier seed and zero-/one-B work.

For `N` denominators shifted by `a_l`, write `S=sum_l a_l` and strip
`i/(16 pi^2 epsilon)`. The needed tensor moments, with unnormalized
symmetrization, are

\[
\begin{array}{ll}
N=3:&\operatorname{Res}(q^\mu q^\nu)=g^{\mu\nu}/4,\\[2pt]
&\operatorname{Res}(q^\mu q^\nu q^\rho)
=(g^{\mu\nu}S^\rho+g^{\mu\rho}S^\nu+g^{\nu\rho}S^\mu)/12,\\[2pt]
N=4:&\operatorname{Res}(q^\mu q^\nu q^\rho q^\sigma)
=\operatorname{sym}(gg)/24,\\[2pt]
&\operatorname{Res}(q^{\mu_1}\cdots q^{\mu_5})
=\operatorname{sym}(ggS)/96.
\end{array}
\]

These follow by shifting the Feynman-parametrized denominator and using
the simplex volume `1/(N-1)!` and first moment `1/N!`. All lower-rank
moments needed here vanish in the UV functional. Four-dimensional
Clifford reduction suffices for these simple residues; any `O(epsilon)`
trace corrections change finite terms only.

As a separate reduction, cancel the `q^2` part of each seed against the
unshifted fermion denominator. This lowers every box to a triangle; the
remaining `-2q.(h+f)` term only requires the rank-four box tensor.
The resulting explicit unweighted residues are as follows. For FF set
`w_F=p+k_i-h`; then

\[
R_{FF}^{ij}=
\frac12(\not K E_jE_i+E_i\not k_iE_j)
-\frac16(\not w_F E_jE_i+E_i\not w_FE_j+E_jE_i\not w_F).
\]

For SS set `u_i=2p-k_i`, `v_j=2p-2k_i-k_j`, and
`w_S=3p-2k_i-k_j-h`; then

\[
R_{SS}^{ij}=\frac13\left[
(\varepsilon_i\cdot\varepsilon_j)\not w_S
+(w_S\cdot\varepsilon_i)E_j+(w_S\cdot\varepsilon_j)E_i\right]
-\frac12\left[(u_i\cdot\varepsilon_i)E_j
 +(v_j\cdot\varepsilon_j)E_i\right].
\]

For FS set `u_j=2p-k_j` and `w_M=2p-k_j-h`; then

\[
R_{FS}^{ij}=\frac12\not k_iE_iE_j
-\frac16(E_jE_i\not w_M+\not w_ME_iE_j)
 +\left[\frac{w_M\cdot\varepsilon_j}{3}
       -\frac{u_j\cdot\varepsilon_j}{2}\right]E_i.
\]

Finally,

\[
R_{SG}=\frac12(2\not p-\not K-\not h)
=\frac12(\not p+\not r).
\]

In the FF simplification one uses
`gamma_mu E gamma^mu=-2E` and
`gamma_mu E_j slash v E_i gamma^mu=-2 E_i slash v E_j` at the
residue level. The scalar factors and ordered gamma chains have been
retained through integration. Substituting the seven individual residues
into their physical weighted sum gives the boxed `R_2` above. All `h`
dependence cancels, with `p=r+h+k_1+k_2` understood.

## Exact standalone reproduction

Run

```
python3 -I -S -B -O verify_two_b.py
python3 -I -S -B -O verify_two_b.py --json
```

The standard-library checker imports no repository code. It expands each
complete routed numerator in the loop components and applies the
Feynman-simplex tensor functional directly in an exact rational Clifford
algebra. Separately, it evaluates the displayed denominator-cancellation
formulas and compares each diagram before comparing the weighted sum to
the chosen action kernel. The checks are explicit failures, and remain
active under `python -O`.

UV power counting establishes that the result is linear in the four
independent external vectors `p,h,k_1,k_2`, and bilinear in the two
polarizations. The checker exhausts all `4*4*4*4=256` coefficients of
this polynomial: one of four momentum slots, one of four components, and
one of four components of each polarization. It also tests generic
rational off-shell vectors, Bose symmetry, and deliberate failures from
omitting or reversing the seagull sign, assigning equal charges to the
two internal lines, omitting a fermion ordering, and omitting the
field-strength connection term.

The basis substitutions evaluate the local UV polynomial already defined
at generic off-shell momenta. They are **not** prescriptions to set
exceptional original massless loop integrals to zero as scaleless, and
make no statement about their IR poles. The generic checks retain
arbitrary components, with no on-shell or transverse-polarization
restriction. The default run passes **2232 exact checks**, including 1792
individual-diagram coefficient comparisons, 256 weighted-sum coefficient
comparisons, 12 generic rational cases, and 60 mutation rejections.
`TWO_B_CHECKS.json` records the complete count breakdown and two explicit
sets of seven diagram residues.

## Meaning of the match

The prior zero- and one-B counterterm action passes an additional,
independently integrated two-B 1PI test at this coupling order. The match
includes the scalar seagull and the covariant-derivative connection in
`Q_sigma`; omitting either fails. It fixes no additional physical beta
function, performs no light-field or gauge-field EOM reduction, and does
not establish a complete operator basis, finite matching, or other
external-field sectors.
