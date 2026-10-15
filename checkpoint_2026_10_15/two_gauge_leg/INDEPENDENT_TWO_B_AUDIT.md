# Independent two-hypercharge-leg UV audit

15 October 2026. Baseline: repository HEAD
`bac2240c12ba59d65b1407605ce2744c30e045f3`, and the literal action and Fourier
conventions in `checkpoint_2026_10_15/operator_map/ACTION_SIGN.md`.
This note and its checker were written outside the repository.

## Result and independence

An independent large-loop-momentum extraction of all seven specified diagrams
gives

\[
R_{BB}(p,h;k_1,k_2;\varepsilon_1,\varepsilon_2)
=(\varepsilon_1\!\cdot\varepsilon_2)\not p
-\frac14\left([\not k_1,\not\varepsilon_1]\not\varepsilon_2
              +[\not k_2,\not\varepsilon_2]\not\varepsilon_1\right).
\]

Consequently, with `F=C_(7,seed)Y_nu*`, `y_L=-1/2`, and the rightmost `P_L`
understood,

\[
i\mathcal M_{BB}^{\rm UV}
=-\frac{iF(g'y_L)^2}{16\pi^2\epsilon}
 \bar u_L R_{BB}P_Lu_N.
\]

The existing positive counterterm action has the opposite vertex pole and
cancels this result. Its predicted two-B kernel is now checked by explicit
loop integration in this selected sector. No additional coefficient is needed
for this particular two-B vertex.

The computation starts from the seven routed numerators, then extracts their
homogeneous UV terms using angular moments. It does **not** define the loop
residue from the target action or from Ward identities. The action comparison
is performed only after each diagram has been independently reduced. The code
imports no repository or producer module; it uses rational Clifford blades
instead of the producer's Feynman-simplex tensor implementation. The previously
audited one-B residue appears only in the longitudinal control.

## Ordered diagrams and phases

Use the (+---) metric and the baseline source orientation. Write

\[
K=k_1+k_2,\qquad t=p-K=r+h,\qquad
S_A(q)=2t^2+(q-h-A)^2,
\]

\[
P_a=(q-a)^2+i0,\qquad
D_a=(q-a)^2-m_H^2+i0.
\]

In an ordered same-line diagram, label 1 is emitted nearest the Yukawa vertex
and label 2 nearest the seed. In a mixed diagram, the first label is on the
fermion. Suppress the integral measure and the common rightmost projector.

| Graph | Numerator | Denominator |
|---|---|---|
| FF12 | `S_K (/q-/K) /eps2 (/q-/k1) /eps1 /q` | `P_0 P_k1 P_K D_p` |
| FF21 | interchange `(k1,eps1)` and `(k2,eps2)` in FF12 | corresponding interchange |
| SS12 | `S_0 /q [eps1.(2p-2q-k1)] [eps2.(2p-2q-2k1-k2)]` | `P_0 D_p D_(p-k1) D_t` |
| SS21 | interchange `(k1,eps1)` and `(k2,eps2)` in SS12 | corresponding interchange |
| FS12 | `S_k1 (/q-/k1) /eps1 /q [eps2.(2p-2q-k2)]` | `P_0 P_k1 D_p D_(p-k2)` |
| FS21 | interchange `(k1,eps1)` and `(k2,eps2)` in FS12 | corresponding interchange |
| SG | `(eps1.eps2) S_0 /q` | `P_0 D_p D_t` |

There is one seagull graph: both identical-field assignments are included in
its vertex factor. The source propagating fields are `chi=L^c` of charge
`-y_L` and `phi=tilde H` of charge `y_L`. The fermion and scalar emission
vertices are respectively

\[
-ig'(-y_L)\gamma^\mu P_R,\qquad
-ig'y_L(s+s')^\mu.
\]

Each ordinary insertion adds one propagator and one gauge vertex, with
relative phase `i(-i)=+1`. Thus FF and SS have positive charge weights and
FS has negative charge weight after stripping `(g'y_L)^2`.

The quadratic scalar kinetic term is `+(g'y_L)^2 phi^dagger phi B_mu B^mu`.
Two B derivatives and multiplication by `i` give the seagull vertex
`+2i(g'y_L)^2 g^{mu nu}`. Its extra propagator therefore gives relative
factor `i(2i)=-2`. This sign precedes any Ward test. The common zero-gauge
source phase is `(i)(i)(-i)(i)=-1`. There is no closed fermion loop.
Accordingly the reduced integral combination is

\[
I_{BB}=I_{FF12}+I_{FF21}+I_{SS12}+I_{SS21}
        -I_{FS12}-I_{FS21}-2I_{SG}.
\]

The derivative seed acts on a neutral complete composite and has no direct
B contact. Radiation from external legs is reducible and is not included in
this 1PI vertex. The symmetric seed's previously fixed 2:1 contraction weights
remain inside `S_A`; there is no further factor of two.

## Homogeneous ultraviolet extraction

For four propagators with shifts `a_j`, the largest numerator degree is five.
Let `N_n` denote the degree-n homogeneous part in q. The UV logarithm is

\[
\operatorname{Res}I_{\rm box}
=\mathcal A[N_4]+\mathcal A[2q\!\cdot(\sum_j a_j)N_5],
\]

where each angular monomial is divided by the matching power of `q^2`.
Equivalently the two terms have denominators `(q^2)^4` and `(q^2)^5`
before stripping their radial logarithm. For the seagull triangle,

\[
\operatorname{Res}I_{\rm SG}
=\mathcal A[N_2]+\mathcal A[2q\!\cdot(p+t)N_3].
\]

The normalized angular operation obeys

\[
\mathcal A[q^\mu q^\nu]=\frac{g^{\mu\nu}}4,
\qquad
\mathcal A[q^\mu q^\nu q^\rho q^\sigma]
=\frac{g^{\mu\nu}g^{\rho\sigma}
       +g^{\mu\rho}g^{\nu\sigma}
       +g^{\mu\sigma}g^{\nu\rho}}{24},
\]

and the rank-six moment has the complete Wick pairing sum divided by 192.
An auxiliary infrared separation is understood: these are UV radial
coefficients, not an assertion that an unseparated scaleless integral has a
nonzero dimensionally regulated value.

Only `q^2-2q.(h+A)` from the seed reaches these degrees. The direct pairing
`2t^2` and the constant `(h+A)^2` are UV finite separately in these diagrams.
Masses first appear at second order in the propagator expansion, beyond the
one order available to produce this logarithm. Hence there is no scalar-mass
dependence; no sterile propagator is present. The residue is jointly linear
in the external momenta and bilinear in the two polarization vectors.

Four-dimensional angular and Clifford factors give the simple pole, since
their O(epsilon) corrections only affect finite terms. This statement does
not determine any finite or evanescent matching coefficient.

## Individual analytic residues

These formulas are independently compared with the raw angular extraction,
before forming the charge-weighted sum. A slash string keeps the displayed
order. For FF and SS let `b,eb` be the first emission and `a,ea` the second,
so `a+b=K`.

\[
\begin{aligned}
R_{FF}(b,a)={}&\frac12\not K\not e_a\not e_b
             +\frac12\not e_b\not b\not e_a\\
 &-\frac16\left(\not d\not e_a\not e_b
                +\not e_b\not d\not e_a
                +\not e_a\not e_b\not d\right),
\qquad d=b+p-h.
\end{aligned}
\]

\[
\begin{aligned}
R_{SS}(b,a)={}&-\frac12(u\!\cdot e_b)\not e_a
              -\frac12(w\!\cdot e_a)\not e_b\\
 &+\frac13\left[(e_a\!\cdot e_b)\not d
                +(d\!\cdot e_a)\not e_b
                +(d\!\cdot e_b)\not e_a\right],\\
u={}&2p-b,\qquad w=2p-2b-a,\qquad d=3p-2b-a-h.
\end{aligned}
\]

For FS, let `a,ea` be the fermion emission and `b,eb` the scalar emission:

\[
\begin{aligned}
R_{FS}(a,b)={}&-\frac12(u\!\cdot e_b)\not e_a
              +\frac12\not a\not e_a\not e_b\\
 &+\frac13\left[(e_a\!\cdot e_b)\not d
                +(e_a\!\cdot d)\not e_b
                -2(e_b\!\cdot d)\not e_a\right],\\
u={}&2p-b,\qquad d=h-2p+b.
\end{aligned}
\]

\[
R_{SG}=\frac12(e_1\!\cdot e_2)(2\not p-\not K-\not h)
       =\frac12(e_1\!\cdot e_2)(\not p+\not r).
\]

Adding the six ordered residues with their physical charge weights, and the
seagull with its relative phase, gives the result at the top. In particular
the separate h dependence cancels. No external EOM, transversality, on-shell
condition, or projection onto one tensor structure is used.

## Reproduction and checks

Run:

```sh
python3 -I -S -B -O verify_two_b_angular.py
```

The completed default run reports **2483 checks**, **24 generic rational
samples**, and **256 complete local-polynomial coefficient points**. The
record is `INDEPENDENT_TWO_B_CHECKS.json`; it includes an exact numerical
witness for all seven individual residues and the final Clifford coefficients.
The 256 points span all 16 independent coordinates of `(p,h,k1,k2)` times
the 16 choices of polarization basis pair. Because the computed residue is
jointly degree one in those momenta and degree one in each polarization,
this is an exhaustive coefficient identity for this local polynomial.
It is not a claim to enumerate a full EFT operator basis.

Checks include:

- Each of the seven angular residues agrees with its individual analytic
  formula, at all 24 rational points and all 256 coefficient points.
- The full unprojected Clifford expression agrees with the action prediction.
- Recomputing the diagrams after swapping the two complete gauge labels gives
  the same result.
- Contracting each gauge leg independently gives
  `R_1(p,h+k_i;k_j,e_j)-R_1(p,h;k_j,e_j)`, for three signed scalar masses.
  Here `R_1=R_H-R_L` is the independently audited one-B residue.
- Double contraction gives `(k1.k2) slash(p)`.
- The literal negative loop pole and positive counterterm vertex cancel.
- Six deliberate mutations fail at all 24 rational points: reversing the
  seagull sign, dropping FF21, dropping SG, giving the two internal fields
  equal signed charges, omitting the `Q_sigma` two-B contact, and reversing
  the counterterm sign.

## Limits

This establishes only the displayed one-loop 1PI vertex at one seed
insertion, one conjugate neutrino Yukawa, and two hypercharge couplings.
The flavor order remains `C_(7,seed)Y_nu*`, and the Hermitian-conjugate
channel is supplied separately. It is not an EOM reduction, a full Green
basis, a complete anomalous-dimension calculation, a physical running
coefficient, or a new observable. It validates the two-B content of the
previously chosen action lift in this restricted sector.
