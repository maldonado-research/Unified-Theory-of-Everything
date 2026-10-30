# Mixed lepton–Higgs pole from the derivative seed

26 October 2026 (Pacific). Input commit: `8999d3de1a0e3dc2576fe1586295426da687fc67`.

This note computes the previously missing mixed two-lepton/two-Higgs **proper** one-loop pole with one original derivative seed and two conjugate neutrino Yukawa vertices. It uses the reference action, retains the sterile Majorana field, and does not compute a finite amplitude, a complete physical operator coefficient, or full EFT running.

## Literal conventions and result

Set `eta=(+---)`, `d=4-2epsilon`, `kappa=1/(16 pi^2 epsilon)`, `phi=tilde H`, and

```text
L_seed=(1/2) C_alpha beta (bar L_alpha phi)
                       (-Box)(phi^T L_beta^c)+h.c.,
L_Y=-bar L_alpha Y_alpha i phi N_Ri+h.c.,
L_M=-(bar N_L M N_R+h.c.)/2,
C=C^T, M=M^T, R=Y* M Y^dagger, K=C R.
```

The derivative acts on the **whole** neutral composite. Weak labels below are explicit: repeated `u,v` are summed. In the declared four-dimensional physical-spinor projection, the displayed seed orientation gives the local divergent action

```text
L_loop,UV = (i kappa/2) K_alpha gamma phi^v
   [ (partial_mu bar L_alpha^u) phi_dagger^v
       - bar L_alpha^u (partial_mu phi_dagger^v) ]
                           gamma^mu L_gamma^u + h.c.       (1)
```

The cancelling counterterm action is the negative of (1). Its flavor matrix transforms as `K' = V^dagger K V` under the established passive active/sterile rotations. The weak contraction pairs the two lepton fields and, separately, the two Higgs fields. No free doublet trace or extra Majorana half is appended. With one active flavor and nonzero `C y*^2 M`, (1) need not vanish.

For independent Fourier sources with all momenta incoming and all waves `exp(-ip.x)`, write `p` for the external barred lepton, `k` for `phi`, `h` for `phi_dagger`, and `r` for the external lepton; `p+k+h+r=0`. The displayed seed contributes the local vertex

```text
i Gamma_UV,S = + i kappa K (p-h)_mu gamma^mu P_L/2,        (2)
```

with the weak pairings and external spinors implicit. Equation (1) reproduces (2) because `i partial` acts as incoming momentum. Add its Hermitian conjugate with its own conjugate derivative placements; it is not generally twice the real part of one placement.

## UV integration and derivative placement

Route internal active momentum `q` from the Yukawa vertex to the seed. The denominators and differentiated composite are

```text
D0=q^2,
D1=(q-p-k)^2-m_H^2,
D2=(q+h)^2-M^2,
P0 at the seed=(q-k)^2.
```

The ordered active numerator is `q_mu`. The anomalous same-chirality sterile propagator supplies the retained mass numerator `M`, yielding `Y* M Y^dagger` after its flavor endpoints are contracted. Its next mass insertion decreases UV degree by two and cannot affect this momentum-linear pole. For a complex symmetric `M` no inverse is needed; the source and `K` transform covariantly in the established coordinates.

Combine the denominators with weights `1-y-z,y,z`, where `y,z>=0`, `y+z<=1`, and shift

```text
q=l+v, v=y(p+k)-zh.
```

After odd terms are discarded, the simple pole of the numerator `q_mu(q-k)^2` is

```text
(1/2)(3v-k)_mu int l^2/(l^2-Delta+i0)^3 |_UV.
```

The factor `1/2` follows from `2/d` evaluated at `d=4` for a simple pole. The radial identity

```text
l^2/(l^2-Delta)^3 = 1/(l^2-Delta)^2 + Delta/(l^2-Delta)^3
```

and Gamma residues give `i kappa` for the first term and no UV pole for the second. The exact normalized simplex moments `2 int 1=1`, `2 int y=2 int z=1/3` therefore give

```text
int q_mu(q-k)^2/(D0 D1 D2) |_UV = i kappa (p-h)_mu/2.    (3)
```

The result is independent of internal squared masses at this degree, but the physical Majorana mass in the anomalous numerator remains. Use generic off-shell external momenta or an independently controlled positive auxiliary mass to separate UV from IR; a scaleless integral is not used to infer a zero.

As a separate algebraic check, the leading `q_mu q^2` numerator cancels `D0`. The resulting bubble has vector residue `(p+k-h)_mu/2`; the `-2q_mu q.k` triangle has residue `-k_mu/2`, while the `k^2 q_mu` term is UV finite. Their sum is (3). This numerator reduction is a second check of the integration, not a second independent Wick calculation.

The seed's internal charged scalar and lepton can instead lie in the **same** neutral composite. Their loop momenta then cancel within `P0`, leaving only external momentum squared. The remaining numerator has one loop momentum over three quadratic denominators, with superficial degree `-1`; these same-composite proper graphs are UV finite. The preceding [topology audit](../../checkpoint_2026_10_23/seed_yukawa/SEED_TOPOLOGY.md) separately exhausts the literal assignments and bridges.

## Wick sign, half and momentum routing

In the chiral reference basis set `E=[[0,1],[-1,0]]`, `L=(l,0)` and `N_L=(n,0)`. The relevant terms and ordered propagators are

```text
S=-(1/2) C phi phi bar l E P0 bar l^T,
B=-Y* phi_dagger l^T E n,
<bar l(x) l(y)> = -i sigma.q exp[-iq.(y-x)]/q^2,
<n n^T>(q)=-i M E/(q^2-M^2),
<phi phi_dagger>(q)=+i/(q^2-m_H^2).
```

The negative sign in the ordered active contraction matters: `q` here flows **into the Yukawa endpoint and out of the seed**. If instead `q_in=-q` is defined incoming at the seed, the ordered active numerator is `+i sigma.q_in`. Reversing the routing changes both this numerator and the vector integral, leaving (1) unchanged. An initially opposite sign obtained using the incoming-seed momentum together with the outgoing-seed numerator was corrected by this Fourier check.

For one labelled opposite-composite contraction take fermion order

```text
[bar l_alpha, bar l_beta, l_delta, n_i, l_gamma, n_j],
```

contract `(bar l_beta,l_delta)` and `(n_i,n_j)`, and leave `(bar l_alpha,l_gamma)` in that order. Moving the contracted fields to the front has odd permutation parity, giving Wick sign `-1`. Direct Pauli-matrix multiplication gives

```text
sum E_rs sigma^mu_{t s} E_tu E_vw E_uw = -bar sigma^mu_{r v}.
```

The seed vertex and the two identical Yukawa vertices including their `1/2!` have phase `+i/4`. The active, anomalous sterile and scalar propagator phases multiply to `(-i)(-i)(+i)=-i`. Including the Wick and spin-tensor signs gives `+1/4` for this labelled assignment. The two opposite-composite seed choices and the two Yukawa-label assignments give four equal terms, hence unit total weight before (3). This is a contraction count; the seed half and expansion factorial must not be removed separately. The separate auditor obtains the same action and Fourier polynomial using the opposite internal momentum routing.

## Hermitian currents and source transport

Define

```text
H=(K+K^dagger)/2, A=(K-K^dagger)/(2i),
rho=phi^dagger phi,
J_phi^mu=i[phi^dagger partial^mu phi-(partial^mu phi_dagger)phi].
```

Both `H` and `A` are Hermitian. Integrating a displayed boundary term by parts rewrites the full (1) as

```text
L_loop,UV = (kappa/2)[ -rho bar L H i<->slash partial L
                         + J_phi,mu bar L H gamma^mu L ]
             + kappa (partial_mu rho) bar L A gamma^mu L,  (4)
```

where `bar L i<->partial L=i bar L partial L-i (partial bar L)L`. The coefficient `Tr(A)=Im Tr(CR)` is the imaginary scalar derivative coefficient found in the preceding four-scalar pole. This equality is an algebraic relation between two primitive contributions; it does not show that physical effects cancel.

For example, at first order in the loop pole the local changes

```text
delta L=i kappa rho A L,
delta phi=i kappa rho Tr(A) phi
```

cancel the last term of (4) and the preceding scalar term `kappa Tr(A)(partial rho).J_phi` through the variations of their kinetic actions. The declared neutrino Yukawa interaction then acquires

```text
delta L_Y=+i kappa rho bar L [A-Tr(A) I] Y phi N_R+h.c.  (5)
```

For one active flavor `A=Tr(A)`, so `phi^dagger L` and the displayed neutrino Yukawa are unchanged by this joint phase. Additional spectator Yukawas and external probes still transform and must be kept in any larger reference action. This is a first-order sourced action-coordinate identity, not the complete quantum effective-action map or a physical cancellation.

## Reproduction and limitations

Run `python -I -B -O verify_mixed_seed_pole.py` in the pinned symbolic environment. The standalone explicit-guard checker performs exact simplex moments, Gamma residue and bubble-reduction controls, four Pauli spin tensors, the representative Wick phase, complex flavor covariance, and canonical-epsilon homogeneity. [Its receipt](MIXED_SEED_RESULTS.json) reports 41 positive checks and 9 rejected mutations. The independent [mixed audit](INDEPENDENT_MIXED_AUDIT.md) checks the opposite routing and the overall action separately.

This completes the three **declared proper matter field-content sectors** at one seed and two pure neutrino-Yukawa vertices when combined with the previously calculated four-lepton and four-scalar poles. It does not collect the reducible earlier three-point source counterterms, the full physical EOM/IBP and evanescent basis, other interactions or operator insertions, quantum source transport, finite matching or observables. No sterile field has been integrated out. Four-dimensional external-spinor identities do not define a complete `d`-dimensional evanescent prescription.
