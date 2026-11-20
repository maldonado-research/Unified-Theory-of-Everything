# Independent mixed lepton–Higgs pole audit

26 October 2026 (Pacific). Input commit: `8999d3de1a0e3dc2576fe1586295426da687fc67`.

This audit independently extracts the UV residue of the one-loop mixed `bar L L phi phi^dagger` triangle with one original derivative seed and two neutrino-Yukawa vertices. It uses an auxiliary-mass large-momentum expansion, explicit Weyl spin contractions and literal Wick phases. No producer code or saved coefficient is imported. The result is a primitive zero-gauge 1PI pole in the declared coordinates; it is not a complete physical Wilson coefficient, finite amplitude or full EFT beta function.

Use `eta=(+---)`, `d=4-2epsilon`, `kappa=1/(16 pi^2 epsilon)`, `phi=tilde H`, `C=C^T`, `M=M^T`, and

```text
L_seed=(1/2) C_ab(bar L_a phi)(-Box)(phi^T L_b^c)+h.c.,
L_Y=−bar L Y phi N_R+h.c.,
L_mass=−(bar N_L M N_R+h.c.)/2,
R=Y* M Y^dagger, K=C R.
```

The derivative acts on the entire neutral composite. The anomalous sterile mass numerator is `M`, in the same convention as the preceding four-lepton calculation. No inverse mass or sterile-field elimination is used.

## Route and independent vector-pole extraction

Treat the four external fields as independent Fourier sources with all momenta incoming and mode `exp(−ip.x)`:

```text
bar L(P),  L(R_ext),  phi(S),  phi^dagger(T),
P+R_ext+S+T=0.
```

In the opposite-composite triangle, let `q` be the momentum entering the seed at its **internal barred-lepton** leg. The remaining propagator momenta can be chosen so that the denominators and seed factor are

```text
D_L=q^2,
D_phi=(q+P+S)^2−m_H^2,
D_N=(q−T)^2−M^2,
P0_seed=(q+S)^2.
```

The leading anomalous sterile propagator supplies the mass matrix. Stripping this matrix and the spin/phase factors, the required vector integral is

```text
J_mu=integral d^d q/(2pi)^d q_mu(q+S)^2/(D_L D_phi D_N).
```

Introduce a positive auxiliary squared mass `a` and put `Q=q^2−a`. At most one external-momentum power can contribute to the UV pole, since the superficial degree is one. The odd leading term integrates to zero. Expanding around `Q` yields the two surviving contributions

```text
2 q_mu(q.S)/Q^3
  −2 q_mu q^2[q.(P+S−T)]/Q^4.
```

Use isotropy `q_mu q_nu -> eta_mu nu q^2/d`, evaluating `1/d` at four only after selecting the simple pole. Gamma-function radial residues are

```text
integral q^2/Q^3 |_UV = i kappa,
integral q^4/Q^4 |_UV = i kappa.
```

They give

```text
J_mu|_UV = i kappa[S_mu/2−(P+S−T)_mu/2]
          = i kappa(T−P)_mu/2.                           (1)
```

Auxiliary and physical squared masses enter only odd terms or UV-convergent higher orders. The physical Majorana matrix still survives in the anomalous numerator. Generic off-shell momenta or the positive auxiliary mass separate this UV residue from possible infrared terms; no scaleless integral is used to infer a UV zero.

If both internal seed legs occupy the *same* neutral composite, their loop momenta cancel before squaring `P0`. The seed factor has no loop momentum and the remaining single-slash triangle has superficial degree `−1`. It is UV finite in this sector.

## Explicit Weyl/Wick factor and routing sign

Use the chiral convention of the earlier action: `L=(l,0)`, `epsilon_spin=[[0,1],[-1,0]]`, and

```text
S_seed=−(C/2) phi phi bar l epsilon_spin (-Box) bar l^T,
B_Y=−Y* phi^dagger l^T epsilon_spin n,
<n n^T>=−i M epsilon_spin/(q^2−M^2), n=N_L.
```

For the ordered fermion sequence

```text
bar l_0, bar l_1, l_1, n_1, l_2, n_2,
```

contract one barred lepton with one Yukawa lepton, and contract the two sterile fields. Bring the two contracted pairs before the remaining ordered `(bar l_external,l_external)` pair. For `(bar position, Yukawa lepton position)=(0,2),(0,4),(1,2),(1,4)`, the resulting Grassmann parities are `+1,−1,−1,+1`. Explicit sums over all epsilon components and each of the four Pauli matrices give the **same** `bar sigma.q` tensor for all four assignments after their respective parity signs. The two opposite seed-composite placements have equal `(q+S)^2` factors, and the two Yukawa assignments have the same `K=C Y* M Y^dagger` after the symmetric mass/flavor contractions.

The seed half, two negative Yukawa signs and second-order expansion factor give `iS(iB_Y)^2/2! = +i/4` before propagators. The scalar and anomalous sterile propagator phases are `+i` and `−i`. For this audit's **incoming seed momentum** `q`, the ordered active contraction is

```text
<bar l(x) l(y)> = −<l(y) bar l(x)>^T
                -> +i sigma.q/q^2
```

when the seed field at `x` is Fourier-expanded with `exp(−iq.x)`. Thus the propagator product is `+i`, giving `−1/4` per assignment and `−1` for the four assignments. Multiplying (1) yields the displayed-seed vertex

```text
i Gamma_S,UV = +i kappa/2 (C R)_alpha gamma
                [bar l_alpha bar sigma.(P−T) l_gamma]
                delta_(bar L,L) delta_(phi,phi^dagger).  (2)
```

The action corresponding to (2), before adding its separately conjugated orientation, is

```text
L_S,UV=(kappa/2) K_alpha gamma {
 rho [i(partial_mu bar L_alpha)] gamma^mu L_gamma
 −phi^a[i(partial_mu phi_dagger^a)]
                    bar L_alpha gamma^mu L_gamma },
rho=phi^dagger phi.                                     (3)
```

The complete real action includes `(3)+h.c.` and the cancelling counterterm has the opposite sign. Derivatives in (3) reproduce `P−T` for the declared Fourier waves. No free weak-doublet trace or extra Majorana half is added. The result is valid for arbitrary symmetric complex `M`; `R` is symmetric and transforms as `R'=V^T R V`, while `K'=V^dagger K V`.

An earlier provisional version of this audit used the ordered propagator `−i sigma.q` while simultaneously defining `q` as entering the seed. That mixes opposite routing conventions and reverses the action sign. The producer instead labels momentum `k=−q` flowing from the seed to the Yukawa vertex, for which the ordered propagator is `−i sigma.k`. Its residue is proportional to `P−T`, exactly as (2). The sign was corrected before the executable result and this note were finalized.

## Reproduction and limits

Run `python -I -B -O independent_mixed_audit.py` in the pinned SymPy environment. The standalone checker uses explicit failures and no project imports; its [receipt](INDEPENDENT_MIXED_RESULTS.json) records exact Gamma residues, all four Wick assignments and sixteen Pauli-component contractions, the Fourier action, complex flavor covariance and rejected momentum/sign mutations. These checks are internal AI-assisted research, not independent human peer review. They do not evaluate finite terms, external gauge vertices, other operator insertions or the complete source/EOM and evanescent map.
