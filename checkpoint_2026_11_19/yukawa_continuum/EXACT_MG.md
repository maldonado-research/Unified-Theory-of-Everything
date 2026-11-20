# MadGraph check of the full three-particle cut

At $s=8m^2$, the regular exchange diagrams increase the magnitude of the $J=0$ three-particle cut by **15.3333%**. MadGraph5_aMC@NLO generates both polarized tree amplitudes entering the unitarity sum. At a fixed kinematic point the full-to-pole spin overlap ratio is $1.257438753877467$. The phase-space integration gives $1.153332814$. The intermediate state is $\mu^-N_1\phi^-$, with $N_1$ and $\phi$ of mass $m$. These results test this intermediate state and the relative sign of its tree diagrams. [Numerical receipt](MADGRAPH.json) · [Standalone driver](madgraph_reference.py)

## Generated model and helicities

The driver writes a minimal UFO with massless charged Dirac leptons, a complex charged scalar of mass one, and self-conjugate fermions of masses $1,4,10,20$. Widths vanish. The muon couples to the first two neutral fermions with Yukawas $(1,i/2)$ and the electron to the last two with $(1,i\sqrt2)$. The UFO vertices use $iY$ with the ordered fields $(N,\ell^-,\phi^+)$ and the Lorentz tensor `ProjP(2,1)`, and $iY^*$ with $(\ell^+,N,\phi^-)$ and `ProjM(2,1)`. This field ordering and chiral assignment also occur in the upstream left-slepton–neutralino vertices. ALOHA generates the wavefunction, propagator and amplitude routines. At the reference point the entire helicity table has two nonzero entries, $(-,-;-,-,0)$ and $(-,-;-,+,0)$ in the incoming-lepton, outgoing-muon, $N_1$, scalar order. [MadGraph model conventions](https://github.com/mg5amcnlo/mg5amcnlo/tree/7eba004757af0646e8e1db9d32e3d86476d9c9d5/models/MSSM_SLHA2) · [ALOHA](https://doi.org/10.1016/j.cpc.2012.05.004)

The process generator finds eight diagrams for $\mu^-\mu^-\to\mu^-N_1\phi^-$ and four for $e^-e^-\to\mu^-N_1\phi^-$. In the muon pole graph selection, the final $\mu^-$ and $N_1$ must be attached at the same charged-scalar vertex, giving four graphs including the two possible internal neutral-fermion species. MadGraph fixes the fermion permutations and relative signs before this selection. All four electron graphs have the scalar pole. Define $I_f=\sum_{h_N}T_\mu^{f}T_e^*$ for $f=\mathrm{full},\mathrm{pole}$. At the symmetric point with $p_\phi=(\sqrt{5}/2,1/2,0,0)$ and all outgoing longitudinal three-momenta zero, the overlaps are $-0.009903039068135592$ and $-0.007875563750201233$. The two additional noncoplanar points in the receipt have complex overlap ratios. [Four-vectors and spin sums](MADGRAPH.json)

## Integrated projection

The electron side is a pure scalar-pole amplitude. After removing the external incident spinor phase its beam orientation enters only through $S_e(x)=F_e(h-vx)+F_e(h+vx)$, where $F_e(q)=-30q/[(100+q)(400+q)]$, $h=(s-1-U)/2$, $v=\sqrt{h^2-U}$ and $x$ is the angle cosine of the outgoing pole scalar with the muon beam. Write $U=(p_\mu+p_N)^2$. The calculation evaluates both generated amplitudes with aligned incident axes and weights their overlap by $\int_{-1}^1 S_e(x')\,dx'/[2S_e(x)]$. This analytically projects the electron beam over directions for $J=0$. At a fixed scalar direction the driver retains both decay angles in the generated muon amplitudes. The implemented logarithmic electron angular average and the unitarity normalization are external inputs to the generator check.

With normalized angular averages, the mass measure including the unitarity factor is

$$
\frac12\,\frac{dU}{2\pi}
\frac{2v}{8\pi s}
\frac{U-1}{8\pi U}.
$$

The identical incident leptons give a symmetric physical integral under reversal of the muon beam, so the scalar beam cosine is integrated over $[0,1]$. The mass map $U=1+[(\sqrt{s}-1)^2-1](1-a^2)$ resolves the upper recoil endpoint. Gauss–Legendre rules integrate $a$, the scalar beam cosine and the lepton decay cosine. Equally spaced midpoint angles integrate the decay azimuth.

| Mass, scalar cosine, decay cosine, azimuth orders | Full $J=0$ cut | Scalar-pole $J=0$ cut |
|---|---:|---:|
| 6, 4, 4, 6 | $-1.457503416\times10^{-6}$ | $-1.263573727\times10^{-6}$ |
| 10, 6, 6, 10 | $-1.457577288\times10^{-6}$ | $-1.263795665\times10^{-6}$ |
| 14, 8, 10, 14 | $-1.457577493\times10^{-6}$ | $-1.263796082\times10^{-6}$ |
| 18, 12, 12, 16 | $-1.457577494\times10^{-6}$ | $-1.263796084\times10^{-6}$ |

## Reproduction

Use Python 3 with NumPy. The recorded run used Python 3.12 and the upstream tag `v3.5.7`, commit `7eba004757af0646e8e1db9d32e3d86476d9c9d5`. The driver creates its UFO in a temporary directory and evaluates MadGraph's Python ALOHA/HELAS implementation. The standard repository runner can read the stored receipt without downloading MadGraph. [Pinned release](https://github.com/mg5amcnlo/mg5amcnlo/tree/7eba004757af0646e8e1db9d32e3d86476d9c9d5) · [MadGraph5_aMC@NLO publication](https://doi.org/10.1007/JHEP07(2014)079)

```bash
curl -fL https://codeload.github.com/mg5amcnlo/mg5amcnlo/tar.gz/refs/tags/v3.5.7 -o mg5.tar.gz
echo 'c3fbdf8ec813270cfbfbb4d7aff10c91c30c32ce3ea0ffc59e7c99ce03134489  mg5.tar.gz' | sha256sum -c -
tar xf mg5.tar.gz
python -I -B madgraph_reference.py --mg5 mg5amcnlo-3.5.7 --output MADGRAPH.json
```
