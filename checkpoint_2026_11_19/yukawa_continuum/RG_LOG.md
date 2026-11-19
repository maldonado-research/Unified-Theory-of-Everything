# Explicit pure-Yukawa scale logarithm in the physical four-lepton amplitude

For the four-singlet benchmark, the coefficient of $L=\ln(\mu/\mu_*)$ in the fractional polarized rate at $s=6m^2$ is $-0.19674991|y_e|^2-0.11752031|g|^2$. In the endpoint ratio $\sigma(10m^2)/\sigma(6m^2)$ the fractional coefficient is $-0.000919096|y_e|^2+0.0000356549|g|^2$. The implicit variation of the leading box with its running parameters cancels these explicit terms. The finite two-loop pure-Yukawa contribution is an independent function of energy and angle. [Numerical results](RG_LOG.json)

## What is held fixed

Use the unbroken gaugeless model and $\mu_L^-\mu_L^-\to e_L^-e_L^-$ in the phase convention $B_{34}^*B_{12}=s$. The scalar pole mass $m$, physical external momenta and helicities are fixed. At a common reference scale $\mu_*$, take positive sterile masses $(10m,20m,m,4m)$ and disjoint electron and muon rows $(y_e,i\sqrt2 y_e,0,0)$ and $(0,0,g,-ig/2)$. Renormalize the sterile masses and Yukawas by full-theory $\overline{\rm MS}$ with all four singlets retained. The result isolates one additional Yukawa pair multiplying the leading amplitude. [Action and phase](../uv_completion/MATCHING.md) · [Leading full box](../uv_completion/DISPERSION.md)

Let $D$ act on the parameters with the one-loop beta functions stripped of $1/(16\pi^2)$. In fixed field coordinates, with $H=Y^\dagger Y$, $P=YY^\dagger$ and $T=\operatorname{tr}H$, they give
$$
 D Y=\tfrac32YH+TY,\qquad D M_N=H^T M_N+M_N H .
$$
Complex conjugates run with the conjugate equations. If $\mathcal M_0$ is the leading one-loop physical amplitude, the RG equation for the two-loop part $\mathcal M_1^{Y}$ in this coupling sector is
$$
 \left.\partial_{\ln\mu}\mathcal M_1^{Y}\right|_{Y,M_N,m,s,c}
       =-\frac{1}{16\pi^2}D\mathcal M_0 .
$$
The reference amplitude is finite with no explicit subtraction-scale logarithm. Integration of the displayed equation yields one explicit power of $L$ at two loops. A convention using $\ln(\mu^2/\mu_*^2)$ gives one half of each listed coefficient. Dimensionless logarithms of kinematics and mass ratios can also occur in the integration constant. [One-loop full-theory equations and coordinate conventions](../radiative_stability/RUNNING.md)

The physical amplitude includes the residues of the external leptons. For $L_{\rm bare}=Z_L^{1/2}L_R$, the one-loop field anomalous dimension is $\gamma_L=P/(32\pi^2)$. The amputated renormalized vertex has a positive field anomalous-dimension contribution from each leg to its total RG derivative. The square root of the propagator pole residue of $L_R$ has the opposite derivative for that leg. Here $P$ is diagonal, so their sum for two electron and two muon legs cancels directly in LSZ reduction. Internal-field renormalizations are already accounted for by the full-theory RG identity. Holding the scalar pole mass fixed also incorporates, in this identity, the conversion of its mass parameter in the leading box. Applying a scalar running-mass derivative as well would count that contribution twice. [Lepton and scalar two-point poles](../../checkpoint_2026_10_16/yukawa_counterterms/YUKAWA_TWO_POINT.md)

## Applying the flow to the entire leading box

Use $E=|y_e|^2$, $G=|g|^2$, and $\gamma=y_e^2(g^*)^2$ at $\mu_*$. At general running parameters write $b_i=M_{\mu i}^2$, $d_j=M_{e j}^2$, $W_i=M_{\mu i}(Y_{\mu i}^*)^2$, and $V_j=M_{e j}Y_{e j}^2$. With $t,u=-s(1\mp c)/2$ the full amplitude is
$$
 \mathcal M_0 = \frac{s}{16\pi^2}\sum_{i,j}W_i V_j\,\mathcal B_{ij},\quad
 \mathcal B_{ij}=D_0(s,t;b_i,d_j)+D_0(s,u;b_i,d_j).
$$
$D_0$ denotes the upper-rim scalar box with two scalar masses $m$. In the propagator convention of the preceding calculation its Feynman-parameter integrand is $[m^2(x_1+x_2)+b_i x_3+d_j x_4-s x_1x_2-t x_3x_4-i0]^{-2}$. The normalization divided out in the tables is the **reference** $\gamma$. Its running inside the full amplitude contributes to the result. [Box normalization and dispersion integral](../uv_completion/DISPERSION.md)

In instantaneous positive Takagi coordinates the Yukawa derivative is $DY+YA$. The anti-Hermitian $A$ solves $(D M_N)_{ij}+M_i A_{ij}-M_jA_{ij}^*=0$ for $i\ne j$. At the benchmark the two entries are $A_{E_1E_2}=i\sqrt2 E/3$ and $A_{N_1N_2}=-3iG/10$, with vanishing diagonal and no mixing between pairs. $D M_i=2M_i H_{ii}$ gives the positive diagonal mass derivatives. The following multipliers include both effects. A residue derivative in the table denotes $(D V)/V$ or $(D W)/W$, also for complex and negative residues. Here $T=3E+5G/4$. [Derivation and checks of the Takagi rotation](../radiative_stability/RUNNING.md)

| State | Squared mass divided by $m^2$ | $D\ln M_i^2$ | Residue divided by $m y_e^2$ or $m(g^*)^2$ | Residue derivative divided by residue |
|---|---:|---:|---:|---:|
| $E_1$ | 100 | $4E$ | 10 | $2T+29E/3$ |
| $E_2$ | 400 | $8E$ | $-40$ | $2T+41E/3$ |
| $N_1$ | 1 | $4G$ | 1 | $2T+109G/20$ |
| $N_2$ | 16 | $G$ | $-1$ | $2T+109G/20$ |

In particular $(W_1+W_2)'=0$ and $(V_1/d_1+V_2/d_2)'=0$ at the tuned point. The denominators at nonzero transfer still change. Let $k_{e j}=(4,8)_j$, $k_{\mu i}=(4,1)_i$, $a_{e j}=(65/3,77/3)_j$ and $a_\mu=209/20$. For each mass pair, differentiation gives
$$
 \frac{D(W_iV_j\mathcal B_{ij})}{W_iV_j}
   = E\,[a_{e j}\mathcal B_{ij}+k_{e j}d_j\partial_{d_j}\mathcal B_{ij}]
    +G\,[a_\mu\mathcal B_{ij}+k_{\mu i}b_i\partial_{b_i}\mathcal B_{ij}] .
$$
This expression determines the explicit coefficient for every scattering angle and energy in the domain of the box, with the same analytic continuation as the leading amplitude. At cut or threshold boundaries derivatives require the associated limiting prescription.

## Physical rates and the endpoint ratio

Define $u=\mathcal M_0/\gamma$ and $C_e,C_\mu$ by $\partial_L\mathcal M_1^Y/\gamma=EC_e+GC_\mu$ at fixed renormalized inputs. For $r=m^2\sigma/|\gamma|^2$, the coefficients $q_f$ in the explicit rate term $L(Eq_e+Gq_\mu)$ are
$$
 r_0=\frac{m^2}{32\pi s}\left\langle |u|^2\right\rangle_c,\qquad
 q_f=\frac{m^2}{16\pi s}\operatorname{Re}\left\langle u^* C_f\right\rangle_c,
 \qquad \langle f\rangle_c=\tfrac12\int_{-1}^1f(c)\,dc .
$$
The rate interference retains the angular dependence and the absorptive and dispersive parts of both amplitudes.

| $s/m^2$ | $10^5 C_e(s,0)$ | $10^5 C_\mu(s,0)$ | $q_e/r_0$ | $q_\mu/r_0$ |
|---:|---:|---:|---:|---:|
| 6 | $1.398233+0.624022\,i$ | $0.844696+0.352198\,i$ | $-0.19674991$ | $-0.11752031$ |
| 8 | $1.718830+1.204873\,i$ | $1.046747+0.687846\,i$ | $-0.19718252$ | $-0.11753398$ |
| 10 | $1.944583+1.816075\,i$ | $1.191318+1.043271\,i$ | $-0.19766901$ | $-0.11748466$ |

The leading endpoint ratio is $R_0=1.7893622161$. Subtracting the fractional slopes at the two settings gives
$$
 \frac{1}{R_0}\left.\partial_L R_1^{Y}\right|_{Y,M_N}
       =-0.0009190956953 E+0.00003565487369 G .
$$
The small coefficient characterizes the cancellation of the explicit common-scale logarithm in this ratio. It places no bound on a finite, energy-dependent pure-Yukawa correction. Converting the sterile masses to pole inputs or using a running scalar mass moves explicit logarithms between the leading parameter conversion and the higher-loop term. A finite redefinition with no explicit scale dependence preserves this first scale derivative at the stated order after parameters are translated. [Absolute coefficients and leading rates](RG_LOG.json)

## Numerical verification

Run `python -I -B -O check_rg.py --output RG_LOG.json`. The first method differentiates the one-parameter upper-rim dispersion formula at fixed Feynman parameter and integrates the mass and external-angle dependence using 42-digit arithmetic. The independent check advances the full complex Yukawa and nondiagonal symmetric mass matrices using their fixed-coordinate beta functions, rediagonalizes the varied matrices, and evaluates the direct Feynman simplex including its imaginary pole contribution at every step. A fourth-order finite difference of the central amplitude and the integrated squared amplitude agrees within $3\times10^{-7}$ relative for each of two linearly independent nonzero coupling choices at each energy. The dispersion quadrature also passes a $3\times10^{-12}$ refinement tolerance. [Checker](check_rg.py) · [Independent box integrator](../uv_completion/check_box_direct.py) · [Numerical receipt](RG_LOG.json)
