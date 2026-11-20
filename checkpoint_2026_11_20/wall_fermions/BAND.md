# Resolving the spectrum through shell 2048

The heavier species has the larger expected particle number at the expanding readout once shells above 32 are included. By the turnaround the lighter species has the larger number. Both orderings survive the maximum contribution from the omitted shells on this interpolated neutral-parent geometry. [Computed high band](BAND.json) · [Low band](MODES.json) · [Tail estimate beyond 2048](TAIL.md)

| Slice | $N_{\le2048}$, $ML=1$ | Upper bound on omitted $N$, $ML=1$ | $N_{\le2048}$, $ML=4$ | Upper bound on omitted $N$, $ML=4$ |
|---|---:|---:|---:|---:|
| Bounce | 0.01288064 | 0.00000658 | 0.04523021 | 0.00010525 |
| Expanding at $t=11.03911$ | 0.11593655 | 0.00022632 | 0.13292642 | 0.00362173 |
| Turnaround | 0.09913113 | 0.00001355 | 0.06617083 | 0.00021689 |

The resolved numbers are numerical expectation values for instantaneous positive-energy projectors on the stated proper-time slices. The omitted-number inequalities apply to $n\ge2049$ for the specified spline. Adding the resolved contributions provides all-shell estimates with a separately bounded ultraviolet remainder. The finite-band solver tolerances below quantify numerical agreement for the resolved terms.

The summed turnaround excitation energies $E_{\le2048}L$ are 0.14747794 and 0.45498880. The corresponding remaining-energy bounds are 0.00427324 and 0.06844023. The heavy species therefore carries more excitation energy at the turnaround. In the expanding readout the last two geometrically spaced shell bands contribute nearly equal excitation energies for each mass. This is consistent with the logarithmic ultraviolet behavior of the instantaneous energy on a slice with $\dot a\ne0$. [Tail analysis](TAIL.md)

## Turnaround spectrum by shell range

| Inclusive shells | $N$, $ML=1$ | $EL$, $ML=1$ | $N$, $ML=4$ | $EL$, $ML=4$ |
|---|---:|---:|---:|---:|
| 0–32 | 0.097367822 | 0.13102701 | 0.045241176 | 0.20270956 |
| 33–128 | 0.001548775 | 0.00870599 | 0.01751732 | 0.12767969 |
| 129–512 | 0.000177157 | 0.00427316 | 0.00281424 | 0.06899822 |
| 513–2048 | 0.000037378 | 0.00347178 | 0.00059810 | 0.05560133 |

Each integer shell is evolved and weighted separately by $2(n+1)(n+2)$. The entries are sums over those shells. At $n=2048$ the two turnaround probabilities are $4.81710\times10^{-16}$ and $7.72149\times10^{-15}$. [Group sums and representative probabilities for all three slices](BAND.json)

## Resolving small mixing amplitudes

Rotate the instantaneous positive and negative eigenvectors into the columns of $R=(u_+,u_-)$. With $Y=R^\dagger u$ the proper-time equation becomes

$$
\dot Y=A(t)Y,\qquad
A=-i\omega\sigma_z+i\dot\theta\sigma_y,\qquad
\dot\theta=-\frac{\mu k\dot a}{2(k^2+\mu^2a^2)},\qquad
Y(0)=(0,1)^T.
$$

Here $\mu=ML$ and $k=n+3/2$. The positive-energy amplitude starts exactly at zero and its probability is obtained directly as $|Y_1|^2$. This representation avoids subtraction of nearly equal fixed-basis amplitudes in the primary ultraviolet calculation.

Within each archived Hermite cell of width $h$, sample the generator at its two Gauss nodes, $t_{\rm mid}\mp h\sqrt3/6$. Call the matrices $A_-,A_+$. With $b_\pm=1/4\pm\sqrt3/6$, the fourth-order commutator-free approximation propagates each state with

$$
e^{h(b_-A_-+b_+A_+)} e^{h(b_+A_-+b_-A_+)} .
$$

The exponentials use the exact sine and cosine formulas for traceless anti-Hermitian two-by-two matrices. All changes of Hermite pieces and readouts are explicit mesh boundaries. Each of the 2,016 high shells is evolved on the archived mesh and after bisecting every cell. The receipt reports the bisected-mesh solution. At each readout, the sum of absolute probability changes times the multiplicity is below $3\times10^{-9}$. Its energy-weighted counterpart is below $2\times10^{-7}$.

At shells 32, 33, 128, 512, 1024 and 2048 the mesh is bisected a second time. A second calculation uses the fixed spinor basis, the original generator $-i(\mu\sigma_z+k\sigma_x/a)$, and explicit endpoint projection. Both the extra refinement and the basis change satisfy $|\Delta P| \le 3\times10^{-5}P_{\rm ref}+2\times10^{-23}$. The overlapping shell 32 also agrees with the separate low-band calculation within $8\times10^{-14}$ absolute. These checks cover oscillatory phase accumulation as well as small occupation amplitudes. [Driver and tests](check_band.py) · [Propagation kernel](band_evolution.cpp) · [Low-band implementation](check_modes.py)

Run `python -I -B -O check_band.py --output BAND.json` in this directory with NumPy, SciPy and a C++17 compiler. The driver builds the kernel in a temporary directory. It uses `g++` by default or the compiler path in `CXX`.
