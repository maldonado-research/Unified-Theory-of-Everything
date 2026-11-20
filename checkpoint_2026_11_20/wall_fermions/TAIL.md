# Ultraviolet bounds at the stationary turnaround

For a free Majorana spectator initialized in the negative instantaneous eigenstate on the closed wall, the shells above the numerical band obey the following bounds. The geometry is the piecewise cubic Hermite curve defined by the archived binary64 knots, radii and proper-time slopes in [check_modes.py](check_modes.py). Both helicities are counted. Quantities in the calculation use the curvature length \(L\) as the unit of length and time.

| Mass \(mL\) | Upper bound on \(N_{n\ge33}\) | Upper bound on \(L E_{n\ge33}\) | Upper bound on \(N_{n\ge2049}\) | Upper bound on \(L E_{n\ge2049}\) |
|---:|---:|---:|---:|---:|
| 1 | 0.32181927 | 1.77517701 | 0.00001355 | 0.00427324 |
| 4 | 6.87931081 | 45.69204931 | 0.00021689 | 0.06844023 |

Each entry bounds a sum of nonnegative instantaneous occupations, weighted for energy by the final one-particle frequency. The bounds permit substantial contributions above shell 32. The calculation addresses finite-time basis occupation. A renormalized stress tensor requires its own state and subtraction prescription.

Repeating the stationary-endpoint calculation on the prefix ending at the bounce, and retaining the leading boundary term on the expanding slice, gives the following independent readout bounds.

| Mass \(mL\) | Bounce \(N_{n\ge33}\) upper | Bounce \(N_{n\ge2049}\) upper | Bounce \(LE_{n\ge2049}\) upper | Expanding slice \(N_{n\ge2049}\) upper |
|---:|---:|---:|---:|---:|
| 1 | 0.00183940 | 0.00000658 | 0.01075080 | 0.00022632 |
| 4 | 0.02945485 | 0.00010525 | 0.17201378 | 0.00362173 |

## Transition amplitude

Use conformal time \(d\eta=dt/a\), \(\kappa=n+3/2\), and dimensionless mass \(m\). The rescaled two-level Hamiltonian is \(\kappa\sigma_x+\mu\sigma_z\) with \(\mu=ma\) and \(\omega=(\kappa^2+\mu^2)^{1/2}\). Primes below denote conformal-time derivatives. With dynamical phases removed, choose signs of the mixing amplitude so that

\[
 \beta'=q\alpha e^{i\phi},\qquad
 \alpha'=-q\beta e^{-i\phi},\qquad
 \phi'=2\omega,\qquad q=\frac{\kappa\mu'}{2\omega^2}.
\]

The initial conditions are \(\alpha=1,\ \beta=0\), and unitarity gives \(|\alpha|,|\beta|\le1\). Throughout the trajectory, \(|\beta|\le b=\min(1,\int |\mu'|\,d\eta/(2\kappa))\). Write

\[
 F=\frac{q}{2\omega}=\frac{\kappa\mu'}{4\omega^3},\qquad
 G=\frac{F'}{2\omega}=\frac{\kappa}{8}
 \left(\frac{\mu''}{\omega^4}-\frac{3\mu(\mu')^2}{\omega^6}\right).
\]

Integration of \(\beta'\) by parts gives a boundary term \([F\alpha e^{i\phi}/i]\), an oscillatory integral of \(F'\alpha\), and an integral bounded by \(b\int |Fq|\,d\eta\). The boundary term vanishes because both endpoint slopes vanish. To bound the remaining oscillatory integral on a block \(B\), take the smaller of

\[
 \int_B |F'|\,d\eta
 \quad\text{and}\quad
 |G(B_\mathrm{left}^{+})|+|G(B_\mathrm{right}^{-})|
 +\sum_{j\,\mathrm{interior}}|G_j^+-G_j^-|
 +\int_B\left(|G'|+b|Gq|\right)d\eta .
\]

The first integration by parts was taken globally. Its \(F\) boundary terms cancel across all interior splits by continuity. Only the \(F'\) remainder is subdivided when choosing a bound on each block. The second inequality follows by another integration by parts on each open cubic cell. At each interior knot \(\mu,\mu'\) are continuous, so the jump is exactly \(\kappa|\Delta\mu''|/(8\omega^4)\). The classical derivative inside a cell is

\[
 G'=\frac{\kappa}{8}\left[
 \frac{\mu'''}{\omega^4}
 -\frac{10\mu\mu'\mu''+3(\mu')^3}{\omega^6}
 +\frac{18\mu^2(\mu')^3}{\omega^8}\right].
\]

All the terms in the second bound, including the separate \(b\int |Fq|\), are \(O(\kappa^{-3})\) or smaller. Thus at the turnaround \(p_n=|\beta_n|^2=O(\kappa^{-6})\). With \(g_n=2(n+1)(n+2)=2(\kappa^2-1/4)\), both \(\sum g_n p_n\) and \(\sum g_n\sqrt{m^2+\kappa^2/a_T^2}\,p_n\) converge.

## Deterministic enclosure

For each cell, use its exact rational values of \(a_L,a_R,\dot a_L,\dot a_R,\Delta t\) after the shared binary64 geometry construction. The Bernstein control values of \(a\) are \(a_L,a_L+\Delta t\dot a_L/3,a_R-\Delta t\dot a_R/3,a_R\). Those of \(\dot a\) are \(\dot a_L,3(a_R-a_L)/\Delta t-\dot a_L-\dot a_R,\dot a_R\). Differencing successive controls gives derivative controls with the usual degree factor divided by \(\Delta t\). Their convex hulls bound the entire cubic and its derivatives. For \(D=a\,d/dt\),

\[
 Da=a\dot a,\quad D^2a=a\dot a^2+a^2\ddot a,\quad
 D^3a=a(\dot a^3+4a\dot a\ddot a+a^2\dddot a).
\]

Use these hulls, \(\Delta\eta\le\Delta t/\min a\) from the positive lower hull, and \(\omega\ge\kappa\) in the amplitude inequalities. One-sided \(D^2a\) values retain their signs when computing jumps. The script sums positive cellwise majorants, rounding upwards to multiples of \(10^{-12}\). It evaluates bounds for shells 33 through 2048 individually, comparing the whole interval with two blocks split at \(t/L=0.01\). This split isolates the short early acceleration pulse.

For \(n\ge2049\), let \(K=2049+3/2\). The whole-interval second bound on \(\kappa^3 |\beta|\), after replacing \(\omega\) by its lower bound \(\kappa\), is nonincreasing in \(\kappa\). Its value \(C\) at \(K\) bounds every remaining shell. Using \(g_n\le2\kappa^2\) and the integral bounds for decreasing powers gives

\[
 N_{n\ge2049}\le \frac{2C^2}{3(K-1)^3},\qquad
 L E_{n\ge2049}\le \frac{C^2\sqrt{a_T^{-2}+(m/K)^2}}{(K-1)^2}.
\]

Run `python check_tail.py` in this directory to reproduce [TAIL.json](TAIL.json). The bound concerns the exact interpolant through the supplied nodes. Mesh refinement of the underlying wall dynamics is a separate source of uncertainty.

## Why the readout time matters

At a nonstationary readout with \(a'(T)\ne0\) and a stationary initial endpoint, the leading boundary term survives. Since \(|\beta|=O(\kappa^{-1})\) uniformly from the sweep bound, \(\alpha(T)=1+O(\kappa^{-2})\) by its evolution equation. The second integration by parts above still bounds the residual by \(O(\kappa^{-3})\). Consequently

\[
 p_n(T)=\frac{m^2a'(T)^2}{16\kappa^4}+O(\kappa^{-5}),\qquad
 g_n\sqrt{m^2+\kappa^2/a_T^2}\,p_n(T)
 =\frac{m^2a'(T)^2}{8a_T\kappa}+O(\kappa^{-2}).
\]

For a massive spectator the gross instantaneous number still converges there, while the gross instantaneous energy diverges logarithmically. For example, on the archived expanding slice at \(t/L=11.03911\), the Hermite curve has \(a'\simeq1.50236\). Finite-band energies remain well-defined. Stationarity is a property of the radius interpolant itself. For the historical PCHIP radius construction, its proper-time derivative at the root of the separately interpolated Hubble rate is approximately \(-8.24\times10^{-6}\). The stationary endpoint proof above uses the explicitly constrained Hermite slope. The bounce radius for the constrained interpolant is evaluated from the native Hermite curve before inserting the zero-slope knot.

To enclose the expanding-slice number, evaluate \(a,a',a''\) at the exact readout inside its cubic cell. Put \(d=m|a'|/4\). The surviving boundary term has absolute value at most \(d/\kappa^2\). On a prefix all positive interior majorants are at most their full-interval totals. Add the new right-end \(|G|\) majorant to the full-interval second bound, which conservatively retains the old endpoint as well. If \(C_T\) denotes the resulting bound on the scaled residual \(\kappa^3|\mathrm{remainder}|\) for \(\kappa\ge K\), the integral test gives

\[
 N_{\kappa\ge K}\le
 2\left[\frac{d^2}{K-1}+\frac{dC_T}{(K-1)^2}
 +\frac{C_T^2}{3(K-1)^3}\right].
\]
