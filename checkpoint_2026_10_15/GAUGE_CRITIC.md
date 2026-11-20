# Independent Feynman-parameter audit of the internal gauge attachments

This audit derives the two specified routed integrals directly by shifting the
Feynman-parameter denominator and contracting the original numerator. It does
not infer a transverse coefficient from a Ward identity. The metric is
`(+---)`, `d=4-2 epsilon`, and every displayed residue has
`i/(16 pi² epsilon)` stripped. Charges and the common interaction prefactor are
not included in `F` and `S`. A signed parameter `m²` appears only through the
scalar denominators. These conventions fix integral residues, not an action
counterterm sign or a physical beta function.

The exact input integrands are

\[
F^\nu=\operatorname{UV}\int_q
\frac{[2(p-k)^2+(q-h-k)^2](\slashed q-\slashed k)
\gamma^\nu\slashed q}
{q^2(q-k)^2[(q-p)^2-m^2]},
\]
\[
S^\nu=\operatorname{UV}\int_q
\frac{[2(p-k)^2+(q-h)^2]\slashed q(2p-2q-k)^\nu}
{q^2[(q-p)^2-m^2][(q-p+k)^2-m^2]}.
\]

## Shift and moment derivation

For the fermion attachment use
`a=x k+y p` and
`Delta=y m²-x k²-y p²+a²`.
For the scalar attachment use
`a=x p+y(p-k)` and
`Delta=(x+y)m²-x p²-y(p-k)²+a²`.
In both cases `q=ell+a`, `x,y>=0`, and `x+y<=1`.
The denominator is `(ell²-Delta)^3`; its Feynman factor `2!` is
included in the normalized simplex moments

\[
\langle x^a y^b\rangle=\frac{2a!b!}{(a+b+2)!}.
\]

At the UV pole,

\[
\operatorname{UV}\int_\ell\frac{1}{(\ell^2-\Delta)^3}=0,
\qquad
\operatorname{UV}\int_\ell\frac{\ell^\mu\ell^\nu}
{(\ell^2-\Delta)^3}=\frac{g^{\mu\nu}}4,
\]
\[
\operatorname{UV}\int_\ell
\frac{\ell^\mu\ell^\nu\ell^\rho\ell^\sigma}
{(\ell^2-\Delta)^3}
=\frac{\Delta}{8}
(g^{\mu\nu}g^{\rho\sigma}+g^{\mu\rho}g^{\nu\sigma}
+g^{\mu\sigma}g^{\nu\rho}).
\]

These follow from the scalar UV residues
`UV integral ell²/(ell²-Delta)^3=1` and
`UV integral (ell²)^2/(ell²-Delta)^3=3 Delta`, with the
`d`-dimensional tensor denominators `d` and `d(d+2)` evaluated at
`d=4` for the simple pole. Terms of order `epsilon` can affect finite terms,
which are outside this calculation. The code applies these tensor moments
directly to every monomial of each shifted matrix-valued numerator.

As a compact intermediate check, define
`A=a-k`, `B=a`, `b=a-h-k`, `c=b²+2(p-k)²` for `F`. Its
parameter integrand after UV extraction is

\[
\slashed A\gamma^\nu\slashed B
+\tfrac12(\slashed A\gamma^\nu\slashed b+
\slashed b\gamma^\nu\slashed B)
-\tfrac12(3\Delta+c)\gamma^\nu.
\]

For `S`, take `b=a-h`, `c=b²+2(p-k)²`, `v=2p-k-2a`; the
corresponding integrand is

\[
\slashed a v^\nu+\tfrac12\slashed b v^\nu-
\slashed a b^\nu-\tfrac12(3\Delta+c)\gamma^\nu.
\]

## Complete off-shell residues

Define the explicitly antisymmetric Clifford tensor

\[
E(a,b)^\nu=\slashed a\gamma^\nu\slashed b
-a^\nu\slashed b-b^\nu\slashed a+(a\cdot b)\gamma^\nu
=a_\mu b_\rho\gamma^{[\mu}\gamma^\nu\gamma^{\rho]}.
\]

With unit-weight antisymmetrization, the result is

\[
\begin{aligned}
F^\nu={}&\left[-\frac{m^2}{2}-k^2-\frac{7p^2}{6}
+\frac52p\cdot k-\frac56 k\cdot h+\frac23p\cdot h
-\frac{h^2}{2}\right]\gamma^\nu\\
&+\frac23p^\nu\slashed p
-\frac12(k^\nu\slashed p+p^\nu\slashed k)
+\frac16(k^\nu\slashed h+h^\nu\slashed k)
-\frac13(p^\nu\slashed h+h^\nu\slashed p)
+\frac12[E(k,h)^\nu-E(k,p)^\nu],
\end{aligned}
\]

\[
\begin{aligned}
S^\nu={}&\left[-m^2-p^2+2p\cdot k-\frac56 k^2
+\frac23p\cdot h-\frac13k\cdot h-\frac{h^2}{2}\right]
\gamma^\nu\\
&-\frac16k^\nu\slashed k
-\frac16(2p-k)^\nu\slashed h
+\frac13h^\nu(2\slashed p-\slashed k).
\end{aligned}
\]

For `Q_L=-Q_H`, the charge-weighted sum is `Q_H(S-F)`. A useful
short form for its stripped tensor is

\[
\begin{aligned}
(S-F)^\nu={}&\left[\frac{p^2+k^2}{6}-\frac{m^2}{2}\right]
\gamma^\nu+h^\nu\slashed p-\frac23p^\nu\slashed p
+\frac12k^\nu\slashed h-\frac16k^\nu\slashed k\\
&+\frac12\slashed k\gamma^\nu(\slashed p-\slashed h).
\end{aligned}
\]

The mass term is essential in the internal sum. It cancels on subtracting the
specified local three-point representative; there is no mass-dependent
transverse remainder for these two routed integrals and this representative.

## Ward anchor computed independently

The two-propagator numerator gives

\[
T(p,h)=\operatorname{UV}\int_q
\frac{[2p^2+(q-h)^2]\slashed q}{q^2[(q-p)^2-m^2]}
=\left[m^2+p^2-\frac23p\cdot h+\frac{h^2}{2}\right]\slashed p
+\left[\frac{p^2}{6}-\frac{m^2}{2}\right]\slashed h.
\]

This uses the separate shift `q=ell+x p`,
`Delta=x m²-x(1-x)p²`, with denominator power two. The UV
moments are `1` at rank zero and `Delta g^{mu nu}/2` at rank two.
Writing `t=p-k` and `I(p)=slashed p/2`, the individual identities are

\[
k_\nu F^\nu=T(t,h)-T(p,h+k)+2(p^2-t^2)I(p),
\]
\[
k_\nu S^\nu=T(t,h)-T(p,h)+2(p^2-t^2)I(p).
\]

Their charge-weighted difference therefore obeys

\[
k_\nu(S-F)^\nu=T(p,h+k)-T(p,h)
=\left[\frac{p^2}{6}-\frac{m^2}{2}\right]\slashed k
+\left[k\cdot h-\frac23p\cdot k+\frac{k^2}{2}\right]\slashed p.
\]

## Transverse subtraction and on-shell projection

For the specified local representative,

\[
G_J^\nu=\left[-\frac23p+h+\frac12k\right]^\nu\slashed p
+\left[\frac{p^2}{6}-\frac{m^2}{2}\right]\gamma^\nu.
\]

Set `r=p-h-k`. Direct subtraction, before any on-shell condition, gives

\[
\boxed{\Delta^\nu=(S-F-G_J)^\nu
=\frac12[r^\nu\slashed k-(r\cdot k)\gamma^\nu+E(k,r)^\nu]
+\frac13[k^\nu\slashed k-k^2\gamma^\nu].}
\]

Both bracketed tensors are separately transverse. The antisymmetric
three-gamma term is present off shell and cannot be discarded by a Ward test.
For a real transverse gauge boson, `k²=k·epsilon=0`,

\[
\Delta\cdot\epsilon=\frac12\slashed k\slashed\epsilon\slashed r.
\]

If the outgoing massless lepton obeys `ubar slashed r=0`, anticommuting
`slashed r` to the left yields

\[
\bar u\,\Delta\cdot\epsilon
=\bar u[(r\cdot\epsilon)\slashed k-(r\cdot k)\slashed\epsilon].
\]

This on-shell projection loses the separate divergence tensor. It is a check
of the displayed residue and does not determine a complete operator basis or
the sourced equations-of-motion descendants.

## Reproduction and mutation checks

Run `python gauge_critic.py` in this directory with SymPy installed and without
Python `-O` or `-OO`; optimized mode is explicitly rejected because assertions
are scientific checks. The code
checks all eight triangle components as exact polynomials in generic symbolic
four-vectors, computes the bubble separately, checks both individual Ward
identities and their charge sum, verifies the subtraction and exact mass
cancellation, and tests the on-shell relation on four exact left-kernel rows
covering two distinct momentum configurations.

The intentionally wrong variants exercise meaningful limitations:

* Equal internal charges, omitted scalar attachment, a common sign flip of the
  triangles relative to the independently anchored bubble, and omission of the
  internal mass term each fail the Ward comparison.
* Dropping the three-gamma term passes the Ward identity but fails both the full
  residue and the on-shell coefficient.
* Dropping the divergence tensor passes the Ward identity and the real-gauge-leg
  on-shell check but fails the generic off-shell residue.

A common redefinition of every action and vertex convention would not be
detected by the internal Ward comparison. Four-dimensional matrix checks also
do not identify evanescent finite terms. Neither issue alters the explicitly
declared simple-pole integral calculation.

Validation receipt, 15 October 2026: `python -u gauge_critic.py` completed with
exit code 0 under SymPy 1.14.0; all symbolic component, bubble, Ward, subtraction,
four-row on-shell and mutation checks passed. `python -O gauge_critic.py`
returned exit code 1 with the intended explicit optimization-rejection error.
