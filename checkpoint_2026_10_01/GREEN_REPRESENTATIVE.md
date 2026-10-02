# A local Green-operator representative of the routed seed pole

## Scope

This calculation starts from the public routed integral in
`checkpoint_2026_09_06/routed_loop/ROUTED_BUBBLE.md` and its separate seed-vertex
review. It gives a local gauge-covariant representative of that **one** pole,
and reduces it with explicitly sourced equations of motion (EOM). It is not a
complete dimension-seven basis, a full counterterm calculation, or a physical
beta function. The additional operators exhibited below demonstrate a concrete
limitation of matching only a three-point Green function with no external gauge
field. No claim of a new fundamental interaction is made.

## 1. Conventions and local operators

Use metric `(+---)`, Fourier waves `exp(-ip.x)` for incoming sterile fields and
`exp(+ir.x)`, `exp(+ih.x)` for the outgoing lepton and Higgs wave functions in
the displayed monomial. Thus `p=r+h`. Use a Majorana sterile field `N=N^c`,
`N_L=P_L N`, `N_R=P_R N`, and a positive diagonal sterile mass matrix `M`.
The outgoing lepton is left chiral. Define

\[
\widetilde H=i\tau^2H^*,\quad
D_\mu=\partial_\mu+i\Gamma_\mu,\quad
[D_\mu,D_\nu]=i\mathcal F_{\mu\nu},\quad
\sigma^{\mu\nu}=\frac i2[\gamma^\mu,\gamma^\nu].
\]

Here `tau^2` is a Pauli matrix; `mathcal F` includes the gauge coupling and the
generators in the representation on which it acts. For the lepton doublet,
`Y_L=-1/2`. For rows, `D_mu bar L = partial_mu bar L - i bar L Gamma_mu`.
The derivative on `N` is ordinary because it is a gauge singlet. Let

\[
A=i\!\not\!\partial,\qquad
\mathcal K\bar X=-i(D_\mu\bar X)\gamma^\mu.
\]

The expressions below specify derivative placement; parentheses are part of
the definition. Flavor indices are suppressed, but every monomial contains
`bar L_alpha ... N_beta`. Define

\[
\begin{aligned}
O_N&=\bar L\widetilde H\, A(-\partial^2)N_L,\\
O_L&=(-D^2\bar L)\widetilde H\, A N_L,\\
O_H&=\bar L(-D^2\widetilde H)\, A N_L,\\
O_R&=(\mathcal K\bar L)\widetilde H\,(-\partial^2)N_L,\\
O_{mN}&=m_H^2\bar L\widetilde H\, A N_L,\\
O_{mL}&=m_H^2(\mathcal K\bar L)\widetilde H\,N_L.
\end{aligned}
\]

The product `bar L tilde H` means contraction of electroweak indices, leaving
the displayed spinor row. In `O_R`, for example, the gamma matrix contracts
the derivative on `bar L` and acts on `N_L`; it does not act on the scalar.

Every monomial has mass dimension seven: the three fields have total
dimension four, and the remaining dimension is either three derivatives or
one derivative plus a mass squared. `m_H^2` is the signed scalar mass
parameter entering the propagator `q^2-m_H^2+i0`. All monomials are Lorentz
and electroweak scalars. The hypercharges of `bar L` and `tilde H` cancel;
their doublet and antidoublet indices contract. Since
`bar L=bar ell P_R` and `P_R gamma^mu=gamma^mu P_L`, the odd gamma chain is
compatible with `N_L` and is not zero. Its conjugate must also be included
in a Hermitian Lagrangian.

At zero gauge field their ordered three-point kernels, stripped of the
common vertex factor `i`, are respectively

\[
\left(p^2\not p,\quad r^2\not p,\quad h^2\not p,\quad
p^2\not r,\quad m_H^2\not p,\quad m_H^2\not r\right)P_L.
\]

In particular, `A` on the incoming sterile field gives `slash p`, whereas
`mathcal K` on the outgoing row gives `slash r`. Both `-partial^2` and
`-D^2` give a positive squared external momentum. Using `+i D_mu bar L`
instead would reverse the two `slash r` kernels.

Consequently the local representative is

\[
\boxed{O_J=\frac56O_N+\frac13O_L+\frac16O_H-\frac16O_R
 +\frac12O_{mN}+\frac12O_{mL}.}
\]

Its three-point kernel equals the public result

\[
\not J=\left(\frac56p^2+\frac13r^2+\frac16h^2+
\frac12m_H^2\right)\not p+
\left(\frac12m_H^2-\frac16p^2\right)\not r.
\]

This is a choice of representative, not a claim that these six monomials
form an independent or complete Green basis. Integration by parts and EOM
relations have deliberately not been imposed in defining the map.

### Pole versus counterterm sign

Let `F=C_(7,seed) Y_nu^*`, of mass dimension minus three. The source defines
the seed by its positive momentum-squared vertex, and gives

\[
i\mathcal M_{\rm loop}^{\rm UV}
=-\frac{i}{16\pi^2\epsilon}F\,\bar u_L\not J P_Lu_N
\]

with its indicated electroweak epsilon tensor understood. The mixed-field
counterterm action

\[
\mathcal L_{\rm ct}=+\frac{1}{16\pi^2\epsilon}
F_{\alpha\beta}(O_J)_{\alpha\beta}+\mathrm{h.c.}
\]

has vertex `+i F J/(16 pi^2 epsilon)` and cancels this pole in the declared
convention. There is no extra identical-field factor `1/2` in this
mixed `bar L-H-N` monomial. This fixes its action-to-vertex map, but it does
not replace the source's warning that a literal `+J^T W_1 Box J/2` seed
has the opposite Fourier sign from a definition using positive `W_1 q^2`.
That separate seed convention must be translated before identifying
`C_(7,seed)` with a coefficient in any other action.

## 2. Sourced EOM, including curvature and the scalar mass

For definiteness the light scalar kinetic and potential terms are

\[
\mathcal L_H=(D_\mu H)^\dagger D^\mu H
-m_H^2H^\dagger H-\lambda(H^\dagger H)^2+\mathcal L_{\rm int}.
\]

Define sources from the *same declared reference action* by

\[
i\!\not D L=S_L,\qquad
(D^2+m_H^2+2\lambda H^\dagger H)H=S_H,
\qquad S_H=\frac{\delta\mathcal L_{\rm int}}{\delta H^\dagger}.
\]

For renormalizable lepton interactions
`-bar L Y_e H e_R - bar L Y_nu tilde H N_R + h.c.`, one has
`S_L=Y_e H e_R+Y_nu tilde H N_R`. The scalar source includes lepton and
quark Yukawa terms. Writing it as a functional derivative retains all of
them without imposing an arbitrary zero-source restriction. If lower
effective operators are present at the power-counting order in question,
their functional derivatives must also be included. The identities here
remain true for the resulting total sources; choosing only the
renormalizable sources does not constitute reduction of a complete
sequential EFT packet.

Dirac adjunction and the squared covariant Dirac operator give

\[
\mathcal K\bar L=\bar S_L,\qquad
(i\!\not D)^2=-D^2-\frac12\sigma^{\mu\nu}\mathcal F_{\mu\nu},
\]

and hence

\[
\boxed{-D^2\bar L=\mathcal K\bar S_L+
\frac12\bar L\mathcal F_{\mu\nu}\sigma^{\mu\nu}.}
\]

The curvature sign follows from
`gamma^mu gamma^nu=g^(mu nu)-i sigma^(mu nu)` and
`[D_mu,D_nu]=i mathcal F_(mu nu)`. Changing the sign convention for `D`
requires changing the definition of `mathcal F` consistently.

With `S_tildeH=i tau^2 S_H^*`, the Higgs equation is

\[
\boxed{-D^2\widetilde H=m_H^2\widetilde H+
2\lambda(H^\dagger H)\widetilde H-S_{\widetilde H}.}
\]

Substituting these equations into the local representative gives

\[
\begin{aligned}
O_J\simeq{}&\frac56\bar L\widetilde H A^3N_L
+\frac23m_H^2\bar L\widetilde H A N_L
+\frac\lambda3(H^\dagger H)\bar L\widetilde H A N_L\\
&+\frac13(\mathcal K\bar S_L)\widetilde H A N_L
+\frac16\bar L\mathcal F_{\mu\nu}\sigma^{\mu\nu}
 \widetilde H A N_L
-\frac16\bar L S_{\widetilde H} A N_L\\
&-\frac16\bar S_L\widetilde H A^2N_L
+\frac12m_H^2\bar S_L\widetilde H N_L.
\end{aligned}
\]

Here `simeq` denotes use of the reference-action EOM, as appropriate for a
field redefinition at first order in the inserted coefficient. It is not
an equality of arbitrary off-shell Green functions. For non-diagonal
flavor tensors the contractions inherited from the original operators
must be kept; the source notation does not license commuting flavor
matrices. The scalar-mass coefficient is `1/2+1/6=2/3`, rather than `1/2`:
replacing `h^2` by zero while retaining a nonzero scalar mass would miss
this descendant. In the genuinely massless free-light limit all sources,
curvature, `lambda`, and `m_H^2` vanish, recovering the source's `5/6`
projection. A symmetric-phase negative `m_H^2` is a signed parameter in
these algebraic identities, not an instruction to use tachyonic external
states as physical asymptotic particles.

### Sterile source tower

Define the two chiral source equations in the fixed mass basis by

\[
A N_L=M N_R+S_R^N,\qquad A N_R=M N_L+S_L^N.
\]

For the renormalizable neutrino Yukawa interaction they are the mutually
charge-conjugate currents built from `Y_nu^dagger tilde H^dagger L` and
its conjugate. Their definition by these equations fixes their signs.
The constant matrix `M` commutes with the spacetime derivative, giving

\[
\begin{aligned}
A^2N_L&=M^2N_L+M S_L^N+A S_R^N,\\
A^3N_L&=M^3N_R+M^2S_R^N+M A S_L^N+A^2S_R^N.
\end{aligned}
\]

Thus the first line of the reduced representative includes a
Yukawa-shaped monomial with matrix
`F[(5/6)M^3+(2/3)m_H^2 M]`, **and also** its correlated source operators.
This matrix records one seed's pole contribution in this explicitly
chosen representative. It is not a beta function. In particular,
`A^3 N_L -> M^3 N_R` discards three nonzero source terms. The two remaining
`A` and `A^2` actions in the preceding displayed reduction must be reduced
with the same identities if a fully expanded source packet is desired.

## 3. An explicit ambiguity invisible to this three-point calculation

Consider the gauge-invariant dimension-seven operator

\[
Q_B=\bar L\widetilde H\gamma^\mu i\partial^\nu N_L B_{\mu\nu}
+\mathrm{h.c.}
\]

The abelian field strength is gauge invariant. This operator has no
`bar L-H-N` vertex at zero external gauge field, so `O_J+a Q_B` has the
same three-point polynomial for every coefficient `a`. A factor of `g'`
may be absorbed into the definition of `a`; that bookkeeping choice does
not remove the null direction of the three-point restriction map.
Its one-gauge-boson kernel is proportional to

\[
\bar u_L(r)\left[\not k\,(p\cdot\varepsilon)
-\not\!\varepsilon\,(p\cdot k)\right]P_Lu_N(p).
\]

It obeys the photon Ward test (`epsilon -> k` makes it zero), but it is
not identically zero on shell. For example, in consistent momentum units
take `p=(2,0,0,0)`, `r=(1,0,0,1)`,
`h=k=(1/2,0,0,-1/2)`, and `epsilon=(0,1,0,0)`.
Then `p=r+h+k`, `p^2=4`, and `r^2=h^2=k^2=0`.
The spin-summed squared stripped kernel is `4` (or `2` averaged over the
two initial sterile spin states) with standard completeness relations.
Thus the ambiguity can change an amplitude. This observation does not
assert that `Q_B` belongs to a particular minimal published basis; any
EOM/IBP reduction must preserve its nonzero amplitude by mapping it to
other operators.

This supplies a precise obstruction to *uniquely reconstructing the
gauge-completed EFT* from the three-point polynomial alone. It does not
prove that `Q_B` contributes to the same coupling-order scalar beta
coefficient as this bubble. That question needs the chosen power
counting and additional diagrams. An actual running coefficient also
requires the complete action and normalization, all counterterms at the
claimed order, field and parameter renormalization, source contacts, and
a subtraction scheme. The public Schur cancellation remains compatible
with any common proposed coefficient and cannot fill these omissions.

## 4. Reproduction and falsification

`verify_green_representative.py` uses only the Python standard library.
It checks the Fourier map by exact rational arithmetic, the two forms of
the routed polynomial, the massless and massive-light on-shell
projections, the sterile source tower as a formal polynomial, and an
explicit Clifford-algebra/spin-trace control. Mutation checks require
wrong row-Fourier signs, a dropped scalar-mass descendant, a dropped
sterile contact, and the wrong curvature sign to fail. Checks use
explicit exceptions and remain enabled with `python -O`.

These tests validate the displayed conventions and algebra; they do not
integrate an additional loop, establish basis completeness, supply
external human review, or measure an observable. The useful next
calculation is matching Green functions with an external gauge field
and the source-contact external legs under one fixed full EFT action.
