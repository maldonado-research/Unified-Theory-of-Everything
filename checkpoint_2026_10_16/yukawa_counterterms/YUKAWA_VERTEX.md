# Proper neutrino-Yukawa vertex and the renormalizable pure-Yukawa benchmark

15 October 2026, Pacific client date; the containing checkpoint uses
16 October UTC. Reference action and charge conjugation are those of
[`ACTION_SIGN.md`](../../checkpoint_2026_10_15/operator_map/ACTION_SIGN.md).
This is a known one-loop **renormalizable pure-neutrino-Yukawa** benchmark.
It is not a new TOE prediction, a seed-operator beta function, or a finite
matching calculation.

The proper vertex is checked here independently by enumerating its actual
allowed contractions and projecting the propagator numerators. The running
coefficients additionally use the separately calculated
[two-point pole inputs](YUKAWA_TWO_POINT.md). The
[independent vertex audit](INDEPENDENT_VERTEX_AUDIT.md) enumerates literal
interaction slots and checks the same result through a separately authored
implementation.

## 1. Action, external field orientation, and proper topology

At the unbroken scalar origin keep

    L0=bar L i slash partial L+(partial phi)^dagger partial phi
       −m_H^2 phi^dagger phi+(1/2)bar N(i slash partial−M)N
       −[bar L Y phi N_R+h.c.],
    phi=tilde H,       Y=Y_nu,       y_phi=y_L=−1/2,
    N=N^c,             N_L=N_R^c.                        (V1)

There is no seed insertion in this benchmark. Set other Yukawas and gauge
couplings to zero for this sector. A scalar quartic or the other Standard
Model couplings have their own contributions outside the stated pure-`Y`
result. Retain the real positive sterile mass when classifying propagators;
it is not silently set to zero to obtain the vertex result.

Call `A` the cubic monomial `bar L phi N_R` and `B` its Hermitian
conjugate `bar N_R phi^dagger L`. In independent left-Weyl variables,
`A` has two barred Weyl fermions and `B` two unbarred ones. The target
dimension-four coefficient has ordered external fields
`(bar L_alpha^a,phi^b,N_Ri)`.

At order `Y^3`, three cubic vertices and three internal lines make one
loop. A connected one-particle-irreducible three-point graph must have
one external field at each vertex. A vertex with two external fields has
only one internal edge, which is a bridge and makes the graph reducible.
With internal degree two at each of three vertices, the only connected
proper topology is the triangle. Before imposing field-type compatibility,
the checker enumerates all connected three-edge cubic multigraphs: one
proper topology remains and twelve connected reducible topology candidates
are excluded. These are labelled topology candidates, not twelve allowed
Yukawa diagrams or additional amplitude symmetry factors. The actual
allowed proper contractions are enumerated separately below.

## 2. The target orientation has no allowed Yukawa triangle

The external `bar L` and `phi` force their vertices to be `A`. The target
external `N_R` coefficient also forces `A` at its third vertex. After
removing those external fields, an `A-A-A` triangle would require

    <phi phi>,       <bar L bar L>,       <N_R N_R>.      (V2)

The first two propagators vanish. The scalar pair has total hypercharge
`−1` and the light-fermion pair `+1`; neither is an allowed propagator in
the unbroken reference action. The light field also has no Majorana mass.
The sterile Majorana anomalous propagator cannot repair either missing
charged contraction.

Exhausting both choices of internal endpoint assignment at all three
vertices therefore gives **zero allowed target-orientation triangles**.
This proves the absence of the local proper `bar L phi N_R` UV pole
without importing a known beta function or canceling UV and IR poles in
a scaleless integral.

The conjugate ordered amplitude obeys the conjugate argument. Reversing
the bookkeeping direction of a Majorana line does not create an extra
independent contraction.

## 3. The other Majorana attachment is allowed but UV finite

For an external four-component Majorana field, its opposite Weyl source
can attach at `B`. Fixing `bar L` and `phi` as above leaves exactly one
inequivalent `A-A-B` triangle, with propagators

    <N_R N_R>,       <phi phi^dagger>,       <L bar L>.   (V3)

This has the opposite external Majorana orientation (`bar N_R`, or the
equivalent `N_L` source), rather than the target scalar Yukawa monomial.
Its UV convergence must still be checked rather than discarding it by
name.

Let `S_N(q)` denote the sterile numerator `slash q+M`, with propagator
denominators displayed separately. Using `C=i gamma^2 gamma^0`, the
same-chirality contraction is

    P_R(slash q+M) C^T P_R^T
       =M P_R C^T P_R^T.                               (V4)

The slash numerator vanishes between the two right endpoints. The other
fermion line is massless and contributes one `slash q`. Thus the stripped
open-chain tensor has the form `M slash q P_L`, with the external source
orientation fixed as above. In explicit Weyl blocks, the checker contracts
the anomalous sterile epsilon, the light-fermion propagator, and the
Majorana source conversion and obtains precisely this one-power momentum
structure. An overall amplitude phase is unnecessary for the UV-zero
claim, and no finite form factor is evaluated here.

There are three quadratic denominators and only one loop momentum in the
numerator. The superficial UV degree is

    omega=4+1−2(3)=−1.                                  (V5)

A normal two-slash Yukawa triangle would instead have degree zero, but it
is forbidden by (V2). Extra mass-squared denominator insertions or external
momentum expansion terms further reduce (V5). The allowed graph has no
logarithmic UV divergence even before using the oddness of its leading
momentum numerator. Its weak contraction is `delta_ab`, with no extra
factor of two from a closed weak trace.

For UV/IR separation take generic nonexceptional spacelike external
momenta, or temporarily give the denominators a common positive auxiliary
mass. At large Euclidean loop momentum the radial behavior is bounded by
`integral dq q^(−2)` in four dimensions; it remains UV convergent near
`d=4−2epsilon`. No prescription that sets the original scaleless graph
to zero is used. Removing an IR regulator can affect finite or infrared
terms, not the ultraviolet power-counting conclusion. Finite Majorana
form factors are not claimed to vanish.

Accordingly the proper dimension-four pure-`Y` vertex pole is

    Gamma^(1),UV_(bar L phi N_R)=0.                      (V6)

External-leg self-energy insertions are reducible. Their field
counterterms must instead be included in the bare-parameter relation
below, once each.

## 4. Two-point inputs and parameter versus action counterterms

Write `ell_loop=1/(16 pi^2)`, `kappa=ell_loop/epsilon`, and

    H=Y^dagger Y,       H_L=Y Y^dagger,
    T=Tr H,            Q_M=Tr(H M^dagger M).             (V7)

The independently calculated two-point poles give

    delta Z_L=−kappa H_L/2,
    delta Z_NR=−kappa H,
    delta Z_NL=−kappa H^T,
    delta Z_phi=−kappa T.                               (V8)

The same calculation finds no proper sterile chirality-flipping mass
pole: its internal light fermion is massless, and its same-oriented
scalar/light contractions are absent. Its scalar two-point loop pole is
`+i kappa[T p^2−2Q_M]`, so the raw scalar mass coefficient in the action
must obey `delta m_coeff^2=−2kappa Q_M`. These integral coefficients are
inputs from `YUKAWA_TWO_POINT.md`, not new integrations in this note.

Define the bare fields and parameters in fixed flavor coordinates by

    L0=(1+delta Z_L/2)L,
    N_R0=(1+delta Z_NR/2)N_R,
    N_L0=(1+delta Z_NR*/2)N_L,
    phi0=(1+delta Z_phi/2)phi,
    Y0=mu^epsilon(Y+delta Y),
    M0=M+delta M,       m_H0^2=m_H^2+delta m_H^2.         (V9)

The **total action coefficient** of the Yukawa counterterm is

    delta Y_action=delta Y+(delta Z_L^dagger)Y/2
                         +Y delta Z_NR/2+delta Z_phi Y/2.

Equation (V6) requires that total coefficient to be zero. It does not
make the bare-parameter pole zero. Substituting (V8) gives

    delta Y=kappa rY,
    rY=H_L Y/4+Y H/2+T Y/2
       =3Y H/4+T Y/2.                                  (V10)

The equality `H_L Y=Y H` is associativity, not a commutation of unrelated
flavor matrices. The scalar field trace is included once.

For the literal Majorana action
`−(1/2)[bar N_L M N_R+h.c.]`, the total mass coefficient is

    delta M_action=delta M+delta Z_NR^T M/2+M delta Z_NR/2.

Its proper loop mass pole is zero, hence

    delta M=kappa rM,
    rM=(H^T M+M H)/2.                                   (V11)

The transpose follows from Majorana completion; replacing it by an
adjoint gives a different result and generally destroys symmetry of the
mass coefficient. There is no additional Majorana half to append.

The scalar potential has the literal sign `−m_H^2 phi^dagger phi`.
Its raw mass action coefficient includes its field normalization:

    delta m_coeff^2=delta m_H^2+m_H^2 delta Z_phi,
    delta m_H^2=kappa r_m2,
    r_m2=m_H^2 T−2Q_M.                                  (V12)

Thus the raw action mass pole and the scalar mass parameter pole differ.
Equations (V11)-(V12) also hold for a general complex symmetric `M`, using
the declared Majorana convention and `M^dagger M` in `Q_M`. Their unitary
covariance is checked without any mass inverse.

## 5. Deriving the canonical-epsilon factors and running

The running parameters here are dimensionful `M,m_H^2` and dimensionless
`Y`; no additional `mu` power is introduced for either mass parameter.
The `d=4−2epsilon` Yukawa has the explicit `mu^epsilon` in (V9).

At leading order `beta_Y^(d)=−epsilon Y`. Both `Y` and `Y*` have that
canonical scaling. The residue `rY` is homogeneous of degree three in
these Yukawa variables, while `rM,r_m2` are degree two. Expanding the
bare equation at first loop order,

    0=mu d/dmu[mu^epsilon(Y+ell_loop rY/epsilon)],
    beta_Y^(d)=−epsilon Y+ell_loop bY,

gives `0=ell_loop[rY+bY−3rY]`, hence `bY=2rY`.
For a mass parameter there is no prefactor term `+r`: its bare equation
instead gives `bM−2rM=0` or `b_m2−2r_m2=0`. Derivatives involving the
one-loop mass running inside a residue enter only at the next loop order.

The complete four-dimensional running **within the declared pure-neutrino-
Yukawa one-loop sector** is therefore

    16 pi^2 beta_Y = (3/2)Y(Y^dagger Y)+Tr(Y^dagger Y)Y,
    16 pi^2 beta_M = (Y^dagger Y)^T M+M(Y^dagger Y),
    16 pi^2 beta_m_H^2
       =2m_H^2 Tr(Y^dagger Y)−4Tr(Y^dagger Y M^dagger M). (V13)

As a separate check, (V8) gives
`gamma_L=H_L/(32pi^2)`, `gamma_NR=H/(16pi^2)`, and
`gamma_phi=T/(16pi^2)`. With a zero proper vertex pole,
`gamma_L Y+Y gamma_NR+gamma_phi Y` agrees with (V13).

For one active and one sterile flavor, these formulas give
`16pi^2 beta_y=(5/2)y^3` for real `y`,
`16pi^2 beta_M=2y^2 M`, and
`16pi^2 beta_m_H^2=2y^2 m_H^2−4y^2 M^2`.
If `M=0`, its running remains zero in this sector. These are known
renormalizable consistency checks, not new predictions of the hypothesis.

## 6. Reproduction and boundaries

Run `python -I -S -B -O verify_yukawa_vertex.py` here. The checker uses
only the standard library and exact Gaussian rational arithmetic. It
enumerates cubic proper/reducible topologies and all internal endpoint
assignments of the literal `A,B` vertices, projects full Majorana
propagator numerators with explicit gamma matrices, contracts the allowed
open-chain Weyl tensors, and checks ultraviolet degrees and weak flow.

The subsequent parameter/RG controls state their two-point inputs
explicitly. They use rectangular noncommuting Yukawa matrices and a
complex symmetric mass, check actual field-factor expansion of the bare
action, verify Yukawa-degree homogeneity, and test active/sterile unitary
covariance, the zero-mass limit, and one-flavor coefficients. Negative
controls omit chiral projectors or field factors, replace a mass numerator
by a second slash, permit charged same-orientation propagators, append a
weak factor of two, miss the canonical-epsilon factor, or replace the
Majorana transpose by an adjoint.

The completed optimized isolated replay passes **80 positive checks**
and rejects **21 mutations**, for **101 checks total**. These are
implementation controls within the stated derivations; the two-point
integrals remain separately verified inputs.

Equation (V13) provides a renormalizable running benchmark for subsequent
calculations. A complete Standard Model/sterile EFT calculation additionally
requires the other declared interaction sectors and the operator
counterterms. In particular it supplies no complete `C7` or selected-seed
beta function, no finite photon amplitude, and no new fundamental discovery.
