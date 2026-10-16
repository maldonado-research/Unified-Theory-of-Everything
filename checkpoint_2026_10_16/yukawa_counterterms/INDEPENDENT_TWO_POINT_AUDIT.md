# Independent pure-Yukawa two-point UV audit

15 October 2026, Pacific date; the directory uses the 16 October UTC label.
Input conventions: candidate.11, `0770c304`, and the literal
`ACTION_SIGN.md` reference action. This calculation independently evaluates
the pure-neutrino-Yukawa renormalizable two-point UV poles. Its large-momentum
expansion and radial Gamma-function route were derived before comparison
with the producer's Feynman-parameter calculation.

The two derivations agree. There is no chirality-flipping fermion loop pole
in this sector. There is, however, a nonzero Majorana **parameter** counterterm
from field renormalization, and the scalar two-point function has both a
kinetic pole and a retained-mass-squared pole. Those distinctions matter
when transporting a seed operator or interpreting a mass counterterm.

## Literal action, amplitudes and IR separation

Use `eta=(+---)`, `d=4-2epsilon`, `S=T exp(i integral L_int)`,
`kappa=1/(16 pi^2 epsilon)`, and propagators with `+i0`. The interactions are

    L_Y=-bar L_alpha^a Y_alpha i phi^a N_Ri+h.c.,
    L_N,mass=-1/2[bar N_L M N_R+bar N_R M^dagger N_L],
    phi=tilde H,   N_R=N_L^c,   M=M^T,
    L_phi,kin+mass=(partial phi)^dagger(partial phi)-m_H^2 phi^dagger phi.

Gauge interactions and other Yukawas are excluded from the computed sector.
The tree scalar mass is the signed parameter `m_H^2`. Start in a real
diagonal sterile-mass basis to evaluate the individual integrals; the
general complex symmetric result is fixed and checked below. No sterile
field is integrated out, and no mass inverse is used.

The formulas below denote the amputated **sum of loop graphs itself** by
`i Pi_loop`. This avoids a convention-dependent choice of whether an author
calls that graph `+i Sigma` or `-i Sigma`. Counterterm vertices are `+i`
times their action kernels and cancel the displayed loop graphs.

Keep the external momentum generic and off shell while extracting the UV
pole. In the asymptotic expansion introduce an independent positive auxiliary
mass squared `r` through `Q=q^2-r+i0`. It regulates the expanded vacuum
integrals at small q. All r-dependent UV residues cancel in the result.
Thus the computation does not infer a UV zero by setting a massless,
zero-external-momentum integral to zero in dimensional regularization.

## Independently evaluated angular and radial UV moments

For integers `a>=0,b>0`, Wick rotation gives the radial integral

    integral d^d q/(2pi)^d (q^2)^a/(q^2-r+i0)^b
      = i(-1)^(a-b)/(4pi)^(d/2)
        r^(d/2+a-b) Gamma(a+d/2) Gamma(b-a-d/2)
        /[Gamma(d/2)Gamma(b)].                              (T1)

The code extracts the Gamma pole at `epsilon=0`; no wavefunction coefficient
is inserted as a target. In units of `i kappa`, useful radial residues are

    integral 1/Q       : r,
    integral 1/Q^2     : 1,
    integral q^2/Q^2   : 2r,
    integral q^2/Q^3   : 1,
    integral q^4/Q^4   : 1.                                (T2)

Angular averaging uses
`q_mu q_nu -> eta_mu nu q^2/d`, with the corresponding higher isotropic
tensors when needed. The d-dependent prefactors are set to d=4 only after
the simple UV pole is identified. Their O(epsilon) corrections affect finite
terms, which are not calculated here.

For the open fermion loop, let

    D_1=q^2-m_N^2+i0,
    D_2=(q-p)^2-m_s^2+i0.

The only logarithmically divergent vector term after angular averaging is

    slash q/(D_1 D_2) -> 2 slash q(q.p)/Q^3,
    Res integral slash q/(D_1 D_2)=i kappa slash p/2.       (T3)

Both physical internal masses and the auxiliary r drop out. As an
independent algebraic cross-check, contract with p and use

    2 p.q=D_1-D_2+p^2+m_N^2-m_s^2.

The two tadpole residues and the bubble residue from (T1) give exactly the
same vector coefficient. This cross-check also uses no Feynman parameters.

For the scalar graph the internal lepton is massless. Its chiral trace is

    tr[P_R(slash q+m_N)P_L(slash q-slash p)]
       =2 q.(q-p).                                      (T4)

The mass term is killed by `P_R P_L=0`. Let `z=q.p`. Expanding both
denominators around Q through second order gives every possible UV term;
higher orders are convergent by radial degree. The nonzero residues are

    leading numerator 2q^2/Q^2                    : 4r,
    numerator 2q^2 times first denominator shift  : 2(-p^2-2r+m_N^2),
    numerator -2z times linear momentum shift    : -4p^2/d,
    numerator 2q^2 times quadratic momentum shift : +8p^2/d.

Their sum at d=4 is

    Res integral 2q.(q-p)/[D_1 (q-p)^2]
       =i kappa(2m_N^2-p^2).                            (T5)

In particular the auxiliary mass cancels before any physical conclusion.
Dropping the leading superficially quadratic term would leave a spurious
regulator dependence and the wrong mass residue. The numerator identity

    2q.(q-p)=D_1+(q-p)^2+m_N^2-p^2

provides a separate tadpole/bubble check of (T5).

## Chiral, weak and loop-sign factors

For an open self-energy chain, the two Yukawa vertices and two propagators
give `(-i)^2 i^2=+1`. The scalar graph has the additional closed-fermion-loop
minus. Consequently (T5) contributes the kernel `+p^2-2m_N^2` after that
minus is included. The factor two in (T4) is the actual chiral Dirac trace;
it must not be replaced by a full Dirac trace of four or halved again.

The weak vertex is `delta_ab`. For an external lepton or external scalar
component, contracting the internal weak index gives `delta_ab`: there is
no summed external-doublet multiplicity. For the sterile two-point
function both weak components circulate, giving `sum_a delta_aa=2`.

The external Majorana two-point function contains the right-chiral channel
and its required charge-conjugate channel. With `H=Y^dagger Y`, the first is
proportional to `H slash p P_R`; charge conjugation gives
`H^T slash p P_L`. They do not supply an additional common factor two.
The scalar loop with one external phi and one external phi-dagger already
has its fixed Yukawa orientation; copying the charge-conjugate description
as a second scalar loop would double-count that contraction.

There is no internal-M numerator contribution to the open L loop because
`P_R M P_L=0`. The open sterile loop has a massless chiral lepton internally,
so it also has no chirality-flipping numerator. A scalar contraction between
two same-charge Yukawa vertices is absent for the complex doublet. Therefore
no additional Majorana mass-type one-particle-irreducible pole is hidden in
these graphs.

## Computed two-point poles

Define `H=Y^dagger Y`, `T_Y=tr H`, and for a general complex mass
`T_M=tr(H M^dagger M)`. The amplitudes are

    i Pi_L,ab^UV=i kappa delta_ab (Y Y^dagger) slash p P_L/2,

    i Pi_N^UV=i kappa slash p[H P_R+H^T P_L],

    i Pi_phi,ab^UV=i kappa delta_ab[T_Y p^2-2 T_M].       (T6)

Their displayed flavor matrices are obtained by summing the actual two
Yukawa endpoints. The scalar result has no extra weak factor two. The
sterile result has the internal weak sum, the chiral conjugate channel and
the Majorana action normalization accounted for separately.

The counterterm field conventions are

    L_B=(1+delta Z_L/2)L,
    N_LB=(1+delta Z_NL/2)N_L,
    N_RB=(1+delta Z_NL*/2)N_R,
    phi_B=(1+delta Z_phi/2)phi.

Cancelling (T6) gives

    delta Z_L=-kappa Y Y^dagger/2,
    delta Z_NL=-kappa H^T,       delta Z_NR=-kappa H,
    delta Z_phi=-kappa T_Y.                                (T7)

The full Majorana kinetic counterterm has its action factor `1/2`; taking
both constrained slots produces the complete kernel in (T6). Equivalently
it is one ordinary left-Weyl kinetic term with coefficient `delta Z_NL`.
The independent exterior-algebra check verifies this equality and rejects
applying the Majorana half twice.

For clarity, distinguish a coefficient in the counterterm **action** from
a counterterm for a bare **parameter**. Write the scalar mass counterterm
action as `-Delta m_action^2 phi^dagger phi`. Its value is

    Delta m_action^2=-2 kappa T_M.                         (T8)

If `m_H,B^2=m_H^2+delta m_parameter^2`, then

    Delta m_action^2=delta m_parameter^2+m_H^2 delta Z_phi,
    delta m_parameter^2=kappa(m_H^2 T_Y-2 T_M).            (T9)

For the Majorana mass, the total mass-action coefficient is

    Delta M_action=delta M_parameter
        +(delta Z_NL^dagger M+M delta Z_NR)/2.

There is no chirality-flipping loop pole, so `Delta M_action=0`. Nonetheless

    delta M_parameter=kappa(H^T M+M H)/2.                 (T10)

It is symmetric, as required by `M=M^T`. Setting this parameter counterterm
to zero would leave the mass variation induced by the field counterterms
uncanceled. Equations (T8)-(T10) are pole counterterms, not finite threshold
corrections or directly measured mass shifts.

## General complex mass and flavor covariance

With the stated four-component mass orientation the inverse free sterile
kernel is `slash q-M P_R-M^dagger P_L`. Expanding its inverse at large q
shows that the mass-squared term in `P_R S_N P_L` contains **`M^dagger M`**.
That is why the scalar trace in (T6) is `tr(H M^dagger M)`. It must not be
replaced by `tr(H M M^dagger)` or by an ordinary complex `M^2`.
The earlier left-chiral source recurrence contains `M M^dagger` instead;
the difference is fixed by the propagator/source chirality, not a choice
to commute flavor matrices.

For `L_old=V L_new`, `N_Lold=U N_Lnew`, with constant unitary V and U,

    Y'=V^dagger Y U*,     M'=U^dagger M U*,
    H'=U^T H U*,
    delta Z_L'=V^dagger delta Z_L V,
    delta Z_NL'=U^dagger delta Z_NL U,
    delta M_parameter'=U^dagger delta M_parameter U*.

Both `T_Y` and `T_M` are invariant. These transformation laws are verified
for generic, rank-deficient and zero complex symmetric masses; no inverse
of M is needed. Real diagonal or real symmetric masses reduce to the
usual `tr(H M^2)` special case.

## Reproduction, controls and limits

Run with SymPy 1.14.0:

```sh
python -I -B -O independent_two_point_audit.py > INDEPENDENT_TWO_POINT_RESULTS.json
```

The standalone code imports no project calculation and uses explicit
failures, so optimization does not disable checks. It computes the radial
Gamma residues, expands the actual integrands, performs angular averages,
and records each asymptotic contribution. It then checks independent
numerator/tadpole reductions, auxiliary-mass cancellation, every chiral
trace coefficient, vertex/propagator phases, weak contractions, the
Majorana conjugate flavor channel, and a genuine exterior-algebra kinetic
normalization witness. Complex flavor data use three active and two sterile
flavors. The mass checks include nonnormal complex symmetric, singular and
zero matrices and actual propagator-series multiplication.

The optimized result is **91 exact controls passed** and **13 rejected
mutations**. The computed normalized integral residues are `1/2` for the
open vector moment and `2m_N^2-p^2` for the closed chiral trace before its
closed-loop minus. Negative controls omit required asymptotic pieces,
closed-loop or Majorana factors, transpose flavor channels, change the
complex-mass product, or confuse mass-action and bare-parameter coefficients.

These are independent computations of the declared renormalizable two-point
subsector, with internal AI-assisted algebra audits. They do not supply a
Yukawa vertex calculation, a complete beta function, additional seed-loop
residues, all gauge/other-Yukawa contributions, finite matching, or an
observable. A field-induced seed counterterm is only one contribution to
the complete seed renormalization and must be labelled accordingly.
