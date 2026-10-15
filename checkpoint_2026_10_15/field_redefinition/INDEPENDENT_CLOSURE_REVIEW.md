# Independent review of the lower-order closure witness

15 October 2026. AI-assisted analytic and component review. No producer code is imported. Scope: the declared classical two-component Weyl action, through inverse-scale order three; this is a concrete omission witness, not a new theorem, loop matching result or completed EFT basis.

Verdict: the displayed field-change expansion is correct within its stated no-order-one-shift scope, with the dimensional and local-normalization premises made explicit below.

For constant Hermitian `Z`, a constant unitary matrix diagonalizes `K=I+t^2 rho Z` at every spacetime point. In that eigenbasis the exact local canonical factor is the scalar function `S_a=(1+t^2 rho z_a)^(-1/2)`. Thus `S_a K_a S_a=1`, and its kinetic connection is `S_a K_a partial S_a-(partial S_a)K_a S_a=0` even for nonconstant `rho`. Expanding this exact map gives the producer's `S=I-t^2 rho Z/2+O(t^4)`. The cancellation follows from this single constant-matrix family; it does not extend automatically to several noncommuting field-dependent kinetic matrices.

For anticommuting left Weyl components, the contracted bilinear `n_i n_j` is symmetric in flavor. Varying its two field slots gives `delta(n^T C n)=(delta n)^T C n+n^T C delta n`, with no new Grassmann minus between these two already ordered slot variations. Inserting `delta n=-(t^2 rho/2) Z n` produces `-(t^2 rho/2)n^T(Z^T C+C Z)n`. Multiplication by the input `t rho` therefore yields the stated `Delta D=-(Z^T C+C Z)/2`. The same slot calculation produces the stated mass and source shifts. A Hermitian adjoint would be wrong for these holomorphic congruences. No mass inverse is introduced, so singular or vanishing `M` is harmless.

The separate executable `independent_grassmann_closure_check.py` reconstructs this slot variation as an exterior algebra of four Weyl Grassmann components rather than expanding the producer's matrix congruence. It checks the Majorana factor `(n_1 n_1)=2 n_1^1 n_1^2`, disappearance of antisymmetric flavor coefficients, a rank-one complex symmetric mass, an incorrect-dagger control, and exact spectral normalization/connection cancellation. All seven checks pass under optimized Python.

A concrete independent fixture is

    Z = [[1, i], [-i, 2]],
    C = [[2, 1+i], [1+i, 3i]],
    M = [[1, i], [i, -1]],          rank(M)=1,
    Delta D = [[-3+i, -3-5i/2], [-3-5i/2, 1-7i]].

The nonzero contact survives this singular mass fixture; its origin is the dimension-five source acted on by the dimension-six field normalization.

Necessary premises to add or make explicit:

- With the declared dimensionless `Y` and canonical `[n]=3/2`, the source must satisfy `[J]=5/2`. It should be a Grassmann-odd, left-Weyl, gauge-singlet light-field composite independent of `n` and with no implicit inverse-scale expansion relevant here. An elementary Weyl source would instead require a dimension-one `Y`.
- This is a formal perturbative normalization near `t=0`, or a finite local normalization in a region where `K` is positive and invertible. Arbitrary Hermitian `Z` does not guarantee a globally positive kinetic matrix at arbitrarily large `rho`. This is a kinetic-domain condition, separate from singular `M`.
- The source and `rho` are held as functions of unchanged light fields. No claim about a quantum Jacobian or physical coefficient invariance is required by this classical omission witness.

Dimensions otherwise agree: `[t]=-1`, `[rho]=2`, each term of the mass matrix has dimension one, `[Z]=[C]=[D]=0`, and the order-two mass-shaped and order-three scalar contacts have the required dimensions. The unitary covariance follows from `U* U^T=I`: `Z'=U^dagger Z U`, `C'=U^T C U`, and `Delta D'=U^T Delta D U`.

Run:

    python3 -I -B -O independent_grassmann_closure_check.py

Receipt: `independent_grassmann_closure_check.json`. This review does not replay any older private manual packet or claim their generic field-change bookkeeping as new theory.

