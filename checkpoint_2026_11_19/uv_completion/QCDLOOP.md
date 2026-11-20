# Independent QCDLoop reference for the finite boxes

QCDLoop 2.1.0 evaluates the scalar integrals in the leading polarized $\mu^-_L\mu^-_L\to e^-_Le^-_L$ amplitude. The calculation uses the tagged upstream source at commit `52d94ce9e628c1b1f9784f66a6ef9de285359ae3`. Stefano Carrazza, R. Keith Ellis and Giulia Zanderighi describe the library in *Computer Physics Communications* **209** (2016), 134–143. [Upstream source](https://github.com/scarrazza/qcdloop/tree/52d94ce9e628c1b1f9784f66a6ef9de285359ae3) · [Journal article](https://doi.org/10.1016/j.cpc.2016.07.033)

The scalar squared mass is one in units of $m^2$. Write $b$ and $c$ for the squared masses of the muon and electron singlets in these units. The QCDLoop squared-mass vector is $(1,b,1,c)$ and its six external invariants are $(0,0,0,0,s,t)$. In the library's modified Cayley matrix $2Y_{13}=2-s$ and $2Y_{24}=b+c-t$. This places the two scalars opposite one another across the physical $s$ cut. In the Euclidean check $(s,t,b,c)=(-1,-2,4,9)$ the scalar integral is positive and agrees with a separate smooth Feynman-simplex quadrature. Both Laurent poles returned by the library vanish for all sample points. [QCDLoop box conventions](https://github.com/scarrazza/qcdloop/blob/52d94ce9e628c1b1f9784f66a6ef9de285359ae3/src/box.cc#L156-L186)

At $M/m=10$ and $y=g=1$, the muon residues for $b=1,16$ are $(1,-1)$ and the electron residues for $c=100,400$ are $(10,-40)$. The driver sums the finite $\epsilon^0$ boxes for $t=-s(1-z)/2$ and $u=-s(1+z)/2$, using $z=\cos\theta$. It multiplies the weighted sum by $s/(16\pi^2)$ for the external phase $B_{34}^*B_{12}=s$. The rate is $m^2\sigma=\int_0^1 |\mathcal M|^2 dz/(32\pi s)$ in these units, including the identical-electron factor. QCDLoop evaluates the boxes and the driver sums their signed contributions in GNU binary128. Normalization and angular quadrature use binary64. The 8- and 20-point Gauss–Legendre rules check convergence independently of the direct Feynman-parameter integration. [Standalone driver](qcdloop_reference.cpp) · [Reference values](QCDLOOP.json)

| $s/m^2$ | $\mathrm{Re}\,\mathcal M(s,0)$ | $\mathrm{Im}\,\mathcal M(s,0)$ | $m^2\sigma$ |
|---:|---:|---:|---:|
| 6 | $-1.4143753842\,10^{-4}$ | $-6.5067129383\,10^{-5}$ | $4.0092537001\,10^{-11}$ |
| 10 | $-1.9338112607\,10^{-4}$ | $-1.8757070705\,10^{-4}$ | $7.1740070858\,10^{-11}$ |

## Reproduce the reference

Use an x86-64 Linux system with Git, Bash, GNU C++ and libquadmath. The run used GCC 14.2.0. Start in a scratch directory containing a copy of `qcdloop_reference.cpp`. The commands compile the upstream box implementation and its dependencies directly. The generated version header uses the substitution from upstream CMake. An independent checkout makes the version check explicit.

```bash
git clone --depth 1 --branch 2.1.0 https://github.com/scarrazza/qcdloop.git
test "$(git -C qcdloop rev-parse HEAD)" = 52d94ce9e628c1b1f9784f66a6ef9de285359ae3 || exit 1
mkdir -p build/qcdloop
sed 's/@VERSION@/2.1.0/g' qcdloop/src/qcdloop/config.h.in > build/qcdloop/config.h
g++ -std=gnu++11 -fext-numeric-literals -O2 -Ibuild -Iqcdloop/src \
  qcdloop_reference.cpp qcdloop/src/{box,cache,topology,tools,types}.cc \
  -lquadmath -o build/qcdloop_reference
./build/qcdloop_reference --output QCDLOOP.json
```

This external evaluation tests the scalar-integral values and their assembly for the specified model. The spin, flavor and identical-particle normalization still come from the benchmark derivation. The JSON records the finite full rate, the rate from the squared imaginary part and its $J=0$ contribution.
