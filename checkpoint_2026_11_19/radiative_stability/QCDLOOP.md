# Independent scalar triangles for the quartic correction

QCDLoop 2.1.0 was built from commit `52d94ce9e628c1b1f9784f66a6ef9de285359ae3` and run separately for the sixteen scalar triangle integrals used as reference values. With invariants $(0,0,s)$ and squared masses $(1,b,1)$, the two unit-mass scalar lines face the $s$ channel. QCDLoop's Minkowski triangle is the negative of the positive-simplex $T_b$ used here. Values at $s=-1$ check that Euclidean sign and normalization. Values at $s=6,8,10$ check both components above threshold for $b=1,16,100,400$, all in units $m=1$. The returned Laurent poles vanish. The driver evaluates each integral in binary128 and prints the finite part as binary64. [Library source](https://github.com/scarrazza/qcdloop/tree/52d94ce9e628c1b1f9784f66a6ef9de285359ae3) · [QCDLoop article](https://doi.org/10.1016/j.cpc.2016.07.033)

From a scratch directory containing the driver, use GNU C++, libquadmath, Bash and the tagged source on Linux.
```bash
git clone --depth 1 --branch 2.1.0 https://github.com/scarrazza/qcdloop.git
test "$(git -C qcdloop rev-parse HEAD)" = 52d94ce9e628c1b1f9784f66a6ef9de285359ae3 || exit 1
mkdir -p build/qcdloop
sed 's/@VERSION@/2.1.0/g' qcdloop/src/qcdloop/config.h.in > build/qcdloop/config.h
g++ -std=gnu++11 -fext-numeric-literals -O2 -Ibuild -Iqcdloop/src \
  qcdloop_triangles.cpp qcdloop/src/{triangle,cache,topology,tools,types}.cc \
  -lquadmath -o build/qcdloop_triangles
./build/qcdloop_triangles --output QCDLOOP_TRIANGLES.json
```
The offline Python replay compares complex triangle values with the retained library output. Tiny imaginary roundoff below threshold is included in the complex error norm. The quartic vertex factor, the pole-mass subtraction and the four-lepton normalization have independent derivations in the research note. [C++ driver](qcdloop_triangles.cpp) · [External values](QCDLOOP_TRIANGLES.json) · [Amplitude calculation](RESCATTERING.md)
