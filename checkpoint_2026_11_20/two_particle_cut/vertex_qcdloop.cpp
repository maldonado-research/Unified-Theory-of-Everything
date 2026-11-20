// Scalar-integral evaluation of the on-shell-production triangle at a=m²=1.
#include <qcdloop/bubble.h>
#include <qcdloop/triangle.h>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

int main() {
    using Number=ql::qdouble;
    using Complex=ql::qcomplex;
    ql::Triangle<Complex,Number,Number> triangle;
    ql::Bubble<Complex,Number,Number> bubble;
    std::cout << std::setprecision(16);
    for (Number b : {1,16,100,400}) {
        for (Number q : {.125,1.,4.,9.}) {
            Number a=1, t=-q;
            std::vector<Complex> tri(3), left(3), right(3);
            triangle.integral(tri,1,{a,b,0},{0,a,t});
            bubble.integral(left,1,{b,0},{a});
            bubble.integral(right,1,{a,0},{t});
            std::vector<Complex> integral(3);
            for (int j=0;j<3;++j)
                integral[j]=(left[j]-right[j]+(2*a-b-t)*tri[j])/(t-a);
            if(cabsq(integral[1])+cabsq(integral[2])>1e-23Q ||
               fabsq(cimagq(integral[0]))>1e-18Q || !finiteq(crealq(integral[0])))
                throw std::runtime_error("unexpected pole, imaginary part or nonfinite value");
            std::cout << static_cast<int>(b) << " " << static_cast<double>(q)
                      << " " << static_cast<double>(crealq(integral[0])) << "\n";
        }
    }
}
