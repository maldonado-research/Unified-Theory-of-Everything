// Positive-simplex T is minus the Minkowski scalar triangle returned by QCDLoop.
#include <qcdloop/triangle.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    if (argc != 3 || std::string(argv[1]) != "--output") {
        std::cerr << "Usage: qcdloop_triangles --output QCDLOOP_TRIANGLES.json\n";
        return 2;
    }
    std::ofstream out(argv[2]);
    out.exceptions(std::ios::failbit | std::ios::badbit);
    out << std::setprecision(16)
        << "{\n  \"qcdloop_commit\": \"52d94ce9e628c1b1f9784f66a6ef9de285359ae3\",\n"
        << "  \"integral\": \"minus QCDLoop C0(0,0,s;1,b,1), m=1\",\n  \"rows\": [\n";
    ql::Triangle<ql::qcomplex, ql::qdouble, ql::qdouble> triangle;
    bool first = true;
    for (int energy : {-1, 6, 8, 10}) {
        for (int mass_squared : {1, 16, 100, 400}) {
            std::vector<ql::qcomplex> result(3);
            const ql::qdouble b = mass_squared;
            const ql::qdouble s = energy;
            triangle.integral(result, 1, {1, b, 1}, {0, 0, s});
            if (cabsq(result[1]) + cabsq(result[2]) != 0)
                throw std::runtime_error("Unexpected Laurent pole");
            if (!first) out << ",\n";
            first = false;
            out << "    {\"s\": " << energy << ", \"b\": " << mass_squared
                << ", \"real\": " << static_cast<double>(-crealq(result[0]))
                << ", \"imag\": " << static_cast<double>(-cimagq(result[0])) << "}";
        }
    }
    out << "\n  ]\n}\n";
}
