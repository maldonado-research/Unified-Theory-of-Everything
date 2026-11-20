// Standalone cross-check with QCDLoop 2.1.0. Build instructions in QCDLOOP.md.
#include <qcdloop/box.h>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using Complex = std::complex<double>;
using Integral = ql::Box<ql::qcomplex, ql::qdouble, ql::qdouble>;
constexpr double pi = 3.14159265358979323846;

struct Node { double x, w; };

// Gauss-Legendre rule on [0,1]; reflection gives the other half of the nodes.
std::vector<Node> quadrature(int n) {
    std::vector<Node> rule(n);
    for (int i = 0; i < n/2; ++i) {
        double root = std::cos(pi*(i+0.75)/(n+0.5)), derivative = 0;
        for (int iteration = 0; iteration < 20; ++iteration) {
            double p0 = 1, p1 = root;
            for (int j = 2; j <= n; ++j) {
                double next = ((2*j-1)*root*p1-(j-1)*p0)/j;
                p0 = p1;
                p1 = next;
            }
            derivative = n*(root*p1-p0)/(root*root-1);
            double step = p1/derivative;
            if (std::abs(step) < 2e-16) break;
            root -= step;
        }
        double weight = 1/((1-root*root)*derivative*derivative);
        rule[i] = {(1-root)/2, weight};
        rule[n-i-1] = {(1+root)/2, weight};
    }
    return rule;
}

ql::qcomplex scalar_box(Integral& evaluator, double s, double t, double b, double c) {
    std::vector<ql::qcomplex> result(3);
    evaluator.integral(result, 1, {1,b,1,c}, {0,0,0,0,s,t});
    if (cabsq(result[1])+cabsq(result[2]) != 0)
        throw std::runtime_error("Unexpected Laurent pole in the massive box");
    return result[0];
}

Complex as_double(ql::qcomplex value) {
    return {static_cast<double>(crealq(value)), static_cast<double>(cimagq(value))};
}

Complex amplitude(Integral& evaluator, double s, double z) {
    ql::qcomplex sum = 0;
    for (int b : {1,16}) {
        int muon = b == 1 ? 1 : -1;
        for (int c : {100,400}) {
            int electron = c == 100 ? 10 : -40;
            auto direct = scalar_box(evaluator, s, -s*(1-z)/2, b, c);
            auto exchange = scalar_box(evaluator, s, -s*(1+z)/2, b, c);
            sum += muon*electron*(direct+exchange);
        }
    }
    return as_double(sum)*s/(16*pi*pi);
}

struct Rates { double full, absorptive, j0; };

Rates integrate(Integral& evaluator, double s, int order) {
    double total = 0, cut = 0, mean_im = 0;
    for (auto point : quadrature(order)) {
        Complex value = amplitude(evaluator, s, point.x);
        total += point.w*std::norm(value);
        cut += point.w*value.imag()*value.imag();
        mean_im += point.w*value.imag();
    }
    double flux = 32*pi*s;
    return {total/flux, cut/flux, mean_im*mean_im/flux};
}

double euclidean_simplex() {
    // s=-1,t=-2, original parameter order (scalar,scalar,muon,electron).
    const auto rule = quadrature(32);
    double value = 0;
    for (auto u : rule) for (auto v : rule) for (auto w : rule) {
        double x1 = u.x, x2 = (1-x1)*v.x, left = (1-x1)*(1-v.x);
        double x3 = left*w.x, x4 = left*(1-w.x);
        double denominator = x1+x2+4*x3+9*x4+x1*x2+2*x3*x4;
        value += u.w*v.w*w.w*(1-x1)*left/(denominator*denominator);
    }
    return value;
}

void require_close(double observed, double expected, double rel) {
    if (!std::isfinite(observed) || std::abs(observed-expected) > rel*std::abs(expected))
        throw std::runtime_error("Numerical reference check failed");
}

void print_rates(std::ostream& out, Rates r) {
    out << "{\"full\": " << r.full << ", \"absorptive\": " << r.absorptive
        << ", \"J0_absorptive\": " << r.j0 << "}";
}

int main(int argc, char** argv) {
    if (argc != 3 || std::string(argv[1]) != "--output") {
        std::cerr << "Usage: qcdloop_reference --output QCDLOOP.json\n";
        return 2;
    }
    Integral evaluator;
    Complex euclidean = as_double(scalar_box(evaluator, -1, -2, 4, 9));
    double independent = euclidean_simplex();
    require_close(euclidean.real(), independent, 5e-13);
    if (!(std::abs(euclidean.imag()) < 1e-28))
        throw std::runtime_error("Unexpected Euclidean imaginary part");
    std::ostringstream out;
    out << std::setprecision(15)
        << "{\n  \"library\": \"QCDLoop 2.1.0\",\n"
        << "  \"git_commit\": \"52d94ce9e628c1b1f9784f66a6ef9de285359ae3\",\n"
        << "  \"publication_doi\": \"10.1016/j.cpc.2016.07.033\",\n"
        << "  \"scalar_integrals\": \"GNU binary128, accumulated before conversion to binary64\",\n"
        << "  \"normalized_rate\": \"m^2 sigma, M/m=10, y=g=1\",\n"
        << "  \"central_cos_theta\": 0,\n"
        << "  \"euclidean_box_s_t_b_c\": [-1,-2,4,9],\n"
        << "  \"euclidean_box_real\": " << euclidean.real() << ",\n"
        << "  \"euclidean_simplex_real\": " << independent << ",\n"
        << "  \"rows\": [\n";
    for (int s : {6,10}) {
        Complex central = amplitude(evaluator, s, 0);
        Rates coarse = integrate(evaluator, s, 8), fine = integrate(evaluator, s, 20);
        require_close(coarse.full, fine.full, 1e-12);
        require_close(coarse.absorptive, fine.absorptive, 1e-12);
        require_close(coarse.j0, fine.j0, 1e-12);
        out << "    {\n      \"s_over_m2\": " << s << ",\n"
            << "      \"central_real\": " << central.real() << ",\n"
            << "      \"central_imaginary\": " << central.imag() << ",\n"
            << "      \"angle_order_8\": ";
        print_rates(out, coarse);
        out << ",\n      \"angle_order_20\": ";
        print_rates(out, fine);
        out << "\n    }" << (s == 6 ? "," : "") << "\n";
    }
    out << "  ]\n}\n";
    std::ofstream receipt(argv[2]);
    receipt.exceptions(std::ios::failbit | std::ios::badbit);
    receipt << out.str();
}
