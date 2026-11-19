// External scalar-continuum benchmark. Build and normalization in EXTERNAL.md.
#include <qcdloop/box.h>
#include <array>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using Complex = std::complex<double>;
using Box = ql::Box<ql::qcomplex, ql::qdouble, ql::qdouble>;
constexpr double pi = 3.14159265358979323846;
constexpr double loop = 16*pi*pi;
struct Node { double x, weight; };
using Rule = std::vector<Node>;

Rule legendre(int n) {
    // Only even orders are used. Reflect each root of P_n onto [0,1].
    Rule result(n);
    if (n < 2 || n % 2) throw std::invalid_argument("expected positive even order");
    for (int i = 0; i < n/2; ++i) {
        double x = std::cos(pi*(i+0.75)/(n+0.5)), derivative = 0;
        for (int step = 0; step < 30; ++step) {
            double previous = 1, current = x;
            for (int k = 2; k <= n; ++k) {
                double next = ((2*k-1)*x*current-(k-1)*previous)/k;
                previous = current;
                current = next;
            }
            derivative = n*(x*current-previous)/(x*x-1);
            double dx = current/derivative;
            if (std::abs(dx)<2e-16) break;
            x -= dx;
        }
        const double w = 1/((1-x*x)*derivative*derivative);
        result[i] = {(1-x)/2,w};
        result[n-i-1] = {(1+x)/2,w};
    }
    return result;
}

Complex amplitude(Box& box, double s, double z, double U) {
    std::vector<ql::qcomplex> result(3);
    ql::qcomplex sum = 0;
    for (int b : {1,16}) for (int d : {100,400}) {
        for (double t : {-s*(1-z)/2,-s*(1+z)/2}) {
            box.integral(result,1,{U,ql::qdouble(b),1,ql::qdouble(d)},{0,0,0,0,s,t});
            if (cabsq(result[1])+cabsq(result[2]) != 0)
                throw std::runtime_error("unexpected Laurent poles");
            sum += (b==1 ? 1 : -1)*(d==100 ? 10 : -40)*result[0];
        }
    }
    Complex value(double(crealq(sum)),double(cimagq(sum)));
    if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
        throw std::runtime_error("nonfinite box sum");
    if (U>std::pow(std::sqrt(s)-1,2) && std::abs(value.imag())>1e-15/U)
        throw std::runtime_error("unexpected imaginary part above the mass endpoint");
    return s/loop*value;
}

// Above the cut the projective simplex has a positive real denominator.
// Integrate the heavy propagator's parameter analytically and the two
// remaining coordinates by Gauss-Legendre. This returns U times a bare box.
double projective_scaled_box(double U, double s, double t, int b, int d, const Rule& rule) {
    const double q=1/U;
    double value=0;
    for (auto h : rule) {
        const double xe=h.x*h.x;
        for (auto k : rule) {
            const double xm=(1-xe)*k.x, xs=(1-xe)*(1-k.x);
            const double L=xs+b*xm+d*xe, C=L-t*xm*xe, B=1+q*(L-s*xs);
            const double D=B*B-4*q*C, root=std::sqrt(D);
            const double J=std::log((B+root)*(B+root)/(4*q*C))/root;
            value+=h.weight*k.weight*2*h.x*(1-xe)*(B/C-2*q*J)/D;
        }
    }
    return value;
}
double high_mass_relative_difference(Box& box, double U) {
    const double s=8;
    const Rule rule=legendre(64);
    double sum=0;
    for (int b : {1,16}) for (int d : {100,400}) {
        sum+=2*(b==1 ? 1 : -1)*(d==100 ? 10 : -40)
             *projective_scaled_box(U,s,-s/2,b,d,rule);
    }
    const double external=U*amplitude(box,s,0,U).real()*loop/s;
    const double difference=std::abs(sum/external-1);
    if (!std::isfinite(difference) || difference>2e-11)
        throw std::runtime_error("high-mass projective box comparison failed");
    return difference;
}

double density(double U, int B) {
    if (U<B) return 0;
    if (B==1) return 1/U;
    const double ratio=(U-B)/(U-1);
    return ratio*ratio/U;
}

struct Coeff { Complex g{}, y{}; };
Coeff continuum(Box& box, double s, double z, const Rule& rule) {
    Coeff value;
    const double end=std::pow(std::sqrt(s)-1,2);
    std::array<double,5> ends={1,end,16,100,400};
    auto add = [&](double U, double weight) {
        const Complex weighted=amplitude(box,s,z,U)*(2*weight/loop);
        value.g += weighted*(density(U,1)+density(U,16)/4);
        value.y += weighted*(density(U,100)+2*density(U,400));
    };
    for (size_t j=1; j<ends.size(); ++j) {
        const double lo=ends[j-1], span=std::log(ends[j]/lo);
        for (auto p : rule) {
            const double v=std::sin(pi*p.x/2);
            const double U=lo*std::exp(span*v*v);
            add(U,p.weight*U*span*pi*std::sin(pi*p.x)/2);
        }
    }
    for (auto p : rule) {
        const double U=400/(p.x*p.x);
        add(U,p.weight*2*U/p.x);
    }
    return value;
}

struct Row {
    Complex leading{}, mean_g{}, mean_y{};
    Coeff central;
    double sigma0=0, slope_g=0, slope_y=0;
};
Row compute(Box& box, double s, int n, int angle) {
    const Rule rule=legendre(n);
    Row row;
    row.leading=amplitude(box,s,0,1);
    row.central=continuum(box,s,0,rule);
    for (auto p : legendre(angle)) {
        const Complex born=amplitude(box,s,p.x,1);
        const Coeff delta=continuum(box,s,p.x,rule);
        row.mean_g+=p.weight*delta.g;
        row.mean_y+=p.weight*delta.y;
        const double flux=p.weight/(32*pi*s);
        row.sigma0+=flux*std::norm(born);
        row.slope_g+=2*flux*std::real(std::conj(born)*delta.g);
        row.slope_y+=2*flux*std::real(std::conj(born)*delta.y);
    }
    return row;
}

void write_complex(std::ostream& out, Complex value) {
    out << "[" << value.real() << ", " << value.imag() << "]";
}
void write_row(std::ostream& out, int s, int n, int angle, const Row& r) {
    out << "    {\"s_over_m2\": " << s << ", \"mass_order\": " << n
        << ", \"angle_order\": " << angle << ", \"leading_central\": ";
    write_complex(out,r.leading);
    out << ", \"delta_g_central\": "; write_complex(out,r.central.g);
    out << ", \"delta_y_central\": "; write_complex(out,r.central.y);
    out << ", \"delta_g_mean\": "; write_complex(out,r.mean_g);
    out << ", \"delta_y_mean\": "; write_complex(out,r.mean_y);
    out << ", \"m2_sigma0\": " << r.sigma0 << ", \"m2_sigma_slope_g\": " << r.slope_g
        << ", \"m2_sigma_slope_y\": " << r.slope_y << "}";
}
int main(int argc, char** argv) {
    if (argc!=3 || std::string(argv[1])!="--output") {
        std::cerr << "Usage: qcdloop_continuum --output reference.json\n";
        return 2;
    }
    Box box;
    std::ofstream out;
    out.exceptions(std::ios::badbit|std::ios::failbit);
    out.open(argv[2]);
    out << std::setprecision(16) << "{\n  \"library\": \"QCDLoop 2.1.0\",\n"
        << "  \"upstream_commit\": \"52d94ce9e628c1b1f9784f66a6ef9de285359ae3\",\n"
        << "  \"doi\": \"10.1016/j.cpc.2016.07.033\",\n"
        << "  \"amplitude_components\": [\"real\", \"imaginary\"],\n"
        << "  \"high_mass_projective_check\": [\n";
    bool initial=true;
    for (double U : {4e4,1e7,2e8}) {
        if (!initial) out << ",\n";
        initial=false;
        out << "    {\"U\": " << U << ", \"s_over_m2\": 8, \"cos_theta\": 0, "
            << "\"relative_difference\": " << high_mass_relative_difference(box,U) << "}";
    }
    out << "\n  ],\n  \"rows\": [\n";
    bool first=true;
    for (auto orders : {std::array<int,2>{10,4}, {18,6}, {28,10}}) {
        for (int s : {6,8,10}) {
            if (!first) out << ",\n";
            first=false;
            write_row(out,s,orders[0],orders[1],compute(box,s,orders[0],orders[1]));
            out.flush();
        }
    }
    out << "\n  ]\n}\n";
}
