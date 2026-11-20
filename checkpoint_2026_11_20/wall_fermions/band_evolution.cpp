// Fourth-order commutator-free unitary evolution for independent Majorana modes.
// The Python driver supplies two Gauss-node samples on each Hermite interval.
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using Complex = std::complex<double>;
using State = std::array<Complex, 2>;
using Step = std::array<double, 5>;  // dt, a_early, adot_early, a_late, adot_late
using Axis = std::array<double, 2>;  // real off-diagonal component and z component
constexpr Complex I{0.0, 1.0};

Axis axis(double k, double mass, double a, double adot, bool rotating) {
    const double p = k / a;
    if (!rotating) return {p, mass};
    const double omega2 = p * p + mass * mass;
    const double dtheta = -mass * p * (adot / a) / (2 * omega2);
    return {dtheta, std::sqrt(omega2)};
}

void advance(State& y, double dt, Axis early, Axis late, double w1, double w2,
             bool rotating) {
    const double transverse = w1 * early[0] + w2 * late[0];
    const double z = w1 * early[1] + w2 * late[1];
    const double omega = std::sqrt(z * z + transverse * transverse);
    const double c = std::cos(omega * dt);
    const double s = std::sin(omega * dt) / omega;
    // In the rotating basis, -iH has real off-diagonal entries dtheta,-dtheta.
    const Complex mix = rotating ? Complex(s * transverse) : -I * (s * transverse);
    const Complex diagonal(c, -s * z);
    y = {diagonal * y[0] + mix * y[1], -std::conj(mix) * y[0] + std::conj(diagonal) * y[1]};
}

double probability(const State& y, double k, double mass, double a, bool rotating) {
    if (rotating) return std::norm(y[0]);
    const double theta = 0.5 * std::atan2(k / a, mass);
    return std::norm(std::cos(theta) * y[0] + std::sin(theta) * y[1]);
}

struct Job {
    int n;
    double mass;
    bool rotating;
    double bounce, first, last, norm_error;
};

void evolve(Job& job, const std::vector<Step>& grid, size_t bounce, size_t split,
            double a0, double ab, double a1, double a2) {
    const double k = job.n + 1.5;
    const double theta = 0.5 * std::atan2(k / a0, job.mass);
    State y = job.rotating ? State{0., 1.} : State{-std::sin(theta), std::cos(theta)};
    const double w1 = 0.25 + std::sqrt(3.) / 6.;
    const double w2 = 0.5 - w1;
    for (size_t j = 0; j < grid.size(); ++j) {
        const Step& row = grid[j];
        const Axis early = axis(k, job.mass, row[1], row[2], job.rotating);
        const Axis late = axis(k, job.mass, row[3], row[4], job.rotating);
        advance(y, row[0], early, late, w1, w2, job.rotating);
        advance(y, row[0], early, late, w2, w1, job.rotating);
        if (j + 1 == bounce)
            job.bounce = probability(y, k, job.mass, ab, job.rotating);
        if (j + 1 == split)
            job.first = probability(y, k, job.mass, a1, job.rotating);
    }
    job.last = probability(y, k, job.mass, a2, job.rotating);
    job.norm_error = std::abs(std::norm(y[0]) + std::norm(y[1]) - 1);
}

int main(int argc, char** argv) {
    if (argc != 9) return 2;
    std::ifstream input(argv[1], std::ios::binary | std::ios::ate);
    if (!input) return 3;
    const auto bytes = static_cast<size_t>(input.tellg());
    if (bytes % sizeof(Step)) return 4;
    std::vector<Step> grid(bytes / sizeof(Step));
    input.seekg(0);
    if (!input.read(reinterpret_cast<char*>(grid.data()), bytes)) return 5;
    const size_t bounce = std::stoul(argv[2]), split = std::stoul(argv[3]);
    const double a0 = std::stod(argv[4]), ab = std::stod(argv[5]);
    const double a1 = std::stod(argv[6]), a2 = std::stod(argv[7]);
    const bool samples_only = std::string(argv[8]) == "samples";
    std::vector<Job> jobs;
    for (double mass : {1., 4.}) {
        for (int n = 32; n <= 2048; ++n) {
            bool spot = n == 32 || n == 33 || n == 128 || n == 512 || n == 1024 || n == 2048;
            if (!samples_only || spot) jobs.push_back({n, mass, true, 0., 0., 0., 0.});
            if (spot) jobs.push_back({n, mass, false, 0., 0., 0., 0.});
        }
    }
    std::atomic<size_t> cursor{0};
    auto worker = [&]() {
        for (size_t j = cursor.fetch_add(1); j < jobs.size(); j = cursor.fetch_add(1))
            evolve(jobs[j], grid, bounce, split, a0, ab, a1, a2);
    };
    std::vector<std::thread> threads;
    const unsigned workers = std::max(1u, std::min(4u, std::thread::hardware_concurrency()));
    for (unsigned i = 1; i < workers; ++i) threads.emplace_back(worker);
    worker();
    for (auto& thread : threads) thread.join();
    std::cout << std::setprecision(17);
    for (const auto& j : jobs)
        std::cout << j.n << ' ' << j.mass << ' ' << j.rotating << ' '
                  << j.bounce << ' ' << j.first << ' ' << j.last << ' ' << j.norm_error << '\n';
}
