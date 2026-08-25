#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "spatial_discretization.h"

namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

double first_derivative_wave_number(int mode, double spacing) {
    const double angle = 2.0 * pi * static_cast<double>(mode) / static_cast<double>(N);
    if constexpr (spatial::order == 2) {
        return std::sin(angle) / spacing;
    } else if constexpr (spatial::order == 4) {
        return (8.0 * std::sin(angle) - std::sin(2.0 * angle)) / (6.0 * spacing);
    } else {
        return (45.0 * std::sin(angle) - 9.0 * std::sin(2.0 * angle)
              + std::sin(3.0 * angle)) / (30.0 * spacing);
    }
}

bool close(double value, double expected, double tolerance = 2.e-11) {
    return std::abs(value - expected) <= tolerance * std::max(1.0, std::abs(expected));
}

} // namespace

int main() {
    constexpr double box_length = 1.4;
    constexpr double spacing = box_length / static_cast<double>(N);
    constexpr int mx = 5;
    constexpr int my = 3;
    constexpr int mz = 2;

    std::vector<double> field(static_cast<std::size_t>(N) * N * N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                const double phase = 2.0 * pi * (mx * i + my * j + mz * k) / N;
                field[spatial::index(i, j, k)] = std::cos(phase);
            }
        }
    }

    const int i = N - 1;
    const int j = 0;
    const int k = N / 3;
    const double value = field[spatial::index(i, j, k)];
    const double phase = 2.0 * pi * (mx * i + my * j + mz * k) / N;
    const double k2 = spatial::effective_momentum_squared(mx, my, mz, spacing);

    const double laplacian = spatial::laplacian(i, j, k, field) / (spacing * spacing);
    if (!close(laplacian, -k2 * value)) {
        std::fprintf(stderr, "order %d: Laplacian eigenvalue mismatch\n", spatial::order);
        return 1;
    }

    const double qx = first_derivative_wave_number(mx, spacing);
    const double qy = first_derivative_wave_number(my, spacing);
    const double dx_field = spatial::first_derivative(0, i, j, k, field, spacing);
    if (!close(dx_field, -qx * std::sin(phase))) {
        std::fprintf(stderr, "order %d: first derivative mismatch\n", spatial::order);
        return 1;
    }

    const double dxx_field = spatial::second_derivative(0, i, j, k, field, spacing);
    const double kx2 = spatial::effective_component_squared(mx, spacing);
    if (!close(dxx_field, -kx2 * value)) {
        std::fprintf(stderr, "order %d: same-axis second derivative mismatch\n", spatial::order);
        return 1;
    }

    const double dxy_field = spatial::mixed_derivative(0, 1, i, j, k, field, spacing);
    if (!close(dxy_field, -qx * qy * value)) {
        std::fprintf(stderr, "order %d: mixed derivative mismatch\n", spatial::order);
        return 1;
    }

    // High-momentum endpoint of Fig. 9 in Appendix C of arXiv:2102.06378,
    // averaged with the same rounded shells and real-FFT multiplicities used
    // by output.cpp for N=128 and L=1.4/m.
    const int bin_count = static_cast<int>(std::sqrt(3.0) * (N / 2)) + 1;
    std::vector<double> sums(bin_count, 0.0);
    std::vector<int> counts(bin_count, 0);
    for (int px_index = 0; px_index < N; ++px_index) {
        const int px = spatial::signed_mode(px_index);
        for (int py_index = 0; py_index < N; ++py_index) {
            const int py = spatial::signed_mode(py_index);
            for (int pz = 0; pz <= N / 2; ++pz) {
                const int multiplicity = (pz == 0 || pz == N / 2) ? 1 : 2;
                const int bin = spatial::shell_index(px, py, pz);
                // Match the production spectrum binning: rounding can place
                // the outermost cube corners just beyond the radial bins.
                if (bin < 0 || bin >= bin_count) continue;
                sums[bin] += multiplicity * std::sqrt(
                    spatial::effective_momentum_squared(px, py, pz, spacing));
                counts[bin] += multiplicity;
            }
        }
    }
    int last_bin = bin_count - 1;
    while (last_bin > 0 && counts[last_bin] == 0) --last_bin;
    const double endpoint = sums[last_bin] / counts[last_bin];
    if constexpr (N == 128) {
        constexpr double expected[] = {316.646339784, 365.611077758, 389.204648973};
        const int slot = spatial::order / 2 - 1;
        if (!close(endpoint, expected[slot], 2.e-8)) {
            std::fprintf(stderr, "order %d: Fig. 9 endpoint mismatch (%.9g)\n",
                         spatial::order, endpoint);
            return 1;
        }
    }

    std::printf("PASS order %d: operators and k_eff endpoint (%.6f)\n",
                spatial::order, endpoint);
    return 0;
}
