// spatial_discretization.h - Finite-difference stencils and effective momenta
//
// All spatial derivatives use centered periodic stencils of order
// SPATIAL_STENCIL_ORDER = 2, 4 or 6 (parameters.h). Along one axis, with
// f_m = f(x + m dx):
//   order 2:  f'' = (f_1 - 2 f_0 + f_-1)/dx^2,
//             f'  = (f_1 - f_-1)/(2 dx)
//   order 4:  f'' = (-f_2 + 16 f_1 - 30 f_0 + 16 f_-1 - f_-2)/(12 dx^2),
//             f'  = (-f_2 + 8 f_1 - 8 f_-1 + f_-2)/(12 dx)
//   order 6:  f'' = (2 f_3 - 27 f_2 + 270 f_1 - 490 f_0 + ...)/(180 dx^2),
//             f'  = (f_3 - 9 f_2 + 45 f_1 - 45 f_-1 + ...)/(60 dx)
// The Laplacian is the sum of the three same-axis second derivatives.
//
// Effective momentum. On a Fourier mode with integer index n, the Laplacian
// stencil has eigenvalue -k_eff^2, with, per axis,
//   k_eff^2 = (4/dx^2) sin^2(pi n/N) [1 + sin^2(pi n/N)/3 (orders 4, 6)
//                                       + 8 sin^4(pi n/N)/45 (order 6)].
// This k_eff (the paper's modified dispersion relation) is used for the
// initial mode frequencies, the momentum labels of all spectra and the GW TT
// projector. The first-derivative stencils have a different Fourier symbol
// (e.g. sin(2 pi n/N)/dx at order 2), so products of first derivatives match
// the Laplacian only away from the lattice cutoff.
//
// Conventions: periodic boundaries; row-major storage with k fastest,
// index(i,j,k) = (i N + j) N + k. laplacian() returns the stencil sum without
// 1/dx^2; first_derivative(), second_derivative() and mixed_derivative()
// include the spacing.

#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

#include "parameters.h"

namespace spatial {

inline constexpr int order = SPATIAL_STENCIL_ORDER;
static_assert(order == 2 || order == 4 || order == 6,
              "SPATIAL_STENCIL_ORDER must be 2, 4, or 6");
static_assert(N >= order + 2,
              "The lattice is too small for the selected periodic stencil");

inline constexpr int radius = order / 2;

// Fast nearest-neighbor helpers used by the default second-order stencil.
// Keeping this path identical to the original implementation avoids adding
// indexing overhead when higher-order stencils are not selected.
inline int increment(int coordinate) {
    return coordinate == N - 1 ? 0 : coordinate + 1;
}

inline int decrement(int coordinate) {
    return coordinate == 0 ? N - 1 : coordinate - 1;
}

inline std::size_t index(int i, int j, int k) {
    return static_cast<std::size_t>(i) * N * N
         + static_cast<std::size_t>(j) * N
         + static_cast<std::size_t>(k);
}

template <int Offset>
inline int shifted(int coordinate) {
    static_assert(Offset >= -3 && Offset <= 3, "Unsupported stencil offset");
    int result = coordinate + Offset;
    if (result < 0) result += N;
    if (result >= N) result -= N;
    return result;
}

inline int shifted(int coordinate, int offset) {
    int result = coordinate + offset;
    if (result < 0) result += N;
    if (result >= N) result -= N;
    return result;
}

// Three-dimensional Laplacian times dx^2 (callers divide by dx^2 and apply
// the scale-factor factors of their equation).
template <typename T>
inline T laplacian(int i, int j, int k, const std::vector<T>& field) {
    if constexpr (order == 2) {
        // Preserve the compact default hot path, including direct interior indexing.
        if (i == 0 || j == 0 || k == 0 || i == N - 1 || j == N - 1 || k == N - 1) {
            return field[index(i, j, increment(k))] + field[index(i, j, decrement(k))]
                 + field[index(i, increment(j), k)] + field[index(i, decrement(j), k)]
                 + field[index(increment(i), j, k)] + field[index(decrement(i), j, k)]
                 - T(6) * field[index(i, j, k)];
        }
        return field[index(i, j, k + 1)] + field[index(i, j, k - 1)]
             + field[index(i, j + 1, k)] + field[index(i, j - 1, k)]
             + field[index(i + 1, j, k)] + field[index(i - 1, j, k)]
             - T(6) * field[index(i, j, k)];
    } else if constexpr (order == 4) {
        const T nearest =
            field[index(shifted<1>(i), j, k)] + field[index(shifted<-1>(i), j, k)]
          + field[index(i, shifted<1>(j), k)] + field[index(i, shifted<-1>(j), k)]
          + field[index(i, j, shifted<1>(k))] + field[index(i, j, shifted<-1>(k))];
        const T next =
            field[index(shifted<2>(i), j, k)] + field[index(shifted<-2>(i), j, k)]
          + field[index(i, shifted<2>(j), k)] + field[index(i, shifted<-2>(j), k)]
          + field[index(i, j, shifted<2>(k))] + field[index(i, j, shifted<-2>(k))];
        return T(4.0 / 3.0) * nearest - T(1.0 / 12.0) * next
             - T(15.0 / 2.0) * field[index(i, j, k)];
    } else {
        const T nearest =
            field[index(shifted<1>(i), j, k)] + field[index(shifted<-1>(i), j, k)]
          + field[index(i, shifted<1>(j), k)] + field[index(i, shifted<-1>(j), k)]
          + field[index(i, j, shifted<1>(k))] + field[index(i, j, shifted<-1>(k))];
        const T next =
            field[index(shifted<2>(i), j, k)] + field[index(shifted<-2>(i), j, k)]
          + field[index(i, shifted<2>(j), k)] + field[index(i, shifted<-2>(j), k)]
          + field[index(i, j, shifted<2>(k))] + field[index(i, j, shifted<-2>(k))];
        const T third =
            field[index(shifted<3>(i), j, k)] + field[index(shifted<-3>(i), j, k)]
          + field[index(i, shifted<3>(j), k)] + field[index(i, shifted<-3>(j), k)]
          + field[index(i, j, shifted<3>(k))] + field[index(i, j, shifted<-3>(k))];
        return T(3.0 / 2.0) * nearest - T(3.0 / 20.0) * next
             + T(1.0 / 90.0) * third - T(49.0 / 6.0) * field[index(i, j, k)];
    }
}

template <typename T>
inline T sample_axis(int dim, int i, int j, int k, int offset, const std::vector<T>& field) {
    if (dim == 0) i = shifted(i, offset);
    else if (dim == 1) j = shifted(j, offset);
    else k = shifted(k, offset);
    return field[index(i, j, k)];
}

// Centered first derivative with the same formal accuracy as the Laplacian.
template <typename T>
inline T first_derivative(int dim, int i, int j, int k, const std::vector<T>& field, double spacing) {
    if constexpr (order == 2) {
        const T half_over_spacing = T(0.5 / spacing);
        if (dim == 0) {
            if (i == 0 || i == N - 1) {
                return (field[index(increment(i), j, k)]
                      - field[index(decrement(i), j, k)]) * half_over_spacing;
            }
            return (field[index(i + 1, j, k)] - field[index(i - 1, j, k)]) * half_over_spacing;
        }
        if (dim == 1) {
            if (j == 0 || j == N - 1) {
                return (field[index(i, increment(j), k)]
                      - field[index(i, decrement(j), k)]) * half_over_spacing;
            }
            return (field[index(i, j + 1, k)] - field[index(i, j - 1, k)]) * half_over_spacing;
        }
        if (k == 0 || k == N - 1) {
            return (field[index(i, j, increment(k))]
                  - field[index(i, j, decrement(k))]) * half_over_spacing;
        }
        return (field[index(i, j, k + 1)] - field[index(i, j, k - 1)]) * half_over_spacing;
    } else if constexpr (order == 4) {
        return (-sample_axis(dim, i, j, k, 2, field)
                + T(8) * sample_axis(dim, i, j, k, 1, field)
                - T(8) * sample_axis(dim, i, j, k, -1, field)
                + sample_axis(dim, i, j, k, -2, field)) / T(12.0 * spacing);
    } else {
        return (sample_axis(dim, i, j, k, 3, field)
                - T(9) * sample_axis(dim, i, j, k, 2, field)
                + T(45) * sample_axis(dim, i, j, k, 1, field)
                - T(45) * sample_axis(dim, i, j, k, -1, field)
                + T(9) * sample_axis(dim, i, j, k, -2, field)
                - sample_axis(dim, i, j, k, -3, field)) / T(60.0 * spacing);
    }
}

// Same-axis second derivative, including the spacing normalization.
template <typename T>
inline T second_derivative(
    int dim, int i, int j, int k, std::size_t center_index,
    const std::vector<T>& field, double spacing)
{
    const T center = field[center_index];
    const T inv_spacing2 = T(1.0 / (spacing * spacing));
    if constexpr (order == 2) {
        if (dim == 0) {
            return (field[index(increment(i), j, k)] - T(2) * center
                  + field[index(decrement(i), j, k)]) * inv_spacing2;
        }
        if (dim == 1) {
            return (field[index(i, increment(j), k)] - T(2) * center
                  + field[index(i, decrement(j), k)]) * inv_spacing2;
        }
        return (field[index(i, j, increment(k))] - T(2) * center
              + field[index(i, j, decrement(k))]) * inv_spacing2;
    } else if constexpr (order == 4) {
        return (-sample_axis(dim, i, j, k, 2, field)
                + T(16) * sample_axis(dim, i, j, k, 1, field) - T(30) * center
                + T(16) * sample_axis(dim, i, j, k, -1, field)
                - sample_axis(dim, i, j, k, -2, field)) * (inv_spacing2 / T(12));
    } else {
        return (T(2) * sample_axis(dim, i, j, k, 3, field)
                - T(27) * sample_axis(dim, i, j, k, 2, field)
                + T(270) * sample_axis(dim, i, j, k, 1, field) - T(490) * center
                + T(270) * sample_axis(dim, i, j, k, -1, field)
                - T(27) * sample_axis(dim, i, j, k, -2, field)
                + T(2) * sample_axis(dim, i, j, k, -3, field)) * (inv_spacing2 / T(180));
    }
}

template <typename T>
inline T second_derivative(
    int dim, int i, int j, int k, const std::vector<T>& field, double spacing)
{
    return second_derivative(dim, i, j, k, index(i, j, k), field, spacing);
}

inline std::size_t offset_two_axes(
    int d1, int offset1, int d2, int offset2, int i, int j, int k)
{
    if (d1 == 0) i = shifted(i, offset1);
    else if (d1 == 1) j = shifted(j, offset1);
    else k = shifted(k, offset1);

    if (d2 == 0) i = shifted(i, offset2);
    else if (d2 == 1) j = shifted(j, offset2);
    else k = shifted(k, offset2);
    return index(i, j, k);
}

// Mixed derivative formed as the tensor product of the selected first-derivative stencil.
template <typename T>
inline T mixed_derivative(
    int d1, int d2, int i, int j, int k, const std::vector<T>& field, double spacing)
{
    static_assert(order == 2 || order == 4 || order == 6, "Unsupported stencil order");
    if constexpr (order == 2) {
        const int ip = (d1 == 0 || d2 == 0) ? increment(i) : i;
        const int im = (d1 == 0 || d2 == 0) ? decrement(i) : i;
        const int jp = (d1 == 1 || d2 == 1) ? increment(j) : j;
        const int jm = (d1 == 1 || d2 == 1) ? decrement(j) : j;
        const int kp = (d1 == 2 || d2 == 2) ? increment(k) : k;
        const int km = (d1 == 2 || d2 == 2) ? decrement(k) : k;
        return (field[index(ip, jp, kp)] - field[index(ip, jm, km)]
              - field[index(im, jp, kp)] + field[index(im, jm, km)])
             / T(4.0 * spacing * spacing);
    } else {
        constexpr int count = order;
        constexpr int shifts6[6] = {-3, -2, -1, 1, 2, 3};
        constexpr double weights4[4] = {1.0 / 12.0, -2.0 / 3.0, 2.0 / 3.0, -1.0 / 12.0};
        constexpr double weights6[6] = {-1.0 / 60.0, 3.0 / 20.0, -3.0 / 4.0,
                                          3.0 / 4.0, -3.0 / 20.0, 1.0 / 60.0};

        double result = 0.0;
        for (int a = 0; a < count; ++a) {
            const int shift_a = shifts6[a + (6 - count) / 2];
            const double weight_a = order == 4 ? weights4[a] : weights6[a];
            for (int b = 0; b < count; ++b) {
                const int shift_b = shifts6[b + (6 - count) / 2];
                const double weight_b = order == 4 ? weights4[b] : weights6[b];
                result += weight_a * weight_b
                        * static_cast<double>(field[offset_two_axes(d1, shift_a, d2, shift_b, i, j, k)]);
            }
        }
        return static_cast<T>(result / (spacing * spacing));
    }
}

// k_eff^2 of one axis for the signed mode index n (see file header).
inline double effective_component_squared(int signed_mode, double spacing) {
    const double sine = std::sin(3.141592653589793238462643383279502884
                               * static_cast<double>(signed_mode) / static_cast<double>(N));
    const double sine2 = sine * sine;
    double correction = 1.0;
    if constexpr (order >= 4) correction += sine2 / 3.0;
    if constexpr (order == 6) correction += 8.0 * sine2 * sine2 / 45.0;
    return 4.0 * sine2 * correction / (spacing * spacing);
}

// |k_eff|^2 of a three-dimensional mode: sum over the three axes.
inline double effective_momentum_squared(int px, int py, int pz, double spacing) {
    return effective_component_squared(px, spacing)
         + effective_component_squared(py, spacing)
         + effective_component_squared(pz, spacing);
}

// Signed k_eff component, used to build the TT projector of the GW spectra.
inline double effective_momentum_component(int signed_mode, double spacing) {
    const double magnitude = std::sqrt(effective_component_squared(signed_mode, spacing));
    return signed_mode < 0 ? -magnitude : magnitude;
}

// FFT index 0..N-1 -> signed mode index in (-N/2, N/2].
inline int signed_mode(int fft_index) {
    return fft_index <= N / 2 ? fft_index : fft_index - N;
}

// Isotropic spectrum bin of a mode: |n| rounded to the nearest integer.
inline int shell_index(int px, int py, int pz) {
    return static_cast<int>(std::lround(std::sqrt(
        static_cast<double>(px * px + py * py + pz * pz))));
}

// Largest per-axis k_eff^2 dx^2 (reached at the Nyquist mode, sin = 1).
inline constexpr double one_dimensional_max_eigenvalue() {
    if constexpr (order == 2) return 4.0;
    if constexpr (order == 4) return 16.0 / 3.0;
    return 272.0 / 45.0;
}

// Leapfrog stability limit for the massless lattice wave equation:
// dt * k_eff,max <= 2, with k_eff,max^2 = 3 * max eigenvalue / dx^2.
inline double courant_dt_over_dx_limit() {
    return 2.0 / std::sqrt(3.0 * one_dimensional_max_eigenvalue());
}

} // namespace spatial
