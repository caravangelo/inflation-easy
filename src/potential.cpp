// potential.cpp — Inflationary potential (analytic or tabulated)
//
// This file defines the inflationary potential and its derivatives.
// For analytic models (e.g. m²ϕ²), customize `analytic_potential` and its derivative.
// For numerical models, potential values are interpolated from tabulated data.

#include "main.h"

#if numerical_potential
namespace {
int find_bracketing_index(double field_value, int start_index) {
    const int n = static_cast<int>(field_numerical.size());
    if (n < 2) {
        std::fprintf(stderr, "Numerical potential table must contain at least two points.\n");
        std::exit(1);
    }

    // Current implementation assumes field_numerical is sorted in descending order.
    if (field_value > field_numerical.front() || field_value < field_numerical.back()) {
        std::fprintf(
            stderr,
            "Field value %.17g is outside the numerical potential range [%.17g, %.17g].\n",
            field_value,
            field_numerical.back(),
            field_numerical.front());
        std::exit(1);
    }

    int l = std::clamp(start_index, 1, n - 1);
    while (l < n && field_numerical[l] >= field_value) ++l;
    while (l > 1 && field_numerical[l - 1] < field_value) --l;

    if (l <= 0 || l >= n) {
        std::fprintf(stderr, "Interpolation index out of bounds for field value %.17g.\n", field_value);
        std::exit(1);
    }

    return l;
}

double interpolate_from_table(
    double field_value,
    int l,
    const std::vector<double>& table_values)
{
    const double x0 = field_numerical[l - 1];
    const double x1 = field_numerical[l];
    const double y0 = table_values[l - 1];
    const double y1 = table_values[l];

    if (x0 == x1) {
        std::fprintf(stderr, "Interpolation failed due to repeated field grid points.\n");
        std::exit(1);
    }

    return y1 + (field_value - x1) * (y0 - y1) / (x0 - x1);
}
} // namespace
#endif

#if !numerical_potential
// -------------------------------------------------------------
// Analytic potential (default: quadratic V(ϕ) = 1/2 m²ϕ²)
// -------------------------------------------------------------

double analytic_potential(double field_value) {
    return V0 * (1.0 - (1.0 - ns) / 4.0 * pw2(field_value)) / pw2(rescale_B);
}

double analytic_potential_derivative(double field_value) {
    return -V0 * (1.0 - ns) / 2.0 * field_value / pw2(rescale_B);
}
#endif

// Nothing to customize in the following functions

// -------------------------------------------------------------
// Return the potential V(ϕ) at a given field value (interpolated or analytic)
// -------------------------------------------------------------
double potential(double field_value) {
#if numerical_potential
    const int l = find_bracketing_index(field_value, 1);
    return interpolate_from_table(field_value, l, potential_numerical) / pw2(rescale_B);
#else
    return analytic_potential(field_value);
#endif
}

// -------------------------------------------------------------
// Return ∂V/∂ϕ at a grid point (interpolated or analytic)
// -------------------------------------------------------------
double potential_derivative(int i, int j, int k) {
#if numerical_potential
    const size_t id = idx(i, j, k);
    const int l = find_bracketing_index(f[id], lstart[id]);
    lstart[id] = std::max(1, l - int_err);
    return interpolate_from_table(f[id], l, potential_derivative_numerical) / pw2(rescale_B);
#else
    return analytic_potential_derivative(f[idx(i,j,k)]);
#endif
}

// -------------------------------------------------------------
// Compute the total potential energy on the grid
// -------------------------------------------------------------
double potential_energy() {
    DECLARE_INDICES
    double pot = 0.0;

#if numerical_potential
    int l;
    LOOP {
        const size_t id = idx(i, j, k);
        l = find_bracketing_index(f[id], lstart[id]);
        pot += interpolate_from_table(f[id], l, potential_numerical) / pw2(rescale_B);
        lstart[id] = std::max(1, l - int_err);
    }
#else
    LOOP {
        pot += potential(f[idx(i,j,k)]);
    }
#endif

    pot /= static_cast<double>(gridsize);
    return pot;
}

// -------------------------------------------------------------
// Compute ∂V/∂ϕ divided by V at a grid point
// Used in δN evolution
// -------------------------------------------------------------
double pot_ratio(int i, int j, int k) {
#if numerical_potential
    const size_t id = idx(i, j, k);
    const int l = find_bracketing_index(f[id], lstart[id]);

    double pot = interpolate_from_table(f[id], l, potential_numerical);
    double pot_deriv = interpolate_from_table(f[id], l, potential_derivative_numerical);

    lstart[id] = std::max(1, l - int_errN);
#else
    double pot = potential(f[idx(i,j,k)]);
    double pot_deriv = potential_derivative(i, j, k);
#endif

    return pot_deriv / pot;
}
