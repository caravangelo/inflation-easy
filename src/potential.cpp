// potential.cpp — Inflationary potential (analytic or tabulated)
//
// This file defines the inflationary potential and its derivatives.
// For analytic models (e.g. m²ϕ²), customize `analytic_potential` and its derivative.
// For numerical models, potential values are interpolated from tabulated data.

#include "main.h"

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
    int l = 0;
    while (field_numerical[l] >= field_value)
        l++;

    return (
        potential_numerical[l] +
        (field_value - field_numerical[l]) *
        (potential_numerical[l - 1] - potential_numerical[l]) /
        (field_numerical[l - 1] - field_numerical[l])
    ) / pw2(rescale_B);
#else
    return analytic_potential(field_value);
#endif
}

// -------------------------------------------------------------
// Return ∂V/∂ϕ at a grid point (interpolated or analytic)
// -------------------------------------------------------------
double potential_derivative(int i, int j, int k) {
#if numerical_potential
    int l = lstart[idx(i,j,k)];
    while (field_numerical[l] >= f[idx(i,j,k)])
        l++;

    return (
        potential_derivative_numerical[l] +
        (f[idx(i,j,k)] - field_numerical[l]) *
        (potential_derivative_numerical[l - 1] - potential_derivative_numerical[l]) /
        (field_numerical[l - 1] - field_numerical[l])
    ) / pw2(rescale_B);
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
        l = lstart[idx(i,j,k)];
        while (field_numerical[l] >= f[idx(i,j,k)])
            l++;

        if (field_numerical[l - 1] < f[idx(i,j,k)]) {
            printf("Interpolation Error\n");
            exit(1);
        }

        pot += (
            potential_numerical[l] +
            (f[idx(i,j,k)] - field_numerical[l]) *
            (potential_numerical[l - 1] - potential_numerical[l]) /
            (field_numerical[l - 1] - field_numerical[l])
        ) / pw2(rescale_B);

        lstart[idx(i,j,k)] = l - int_err;
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
    int l = lstart[idx(i,j,k)];
    while (field_numerical[l] >= f[idx(i,j,k)])
        l++;

    if (field_numerical[l - 1] < f[idx(i,j,k)]) {
        printf("Interpolation Error\n");
        exit(1);
    }

    double pot = (
        potential_numerical[l] +
        (f[idx(i,j,k)] - field_numerical[l]) *
        (potential_numerical[l - 1] - potential_numerical[l]) /
        (field_numerical[l - 1] - field_numerical[l])
    );

    double pot_deriv = (
        potential_derivative_numerical[l] +
        (f[idx(i,j,k)] - field_numerical[l]) *
        (potential_derivative_numerical[l - 1] - potential_derivative_numerical[l]) /
        (field_numerical[l - 1] - field_numerical[l])
    );

    lstart[idx(i,j,k)] = l - int_errN;
#else
    double pot = potential(f[idx(i,j,k)]);
    double pot_deriv = potential_derivative(i, j, k);
#endif

    return pot_deriv / pot;
}
