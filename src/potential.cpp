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

float analytic_potential(float field_value) {
    return 0.5 * m2 * pw2(field_value) / pw2(rescale_B);
}

float analytic_potential_derivative(float field_value) {
    return field_value * m2 / pw2(rescale_B);
}
#endif

// Nothing to customize in the following functions

// -------------------------------------------------------------
// Return the potential V(ϕ) at a given field value (interpolated or analytic)
// -------------------------------------------------------------
float potential(float field_value) {
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
float potential_derivative(int i, int j, int k) {
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
float potential_energy() {
    DECLARE_INDICES
    float pot = 0.0;

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

    pot /= gridsize;
    return pot;
}

// -------------------------------------------------------------
// Compute ∂V/∂ϕ divided by V at a grid point
// Used in δN evolution
// -------------------------------------------------------------
float pot_ratio(int i, int j, int k) {
#if numerical_potential
    int l = lstart[idx(i,j,k)];
    while (field_numerical[l] >= f[idx(i,j,k)])
        l++;

    if (field_numerical[l - 1] < f[idx(i,j,k)]) {
        printf("Interpolation Error\n");
        exit(1);
    }

    float pot = (
        potential_numerical[l] +
        (f[idx(i,j,k)] - field_numerical[l]) *
        (potential_numerical[l - 1] - potential_numerical[l]) /
        (field_numerical[l - 1] - field_numerical[l])
    );

    float pot_deriv = (
        potential_derivative_numerical[l] +
        (f[idx(i,j,k)] - field_numerical[l]) *
        (potential_derivative_numerical[l - 1] - potential_derivative_numerical[l]) /
        (field_numerical[l - 1] - field_numerical[l])
    );

    lstart[idx(i,j,k)] = l - int_errN;
#else
    float pot = potential(f[idx(i,j,k)]);
    float pot_deriv = potential_derivative(i, j, k);
#endif

    return pot_deriv / pot;
}
