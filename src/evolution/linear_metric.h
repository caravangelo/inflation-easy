// linear_metric.h - Linear scalar metric correction for inflaton evolution

#pragma once

#include <vector>

// Quantities needed to add the linear metric mass shift to a lattice RHS.
// rhs_coefficient is already expressed in the active code-time convention.
struct LinearMetricCorrection {
    double field_mean = 0.0;
    double rhs_coefficient = 0.0;
};

// Compute the spatially averaged background and the associated metric-induced
// RHS coefficient. Invalid or degenerate inputs return a zero correction.
LinearMetricCorrection compute_linear_metric_correction(
    const std::vector<double>& field,
    const std::vector<double>& field_derivative,
    double scale_factor,
    double scale_factor_derivative,
    double scale_factor_second_derivative);
