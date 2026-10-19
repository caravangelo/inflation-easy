// linear_metric.h - Optional linear metric correction to the inflaton equation
//
// The paper adds a^2 Dm^2_metric (phi - <phi>) to the left-hand side of the
// Klein-Gordon equation. In code units this is the right-hand-side term
// C (phi - <phi>) of the phi'' equation, with C = -a^(2-2s) Dm^2_metric / B^2.

#pragma once

#include <vector>

// field_mean = <phi>; rhs_coefficient = C. Callers add C*(phi - <phi>) to phi''.
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
