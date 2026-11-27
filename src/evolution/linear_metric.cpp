// linear_metric.cpp - Optional linear metric correction to the inflaton equation
//
// Linear scalar metric perturbations in spatially flat gauge shift the
// effective mass of the inflaton fluctuations by (cosmic time t, M_Pl = 1)
//     Dm^2_metric = -(1/a^3) d/dt [ a^3 phidot_bar^2 / H ],
// with phi_bar identified with the lattice average <phi>. The lapse and shift
// are not evolved; the correction assumes the metric perturbations stay linear.

#include "linear_metric.h"

#include "main.h"

namespace {

// Homogeneous quantities reconstructed from the current nonlinear lattice state.
struct LinearMetricBackground {
    double field_mean = 0.0;
    double deriv_mean = 0.0;
    double deriv2_mean = 0.0;
};

// Background quantities <phi>, <phi'> and <phi''> from the lattice state.
// <phi''> follows from averaging the scalar equation over the box:
//     <phi''> = -(2+s)(a'/a) <phi'> - a^(2-2s) <dV~/dphi>
// (the Laplacian averages to zero; the force is averaged over sites rather
// than evaluated at <phi>).
LinearMetricBackground compute_background(
    const std::vector<double>& field,
    const std::vector<double>& field_derivative,
    double scale_factor,
    double scale_factor_derivative,
    double velocity_time_offset)
{
    LinearMetricBackground background;
    double field_sum = 0.0;
    double deriv_sum = 0.0;
    double potential_deriv_sum = 0.0;

#if parallel_calculation
#pragma omp parallel for reduction(+:field_sum,deriv_sum,potential_deriv_sum)
#endif
    for (long long raw_id = 0; raw_id < static_cast<long long>(field.size()); ++raw_id) {
        const size_t id = static_cast<size_t>(raw_id);
        const double field_here = field[id];
        double pot_unused = 0.0;
        double pot_deriv_here = 0.0;
#if numerical_potential
        evaluate_potential_from_value(
            field_here, lstart[id], int_err, nullptr, pot_unused, pot_deriv_here);
#else
        evaluate_potential_from_value(
            field_here, 1, 1, nullptr, pot_unused, pot_deriv_here);
#endif
        field_sum += field_here;
        deriv_sum += field_derivative[id];
        potential_deriv_sum += pot_deriv_here;
    }

    const double inv_size = 1.0 / static_cast<double>(field.size());
    background.field_mean = field_sum * inv_size;
    background.deriv_mean = deriv_sum * inv_size;

    const double code_hubble = scale_factor_derivative / scale_factor;
    const double potential_deriv_mean = potential_deriv_sum * inv_size;
    const double potential_force_coefficient =
        std::pow(scale_factor, 2.0 - 2.0 * rescale_s);
    const double friction = (2.0 + rescale_s) * code_hubble;
    const double force = -potential_force_coefficient * potential_deriv_mean;
    // The periodic Laplacian and the metric fluctuation force both average
    // to zero. Solve the mean kick for the velocity at the field time.
    background.deriv_mean = (background.deriv_mean + velocity_time_offset * force)
                          / (1.0 + velocity_time_offset * friction);
    background.deriv2_mean = force - friction * background.deriv_mean;
    return background;
}

} // namespace

// Returns <phi> and C = -a^(2-2s) Dm^2_metric / B^2. Writing Dm^2_metric in
// code units (phidot = B a^(s-1) phi', H = B a^(s-2) a') gives
//     C = (s+3) <phi'>^2 + 2 <phi'><phi''> a/a' - <phi'>^2 a a''/a'^2.
// RK stage data are synchronized. Leapfrog supplies the half-step offset.
// Degenerate input (empty lattice, a <= 0 or a' = 0) gives C = 0.
LinearMetricCorrection compute_linear_metric_correction(
    const std::vector<double>& field,
    const std::vector<double>& field_derivative,
    double scale_factor,
    double scale_factor_derivative,
    double scale_factor_second_derivative,
    double velocity_time_offset)
{
    LinearMetricCorrection correction;
    if (field.empty() || field.size() != field_derivative.size()
        || !(scale_factor > 0.0) || std::abs(scale_factor_derivative) < 1e-30) {
        return correction;
    }

    const LinearMetricBackground background = compute_background(
        field, field_derivative, scale_factor, scale_factor_derivative, velocity_time_offset);

    correction.field_mean = background.field_mean;
    const double code_hubble = scale_factor_derivative / scale_factor;
    correction.rhs_coefficient = (rescale_s + 3.0) * pw2(background.deriv_mean)
        + 2.0 * background.deriv_mean * background.deriv2_mean / code_hubble
        - pw2(background.deriv_mean) * scale_factor * scale_factor_second_derivative
            / pw2(scale_factor_derivative);
    return correction;
}
