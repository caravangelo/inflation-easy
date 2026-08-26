// linear_metric.cpp - Linear scalar metric correction for inflaton evolution

#include "linear_metric.h"

#include "main.h"

namespace {

// Homogeneous quantities reconstructed from the current nonlinear lattice state.
struct LinearMetricBackground {
    double field_mean = 0.0;
    double deriv_mean = 0.0;
    double deriv2_mean = 0.0;
};

// Construct the homogeneous quantities required by the metric mass shift from
// spatial averages of the nonlinear lattice state. The mean acceleration is
// evaluated with the homogeneous scalar equation in code-time variables.
LinearMetricBackground compute_background(
    const std::vector<double>& field,
    const std::vector<double>& field_derivative,
    double scale_factor,
    double scale_factor_derivative)
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

    const double hubble_code = scale_factor_derivative / scale_factor;
    const double potential_deriv_mean = potential_deriv_sum * inv_size;
    const double potnorm = std::pow(scale_factor, 2.0 - 2.0 * rescale_s);
    background.deriv2_mean = -(2.0 + rescale_s) * hubble_code * background.deriv_mean
                           - potnorm * potential_deriv_mean;
    return background;
}

} // namespace

// Leading scalar metric correction in spatially flat gauge. In physical
// variables, Delta m^2 = -(1/a^3) d/dt [a^3 phidot_bar^2 / H], with M_Pl = 1.
// The returned coefficient carries the sign and scale-factor conversion needed
// for direct addition to the code-time field-acceleration RHS.
LinearMetricCorrection compute_linear_metric_correction(
    const std::vector<double>& field,
    const std::vector<double>& field_derivative,
    double scale_factor,
    double scale_factor_derivative,
    double scale_factor_second_derivative)
{
    LinearMetricCorrection correction;
    if (field.empty() || field.size() != field_derivative.size()
        || !(scale_factor > 0.0) || std::abs(scale_factor_derivative) < 1e-30) {
        return correction;
    }

    const LinearMetricBackground background = compute_background(
        field, field_derivative, scale_factor, scale_factor_derivative);

    correction.field_mean = background.field_mean;
    const double hubble_code = scale_factor_derivative / scale_factor;
    correction.rhs_coefficient = (rescale_s + 3.0) * pw2(background.deriv_mean)
        + 2.0 * background.deriv_mean * background.deriv2_mean / hubble_code
        - pw2(background.deriv_mean) * scale_factor * scale_factor_second_derivative
            / pw2(scale_factor_derivative);
    return correction;
}
