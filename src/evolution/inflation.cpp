// inflation.cpp - Inflationary equations, energy diagnostics, and leapfrog updates

#include "main.h"
#include "evolution_internal.h"
#include "linear_metric.h"

#if calculate_SIGW
// -------------------- Stress-energy tensor (inflaton) --------------------

namespace {

// Scalar-field gradient at one lattice site, evaluated in double precision.
struct ScalarGradient {
    double component[3]; // ∂_x φ, ∂_y φ, ∂_z φ
};

// Cache all scalar first derivatives needed by the six tensor-source components.
INFLATIONEASY_NOINLINE void compute_scalar_gradient(
    int i, int j, int k, ScalarGradient& gradient)
{
    gradient.component[0] = evolution::first_spatial_derivative<double>(0, i, j, k, f);
    gradient.component[1] = evolution::first_spatial_derivative<double>(1, i, j, k, f);
    gradient.component[2] = evolution::first_spatial_derivative<double>(2, i, j, k, f);
}

} // namespace

#endif

// -------------------- Energy Calculations --------------------

// Compute gradient energy density (averaged)
double gradient_energy() {
    DECLARE_INDICES
    double gradient = 0.0;
    const double norm = pw2(1.0 / (a * dx));
    LOOP gradient -= f[idx(i,j,k)] * spatial::laplacian(i, j, k, f);
    return 0.5 * gradient * norm / static_cast<double>(gridsize);
}

// Compute kinetic energy density (averaged)
double kin_energy() {
    DECLARE_INDICES
    double deriv_energy = 0.0;
    LOOP deriv_energy += pw2(fd[idx(i,j,k)]);
    deriv_energy /= static_cast<double>(gridsize);
    return 0.5 * std::pow(a, 2.0 * rescale_s - 2.0) * deriv_energy;
}

namespace {

// Apply one scalar leapfrog kick. Compile-time specialization keeps the metric
// arithmetic entirely outside the per-site loop when the run-time option is off.
template <bool IncludeLinearMetric>
void apply_scalar_leapfrog_kick(
    double step_size,
    double laplacian_coefficient,
    double friction,
    double potential_force_coefficient,
    double metric_coefficient = 0.0,
    double metric_field_mean = 0.0)
{
    DECLARE_INDICES
#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
    LOOP {
        const size_t id = idx(i,j,k);
        double rhs = laplacian_coefficient * spatial::laplacian(i, j, k, f)
                   - friction * fd[id]
                   - potential_force_coefficient * potential_derivative(i, j, k);
        if constexpr (IncludeLinearMetric) {
            rhs += metric_coefficient * (f[id] - metric_field_mean);
        }
        fd[id] += step_size * rhs;
    }
}

} // namespace

// -------------------- Main Field Evolution --------------------

// Apply a full leapfrog kick to the inflationary scalar, scale factor, and
// optional tensor fields while preserving their staggered-time convention.
void evolution::apply_inflation_leapfrog_kick(double step_size) {
    const double laplacian_coefficient = std::pow(a, -2.0 * rescale_s) / pw2(dx);
    const double one_plus_rescale_s = rescale_s + 1.0;
    const double scale_factor_power = -2.0 * rescale_s + 2.0;

    // Update second derivative of scale factor (ad2)
    ad2 = (-2.0 * ad - 2.0 * a / step_size / one_plus_rescale_s * (
    1.0 - std::sqrt(1.0 + 2.0 * step_size * one_plus_rescale_s * ad / a +
    pw2(step_size) * one_plus_rescale_s * std::pow(a, scale_factor_power) *
    (2.0 * gradient_energy() / 3.0 + potential_energy()))
    )) / step_size;

    ad += 0.5 * step_size * ad2;

    // Scalar field evolution
    const double friction = (2.0 + rescale_s) * ad / a;
    const double potential_force_coefficient = std::pow(a, 2.0 - 2.0 * rescale_s);
    if (linear_metric_perturbations) {
        const LinearMetricCorrection metric_correction =
            compute_linear_metric_correction(f, fd, a, ad, ad2);
        apply_scalar_leapfrog_kick<true>(
            step_size, laplacian_coefficient, friction, potential_force_coefficient,
            metric_correction.rhs_coefficient,
            metric_correction.field_mean);
    } else {
        apply_scalar_leapfrog_kick<false>(
            step_size, laplacian_coefficient, friction, potential_force_coefficient);
    }

#if calculate_SIGW
    // Gravitational wave tensor evolution (store as float, compute RHS in double)
    const double tensor_source_prefactor = 2.0 * std::pow(a, -2.0 * rescale_s);

#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
    for (int i=0; i<N; ++i)
    for (int j=0; j<N; ++j)
    for (int k=0; k<N; ++k) {
        const size_t id = idx(i,j,k);

        // Build gradients once and reuse them for all tensor components.
        ScalarGradient gradient;
        compute_scalar_gradient(i, j, k, gradient);

        // Form the six stored components of the symmetric source in double precision.
        const double gx = gradient.component[0];
        const double gy = gradient.component[1];
        const double gz = gradient.component[2];

        const double T_xx = gx * gx;
        const double T_yy = gy * gy;
        const double T_zz = gz * gz;
        const double T_xy = gx * gy;
        const double T_xz = gx * gz;
        const double T_yz = gy * gz;

        // Evaluate each RHS in double precision and cast once on storage.
        {
            const double lap_h = static_cast<double>(spatial::laplacian(i,j,k,hij[0]));
            const double rhs = step_size * ( laplacian_coefficient * lap_h
            - friction * static_cast<double>(hijd[0][id])
            + tensor_source_prefactor * T_xx );
            hijd[0][id] += static_cast<float>(rhs);
        }
        {
            const double lap_h = static_cast<double>(spatial::laplacian(i,j,k,hij[1]));
            const double rhs = step_size * ( laplacian_coefficient * lap_h
            - friction * static_cast<double>(hijd[1][id])
            + tensor_source_prefactor * T_yy );
            hijd[1][id] += static_cast<float>(rhs);
        }
        {
            const double lap_h = static_cast<double>(spatial::laplacian(i,j,k,hij[2]));
            const double rhs = step_size * ( laplacian_coefficient * lap_h
            - friction * static_cast<double>(hijd[2][id])
            + tensor_source_prefactor * T_zz );
            hijd[2][id] += static_cast<float>(rhs);
        }
        {
            const double lap_h = static_cast<double>(spatial::laplacian(i,j,k,hij[3]));
            const double rhs = step_size * ( laplacian_coefficient * lap_h
            - friction * static_cast<double>(hijd[3][id])
            + tensor_source_prefactor * T_xy );
            hijd[3][id] += static_cast<float>(rhs);
        }
        {
            const double lap_h = static_cast<double>(spatial::laplacian(i,j,k,hij[4]));
            const double rhs = step_size * ( laplacian_coefficient * lap_h
            - friction * static_cast<double>(hijd[4][id])
            + tensor_source_prefactor * T_xz );
            hijd[4][id] += static_cast<float>(rhs);
        }
        {
            const double lap_h = static_cast<double>(spatial::laplacian(i,j,k,hij[5]));
            const double rhs = step_size * ( laplacian_coefficient * lap_h
            - friction * static_cast<double>(hijd[5][id])
            + tensor_source_prefactor * T_yz );
            hijd[5][id] += static_cast<float>(rhs);
        }
    }
#endif

    ad += 0.5 * step_size * ad2;
}

namespace evolution {

// Evaluate the inflationary scalar, tensor, and background right-hand sides.
void compute_inflation_rhs(
    const std::vector<double>& f_state,
    const std::vector<double>& fd_state,
#if calculate_SIGW
    const std::vector<float> (&h_state)[6],
    const std::vector<float> (&hd_state)[6],
#endif
    std::vector<double>& dfdt,
    std::vector<double>& dfddt,
#if calculate_SIGW
    std::vector<float> (&dhdt)[6],
    std::vector<float> (&dhddt)[6],
#endif
    double a_state,
    double ad_state,
    double& dadt,
    double& daddt)
{
    DECLARE_INDICES

    const double laplacian_coefficient =
        std::pow(a_state, -2.0 * rescale_s) / pw2(dx);
    const double friction = (2.0 + rescale_s) * ad_state / a_state;
    const double potential_force_coefficient =
        std::pow(a_state, 2.0 - 2.0 * rescale_s);
#if calculate_SIGW
    const double tensor_source_prefactor =
        2.0 * std::pow(a_state, -2.0 * rescale_s);
#endif

    double gradient_sum = 0.0;
    double potential_sum = 0.0;

#if parallel_calculation
#pragma omp parallel for collapse(3) reduction(+:gradient_sum,potential_sum)
#endif
    LOOP {
        const size_t id = idx(i, j, k);
        const double field_here = f_state[id];
        const double lap_f = spatial::laplacian(i, j, k, f_state);
        double pot_here = 0.0;
        double pot_deriv_here = 0.0;
#if numerical_potential
        int next_hint = 1;
        evaluate_potential_from_value(field_here, lstart[id], int_err, &next_hint, pot_here, pot_deriv_here);
        lstart[id] = next_hint;
#else
        evaluate_potential_from_value(field_here, 1, 1, nullptr, pot_here, pot_deriv_here);
#endif

        dfdt[id] = fd_state[id];
        dfddt[id] = laplacian_coefficient * lap_f
                  - friction * fd_state[id]
                  - potential_force_coefficient * pot_deriv_here;

        gradient_sum -= field_here * lap_f;
        potential_sum += pot_here;

#if calculate_SIGW
        const double gx = first_spatial_derivative<double>(0, i, j, k, f_state);
        const double gy = first_spatial_derivative<double>(1, i, j, k, f_state);
        const double gz = first_spatial_derivative<double>(2, i, j, k, f_state);

        const double T_xx = gx * gx;
        const double T_yy = gy * gy;
        const double T_zz = gz * gz;
        const double T_xy = gx * gy;
        const double T_xz = gx * gz;
        const double T_yz = gy * gz;

        dhdt[0][id] = hd_state[0][id];
        dhdt[1][id] = hd_state[1][id];
        dhdt[2][id] = hd_state[2][id];
        dhdt[3][id] = hd_state[3][id];
        dhdt[4][id] = hd_state[4][id];
        dhdt[5][id] = hd_state[5][id];

        dhddt[0][id] = static_cast<float>(
            laplacian_coefficient
                * static_cast<double>(spatial::laplacian(i, j, k, h_state[0]))
            - friction * static_cast<double>(hd_state[0][id])
            + tensor_source_prefactor * T_xx);
        dhddt[1][id] = static_cast<float>(
            laplacian_coefficient
                * static_cast<double>(spatial::laplacian(i, j, k, h_state[1]))
            - friction * static_cast<double>(hd_state[1][id])
            + tensor_source_prefactor * T_yy);
        dhddt[2][id] = static_cast<float>(
            laplacian_coefficient
                * static_cast<double>(spatial::laplacian(i, j, k, h_state[2]))
            - friction * static_cast<double>(hd_state[2][id])
            + tensor_source_prefactor * T_zz);
        dhddt[3][id] = static_cast<float>(
            laplacian_coefficient
                * static_cast<double>(spatial::laplacian(i, j, k, h_state[3]))
            - friction * static_cast<double>(hd_state[3][id])
            + tensor_source_prefactor * T_xy);
        dhddt[4][id] = static_cast<float>(
            laplacian_coefficient
                * static_cast<double>(spatial::laplacian(i, j, k, h_state[4]))
            - friction * static_cast<double>(hd_state[4][id])
            + tensor_source_prefactor * T_xz);
        dhddt[5][id] = static_cast<float>(
            laplacian_coefficient
                * static_cast<double>(spatial::laplacian(i, j, k, h_state[5]))
            - friction * static_cast<double>(hd_state[5][id])
            + tensor_source_prefactor * T_yz);
#endif
    }

    const double gradient = 0.5 * gradient_sum * pw2(1.0 / (a_state * dx))
        / static_cast<double>(gridsize);
    const double potential_avg = potential_sum / static_cast<double>(gridsize);
    const double source = 2.0 * gradient / 3.0 + potential_avg;

    dadt = ad_state;
    daddt = std::pow(a_state, 3.0 - 2.0 * rescale_s) * source
        - (rescale_s + 1.0) * pw2(ad_state) / a_state;

    if (linear_metric_perturbations) {
        const LinearMetricCorrection metric_correction =
            compute_linear_metric_correction(
                f_state, fd_state, a_state, ad_state, daddt);
#if parallel_calculation
#pragma omp parallel for
#endif
        for (long long raw_id = 0; raw_id < static_cast<long long>(f_state.size()); ++raw_id) {
            const size_t id = static_cast<size_t>(raw_id);
            dfddt[id] += metric_correction.rhs_coefficient
                       * (f_state[id] - metric_correction.field_mean);
        }
    }
    return;
}

} // namespace evolution
