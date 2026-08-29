// post_inflation.cpp - Post-inflationary scalar and tensor evolution

#include "main.h"
#include "evolution_internal.h"

#if post_inflation
// -------------------- Post-inflation Evolution --------------------

namespace {

#if calculate_SIGW
// Coefficient multiplying the Newtonian-gauge source in the tensor
// normalization used by InflationEasy.
constexpr double POST_INFLATION_TENSOR_SOURCE_COEFFICIENT = -2.0;
#endif

// Newtonian-potential derivatives reused by all tensor-source components.
struct TensorSourceDerivatives {
    double potential_gradient[3];
    double velocity_gradient[3];
    double potential_hessian[6]; // [xx, yy, zz, xy, xz, yz]
    double newtonian_potential;
};

// Same-axis derivative for the active post-inflationary field state.
inline double second_derivative_along_axis(
    int axis, int i, int j, int k, size_t id, const std::vector<double>& field)
{
    if constexpr (spatial::order == 2) {
        const int ip = axis == 0 ? spatial::increment(i) : i;
        const int im = axis == 0 ? spatial::decrement(i) : i;
        const int jp = axis == 1 ? spatial::increment(j) : j;
        const int jm = axis == 1 ? spatial::decrement(j) : j;
        const int kp = axis == 2 ? spatial::increment(k) : k;
        const int km = axis == 2 ? spatial::decrement(k) : k;
        const double inv_dx2 = 1.0 / (dx * dx);
        if (axis == 0) {
            return (field[idx(ip, j, k)] - 2.0 * field[id]
                    + field[idx(im, j, k)]) * inv_dx2;
        }
        if (axis == 1) {
            return (field[idx(i, jp, k)] - 2.0 * field[id]
                    + field[idx(i, jm, k)]) * inv_dx2;
        }
        return (field[idx(i, j, kp)] - 2.0 * field[id]
                + field[idx(i, j, km)]) * inv_dx2;
    } else {
        return spatial::second_derivative(axis, i, j, k, id, field, dx);
    }
}

// Mixed second derivative for the active post-inflationary field state.
inline double mixed_second_derivative(
    int first_axis,
    int second_axis,
    int i,
    int j,
    int k,
    const std::vector<double>& field)
{
    if constexpr (spatial::order == 2) {
        const int ip = (first_axis == 0 || second_axis == 0)
            ? spatial::increment(i) : i;
        const int im = (first_axis == 0 || second_axis == 0)
            ? spatial::decrement(i) : i;
        const int jp = (first_axis == 1 || second_axis == 1)
            ? spatial::increment(j) : j;
        const int jm = (first_axis == 1 || second_axis == 1)
            ? spatial::decrement(j) : j;
        const int kp = (first_axis == 2 || second_axis == 2)
            ? spatial::increment(k) : k;
        const int km = (first_axis == 2 || second_axis == 2)
            ? spatial::decrement(k) : k;
        return (field[idx(ip, jp, kp)] - field[idx(ip, jm, km)]
              - field[idx(im, jp, kp)] + field[idx(im, jm, km)]) / (4.0 * dx * dx);
    } else {
        return spatial::mixed_derivative(
            first_axis, second_axis, i, j, k, field, dx);
    }
}

// Cache the gradient and Hessian entries reused by the tensor-source terms.
INFLATIONEASY_NOINLINE void compute_tensor_source_derivatives(
    int i, int j, int k, TensorSourceDerivatives& derivatives)
{
    // First derivatives.
    derivatives.potential_gradient[0] =
        evolution::first_spatial_derivative<double>(0, i, j, k, f);
    derivatives.potential_gradient[1] =
        evolution::first_spatial_derivative<double>(1, i, j, k, f);
    derivatives.potential_gradient[2] =
        evolution::first_spatial_derivative<double>(2, i, j, k, f);
    derivatives.velocity_gradient[0] =
        evolution::first_spatial_derivative<double>(0, i, j, k, fd);
    derivatives.velocity_gradient[1] =
        evolution::first_spatial_derivative<double>(1, i, j, k, fd);
    derivatives.velocity_gradient[2] =
        evolution::first_spatial_derivative<double>(2, i, j, k, fd);

    // Hessian entries in sym_idx packing.
    const size_t id = idx(i, j, k);
    derivatives.potential_hessian[sym_idx(0,0)] =
        second_derivative_along_axis(0, i, j, k, id, f); // xx
    derivatives.potential_hessian[sym_idx(1,1)] =
        second_derivative_along_axis(1, i, j, k, id, f); // yy
    derivatives.potential_hessian[sym_idx(2,2)] =
        second_derivative_along_axis(2, i, j, k, id, f); // zz
    derivatives.potential_hessian[sym_idx(0,1)] =
        mixed_second_derivative(0, 1, i, j, k, f); // xy
    derivatives.potential_hessian[sym_idx(0,2)] =
        mixed_second_derivative(0, 2, i, j, k, f); // xz
    derivatives.potential_hessian[sym_idx(1,2)] =
        mixed_second_derivative(1, 2, i, j, k, f); // yz

    derivatives.newtonian_potential = f[id];
}

// Apply one leapfrog kick to the post-inflationary scalar, background, and
// optional tensor fields.
void apply_post_inflation_leapfrog_kick(double step_size) {
    DECLARE_INDICES

    const double laplacian_coefficient = std::pow(a, -2.0 * rescale_s) / pw2(dx);

    // Update second derivative of scale factor (ad2)
    ad2 = - ( rescale_s - 0.5*(1.0 - 3.0*omega)) * pw2(ad) / a;

    ad += 0.5 * step_size * ad2;

    // Scalar field evolution
#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
    LOOP {
        fd[idx(i,j,k)] += step_size * (
        omega * laplacian_coefficient * spatial::laplacian(i, j, k, f)
        - (3.0 * (1.0 + omega) + rescale_s) * ad * fd[idx(i,j,k)] / a
        );
    }

#if calculate_SIGW
    // Precompute time-dependent factors once per call (hot loop optimization)
    const double code_hubble = ad / a;
    const double inverse_code_hubble = 1.0 / code_hubble;
    const double velocity_source_coefficient = 4.0 / (3.0 * (1.0 + omega));

    // Source prefactor in the tensor equation normalization used by this code.
    const double tensor_source_prefactor =
        POST_INFLATION_TENSOR_SOURCE_COEFFICIENT
        * std::pow(a, -2.0 * rescale_s);

#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
    for (int i=0; i<N; ++i)
    for (int j=0; j<N; ++j)
    for (int k=0; k<N; ++k) {
        const size_t id = idx(i,j,k);

        // Build derivatives once per site and reuse them for all components.
        TensorSourceDerivatives derivatives;
        compute_tensor_source_derivatives(i, j, k, derivatives);

        // Construct the unprojected real-space sources.
        const double newtonian_potential = derivatives.newtonian_potential;

        const double gx = derivatives.potential_gradient[0];
        const double gy = derivatives.potential_gradient[1];
        const double gz = derivatives.potential_gradient[2];

        // Gradient of U = Phi'/H_code + Phi in the tensor source.
        const double grad_u_x =
            derivatives.velocity_gradient[0] * inverse_code_hubble + gx;
        const double grad_u_y =
            derivatives.velocity_gradient[1] * inverse_code_hubble + gy;
        const double grad_u_z =
            derivatives.velocity_gradient[2] * inverse_code_hubble + gz;

        // Hessian packing is [xx, yy, zz, xy, xz, yz]
        const double hessian_xx = derivatives.potential_hessian[0];
        const double hessian_yy = derivatives.potential_hessian[1];
        const double hessian_zz = derivatives.potential_hessian[2];
        const double hessian_xy = derivatives.potential_hessian[3];
        const double hessian_xz = derivatives.potential_hessian[4];
        const double hessian_yz = derivatives.potential_hessian[5];

        const double source_xx = 4.0 * newtonian_potential * hessian_xx
            + 2.0 * gx * gx
            - velocity_source_coefficient * (grad_u_x * grad_u_x);
        const double source_yy = 4.0 * newtonian_potential * hessian_yy
            + 2.0 * gy * gy
            - velocity_source_coefficient * (grad_u_y * grad_u_y);
        const double source_zz = 4.0 * newtonian_potential * hessian_zz
            + 2.0 * gz * gz
            - velocity_source_coefficient * (grad_u_z * grad_u_z);
        const double source_xy = 4.0 * newtonian_potential * hessian_xy
            + 2.0 * gx * gy
            - velocity_source_coefficient * (grad_u_x * grad_u_y);
        const double source_xz = 4.0 * newtonian_potential * hessian_xz
            + 2.0 * gx * gz
            - velocity_source_coefficient * (grad_u_x * grad_u_z);
        const double source_yz = 4.0 * newtonian_potential * hessian_yz
            + 2.0 * gy * gz
            - velocity_source_coefficient * (grad_u_y * grad_u_z);

        // Update the six stored components of the symmetric tensor.
        {
            const double lap_h = static_cast<double>(spatial::laplacian(i, j, k, hij[0]));
            const double rhs = step_size * (
            laplacian_coefficient * lap_h
            - (2.0 + rescale_s) * ad * static_cast<double>(hijd[0][id]) / a
            + tensor_source_prefactor * source_xx
            );
            hijd[0][id] += static_cast<float>(rhs);
        }
        {
            const double lap_h = static_cast<double>(spatial::laplacian(i, j, k, hij[1]));
            const double rhs = step_size * (
            laplacian_coefficient * lap_h
            - (2.0 + rescale_s) * ad * static_cast<double>(hijd[1][id]) / a
            + tensor_source_prefactor * source_yy
            );
            hijd[1][id] += static_cast<float>(rhs);
        }
        {
            const double lap_h = static_cast<double>(spatial::laplacian(i, j, k, hij[2]));
            const double rhs = step_size * (
            laplacian_coefficient * lap_h
            - (2.0 + rescale_s) * ad * static_cast<double>(hijd[2][id]) / a
            + tensor_source_prefactor * source_zz
            );
            hijd[2][id] += static_cast<float>(rhs);
        }
        {
            const double lap_h = static_cast<double>(spatial::laplacian(i, j, k, hij[3]));
            const double rhs = step_size * (
            laplacian_coefficient * lap_h
            - (2.0 + rescale_s) * ad * static_cast<double>(hijd[3][id]) / a
            + tensor_source_prefactor * source_xy
            );
            hijd[3][id] += static_cast<float>(rhs);
        }
        {
            const double lap_h = static_cast<double>(spatial::laplacian(i, j, k, hij[4]));
            const double rhs = step_size * (
            laplacian_coefficient * lap_h
            - (2.0 + rescale_s) * ad * static_cast<double>(hijd[4][id]) / a
            + tensor_source_prefactor * source_xz
            );
            hijd[4][id] += static_cast<float>(rhs);
        }
        {
            const double lap_h = static_cast<double>(spatial::laplacian(i, j, k, hij[5]));
            const double rhs = step_size * (
            laplacian_coefficient * lap_h
            - (2.0 + rescale_s) * ad * static_cast<double>(hijd[5][id]) / a
            + tensor_source_prefactor * source_yz
            );
            hijd[5][id] += static_cast<float>(rhs);
        }
    }
#endif

    ad += 0.5 * step_size * ad2;
}

} // namespace

// -------------------- Main Post-Inflation Evolution Loop --------------------

// Run the complete post-inflationary stage. RK modes reuse the shared field
// integrator with the post-inflationary right-hand side selected explicitly.
void run_post_inflation_loop(FILE* output_log) {
    initialize_post_inflation();

    int numsteps = 0;

    evolution::FieldRKWorkspace rk_workspace;

    const auto report_step = [&](int step_index) {
        if (step_index % output_freq == 0) {
            save_post_inflation((step_index % output_infrequent_freq == 0) ? 1 : 0);
        }

        if (screen_updates && step_index % output_freq == 0) {
            printf("scale factor a = %f\n", a);
            printf("numsteps %i\n\n", step_index);
        }

        fprintf(output_log, "scale factor a = %f\n", a);
        fprintf(output_log, "numsteps %i\n\n", step_index);
        if (step_index % output_freq == 0) {
            fflush(output_log);
        }
    };

    switch (post_inflation_integrator) {
        case INTEGRATOR_LEAPFROG: {
            apply_leapfrog_drift(0.5 * dt_post_inflation);

            while (a <= af_post_inflation) {
                apply_post_inflation_leapfrog_kick(dt_post_inflation);
                apply_leapfrog_drift(dt_post_inflation);

                numsteps++;
                report_step(numsteps);
            }
            break;
        }
        case INTEGRATOR_RK4: {
            while (a <= af_post_inflation) {
                evolution::rk4_step_fields(
                    dt_post_inflation,
                    rk_workspace,
                    evolution::SimulationPhase::PostInflation);
                numsteps++;
                report_step(numsteps);
            }
            break;
        }
        case INTEGRATOR_RK45:
        default: {
            double h = evolution::clamp_rk45_step(
                dt_post_inflation, rk45_max_dt);

            while (a <= af_post_inflation) {
                const double hmax = evolution::rk45_hmax_from_base_step(
                    dt_post_inflation);
                h = evolution::clamp_rk45_step(h, hmax);

                const evolution::RK45StepResult step_result =
                    evolution::rk45_accept_step(
                    h,
                    hmax,
                    [&](double h_trial, double h_limit, double& h_next) {
                        return evolution::rk45_step_fields(
                            h_trial,
                            h_limit,
                            h_next,
                            rk_workspace,
                            evolution::SimulationPhase::PostInflation);
                    },
                    [&](int attempts, double failed_h) {
                        std::fprintf(stderr,
                            "RK45 failed to converge in post-inflation at a=%e, t=%e (attempts=%d, h=%e)\n",
                            a, t, attempts, failed_h);
                        std::exit(1);
                    }
                );
                h = step_result.next_step_size;

                numsteps++;
                report_step(numsteps);
            }
            break;
        }
    }

    printf("Saving final inflaton data\n");
    fflush(output_log);
    save_post_inflation(1);
    // Note: save_post_inflation() performs a temporary sync/desync for output.
    // Leave the integrator state unchanged here.
}

namespace evolution {

// Evaluate the post-inflationary scalar, tensor, and background right-hand sides.
void compute_post_inflation_rhs(
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
    const double scalar_friction = (3.0 * (1.0 + omega) + rescale_s) * ad_state / a_state;
#if calculate_SIGW
    const double tensor_friction = (2.0 + rescale_s) * ad_state / a_state;
    const double code_hubble = ad_state / a_state;
    const double inverse_code_hubble = 1.0 / code_hubble;
    const double velocity_source_coefficient = 4.0 / (3.0 * (1.0 + omega));
    const double tensor_source_prefactor =
        POST_INFLATION_TENSOR_SOURCE_COEFFICIENT
        * std::pow(a_state, -2.0 * rescale_s);
#endif

#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
    LOOP {
        const size_t id = idx(i, j, k);
        const double lap_f = spatial::laplacian(i, j, k, f_state);

        dfdt[id] = fd_state[id];
        dfddt[id] = omega * laplacian_coefficient * lap_f
                   - scalar_friction * fd_state[id];

#if calculate_SIGW
        const double gx = first_spatial_derivative<double>(0, i, j, k, f_state);
        const double gy = first_spatial_derivative<double>(1, i, j, k, f_state);
        const double gz = first_spatial_derivative<double>(2, i, j, k, f_state);
        const double gfdx = first_spatial_derivative<double>(0, i, j, k, fd_state);
        const double gfdy = first_spatial_derivative<double>(1, i, j, k, fd_state);
        const double gfdz = first_spatial_derivative<double>(2, i, j, k, fd_state);

        const double hessian_xx =
            second_derivative_along_axis(0, i, j, k, id, f_state);
        const double hessian_yy =
            second_derivative_along_axis(1, i, j, k, id, f_state);
        const double hessian_zz =
            second_derivative_along_axis(2, i, j, k, id, f_state);
        const double hessian_xy =
            mixed_second_derivative(0, 1, i, j, k, f_state);
        const double hessian_xz =
            mixed_second_derivative(0, 2, i, j, k, f_state);
        const double hessian_yz =
            mixed_second_derivative(1, 2, i, j, k, f_state);

        // Gradient of U = Phi'/H_code + Phi in the tensor source.
        const double grad_u_x = gfdx * inverse_code_hubble + gx;
        const double grad_u_y = gfdy * inverse_code_hubble + gy;
        const double grad_u_z = gfdz * inverse_code_hubble + gz;
        const double newtonian_potential = f_state[id];

        const double source_xx = 4.0 * newtonian_potential * hessian_xx
            + 2.0 * gx * gx
            - velocity_source_coefficient * (grad_u_x * grad_u_x);
        const double source_yy = 4.0 * newtonian_potential * hessian_yy
            + 2.0 * gy * gy
            - velocity_source_coefficient * (grad_u_y * grad_u_y);
        const double source_zz = 4.0 * newtonian_potential * hessian_zz
            + 2.0 * gz * gz
            - velocity_source_coefficient * (grad_u_z * grad_u_z);
        const double source_xy = 4.0 * newtonian_potential * hessian_xy
            + 2.0 * gx * gy
            - velocity_source_coefficient * (grad_u_x * grad_u_y);
        const double source_xz = 4.0 * newtonian_potential * hessian_xz
            + 2.0 * gx * gz
            - velocity_source_coefficient * (grad_u_x * grad_u_z);
        const double source_yz = 4.0 * newtonian_potential * hessian_yz
            + 2.0 * gy * gz
            - velocity_source_coefficient * (grad_u_y * grad_u_z);

        dhdt[0][id] = hd_state[0][id];
        dhdt[1][id] = hd_state[1][id];
        dhdt[2][id] = hd_state[2][id];
        dhdt[3][id] = hd_state[3][id];
        dhdt[4][id] = hd_state[4][id];
        dhdt[5][id] = hd_state[5][id];

        dhddt[0][id] = static_cast<float>(
            laplacian_coefficient
                * static_cast<double>(spatial::laplacian(i, j, k, h_state[0]))
            - tensor_friction * static_cast<double>(hd_state[0][id])
            + tensor_source_prefactor * source_xx);
        dhddt[1][id] = static_cast<float>(
            laplacian_coefficient
                * static_cast<double>(spatial::laplacian(i, j, k, h_state[1]))
            - tensor_friction * static_cast<double>(hd_state[1][id])
            + tensor_source_prefactor * source_yy);
        dhddt[2][id] = static_cast<float>(
            laplacian_coefficient
                * static_cast<double>(spatial::laplacian(i, j, k, h_state[2]))
            - tensor_friction * static_cast<double>(hd_state[2][id])
            + tensor_source_prefactor * source_zz);
        dhddt[3][id] = static_cast<float>(
            laplacian_coefficient
                * static_cast<double>(spatial::laplacian(i, j, k, h_state[3]))
            - tensor_friction * static_cast<double>(hd_state[3][id])
            + tensor_source_prefactor * source_xy);
        dhddt[4][id] = static_cast<float>(
            laplacian_coefficient
                * static_cast<double>(spatial::laplacian(i, j, k, h_state[4]))
            - tensor_friction * static_cast<double>(hd_state[4][id])
            + tensor_source_prefactor * source_xz);
        dhddt[5][id] = static_cast<float>(
            laplacian_coefficient
                * static_cast<double>(spatial::laplacian(i, j, k, h_state[5]))
            - tensor_friction * static_cast<double>(hd_state[5][id])
            + tensor_source_prefactor * source_yz);
#endif
    }

    dadt = ad_state;
    daddt = - (rescale_s - 0.5 * (1.0 - 3.0 * omega)) * pw2(ad_state) / a_state;
}

} // namespace evolution

#endif
