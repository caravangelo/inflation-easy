// evolution.cpp - Inflationary evolution driver
//
// The inflationary equations and leapfrog kernels are implemented in
// inflation.cpp. Shared RK4/RK45 machinery lives in integrators.cpp.

#include "main.h"
#include "evolution_internal.h"

namespace {

// Convert the configured inflationary step to the active code-time variable.
inline double inflation_code_time_step(double scale_factor_at_step_start) {
    return dt * std::pow(scale_factor_at_step_start, rescale_s - 1.0);
}

} // namespace

// -------------------- Inflationary Evolution Loop --------------------

void run_inflation_loop(FILE* output_log) {
    // Main inflation driver:
    // - selects integrator backend
    // - performs periodic outputs/logging
    // - guarantees final saved state is synchronized for output
    int numsteps = 0;

    evolution::FieldRKWorkspace rk_workspace;

    const auto report_step = [&](int step_index) {
        if (step_index % output_freq == 0 && a < af) {
            save((step_index % output_infrequent_freq == 0) ? 1 : 0);
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

    switch (integrator) {
        case INTEGRATOR_LEAPFROG: {
            apply_leapfrog_drift(0.5 * dt); // First leapfrog half-step

            while (a <= af) {
                const double code_time_step = inflation_code_time_step(astep);
                evolution::apply_inflation_leapfrog_kick(code_time_step);
                apply_leapfrog_drift(code_time_step);

                numsteps++;
                report_step(numsteps);
                astep = a;
            }
            break;
        }
        case INTEGRATOR_RK4: {
            while (a <= af) {
                const double code_time_step = inflation_code_time_step(astep);
                evolution::rk4_step_fields(
                    code_time_step,
                    rk_workspace,
                    evolution::SimulationPhase::Inflation);

                numsteps++;
                report_step(numsteps);
                astep = a;
            }
            break;
        }
        case INTEGRATOR_RK45:
        default: {
            double h = evolution::clamp_rk45_step(
                inflation_code_time_step(astep), rk45_max_dt);

            while (a <= af) {
                const double base_step = inflation_code_time_step(astep);
                const double hmax = evolution::rk45_hmax_from_base_step(base_step);
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
                            evolution::SimulationPhase::Inflation);
                    },
                    [&](int attempts, double failed_h) {
                        std::fprintf(stderr,
                            "RK45 failed to converge at a=%e, t=%e (attempts=%d, h=%e)\n",
                            a, t, attempts, failed_h);
                        std::exit(1);
                    }
                );
                h = step_result.next_step_size;

                numsteps++;
                report_step(numsteps);
                astep = a;
            }
            break;
        }
    }

    printf("Saving final inflaton data\n");
    fflush(output_log);
    save(1);
    if (inflation_uses_staggered_derivatives()) {
        // Temporarily align the staggered field and velocity states for output.
        apply_leapfrog_drift(-0.5 * inflation_code_time_step(astep));
    }
    save_last();
}
