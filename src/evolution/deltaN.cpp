// deltaN.cpp - Separate-universe evolution and curvature construction

#include "main.h"
#include "evolution_internal.h"

#include <limits>

#if perform_deltaN

using evolution::clamp_rk45_step;
using evolution::DP_A;
using evolution::DP_B;
using evolution::DP_E;
using evolution::rk45_accept_step;
using evolution::rk45_hmax_from_base_step;
using evolution::RK45StepResult;
using evolution::RK4_A;
using evolution::RK4_B;

// -------------------- DeltaN Evolution --------------------

namespace {
// Potential on the selected final hypersurface, used by the generic stopping criterion.
double reference_field_potential = 0.0;

bool deltaN_patch_is_active_with_potential(
    double field_value, [[maybe_unused]] double potential_value)
{
    const bool forward = dN > 0.0;
#if monotonic_potential
    return forward
        ? std::abs(field_value) > std::abs(phiref)
        : std::abs(field_value) < std::abs(phiref);
#elif antimonotonic_potential
    return forward
        ? std::abs(field_value) < std::abs(phiref)
        : std::abs(field_value) > std::abs(phiref);
#else
    return forward
        ? potential_value > reference_field_potential
        : potential_value < reference_field_potential;
#endif
}

bool deltaN_patch_is_beyond_reference(double field_value) {
    const bool forward = dN > 0.0;
#if monotonic_potential
    return forward
        ? std::abs(field_value) < std::abs(phiref)
        : std::abs(field_value) > std::abs(phiref);
#elif antimonotonic_potential
    return forward
        ? std::abs(field_value) > std::abs(phiref)
        : std::abs(field_value) < std::abs(phiref);
#else
    const double potential_value = potential(field_value);
    return forward
        ? potential_value < reference_field_potential
        : potential_value > reference_field_potential;
#endif
}
} // namespace

// Report whether a patch still has to reach the selected deltaN hypersurface.
bool deltaN_patch_is_active(double field_value) {
#if monotonic_potential || antimonotonic_potential
    return deltaN_patch_is_active_with_potential(field_value, 0.0);
#else
    return deltaN_patch_is_active_with_potential(field_value, potential(field_value));
#endif
}

namespace {

// Apply one leapfrog kick to the separate-universe velocity in e-fold time.
void apply_deltaN_leapfrog_kick(double step_size) {
    DECLARE_INDICES
#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
    LOOP {
        fd[idx(i,j,k)] += step_size
            * (-(3.0 - 0.5 * pw2(fd[idx(i,j,k)]))
            * (fd[idx(i,j,k)] + pot_ratio(i, j, k)));
    }
}

// Drift active separate-universe patches and accumulate their local expansion.
// Sites freeze after crossing the selected final hypersurface.
void apply_deltaN_leapfrog_drift(double step_size) {
    DECLARE_INDICES
    t += step_size;

#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
    LOOP {
        if (deltaN_patch_is_active(f[idx(i,j,k)])) {
            deltaN[idx(i,j,k)] += step_size;
            f[idx(i,j,k)] += step_size * fd[idx(i,j,k)];
        }
    }

    a += step_size * ad;
}

} // namespace

namespace {
// Reusable RK stage storage for the separate-universe system.
struct DeltaNRKWorkspace {
    std::vector<double> ftmp;
    std::vector<double> fdtmp;
    std::vector<double> dntmp;
    std::vector<double> kf[7];
    std::vector<double> kfd[7];
    std::vector<double> kdn[7];

    void ensure_size(size_t site_count) {
        if (ftmp.size() == site_count) return;
        ftmp.resize(site_count);
        fdtmp.resize(site_count);
        dntmp.resize(site_count);
        for (int s = 0; s < 7; ++s) {
            kf[s].resize(site_count);
            kfd[s].resize(site_count);
            kdn[s].resize(site_count);
        }
    }
};

// Evaluate the separate-universe RHS in e-fold-time coordinates. ddNdt is one
// before a patch reaches the final hypersurface and zero after it freezes.
void compute_deltaN_rhs(
    const std::vector<double>& f_state,
    const std::vector<double>& fd_state,
    std::vector<double>& dfdt,
    std::vector<double>& dfddt,
    std::vector<double>& ddNdt,
    double& dadt)
{
    DECLARE_INDICES

#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
    LOOP {
        const size_t id = idx(i, j, k);
        const double field_here = f_state[id];
        const double deriv_here = fd_state[id];
        double pot_here = 0.0;
        double pot_deriv_here = 0.0;
#if numerical_potential
        int next_hint = 1;
        evaluate_potential_from_value(field_here, lstart[id], int_errN, &next_hint, pot_here, pot_deriv_here);
        lstart[id] = next_hint;
#else
        evaluate_potential_from_value(field_here, 1, 1, nullptr, pot_here, pot_deriv_here);
#endif
        const double pot_ratio_here = pot_deriv_here / pot_here;
        const bool is_active = deltaN_patch_is_active_with_potential(field_here, pot_here);

        dfddt[id] = -(3.0 - 0.5 * pw2(deriv_here)) * (deriv_here + pot_ratio_here);

        if (is_active) {
            dfdt[id] = deriv_here;
            ddNdt[id] = 1.0;
        } else {
            dfdt[id] = 0.0;
            ddNdt[id] = 0.0;
        }
    }

    dadt = ad;
}

// Assemble the three per-site variables for one deltaN RK stage.
template <int StageCount>
void prepare_deltaN_stage_state(
    double h,
    int stage,
    const double (&A)[StageCount][StageCount],
    DeltaNRKWorkspace& workspace,
    size_t site_count)
{
    for (size_t id = 0; id < site_count; ++id) {
        double sum_f = 0.0;
        double sum_fd = 0.0;
        double sum_dn = 0.0;
        for (int s = 0; s < stage; ++s) {
            const double c = A[stage][s];
            if (c == 0.0) continue;
            sum_f += c * workspace.kf[s][id];
            sum_fd += c * workspace.kfd[s][id];
            sum_dn += c * workspace.kdn[s][id];
        }
        workspace.ftmp[id] = f[id] + h * sum_f;
        workspace.fdtmp[id] = fd[id] + h * sum_fd;
        workspace.dntmp[id] = deltaN[id] + h * sum_dn;
    }
}

// Build and evaluate one complete deltaN RK stage.
template <int StageCount>
void evaluate_deltaN_stage(
    double h,
    int stage,
    const double (&A)[StageCount][StageCount],
    DeltaNRKWorkspace& workspace,
    size_t site_count,
    double (&ka)[StageCount]) {
    prepare_deltaN_stage_state(h, stage, A, workspace, site_count);
    compute_deltaN_rhs(
        workspace.ftmp,
        workspace.fdtmp,
        workspace.kf[stage],
        workspace.kfd[stage],
        workspace.kdn[stage],
        ka[stage]);
}

// Commit a weighted combination of deltaN RK slopes.
template <int StageCount>
void apply_deltaN_weighted_update(
    double h,
    const double (&B)[StageCount],
    DeltaNRKWorkspace& workspace,
    size_t site_count,
    const double (&ka)[StageCount]) {
    for (size_t id = 0; id < site_count; ++id) {
        double sum_f = 0.0;
        double sum_fd = 0.0;
        double sum_dn = 0.0;
        for (int s = 0; s < StageCount; ++s) {
            const double bs = B[s];
            if (bs == 0.0) continue;
            sum_f += bs * workspace.kf[s][id];
            sum_fd += bs * workspace.kfd[s][id];
            sum_dn += bs * workspace.kdn[s][id];
        }
        f[id] += h * sum_f;
        fd[id] += h * sum_fd;
        deltaN[id] += h * sum_dn;
    }
    for (int s = 0; s < StageCount; ++s) {
        const double bs = B[s];
        if (bs == 0.0) continue;
        a += h * bs * ka[s];
    }
    t += h;
}

// Advance all separate-universe patches by one fixed classical RK4 step.
void rk4_step_deltaN(double h, DeltaNRKWorkspace& workspace) {
    const size_t site_count = f.size();
    workspace.ensure_size(site_count);

    double ka[4] = {0.0, 0.0, 0.0, 0.0};
    compute_deltaN_rhs(
        f, fd, workspace.kf[0], workspace.kfd[0], workspace.kdn[0], ka[0]);

    for (int stage = 1; stage < 4; ++stage) {
        evaluate_deltaN_stage(h, stage, RK4_A, workspace, site_count, ka);
    }

    apply_deltaN_weighted_update(h, RK4_B, workspace, site_count, ka);
}

// Attempt one adaptive Dormand-Prince step with the shared RK45 error policy.
bool rk45_step_deltaN(
    double h, double hmax, double& h_next, DeltaNRKWorkspace& workspace)
{
    const size_t site_count = f.size();
    workspace.ensure_size(site_count);

    double ka[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    compute_deltaN_rhs(
        f, fd, workspace.kf[0], workspace.kfd[0], workspace.kdn[0], ka[0]);

    for (int stage = 1; stage < 7; ++stage) {
        evaluate_deltaN_stage(h, stage, DP_A, workspace, site_count, ka);
    }

    double a5 = a;
    double erra = 0.0;
    for (int s = 0; s < 7; ++s) {
        const double bs = DP_B[s];
        const double es = DP_E[s];
        if (bs != 0.0) a5 += h * bs * ka[s];
        if (es != 0.0) erra += h * es * ka[s];
    }

    double err_acc = 0.0;
    std::size_t nvars = 3 * site_count + 1;

    for (size_t id = 0; id < site_count; ++id) {
        double y5f = f[id];
        double y5fd = fd[id];
        double y5dn = deltaN[id];
        double errf = 0.0;
        double errfd = 0.0;
        double errdn = 0.0;
        for (int s = 0; s < 7; ++s) {
            const double bs = DP_B[s];
            const double es = DP_E[s];
            if (bs != 0.0) {
                y5f += h * bs * workspace.kf[s][id];
                y5fd += h * bs * workspace.kfd[s][id];
                y5dn += h * bs * workspace.kdn[s][id];
            }
            if (es != 0.0) {
                errf += h * es * workspace.kf[s][id];
                errfd += h * es * workspace.kfd[s][id];
                errdn += h * es * workspace.kdn[s][id];
            }
        }

        workspace.ftmp[id] = y5f;
        workspace.fdtmp[id] = y5fd;
        workspace.dntmp[id] = y5dn;

        const double sf = rk45_abs_tol + rk45_rel_tol * std::max(std::abs(f[id]), std::abs(y5f));
        const double sfd = rk45_abs_tol + rk45_rel_tol * std::max(std::abs(fd[id]), std::abs(y5fd));
        const double sdn = rk45_abs_tol + rk45_rel_tol * std::max(std::abs(deltaN[id]), std::abs(y5dn));
        err_acc += pw2(errf / sf) + pw2(errfd / sfd) + pw2(errdn / sdn);
    }

    const double sa = rk45_abs_tol + rk45_rel_tol * std::max(std::abs(a), std::abs(a5));
    err_acc += pw2(erra / sa);

    const double err_norm = std::sqrt(err_acc / static_cast<double>(nvars));
    const double safe_err = std::max(err_norm, 1e-16);
    double factor = rk45_safety * std::pow(safe_err, -0.2);
    factor = std::clamp(factor, 0.1, 5.0);

    if (err_norm <= 1.0) {
        f.swap(workspace.ftmp);
        fd.swap(workspace.fdtmp);
        deltaN.swap(workspace.dntmp);
        a = a5;
        t += h;
        h_next = clamp_rk45_step(h * factor, hmax);
        return true;
    }

    h_next = clamp_rk45_step(h * std::max(0.1, std::min(0.5, factor)), hmax);
    return false;
}
} // namespace

namespace {

// Select the automatically determined field value defining the final deltaN slice.
// The serial scan is intentional because tie-breaking must remain deterministic.
double select_automatic_reference_field() {
    DECLARE_INDICES
    double fref = f[idx(0,0,0)];

    // Keep this reduction deterministic: fref selection depends on ordering.
    LOOP {
#if monotonic_potential
        if (std::abs(f[idx(i,j,k)]) < std::abs(fref))
#elif antimonotonic_potential
        if (std::abs(f[idx(i,j,k)]) > std::abs(fref))
#else
        if (potential(f[idx(i,j,k)]) < potential(fref))
#endif
        fref = f[idx(i,j,k)];
    }

    return fref;
}

} // namespace

// Run the complete separate-universe stage with the selected integrator.
void run_deltaN_loop(FILE* output_log) {
    printf("Starting deltaN calculation\n");
    fprintf(output_log, "Starting deltaN calculation\n");

    int numsteps = 0;
    Ne = 0.0;

    initializeN();
    phiref = use_phiref_manual
        ? phiref_manual_value
        : select_automatic_reference_field();
    reference_field_potential = potential(phiref);

    for (double field_value : f) {
        if (deltaN_patch_is_beyond_reference(field_value)) {
            std::fprintf(stderr,
                "The selected phiref lies in the opposite direction from dN for part of the lattice.\n");
            std::exit(1);
        }
    }

    const double budget_tolerance =
        32.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, Nend);
    const auto fixed_step_fits_budget = [&]() {
        return std::abs(dN) <= Nend - std::abs(Ne) + budget_tolerance;
    };

    const bool deltaN_state_is_staggered =
        deltaN_uses_staggered_derivatives() && fixed_step_fits_budget();
    if (deltaN_state_is_staggered) {
        apply_deltaN_leapfrog_drift(0.5 * dN);
    }

    DeltaNRKWorkspace rk_workspace;

    const auto report_step = [&](int step_index) {
        if (screen_updates && step_index % output_freq == 0) {
            printf("N = %f\n\n", Ne);
        }

        fprintf(output_log, "N = %f\n\n", Ne);
        if (step_index % output_freq == 0) {
            fflush(output_log);
        }
    };

    switch (deltaN_integrator) {
        case INTEGRATOR_LEAPFROG: {
            while (fixed_step_fits_budget()) {
                apply_deltaN_leapfrog_kick(dN);
                apply_deltaN_leapfrog_drift(dN);
                Ne += dN;
                numsteps++;
                report_step(numsteps);
            }
            break;
        }
        case INTEGRATOR_RK4: {
            while (fixed_step_fits_budget()) {
                rk4_step_deltaN(dN, rk_workspace);
                Ne += dN;
                numsteps++;
                report_step(numsteps);
            }
            break;
        }
        case INTEGRATOR_RK45:
        default: {
            double h = clamp_rk45_step(dN, rk45_max_dt);

            while (true) {
                const double remaining_budget = Nend - std::abs(Ne);
                if (remaining_budget <= budget_tolerance
                    || remaining_budget < rk45_min_dt) break;

                const double hmax = std::min(
                    rk45_hmax_from_base_step(dN), remaining_budget);
                h = std::copysign(std::min(std::abs(h), hmax), dN);
                h = clamp_rk45_step(h, hmax);

                const RK45StepResult step_result = rk45_accept_step(
                    h,
                    hmax,
                    [&](double h_trial, double h_limit, double& h_next) {
                        return rk45_step_deltaN(
                            h_trial, h_limit, h_next, rk_workspace);
                    },
                    [&](int attempts, double failed_h) {
                        std::fprintf(stderr,
                            "RK45 failed to converge in deltaN loop at N=%e, t=%e (attempts=%d, h=%e)\n",
                            Ne, t, attempts, failed_h);
                        std::exit(1);
                    }
                );

                Ne += step_result.accepted_step_size;
                h = step_result.next_step_size;
                numsteps++;
                report_step(numsteps);
            }
            break;
        }
    }

    if (deltaN_state_is_staggered) {
        apply_deltaN_leapfrog_drift(-0.5 * dN);
    }
    saveN(output_log);
    fflush(output_log);
}
#endif
