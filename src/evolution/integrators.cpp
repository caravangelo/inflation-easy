// integrators.cpp - Shared leapfrog drift and RK4/RK45 stepping machinery

#include "main.h"
#include "evolution_internal.h"

// Apply the leapfrog drift to all active fields and advance code time.
// Inflationary and post-inflationary leapfrog loops share this operation.
void apply_leapfrog_drift(double step_size) {
    DECLARE_INDICES
    t += step_size;

#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
    LOOP {
        f[idx(i,j,k)] += step_size * fd[idx(i,j,k)];
    }

#if calculate_SIGW
#if parallel_calculation
#pragma omp parallel for collapse(4)
#endif
    for (int comp = 0; comp < 6; ++comp)
    LOOP {
        hij[comp][idx(i,j,k)] += static_cast<float>(
            step_size * hijd[comp][idx(i,j,k)]);
    }
#endif

    a += step_size * ad;
}

namespace evolution {
// Reused stage/state buffers for RK methods; ensure_size() allocates only when
// the lattice size changes, avoiding allocations inside time-stepping loops.
void FieldRKWorkspace::ensure_size(size_t site_count) {
    if (ftmp.size() == site_count) return;
    ftmp.resize(site_count);
    fdtmp.resize(site_count);
    for (int s = 0; s < 7; ++s) {
        kf[s].resize(site_count);
        kfd[s].resize(site_count);
    }
#if calculate_SIGW
    for (int c = 0; c < 6; ++c) {
        htmp[c].resize(site_count);
        hdtmp[c].resize(site_count);
    }
    for (int s = 0; s < 7; ++s) {
        for (int c = 0; c < 6; ++c) {
            kh[s][c].resize(site_count);
            khd[s][c].resize(site_count);
        }
    }
#endif
}

// Restrict a proposed adaptive step to the configured and stage-specific bounds.
double clamp_rk45_step(double h, double hmax) {
    const double hmin = std::max(1e-16, rk45_min_dt);
    const double hhi = std::max(hmin, hmax);
    const double direction = std::signbit(h) ? -1.0 : 1.0;
    return direction * std::clamp(std::abs(h), hmin, hhi);
}

// Allow adaptive growth while keeping an accepted step near the requested base step.
double rk45_hmax_from_base_step(double base_step) {
    return std::min(rk45_max_dt, std::max(rk45_min_dt, 2.0 * std::abs(base_step)));
}

namespace {

// Dispatch one RK stage to the equations for the selected simulation phase.
void compute_field_rhs(
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
    double& daddt,
    SimulationPhase phase)
{
    if (phase == SimulationPhase::Inflation) {
#if calculate_SIGW
        compute_inflation_rhs(
            f_state, fd_state, h_state, hd_state, dfdt, dfddt, dhdt, dhddt,
            a_state, ad_state, dadt, daddt);
#else
        compute_inflation_rhs(
            f_state, fd_state, dfdt, dfddt, a_state, ad_state, dadt, daddt);
#endif
        return;
    }

#if post_inflation
#if calculate_SIGW
    compute_post_inflation_rhs(
        f_state, fd_state, h_state, hd_state, dfdt, dfddt, dhdt, dhddt,
        a_state, ad_state, dadt, daddt);
#else
    compute_post_inflation_rhs(
        f_state, fd_state, dfdt, dfddt, a_state, ad_state, dadt, daddt);
#endif
#endif
}

} // namespace

namespace {

// Assemble the lattice fields for one explicit RK stage from prior stage slopes.
template <int StageCount>
void prepare_field_stage_state(
    double h,
    int stage,
    const double (&A)[StageCount][StageCount],
    FieldRKWorkspace& workspace,
    size_t site_count)
{
    for (size_t id = 0; id < site_count; ++id) {
        double sum_f = 0.0;
        double sum_fd = 0.0;
        for (int s = 0; s < stage; ++s) {
            const double c = A[stage][s];
            if (c == 0.0) continue;
            sum_f += c * workspace.kf[s][id];
            sum_fd += c * workspace.kfd[s][id];
        }
        workspace.ftmp[id] = f[id] + h * sum_f;
        workspace.fdtmp[id] = fd[id] + h * sum_fd;
    }
#if calculate_SIGW
    for (int c = 0; c < 6; ++c) {
        for (size_t id = 0; id < site_count; ++id) {
            double sum_h = 0.0;
            double sum_hd = 0.0;
            for (int s = 0; s < stage; ++s) {
                const double cs = A[stage][s];
                if (cs == 0.0) continue;
                sum_h += cs * static_cast<double>(workspace.kh[s][c][id]);
                sum_hd += cs * static_cast<double>(workspace.khd[s][c][id]);
            }
            workspace.htmp[c][id] = static_cast<float>(
                static_cast<double>(hij[c][id]) + h * sum_h);
            workspace.hdtmp[c][id] = static_cast<float>(
                static_cast<double>(hijd[c][id]) + h * sum_hd);
        }
    }
#endif
}

// Assemble the scale factor and its derivative for one explicit RK stage.
template <int StageCount>
void prepare_field_background_state(
    double h,
    int stage,
    const double (&A)[StageCount][StageCount],
    const double (&ka)[StageCount],
    const double (&kad)[StageCount],
    double& a_stage,
    double& ad_stage) {
    a_stage = a;
    ad_stage = ad;
    for (int s = 0; s < stage; ++s) {
        const double c = A[stage][s];
        if (c == 0.0) continue;
        a_stage += h * c * ka[s];
        ad_stage += h * c * kad[s];
    }
}

// Build and evaluate one complete RK stage for the selected physical system.
template <int StageCount>
void evaluate_field_stage(
    double h,
    int stage,
    const double (&A)[StageCount][StageCount],
    FieldRKWorkspace& workspace,
    size_t site_count,
    double (&ka)[StageCount],
    double (&kad)[StageCount],
    SimulationPhase phase) {
    prepare_field_stage_state(h, stage, A, workspace, site_count);
    double a_stage = a;
    double ad_stage = ad;
    prepare_field_background_state(h, stage, A, ka, kad, a_stage, ad_stage);
#if calculate_SIGW
    compute_field_rhs(
        workspace.ftmp, workspace.fdtmp, workspace.htmp, workspace.hdtmp,
        workspace.kf[stage], workspace.kfd[stage],
        workspace.kh[stage], workspace.khd[stage],
        a_stage, ad_stage, ka[stage], kad[stage], phase);
#else
    compute_field_rhs(
        workspace.ftmp, workspace.fdtmp,
        workspace.kf[stage], workspace.kfd[stage],
        a_stage, ad_stage, ka[stage], kad[stage], phase);
#endif
}

// Commit a weighted combination of RK slopes to every active state variable.
template <int StageCount>
void apply_weighted_state_update(
    double h,
    const double (&B)[StageCount],
    FieldRKWorkspace& workspace,
    size_t site_count,
    const double (&ka)[StageCount],
    const double (&kad)[StageCount]) {
    for (size_t id = 0; id < site_count; ++id) {
        double sum_f = 0.0;
        double sum_fd = 0.0;
        for (int s = 0; s < StageCount; ++s) {
            const double bs = B[s];
            if (bs == 0.0) continue;
            sum_f += bs * workspace.kf[s][id];
            sum_fd += bs * workspace.kfd[s][id];
        }
        f[id] += h * sum_f;
        fd[id] += h * sum_fd;
    }
#if calculate_SIGW
    for (int c = 0; c < 6; ++c) {
        for (size_t id = 0; id < site_count; ++id) {
            double sum_h = 0.0;
            double sum_hd = 0.0;
            for (int s = 0; s < StageCount; ++s) {
                const double bs = B[s];
                if (bs == 0.0) continue;
                sum_h += bs * static_cast<double>(workspace.kh[s][c][id]);
                sum_hd += bs * static_cast<double>(workspace.khd[s][c][id]);
            }
            hij[c][id] = static_cast<float>(static_cast<double>(hij[c][id]) + h * sum_h);
            hijd[c][id] = static_cast<float>(static_cast<double>(hijd[c][id]) + h * sum_hd);
        }
    }
#endif
    for (int s = 0; s < StageCount; ++s) {
        const double bs = B[s];
        if (bs == 0.0) continue;
        a += h * bs * ka[s];
        ad += h * bs * kad[s];
    }
    t += h;
}

} // namespace

// Advance the selected coupled field/background system by one classical RK4 step.
void rk4_step_fields(
    double h,
    FieldRKWorkspace& workspace,
    SimulationPhase phase) {
    const size_t site_count = f.size();
    workspace.ensure_size(site_count);

    double ka[4] = {0.0, 0.0, 0.0, 0.0};
    double kad[4] = {0.0, 0.0, 0.0, 0.0};

#if calculate_SIGW
    compute_field_rhs(
        f, fd, hij, hijd,
        workspace.kf[0], workspace.kfd[0], workspace.kh[0], workspace.khd[0],
        a, ad, ka[0], kad[0], phase);
#else
    compute_field_rhs(
        f, fd, workspace.kf[0], workspace.kfd[0],
        a, ad, ka[0], kad[0], phase);
#endif

    for (int stage = 1; stage < 4; ++stage) {
        evaluate_field_stage(
            h, stage, RK4_A, workspace, site_count, ka, kad, phase);
    }

    apply_weighted_state_update(h, RK4_B, workspace, site_count, ka, kad);
}

// Attempt one Dormand-Prince 5(4) step. Accepted candidates replace the global
// state; rejected candidates leave it unchanged and return a smaller proposal.
bool rk45_step_fields(
    double h,
    double hmax,
    double& h_next,
    FieldRKWorkspace& workspace,
    SimulationPhase phase) {
    const size_t site_count = f.size();
    workspace.ensure_size(site_count);

    double ka[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double kad[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

#if calculate_SIGW
    compute_field_rhs(
        f, fd, hij, hijd,
        workspace.kf[0], workspace.kfd[0], workspace.kh[0], workspace.khd[0],
        a, ad, ka[0], kad[0], phase);
#else
    compute_field_rhs(
        f, fd, workspace.kf[0], workspace.kfd[0],
        a, ad, ka[0], kad[0], phase);
#endif

    for (int stage = 1; stage < 7; ++stage) {
        evaluate_field_stage(
            h, stage, DP_A, workspace, site_count, ka, kad, phase);
    }

    double a5 = a;
    double ad5 = ad;
    double erra = 0.0;
    double errad = 0.0;
    for (int s = 0; s < 7; ++s) {
        const double bs = DP_B[s];
        const double es = DP_E[s];
        if (bs != 0.0) {
            a5 += h * bs * ka[s];
            ad5 += h * bs * kad[s];
        }
        if (es != 0.0) {
            erra += h * es * ka[s];
            errad += h * es * kad[s];
        }
    }

    double err_acc = 0.0;
    std::size_t nvars = 2 * site_count + 2;

    for (size_t id = 0; id < site_count; ++id) {
        double y5f = f[id];
        double y5fd = fd[id];
        double errf = 0.0;
        double errfd = 0.0;
        for (int s = 0; s < 7; ++s) {
            const double bs = DP_B[s];
            const double es = DP_E[s];
            if (bs != 0.0) {
                y5f += h * bs * workspace.kf[s][id];
                y5fd += h * bs * workspace.kfd[s][id];
            }
            if (es != 0.0) {
                errf += h * es * workspace.kf[s][id];
                errfd += h * es * workspace.kfd[s][id];
            }
        }

        workspace.ftmp[id] = y5f;
        workspace.fdtmp[id] = y5fd;

        const double sf = rk45_abs_tol + rk45_rel_tol * std::max(std::abs(f[id]), std::abs(y5f));
        const double sfd = rk45_abs_tol + rk45_rel_tol * std::max(std::abs(fd[id]), std::abs(y5fd));
        err_acc += pw2(errf / sf) + pw2(errfd / sfd);
    }

#if calculate_SIGW
    nvars += 12 * site_count;
    for (int c = 0; c < 6; ++c) {
        for (size_t id = 0; id < site_count; ++id) {
            double y5h = static_cast<double>(hij[c][id]);
            double y5hd = static_cast<double>(hijd[c][id]);
            double errh = 0.0;
            double errhd = 0.0;
            for (int s = 0; s < 7; ++s) {
                const double bs = DP_B[s];
                const double es = DP_E[s];
                if (bs != 0.0) {
                    y5h += h * bs * static_cast<double>(workspace.kh[s][c][id]);
                    y5hd += h * bs * static_cast<double>(workspace.khd[s][c][id]);
                }
                if (es != 0.0) {
                    errh += h * es * static_cast<double>(workspace.kh[s][c][id]);
                    errhd += h * es * static_cast<double>(workspace.khd[s][c][id]);
                }
            }

            workspace.htmp[c][id] = static_cast<float>(y5h);
            workspace.hdtmp[c][id] = static_cast<float>(y5hd);

            const double sh = rk45_abs_tol + rk45_rel_tol * std::max(
                std::abs(static_cast<double>(hij[c][id])), std::abs(y5h));
            const double shd = rk45_abs_tol + rk45_rel_tol * std::max(
                std::abs(static_cast<double>(hijd[c][id])), std::abs(y5hd));
            err_acc += pw2(errh / sh) + pw2(errhd / shd);
        }
    }
#endif

    const double sa = rk45_abs_tol + rk45_rel_tol * std::max(std::abs(a), std::abs(a5));
    const double sad = rk45_abs_tol + rk45_rel_tol * std::max(std::abs(ad), std::abs(ad5));
    err_acc += pw2(erra / sa) + pw2(errad / sad);

    const double err_norm = std::sqrt(err_acc / static_cast<double>(nvars));
    const double safe_err = std::max(err_norm, 1e-16);
    double factor = rk45_safety * std::pow(safe_err, -0.2);
    factor = std::clamp(factor, 0.1, 5.0);

    if (err_norm <= 1.0) {
        f.swap(workspace.ftmp);
        fd.swap(workspace.fdtmp);
#if calculate_SIGW
        for (int c = 0; c < 6; ++c) {
            hij[c].swap(workspace.htmp[c]);
            hijd[c].swap(workspace.hdtmp[c]);
        }
#endif
        a = a5;
        ad = ad5;
        t += h;
        h_next = clamp_rk45_step(h * factor, hmax);
        return true;
    }

    h_next = clamp_rk45_step(h * std::max(0.1, std::min(0.5, factor)), hmax);
    return false;
}
} // namespace evolution
