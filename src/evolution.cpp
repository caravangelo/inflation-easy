// evolution.cpp - Time evolution of fields and background
//
// This file implements the time integration of the scalar field,
// optional auxiliary systems (deltaN, tensor perturbations), and
// the background scale factor using the chosen update scheme.


#include "main.h"

// High-level flow in this module:
// 1) Inflation loop (integrator selectable: leapfrog, RK4, RK45).
// 2) Optional deltaN loop (own selectable integrator, same RK back-end pattern).
// 3) Optional post-inflation loop (reuses inflation RK back-end with alternate RHS).
//
// Numerical invariants to preserve when editing:
// - Periodic finite-difference stencils on all spatial derivatives.
// - Leapfrog staggering semantics (state/derivative half-step offsets).
// - RK45 acceptance based on weighted RMS norm with abs/rel tolerances.
// - Output formats and file names consumed by downstream analysis scripts.

namespace {
// Precomputed periodic neighbors (i+1 mod N, i-1 mod N), reused in hot stencils.
std::vector<int> g_periodic_inc;
std::vector<int> g_periodic_dec;

inline void ensure_periodic_index_cache() {
    if (static_cast<int>(g_periodic_inc.size()) == N &&
        static_cast<int>(g_periodic_dec.size()) == N) {
        return;
    }

    g_periodic_inc.resize(N);
    g_periodic_dec.resize(N);
    for (int i = 0; i < N; ++i) {
        g_periodic_inc[i] = (i == N - 1) ? 0 : i + 1;
        g_periodic_dec[i] = (i == 0) ? N - 1 : i - 1;
    }
}
} // namespace

// -------------------- Laplacians --------------------

// Helper for periodic indexing
inline int INCREMENT(int i) {
    if (static_cast<int>(g_periodic_inc.size()) != N) {
        ensure_periodic_index_cache();
    }
    return g_periodic_inc[i];
}

// Decrement index with periodic wrapping (i → i-1 mod N)
inline int DECREMENT(int i) {
    if (static_cast<int>(g_periodic_dec.size()) != N) {
        ensure_periodic_index_cache();
    }
    return g_periodic_dec[i];
}

// Laplacian for real-space arrays (works for double or float vectors)
template <typename T>
inline T lapl(int i, int j, int k, const std::vector<T>& field) {
    if (i == 0 || j == 0 || k == 0 || i == N - 1 || j == N - 1 || k == N - 1) {
        return (
        field[idx(i,j,INCREMENT(k))] + field[idx(i,j,DECREMENT(k))] +
        field[idx(i,INCREMENT(j),k)] + field[idx(i,DECREMENT(j),k)] +
        field[idx(INCREMENT(i),j,k)] + field[idx(DECREMENT(i),j,k)] -
        T(6) * field[idx(i,j,k)]
        );
    } else {
        return (
        field[idx(i,j,k+1)] + field[idx(i,j,k-1)] +
        field[idx(i,j+1,k)] + field[idx(i,j-1,k)] +
        field[idx(i+1,j,k)] + field[idx(i-1,j,k)] -
        T(6) * field[idx(i,j,k)]
        );
    }
}

#if calculate_SIGW || post_inflation
// Central difference for spatial derivative of a scalar field (double or float)
template <typename T>
inline T dfdx(int dim, int i, int j, int k, const std::vector<T>& field) {
    const double half_over_dx = 0.5 / dx;
    if (dim == 0) {
        if (i == 0 || i == N - 1) {
            return T( (field[idx(INCREMENT(i), j, k)] - field[idx(DECREMENT(i), j, k)]) * half_over_dx );
        } else {
            return T( (field[idx(i+1, j, k)] - field[idx(i-1, j, k)]) * half_over_dx );
        }
    } else if (dim == 1) {
        if (j == 0 || j == N - 1) {
            return T( (field[idx(i, INCREMENT(j), k)] - field[idx(i, DECREMENT(j), k)]) * half_over_dx );
        } else {
            return T( (field[idx(i, j+1, k)] - field[idx(i, j-1, k)]) * half_over_dx );
        }
    } else { // dim == 2
        if (k == 0 || k == N - 1) {
            return T( (field[idx(i, j, INCREMENT(k))] - field[idx(i, j, DECREMENT(k))]) * half_over_dx );
        } else {
            return T( (field[idx(i, j, k+1)] - field[idx(i, j, k-1)]) * half_over_dx );
        }
    }
}

inline double d2_same_state(int dim, int i, int j, int k, size_t id, const std::vector<double>& field) {
    const int ip = (dim == 0) ? INCREMENT(i) : i;
    const int im = (dim == 0) ? DECREMENT(i) : i;
    const int jp = (dim == 1) ? INCREMENT(j) : j;
    const int jm = (dim == 1) ? DECREMENT(j) : j;
    const int kp = (dim == 2) ? INCREMENT(k) : k;
    const int km = (dim == 2) ? DECREMENT(k) : k;
    const double inv_dx2 = 1.0 / (dx * dx);

    if (dim == 0) return (field[idx(ip, j,  k)] - 2.0 * field[id] + field[idx(im, j,  k)]) * inv_dx2;
    if (dim == 1) return (field[idx(i,  jp, k)] - 2.0 * field[id] + field[idx(i,  jm, k)]) * inv_dx2;
    return             (field[idx(i,  j,  kp)] - 2.0 * field[id] + field[idx(i,  j,  km)]) * inv_dx2;
}

inline double d2_cross_state(int d1, int d2, int i, int j, int k, const std::vector<double>& field) {
    const int ip = (d1 == 0 || d2 == 0) ? INCREMENT(i) : i;
    const int im = (d1 == 0 || d2 == 0) ? DECREMENT(i) : i;
    const int jp = (d1 == 1 || d2 == 1) ? INCREMENT(j) : j;
    const int jm = (d1 == 1 || d2 == 1) ? DECREMENT(j) : j;
    const int kp = (d1 == 2 || d2 == 2) ? INCREMENT(k) : k;
    const int km = (d1 == 2 || d2 == 2) ? DECREMENT(k) : k;

    const double fpp = field[idx(ip, jp, kp)];
    const double fpm = field[idx(ip, jm, km)];
    const double fmp = field[idx(im, jp, kp)];
    const double fmm = field[idx(im, jm, km)];
    return (fpp - fpm - fmp + fmm) / (4.0 * dx * dx);
}

// -------------------- Stress-energy tensor (inflaton) --------------------

// Pack of first derivatives for one site (double precision for the scalar field)
struct GradPack {
    double g[3];   // ∂_x φ, ∂_y φ, ∂_z φ
} ;

inline void build_grad_pack(int i, int j, int k, GradPack& P) {
    P.g[0] = dfdx<double>(0, i, j, k, f);
    P.g[1] = dfdx<double>(1, i, j, k, f);
    P.g[2] = dfdx<double>(2, i, j, k, f);
}

// Stress-energy tensor component T_{lm} at site (i,j,k)
inline double stress_energy_fast(int l, int m, const GradPack& P) {
    return P.g[l] * P.g[m];
}
#endif

// -------------------- Energy Calculations --------------------

// Compute gradient energy density (averaged)
double gradient_energy() {
    DECLARE_INDICES
    double gradient = 0.0;
    const double norm = pw2(1.0 / (a * dx));
    LOOP gradient -= f[idx(i,j,k)] * lapl<double>(i, j, k, f);
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

// -------------------- Main Field Evolution --------------------

// Leapfrog update of field derivatives
void evolve_derivs(double d) {
    DECLARE_INDICES

    const double laplnorm = std::pow(a, -2.0 * rescale_s) / pw2(dx);
    const double sfev1    = rescale_s + 1.0;
    const double sfev2    = -2.0 * rescale_s + 2.0;

    // Update second derivative of scale factor (ad2)
    ad2 = (-2.0 * ad - 2.0 * a / d / sfev1 * (
    1.0 - std::sqrt(1.0 + 2.0 * d * sfev1 * ad / a +
    pw2(d) * sfev1 * std::pow(a, sfev2) *
    (2.0 * gradient_energy() / 3.0 + potential_energy()))
    )) / d;

    ad += 0.5 * d * ad2;

    // Scalar field evolution
    const double friction = (2.0 + rescale_s) * ad / a;
    const double potnorm  = std::pow(a, 2.0 - 2.0 * rescale_s);

#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
    LOOP {
        const size_t id = idx(i,j,k);

        fd[id] += d * (
        laplnorm * lapl<double>(i, j, k, f)
        - friction * fd[id]
        - potnorm * potential_derivative(i, j, k)
        );
    }

#if calculate_SIGW
    // Gravitational wave tensor evolution (store as float, compute RHS in double)
    const double srcAmp = 2.0 * std::pow(a, -2.0 * rescale_s);

#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
    for (int i=0; i<N; ++i)
    for (int j=0; j<N; ++j)
    for (int k=0; k<N; ++k) {
        const size_t id = idx(i,j,k);

        // build gradients once
        GradPack G; build_grad_pack(i,j,k,G);

        // explicit six components (double)
        const double gx = G.g[0];
        const double gy = G.g[1];
        const double gz = G.g[2];

        const double T_xx = gx * gx;
        const double T_yy = gy * gy;
        const double T_zz = gz * gz;
        const double T_xy = gx * gy;
        const double T_xz = gx * gz;
        const double T_yz = gy * gz;

        // compute each RHS entirely in double; cast once on store
        {
            const double lap_h = static_cast<double>(lapl<float>(i,j,k,hij[0]));
            const double rhs = d * ( laplnorm * lap_h
            - friction * static_cast<double>(hijd[0][id])
            + srcAmp * T_xx );
            hijd[0][id] += static_cast<float>(rhs);
        }
        {
            const double lap_h = static_cast<double>(lapl<float>(i,j,k,hij[1]));
            const double rhs = d * ( laplnorm * lap_h
            - friction * static_cast<double>(hijd[1][id])
            + srcAmp * T_yy );
            hijd[1][id] += static_cast<float>(rhs);
        }
        {
            const double lap_h = static_cast<double>(lapl<float>(i,j,k,hij[2]));
            const double rhs = d * ( laplnorm * lap_h
            - friction * static_cast<double>(hijd[2][id])
            + srcAmp * T_zz );
            hijd[2][id] += static_cast<float>(rhs);
        }
        {
            const double lap_h = static_cast<double>(lapl<float>(i,j,k,hij[3]));
            const double rhs = d * ( laplnorm * lap_h
            - friction * static_cast<double>(hijd[3][id])
            + srcAmp * T_xy );
            hijd[3][id] += static_cast<float>(rhs);
        }
        {
            const double lap_h = static_cast<double>(lapl<float>(i,j,k,hij[4]));
            const double rhs = d * ( laplnorm * lap_h
            - friction * static_cast<double>(hijd[4][id])
            + srcAmp * T_xz );
            hijd[4][id] += static_cast<float>(rhs);
        }
        {
            const double lap_h = static_cast<double>(lapl<float>(i,j,k,hij[5]));
            const double rhs = d * ( laplnorm * lap_h
            - friction * static_cast<double>(hijd[5][id])
            + srcAmp * T_yz );
            hijd[5][id] += static_cast<float>(rhs);
        }
    }
#endif

    ad += 0.5 * d * ad2;
}

// Update fields using current derivatives
void evolve_fields(double d) {
    DECLARE_INDICES
    t += d;

    // Scalar field update
#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
    LOOP {
        f[idx(i,j,k)] += d * fd[idx(i,j,k)];
    }

#if calculate_SIGW
    // Gravitational wave tensor update
#if parallel_calculation
#pragma omp parallel for collapse(4)
#endif
    for (int comp = 0; comp < 6; ++comp)
    LOOP {
        hij[comp][idx(i,j,k)] += static_cast<float>( d * hijd[comp][idx(i,j,k)] );
    }
#endif

    a += d * ad;
}

namespace {
// Reused stage/state buffers for RK methods to avoid per-step allocations.
struct InflationRKScratch {
    std::vector<double> ftmp;
    std::vector<double> fdtmp;
    std::vector<double> kf[7];
    std::vector<double> kfd[7];
#if calculate_SIGW
    std::vector<float> htmp[6];
    std::vector<float> hdtmp[6];
    std::vector<float> kh[7][6];
    std::vector<float> khd[7][6];
#endif

    void ensure_size(size_t gs) {
        if (ftmp.size() == gs) return;
        ftmp.resize(gs);
        fdtmp.resize(gs);
        for (int s = 0; s < 7; ++s) {
            kf[s].resize(gs);
            kfd[s].resize(gs);
        }
#if calculate_SIGW
        for (int c = 0; c < 6; ++c) {
            htmp[c].resize(gs);
            hdtmp[c].resize(gs);
        }
        for (int s = 0; s < 7; ++s) {
            for (int c = 0; c < 6; ++c) {
                kh[s][c].resize(gs);
                khd[s][c].resize(gs);
            }
        }
#endif
    }
};

// Selected at runtime via a scoped guard in post-inflation loops.
// false: inflation equations, true: post-inflation equations.
bool g_use_post_inflation_rhs = false;

double clamp_rk45_step(double h, double hmax) {
    const double hmin = std::max(1e-16, rk45_min_dt);
    const double hhi = std::max(hmin, hmax);
    if (h < hmin) return hmin;
    if (h > hhi) return hhi;
    return h;
}

constexpr int RK45_MAX_ATTEMPTS = 25;
constexpr double RK45_MIN_STEP_GUARD = 1.0 + 1e-12;

inline double inflation_rescaled_step(double astep_value) {
    return dt * std::pow(astep_value, rescale_s - 1.0);
}

inline double rk45_hmax_from_base_step(double base_step) {
    return std::min(rk45_max_dt, std::max(rk45_min_dt, 2.0 * base_step));
}

template <typename StepFunction, typename FailureHandler>
double rk45_accept_step(double h, double hmax, StepFunction&& step_function, FailureHandler&& failure_handler) {
    // Generic adaptive-step accept/reject driver shared by inflation/deltaN/post-inflation.
    bool accepted = false;
    int attempts = 0;
    while (!accepted) {
        double h_suggested = h;
        accepted = step_function(h, hmax, h_suggested);
        h = h_suggested;
        ++attempts;

        if (!accepted && (attempts > RK45_MAX_ATTEMPTS || h <= rk45_min_dt * RK45_MIN_STEP_GUARD)) {
            failure_handler(attempts, h);
        }
    }
    return h;
}

struct ScopedPostInflationRhsMode {
    explicit ScopedPostInflationRhsMode(bool enabled) {
        g_use_post_inflation_rhs = enabled;
    }
    ~ScopedPostInflationRhsMode() {
        g_use_post_inflation_rhs = false;
    }
};

struct DormandPrince45Coefficients {
    static constexpr double a21 = 1.0 / 5.0;
    static constexpr double a31 = 3.0 / 40.0;
    static constexpr double a32 = 9.0 / 40.0;
    static constexpr double a41 = 44.0 / 45.0;
    static constexpr double a42 = -56.0 / 15.0;
    static constexpr double a43 = 32.0 / 9.0;
    static constexpr double a51 = 19372.0 / 6561.0;
    static constexpr double a52 = -25360.0 / 2187.0;
    static constexpr double a53 = 64448.0 / 6561.0;
    static constexpr double a54 = -212.0 / 729.0;
    static constexpr double a61 = 9017.0 / 3168.0;
    static constexpr double a62 = -355.0 / 33.0;
    static constexpr double a63 = 46732.0 / 5247.0;
    static constexpr double a64 = 49.0 / 176.0;
    static constexpr double a65 = -5103.0 / 18656.0;
    static constexpr double a71 = 35.0 / 384.0;
    static constexpr double a73 = 500.0 / 1113.0;
    static constexpr double a74 = 125.0 / 192.0;
    static constexpr double a75 = -2187.0 / 6784.0;
    static constexpr double a76 = 11.0 / 84.0;

    static constexpr double b1 = 35.0 / 384.0;
    static constexpr double b3 = 500.0 / 1113.0;
    static constexpr double b4 = 125.0 / 192.0;
    static constexpr double b5 = -2187.0 / 6784.0;
    static constexpr double b6 = 11.0 / 84.0;

    static constexpr double bs1 = 5179.0 / 57600.0;
    static constexpr double bs3 = 7571.0 / 16695.0;
    static constexpr double bs4 = 393.0 / 640.0;
    static constexpr double bs5 = -92097.0 / 339200.0;
    static constexpr double bs6 = 187.0 / 2100.0;
    static constexpr double bs7 = 1.0 / 40.0;

    static constexpr double e1 = b1 - bs1;
    static constexpr double e3 = b3 - bs3;
    static constexpr double e4 = b4 - bs4;
    static constexpr double e5 = b5 - bs5;
    static constexpr double e6 = b6 - bs6;
    static constexpr double e7 = -bs7;
};

// Compute RHS for the coupled system (field, optional tensors, background).
// Inputs are read-only state snapshots; outputs are time derivatives at that state.
// Note: in numerical-potential mode, lstart[] is updated as a search hint cache.
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

    if (!g_use_post_inflation_rhs) {
        const double laplnorm = std::pow(a_state, -2.0 * rescale_s) / pw2(dx);
        const double friction = (2.0 + rescale_s) * ad_state / a_state;
        const double potnorm  = std::pow(a_state, 2.0 - 2.0 * rescale_s);
#if calculate_SIGW
        const double srcAmp = 2.0 * std::pow(a_state, -2.0 * rescale_s);
#endif

        double gradient_sum = 0.0;
        double potential_sum = 0.0;

#if parallel_calculation
#pragma omp parallel for collapse(3) reduction(+:gradient_sum,potential_sum)
#endif
        LOOP {
            const size_t id = idx(i, j, k);
            const double field_here = f_state[id];
            const double lap_f = lapl<double>(i, j, k, f_state);
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
            dfddt[id] = laplnorm * lap_f
                      - friction * fd_state[id]
                      - potnorm * pot_deriv_here;

            gradient_sum -= field_here * lap_f;
            potential_sum += pot_here;

#if calculate_SIGW
            const double gx = dfdx<double>(0, i, j, k, f_state);
            const double gy = dfdx<double>(1, i, j, k, f_state);
            const double gz = dfdx<double>(2, i, j, k, f_state);

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

            dhddt[0][id] = static_cast<float>(laplnorm * static_cast<double>(lapl<float>(i, j, k, h_state[0]))
                        - friction * static_cast<double>(hd_state[0][id]) + srcAmp * T_xx);
            dhddt[1][id] = static_cast<float>(laplnorm * static_cast<double>(lapl<float>(i, j, k, h_state[1]))
                        - friction * static_cast<double>(hd_state[1][id]) + srcAmp * T_yy);
            dhddt[2][id] = static_cast<float>(laplnorm * static_cast<double>(lapl<float>(i, j, k, h_state[2]))
                        - friction * static_cast<double>(hd_state[2][id]) + srcAmp * T_zz);
            dhddt[3][id] = static_cast<float>(laplnorm * static_cast<double>(lapl<float>(i, j, k, h_state[3]))
                        - friction * static_cast<double>(hd_state[3][id]) + srcAmp * T_xy);
            dhddt[4][id] = static_cast<float>(laplnorm * static_cast<double>(lapl<float>(i, j, k, h_state[4]))
                        - friction * static_cast<double>(hd_state[4][id]) + srcAmp * T_xz);
            dhddt[5][id] = static_cast<float>(laplnorm * static_cast<double>(lapl<float>(i, j, k, h_state[5]))
                        - friction * static_cast<double>(hd_state[5][id]) + srcAmp * T_yz);
#endif
        }

        const double gradient = 0.5 * gradient_sum * pw2(1.0 / (a_state * dx))
            / static_cast<double>(gridsize);
        const double potential_avg = potential_sum / static_cast<double>(gridsize);
        const double source = 2.0 * gradient / 3.0 + potential_avg;

        dadt = ad_state;
        daddt = std::pow(a_state, 3.0 - 2.0 * rescale_s) * source
            - (rescale_s + 1.0) * pw2(ad_state) / a_state;
        return;
    }

#if post_inflation
    const double laplnorm = std::pow(a_state, -2.0 * rescale_s) / pw2(dx);
    const double scalar_friction = (3.0 * (1.0 + omega) + rescale_s) * ad_state / a_state;
#if calculate_SIGW
    const double tensor_friction = (2.0 + rescale_s) * ad_state / a_state;
    const double Htilde = ad_state / a_state;
    const double invH   = 1.0 / Htilde;
    const double coeffU = 4.0 / (3.0 * (1.0 + omega));
    const double srcAmp = 2.0 * std::pow(a_state, -2.0 * rescale_s);
#endif

#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
    LOOP {
        const size_t id = idx(i, j, k);
        const double lap_f = lapl<double>(i, j, k, f_state);

        dfdt[id] = fd_state[id];
        dfddt[id] = omega * laplnorm * lap_f - scalar_friction * fd_state[id];

#if calculate_SIGW
        const double gx = dfdx<double>(0, i, j, k, f_state);
        const double gy = dfdx<double>(1, i, j, k, f_state);
        const double gz = dfdx<double>(2, i, j, k, f_state);
        const double gfdx = dfdx<double>(0, i, j, k, fd_state);
        const double gfdy = dfdx<double>(1, i, j, k, fd_state);
        const double gfdz = dfdx<double>(2, i, j, k, fd_state);

        const double Hxx = d2_same_state(0, i, j, k, id, f_state);
        const double Hyy = d2_same_state(1, i, j, k, id, f_state);
        const double Hzz = d2_same_state(2, i, j, k, id, f_state);
        const double Hxy = d2_cross_state(0, 1, i, j, k, f_state);
        const double Hxz = d2_cross_state(0, 2, i, j, k, f_state);
        const double Hyz = d2_cross_state(1, 2, i, j, k, f_state);

        const double dux = gfdx * invH + gx;
        const double duy = gfdy * invH + gy;
        const double duz = gfdz * invH + gz;
        const double phi_here = f_state[id];

        const double S_xx = 4.0 * phi_here * Hxx + 2.0 * gx * gx - coeffU * (dux * dux);
        const double S_yy = 4.0 * phi_here * Hyy + 2.0 * gy * gy - coeffU * (duy * duy);
        const double S_zz = 4.0 * phi_here * Hzz + 2.0 * gz * gz - coeffU * (duz * duz);
        const double S_xy = 4.0 * phi_here * Hxy + 2.0 * gx * gy - coeffU * (dux * duy);
        const double S_xz = 4.0 * phi_here * Hxz + 2.0 * gx * gz - coeffU * (dux * duz);
        const double S_yz = 4.0 * phi_here * Hyz + 2.0 * gy * gz - coeffU * (duy * duz);

        dhdt[0][id] = hd_state[0][id];
        dhdt[1][id] = hd_state[1][id];
        dhdt[2][id] = hd_state[2][id];
        dhdt[3][id] = hd_state[3][id];
        dhdt[4][id] = hd_state[4][id];
        dhdt[5][id] = hd_state[5][id];

        dhddt[0][id] = static_cast<float>(laplnorm * static_cast<double>(lapl<float>(i, j, k, h_state[0]))
                    - tensor_friction * static_cast<double>(hd_state[0][id]) + srcAmp * S_xx);
        dhddt[1][id] = static_cast<float>(laplnorm * static_cast<double>(lapl<float>(i, j, k, h_state[1]))
                    - tensor_friction * static_cast<double>(hd_state[1][id]) + srcAmp * S_yy);
        dhddt[2][id] = static_cast<float>(laplnorm * static_cast<double>(lapl<float>(i, j, k, h_state[2]))
                    - tensor_friction * static_cast<double>(hd_state[2][id]) + srcAmp * S_zz);
        dhddt[3][id] = static_cast<float>(laplnorm * static_cast<double>(lapl<float>(i, j, k, h_state[3]))
                    - tensor_friction * static_cast<double>(hd_state[3][id]) + srcAmp * S_xy);
        dhddt[4][id] = static_cast<float>(laplnorm * static_cast<double>(lapl<float>(i, j, k, h_state[4]))
                    - tensor_friction * static_cast<double>(hd_state[4][id]) + srcAmp * S_xz);
        dhddt[5][id] = static_cast<float>(laplnorm * static_cast<double>(lapl<float>(i, j, k, h_state[5]))
                    - tensor_friction * static_cast<double>(hd_state[5][id]) + srcAmp * S_yz);
#endif
    }

    dadt = ad_state;
    daddt = - (rescale_s - 0.5 * (1.0 - 3.0 * omega)) * pw2(ad_state) / a_state;
#endif
}

static constexpr double RK4_A[4][4] = {
    {0.0, 0.0, 0.0, 0.0},
    {0.5, 0.0, 0.0, 0.0},
    {0.0, 0.5, 0.0, 0.0},
    {0.0, 0.0, 1.0, 0.0}
};

static constexpr double RK4_B[4] = {
    1.0 / 6.0,
    1.0 / 3.0,
    1.0 / 3.0,
    1.0 / 6.0
};

static constexpr double DP_A[7][7] = {
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {DormandPrince45Coefficients::a21, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {DormandPrince45Coefficients::a31, DormandPrince45Coefficients::a32, 0.0, 0.0, 0.0, 0.0, 0.0},
    {DormandPrince45Coefficients::a41, DormandPrince45Coefficients::a42, DormandPrince45Coefficients::a43, 0.0, 0.0, 0.0, 0.0},
    {DormandPrince45Coefficients::a51, DormandPrince45Coefficients::a52, DormandPrince45Coefficients::a53, DormandPrince45Coefficients::a54, 0.0, 0.0, 0.0},
    {DormandPrince45Coefficients::a61, DormandPrince45Coefficients::a62, DormandPrince45Coefficients::a63, DormandPrince45Coefficients::a64, DormandPrince45Coefficients::a65, 0.0, 0.0},
    {DormandPrince45Coefficients::a71, 0.0, DormandPrince45Coefficients::a73, DormandPrince45Coefficients::a74, DormandPrince45Coefficients::a75, DormandPrince45Coefficients::a76, 0.0}
};

static constexpr double DP_B[7] = {
    DormandPrince45Coefficients::b1,
    0.0,
    DormandPrince45Coefficients::b3,
    DormandPrince45Coefficients::b4,
    DormandPrince45Coefficients::b5,
    DormandPrince45Coefficients::b6,
    0.0
};

static constexpr double DP_E[7] = {
    DormandPrince45Coefficients::e1,
    0.0,
    DormandPrince45Coefficients::e3,
    DormandPrince45Coefficients::e4,
    DormandPrince45Coefficients::e5,
    DormandPrince45Coefficients::e6,
    DormandPrince45Coefficients::e7
};

template <int StageCount>
void prepare_inflation_stage_state(double h, int stage, const double (&A)[StageCount][StageCount], InflationRKScratch& scratch, size_t gs) {
    for (size_t id = 0; id < gs; ++id) {
        double sum_f = 0.0;
        double sum_fd = 0.0;
        for (int s = 0; s < stage; ++s) {
            const double c = A[stage][s];
            if (c == 0.0) continue;
            sum_f += c * scratch.kf[s][id];
            sum_fd += c * scratch.kfd[s][id];
        }
        scratch.ftmp[id] = f[id] + h * sum_f;
        scratch.fdtmp[id] = fd[id] + h * sum_fd;
    }
#if calculate_SIGW
    for (int c = 0; c < 6; ++c) {
        for (size_t id = 0; id < gs; ++id) {
            double sum_h = 0.0;
            double sum_hd = 0.0;
            for (int s = 0; s < stage; ++s) {
                const double cs = A[stage][s];
                if (cs == 0.0) continue;
                sum_h += cs * static_cast<double>(scratch.kh[s][c][id]);
                sum_hd += cs * static_cast<double>(scratch.khd[s][c][id]);
            }
            scratch.htmp[c][id] = static_cast<float>(static_cast<double>(hij[c][id]) + h * sum_h);
            scratch.hdtmp[c][id] = static_cast<float>(static_cast<double>(hijd[c][id]) + h * sum_hd);
        }
    }
#endif
}

template <int StageCount>
void prepare_inflation_background_state(
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

template <int StageCount>
void evaluate_inflation_stage(
    double h,
    int stage,
    const double (&A)[StageCount][StageCount],
    InflationRKScratch& scratch,
    size_t gs,
    double (&ka)[StageCount],
    double (&kad)[StageCount]) {
    prepare_inflation_stage_state(h, stage, A, scratch, gs);
    double a_stage = a;
    double ad_stage = ad;
    prepare_inflation_background_state(h, stage, A, ka, kad, a_stage, ad_stage);
#if calculate_SIGW
    compute_inflation_rhs(
        scratch.ftmp, scratch.fdtmp, scratch.htmp, scratch.hdtmp,
        scratch.kf[stage], scratch.kfd[stage], scratch.kh[stage], scratch.khd[stage],
        a_stage, ad_stage, ka[stage], kad[stage]);
#else
    compute_inflation_rhs(
        scratch.ftmp, scratch.fdtmp, scratch.kf[stage], scratch.kfd[stage],
        a_stage, ad_stage, ka[stage], kad[stage]);
#endif
}

template <int StageCount>
void apply_inflation_weighted_update(
    double h,
    const double (&B)[StageCount],
    InflationRKScratch& scratch,
    size_t gs,
    const double (&ka)[StageCount],
    const double (&kad)[StageCount]) {
    for (size_t id = 0; id < gs; ++id) {
        double sum_f = 0.0;
        double sum_fd = 0.0;
        for (int s = 0; s < StageCount; ++s) {
            const double bs = B[s];
            if (bs == 0.0) continue;
            sum_f += bs * scratch.kf[s][id];
            sum_fd += bs * scratch.kfd[s][id];
        }
        f[id] += h * sum_f;
        fd[id] += h * sum_fd;
    }
#if calculate_SIGW
    for (int c = 0; c < 6; ++c) {
        for (size_t id = 0; id < gs; ++id) {
            double sum_h = 0.0;
            double sum_hd = 0.0;
            for (int s = 0; s < StageCount; ++s) {
                const double bs = B[s];
                if (bs == 0.0) continue;
                sum_h += bs * static_cast<double>(scratch.kh[s][c][id]);
                sum_hd += bs * static_cast<double>(scratch.khd[s][c][id]);
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

void rk4_step_inflation(double h, InflationRKScratch& scratch) {
    // Classical RK4 single accepted step with fixed h.
    const size_t gs = f.size();
    scratch.ensure_size(gs);

    double ka[4] = {0.0, 0.0, 0.0, 0.0};
    double kad[4] = {0.0, 0.0, 0.0, 0.0};

#if calculate_SIGW
    compute_inflation_rhs(f, fd, hij, hijd, scratch.kf[0], scratch.kfd[0], scratch.kh[0], scratch.khd[0], a, ad, ka[0], kad[0]);
#else
    compute_inflation_rhs(f, fd, scratch.kf[0], scratch.kfd[0], a, ad, ka[0], kad[0]);
#endif

    for (int stage = 1; stage < 4; ++stage) {
        evaluate_inflation_stage(h, stage, RK4_A, scratch, gs, ka, kad);
    }

    apply_inflation_weighted_update(h, RK4_B, scratch, gs, ka, kad);
}

bool rk45_step_inflation(double h, double hmax, double& h_next, InflationRKScratch& scratch) {
    // Dormand-Prince 5(4): computes candidate 5th-order state + embedded error estimate.
    // On success: commits state and proposes next h. On failure: keeps state, shrinks h.
    const size_t gs = f.size();
    scratch.ensure_size(gs);

    double ka[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double kad[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

#if calculate_SIGW
    compute_inflation_rhs(f, fd, hij, hijd, scratch.kf[0], scratch.kfd[0], scratch.kh[0], scratch.khd[0], a, ad, ka[0], kad[0]);
#else
    compute_inflation_rhs(f, fd, scratch.kf[0], scratch.kfd[0], a, ad, ka[0], kad[0]);
#endif

    for (int stage = 1; stage < 7; ++stage) {
        evaluate_inflation_stage(h, stage, DP_A, scratch, gs, ka, kad);
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
    std::size_t nvars = 2 * gs + 2;

    for (size_t id = 0; id < gs; ++id) {
        double y5f = f[id];
        double y5fd = fd[id];
        double errf = 0.0;
        double errfd = 0.0;
        for (int s = 0; s < 7; ++s) {
            const double bs = DP_B[s];
            const double es = DP_E[s];
            if (bs != 0.0) {
                y5f += h * bs * scratch.kf[s][id];
                y5fd += h * bs * scratch.kfd[s][id];
            }
            if (es != 0.0) {
                errf += h * es * scratch.kf[s][id];
                errfd += h * es * scratch.kfd[s][id];
            }
        }

        scratch.ftmp[id] = y5f;
        scratch.fdtmp[id] = y5fd;

        const double sf = rk45_abs_tol + rk45_rel_tol * std::max(std::abs(f[id]), std::abs(y5f));
        const double sfd = rk45_abs_tol + rk45_rel_tol * std::max(std::abs(fd[id]), std::abs(y5fd));
        err_acc += pw2(errf / sf) + pw2(errfd / sfd);
    }

#if calculate_SIGW
    nvars += 12 * gs;
    for (int c = 0; c < 6; ++c) {
        for (size_t id = 0; id < gs; ++id) {
            double y5h = static_cast<double>(hij[c][id]);
            double y5hd = static_cast<double>(hijd[c][id]);
            double errh = 0.0;
            double errhd = 0.0;
            for (int s = 0; s < 7; ++s) {
                const double bs = DP_B[s];
                const double es = DP_E[s];
                if (bs != 0.0) {
                    y5h += h * bs * static_cast<double>(scratch.kh[s][c][id]);
                    y5hd += h * bs * static_cast<double>(scratch.khd[s][c][id]);
                }
                if (es != 0.0) {
                    errh += h * es * static_cast<double>(scratch.kh[s][c][id]);
                    errhd += h * es * static_cast<double>(scratch.khd[s][c][id]);
                }
            }

            scratch.htmp[c][id] = static_cast<float>(y5h);
            scratch.hdtmp[c][id] = static_cast<float>(y5hd);

            const double sh = rk45_abs_tol + rk45_rel_tol * std::max(std::abs(static_cast<double>(hij[c][id])), std::abs(y5h));
            const double shd = rk45_abs_tol + rk45_rel_tol * std::max(std::abs(static_cast<double>(hijd[c][id])), std::abs(y5hd));
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
        f.swap(scratch.ftmp);
        fd.swap(scratch.fdtmp);
#if calculate_SIGW
        for (int c = 0; c < 6; ++c) {
            hij[c].swap(scratch.htmp[c]);
            hijd[c].swap(scratch.hdtmp[c]);
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
} // namespace

// -------------------- Main Evolution Loop --------------------

void run_evolution_loop(FILE* output_) {
    // Main inflation driver:
    // - selects integrator backend
    // - performs periodic outputs/logging
    // - guarantees final saved state is synchronized for output
    ensure_periodic_index_cache();
    int numsteps = 0;

    InflationRKScratch rk_scratch;
    rk_scratch.ensure_size(f.size());

    const auto report_step = [&](int step_index) {
        if (step_index % output_freq == 0 && a < af) {
            save((step_index % output_infrequent_freq == 0) ? 1 : 0);
        }
        if (screen_updates && step_index % output_freq == 0) {
            printf("scale factor a = %f\n", a);
            printf("numsteps %i\n\n", step_index);
        }
        fprintf(output_, "scale factor a = %f\n", a);
        fprintf(output_, "numsteps %i\n\n", step_index);
        if (step_index % output_freq == 0) {
            fflush(output_);
        }
    };

    switch (integrator) {
        case INTEGRATOR_LEAPFROG: {
            evolve_fields(0.5 * dt); // First leapfrog step

            while (a <= af) {
                const double dt_rescaled = inflation_rescaled_step(astep);
                evolve_derivs(dt_rescaled);
                evolve_fields(dt_rescaled);

                numsteps++;
                report_step(numsteps);
                astep = a;
            }
            break;
        }
        case INTEGRATOR_RK4: {
            while (a <= af) {
                const double dt_rescaled = inflation_rescaled_step(astep);
                rk4_step_inflation(dt_rescaled, rk_scratch);

                numsteps++;
                report_step(numsteps);
                astep = a;
            }
            break;
        }
        case INTEGRATOR_RK45:
        default: {
            double h = clamp_rk45_step(inflation_rescaled_step(astep), rk45_max_dt);

            while (a <= af) {
                const double base_step = inflation_rescaled_step(astep);
                const double hmax = rk45_hmax_from_base_step(base_step);
                h = clamp_rk45_step(h, hmax);

                h = rk45_accept_step(
                    h,
                    hmax,
                    [&](double h_trial, double h_limit, double& h_next) {
                        return rk45_step_inflation(h_trial, h_limit, h_next, rk_scratch);
                    },
                    [&](int attempts, double failed_h) {
                        std::fprintf(stderr,
                            "RK45 failed to converge at a=%e, t=%e (attempts=%d, h=%e)\n",
                            a, t, attempts, failed_h);
                        std::exit(1);
                    }
                );

                numsteps++;
                report_step(numsteps);
                astep = a;
            }
            break;
        }
    }

    printf("Saving final inflaton data\n");
    fflush(output_);
    save(1);
    if (inflation_uses_staggered_derivatives()) {
        evolve_fields(-0.5 * inflation_rescaled_step(astep)); // sync field values with velocities
    }
    save_last();
}

#if perform_deltaN

// -------------------- DeltaN Evolution --------------------

// Update fd in e-folding time coordinates
void evolve_derivsN(double d) {
    DECLARE_INDICES
#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
    LOOP {
        fd[idx(i,j,k)] += d * (-(3.0 - 0.5 * pw2(fd[idx(i,j,k)])) * (fd[idx(i,j,k)] + pot_ratio(i, j, k)));
    }
}

// Update f and deltaN grid
void evolve_fieldsN(double d) {
    DECLARE_INDICES
    t += d;

#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
    LOOP {
#if monotonic_potential
        if (std::abs(f[idx(i,j,k)]) > std::abs(phiref))
#elif antimonotonic_potential
        if (std::abs(f[idx(i,j,k)]) < std::abs(phiref))
#else
        if (potential(f[idx(i,j,k)]) > potential(phiref))
#endif
        {
            deltaN[idx(i,j,k)] += d;
            f[idx(i,j,k)] += d * fd[idx(i,j,k)];
        }
    }

    a += d * ad;
}

namespace {
struct DeltaNRKScratch {
    std::vector<double> ftmp;
    std::vector<double> fdtmp;
    std::vector<double> dntmp;
    std::vector<double> kf[7];
    std::vector<double> kfd[7];
    std::vector<double> kdn[7];

    void ensure_size(size_t gs) {
        if (ftmp.size() == gs) return;
        ftmp.resize(gs);
        fdtmp.resize(gs);
        dntmp.resize(gs);
        for (int s = 0; s < 7; ++s) {
            kf[s].resize(gs);
            kfd[s].resize(gs);
            kdn[s].resize(gs);
        }
    }
};

double g_deltaN_phiref_potential = 0.0;

// deltaN RHS in e-folding-time coordinates.
// Stopping criterion is encoded in ddNdt (active sites evolve N, inactive sites freeze N).
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
        bool is_active = false;
#if monotonic_potential
        is_active = std::abs(field_here) > std::abs(phiref);
#elif antimonotonic_potential
        is_active = std::abs(field_here) < std::abs(phiref);
#else
        is_active = pot_here > g_deltaN_phiref_potential;
#endif

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

template <int StageCount>
void prepare_deltaN_stage_state(double h, int stage, const double (&A)[StageCount][StageCount], DeltaNRKScratch& scratch, size_t gs) {
    for (size_t id = 0; id < gs; ++id) {
        double sum_f = 0.0;
        double sum_fd = 0.0;
        double sum_dn = 0.0;
        for (int s = 0; s < stage; ++s) {
            const double c = A[stage][s];
            if (c == 0.0) continue;
            sum_f += c * scratch.kf[s][id];
            sum_fd += c * scratch.kfd[s][id];
            sum_dn += c * scratch.kdn[s][id];
        }
        scratch.ftmp[id] = f[id] + h * sum_f;
        scratch.fdtmp[id] = fd[id] + h * sum_fd;
        scratch.dntmp[id] = deltaN[id] + h * sum_dn;
    }
}

template <int StageCount>
void evaluate_deltaN_stage(
    double h,
    int stage,
    const double (&A)[StageCount][StageCount],
    DeltaNRKScratch& scratch,
    size_t gs,
    double (&ka)[StageCount]) {
    prepare_deltaN_stage_state(h, stage, A, scratch, gs);
    compute_deltaN_rhs(scratch.ftmp, scratch.fdtmp, scratch.kf[stage], scratch.kfd[stage], scratch.kdn[stage], ka[stage]);
}

template <int StageCount>
void apply_deltaN_weighted_update(
    double h,
    const double (&B)[StageCount],
    DeltaNRKScratch& scratch,
    size_t gs,
    const double (&ka)[StageCount]) {
    for (size_t id = 0; id < gs; ++id) {
        double sum_f = 0.0;
        double sum_fd = 0.0;
        double sum_dn = 0.0;
        for (int s = 0; s < StageCount; ++s) {
            const double bs = B[s];
            if (bs == 0.0) continue;
            sum_f += bs * scratch.kf[s][id];
            sum_fd += bs * scratch.kfd[s][id];
            sum_dn += bs * scratch.kdn[s][id];
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

void rk4_step_deltaN(double h, DeltaNRKScratch& scratch) {
    // Classical RK4 single accepted step with fixed h.
    const size_t gs = f.size();
    scratch.ensure_size(gs);

    double ka[4] = {0.0, 0.0, 0.0, 0.0};
    compute_deltaN_rhs(f, fd, scratch.kf[0], scratch.kfd[0], scratch.kdn[0], ka[0]);

    for (int stage = 1; stage < 4; ++stage) {
        evaluate_deltaN_stage(h, stage, RK4_A, scratch, gs, ka);
    }

    apply_deltaN_weighted_update(h, RK4_B, scratch, gs, ka);
}

bool rk45_step_deltaN(double h, double hmax, double& h_next, DeltaNRKScratch& scratch) {
    // Dormand-Prince 5(4) with the same error-controller policy as inflation RK45.
    const size_t gs = f.size();
    scratch.ensure_size(gs);

    double ka[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    compute_deltaN_rhs(f, fd, scratch.kf[0], scratch.kfd[0], scratch.kdn[0], ka[0]);

    for (int stage = 1; stage < 7; ++stage) {
        evaluate_deltaN_stage(h, stage, DP_A, scratch, gs, ka);
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
    std::size_t nvars = 3 * gs + 1;

    for (size_t id = 0; id < gs; ++id) {
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
                y5f += h * bs * scratch.kf[s][id];
                y5fd += h * bs * scratch.kfd[s][id];
                y5dn += h * bs * scratch.kdn[s][id];
            }
            if (es != 0.0) {
                errf += h * es * scratch.kf[s][id];
                errfd += h * es * scratch.kfd[s][id];
                errdn += h * es * scratch.kdn[s][id];
            }
        }

        scratch.ftmp[id] = y5f;
        scratch.fdtmp[id] = y5fd;
        scratch.dntmp[id] = y5dn;

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
        f.swap(scratch.ftmp);
        fd.swap(scratch.fdtmp);
        deltaN.swap(scratch.dntmp);
        a = a5;
        t += h;
        h_next = clamp_rk45_step(h * factor, hmax);
        return true;
    }

    h_next = clamp_rk45_step(h * std::max(0.1, std::min(0.5, factor)), hmax);
    return false;
}
} // namespace

// Determine reference φ value for ending deltaN integration
double get_phiref() {
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

// Run main deltaN integration loop
void run_deltaN_loop(FILE* output_) {
    // Main deltaN driver with per-mode integrator dispatch.
    ensure_periodic_index_cache();
    printf("Starting deltaN calculation\n");
    fprintf(output_, "Starting deltaN calculation\n");

    int numsteps = 0;
    Ne = 0.0;

    initializeN();
    if (deltaN_uses_staggered_derivatives()) {
        evolve_fieldsN(0.5 * dN);
    }

    phiref = use_phiref_manual ? phiref_manual_value : get_phiref();
    g_deltaN_phiref_potential = potential(phiref);

    DeltaNRKScratch rk_scratch;
    rk_scratch.ensure_size(f.size());

    const auto report_step = [&](int step_index) {
        if (screen_updates && step_index % output_freq == 0) {
            printf("N = %f\n\n", Ne);
        }

        fprintf(output_, "N = %f\n\n", Ne);
        if (step_index % output_freq == 0) {
            fflush(output_);
        }
    };

    switch (deltaN_integrator) {
        case INTEGRATOR_LEAPFROG: {
            while (Ne <= Nend) {
                evolve_derivsN(dN);
                evolve_fieldsN(dN);
                Ne += dN;
                numsteps++;
                report_step(numsteps);
            }
            break;
        }
        case INTEGRATOR_RK4: {
            while (Ne <= Nend) {
                rk4_step_deltaN(dN, rk_scratch);
                Ne += dN;
                numsteps++;
                report_step(numsteps);
            }
            break;
        }
        case INTEGRATOR_RK45:
        default: {
            double h = clamp_rk45_step(dN, rk45_max_dt);

            while (Ne <= Nend) {
                const double hmax = rk45_hmax_from_base_step(dN);
                h = clamp_rk45_step(h, hmax);

                h = rk45_accept_step(
                    h,
                    hmax,
                    [&](double h_trial, double h_limit, double& h_next) {
                        return rk45_step_deltaN(h_trial, h_limit, h_next, rk_scratch);
                    },
                    [&](int attempts, double failed_h) {
                        std::fprintf(stderr,
                            "RK45 failed to converge in deltaN loop at N=%e, t=%e (attempts=%d, h=%e)\n",
                            Ne, t, attempts, failed_h);
                        std::exit(1);
                    }
                );

                Ne += h;
                numsteps++;
                report_step(numsteps);
            }
            break;
        }
    }

    saveN();
    fflush(output_);
    if (deltaN_uses_staggered_derivatives()) {
        evolve_fieldsN(-0.5 * dN);
    }
}
#endif

#if post_inflation
// -------------------- Post-inflation Evolution --------------------

// ===== derivative pack (uses existing dfdx, idx, INCREMENT/DECREMENT, N, dx) =====
struct DerivPack {
    double gf[3], gfd[3];  // ∂_i f, ∂_i fd
    double Hf[6];          // Hessian(f): [xx, yy, zz, xy, xz, yz] in your sym_idx order
    double f_here;
} ;

// second derivatives on a double field
inline double d2_same(int dim, int i, int j, int k, const std::vector<double>& A) {
    const size_t id = idx(i, j, k);
    return d2_same_state(dim, i, j, k, id, A);
}

inline double d2_cross(int d1, int d2, int i, int j, int k, const std::vector<double>& A) {
    return d2_cross_state(d1, d2, i, j, k, A);
}

inline void build_deriv_pack(int i, int j, int k, DerivPack& P) {
    // first derivatives
    P.gf[0]  = dfdx<double>(0,i,j,k,f);
    P.gf[1]  = dfdx<double>(1,i,j,k,f);
    P.gf[2]  = dfdx<double>(2,i,j,k,f);
    P.gfd[0] = dfdx<double>(0,i,j,k,fd);
    P.gfd[1] = dfdx<double>(1,i,j,k,fd);
    P.gfd[2] = dfdx<double>(2,i,j,k,fd);

    // Hessian entries in your sym_idx packing
    P.Hf[sym_idx(0,0)] = d2_same(0,i,j,k,f);   // xx
    P.Hf[sym_idx(1,1)] = d2_same(1,i,j,k,f);   // yy
    P.Hf[sym_idx(2,2)] = d2_same(2,i,j,k,f);   // zz
    P.Hf[sym_idx(0,1)] = d2_cross(0,1,i,j,k,f);// xy
    P.Hf[sym_idx(0,2)] = d2_cross(0,2,i,j,k,f);// xz
    P.Hf[sym_idx(1,2)] = d2_cross(1,2,i,j,k,f);// yz

    P.f_here = f[idx(i,j,k)];
}


// ===== fast source using sym_idx() and comp_to_indices() =====
inline double stress_energy_post_inflation_fast(int l, int m, const DerivPack& P) {
    const double Htilde = ad / a;             // \tilde{\mathcal H}
    const double invH   = 1.0 / Htilde;
    const double coeffU = 4.0 / (3.0 * (1.0 + omega));

    const double dPhi_l   = P.gf[l];
    const double dPhi_m   = P.gf[m];
    const double d2Phi_lm = P.Hf[sym_idx(l,m)];
    const double dU_l     = P.gfd[l] * invH + dPhi_l;
    const double dU_m     = P.gfd[m] * invH + dPhi_m;

    return 4.0 * P.f_here * d2Phi_lm
    + 2.0 * dPhi_l * dPhi_m
    - coeffU * (dU_l * dU_m);
}

// Leapfrog update of field derivatives
void evolve_derivs_post_inflation(double d) {
    DECLARE_INDICES

    const double laplnorm = std::pow(a, -2.0 * rescale_s) / pw2(dx);

    // Update second derivative of scale factor (ad2)
    ad2 = - ( rescale_s - 0.5*(1.0 - 3.0*omega)) * pw2(ad) / a;

    ad += 0.5 * d * ad2;

    // Scalar field evolution
#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
    LOOP {
        fd[idx(i,j,k)] += d * (
        omega * laplnorm * lapl<double>(i, j, k, f)
        - (3.0 * (1.0 + omega) + rescale_s) * ad * fd[idx(i,j,k)] / a
        );
    }

#if calculate_SIGW
    // Precompute time-dependent factors once per call (hot loop optimization)
    const double Htilde = ad / a;
    const double invH   = 1.0 / Htilde;
    const double coeffU = 4.0 / (3.0 * (1.0 + omega));

    // RHS amplitude (C=2 convention)
    const double srcAmp = 2.0 * std::pow(a, -2.0 * rescale_s);

#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
    for (int i=0; i<N; ++i)
    for (int j=0; j<N; ++j)
    for (int k=0; k<N; ++k) {
        const size_t id = idx(i,j,k);

        // build derivatives once per site
        DerivPack P;
        build_deriv_pack(i,j,k,P);

        // bare sources (not TT-projected), explicit pairs
        const double Phi = P.f_here;

        const double gx = P.gf[0];
        const double gy = P.gf[1];
        const double gz = P.gf[2];

        const double dux = P.gfd[0] * invH + gx;
        const double duy = P.gfd[1] * invH + gy;
        const double duz = P.gfd[2] * invH + gz;

        // Hessian packing is [xx, yy, zz, xy, xz, yz]
        const double Hxx = P.Hf[0];
        const double Hyy = P.Hf[1];
        const double Hzz = P.Hf[2];
        const double Hxy = P.Hf[3];
        const double Hxz = P.Hf[4];
        const double Hyz = P.Hf[5];

        const double S_xx = 4.0 * Phi * Hxx + 2.0 * gx * gx - coeffU * (dux * dux);
        const double S_yy = 4.0 * Phi * Hyy + 2.0 * gy * gy - coeffU * (duy * duy);
        const double S_zz = 4.0 * Phi * Hzz + 2.0 * gz * gz - coeffU * (duz * duz);
        const double S_xy = 4.0 * Phi * Hxy + 2.0 * gx * gy - coeffU * (dux * duy);
        const double S_xz = 4.0 * Phi * Hxz + 2.0 * gx * gz - coeffU * (dux * duz);
        const double S_yz = 4.0 * Phi * Hyz + 2.0 * gy * gz - coeffU * (duy * duz);

        // ---- explicit 6 component updates ----
        {
            const double lap_h = static_cast<double>(lapl<float>(i, j, k, hij[0]));
            const double rhs = d * (
            laplnorm * lap_h
            - (2.0 + rescale_s) * ad * static_cast<double>(hijd[0][id]) / a
            + srcAmp * S_xx
            );
            hijd[0][id] += static_cast<float>(rhs);
        }
        {
            const double lap_h = static_cast<double>(lapl<float>(i, j, k, hij[1]));
            const double rhs = d * (
            laplnorm * lap_h
            - (2.0 + rescale_s) * ad * static_cast<double>(hijd[1][id]) / a
            + srcAmp * S_yy
            );
            hijd[1][id] += static_cast<float>(rhs);
        }
        {
            const double lap_h = static_cast<double>(lapl<float>(i, j, k, hij[2]));
            const double rhs = d * (
            laplnorm * lap_h
            - (2.0 + rescale_s) * ad * static_cast<double>(hijd[2][id]) / a
            + srcAmp * S_zz
            );
            hijd[2][id] += static_cast<float>(rhs);
        }
        {
            const double lap_h = static_cast<double>(lapl<float>(i, j, k, hij[3]));
            const double rhs = d * (
            laplnorm * lap_h
            - (2.0 + rescale_s) * ad * static_cast<double>(hijd[3][id]) / a
            + srcAmp * S_xy
            );
            hijd[3][id] += static_cast<float>(rhs);
        }
        {
            const double lap_h = static_cast<double>(lapl<float>(i, j, k, hij[4]));
            const double rhs = d * (
            laplnorm * lap_h
            - (2.0 + rescale_s) * ad * static_cast<double>(hijd[4][id]) / a
            + srcAmp * S_xz
            );
            hijd[4][id] += static_cast<float>(rhs);
        }
        {
            const double lap_h = static_cast<double>(lapl<float>(i, j, k, hij[5]));
            const double rhs = d * (
            laplnorm * lap_h
            - (2.0 + rescale_s) * ad * static_cast<double>(hijd[5][id]) / a
            + srcAmp * S_yz
            );
            hijd[5][id] += static_cast<float>(rhs);
        }
    }
#endif

    ad += 0.5 * d * ad2;
}

// -------------------- Main Post-Inflation Evolution Loop --------------------

void run_post_inflation_loop(FILE* output_) {
    // Main post-inflation driver:
    // reuses inflation RK machinery while switching RHS via ScopedPostInflationRhsMode.
    ensure_periodic_index_cache();

    initialize_post_inflation();

    int numsteps = 0;

    InflationRKScratch rk_scratch;
    rk_scratch.ensure_size(f.size());

    const auto report_step = [&](int step_index) {
        if (step_index % output_freq == 0) {
            save_post_inflation((step_index % output_infrequent_freq == 0) ? 1 : 0);
        }

        if (screen_updates && step_index % output_freq == 0) {
            printf("scale factor a = %f\n", a);
            printf("numsteps %i\n\n", step_index);
        }

        fprintf(output_, "scale factor a = %f\n", a);
        fprintf(output_, "numsteps %i\n\n", step_index);
        if (step_index % output_freq == 0) {
            fflush(output_);
        }
    };

    switch (post_inflation_integrator) {
        case INTEGRATOR_LEAPFROG: {
            evolve_fields(0.5 * dt_post_inflation);

            while (a <= af_post_inflation) {
                evolve_derivs_post_inflation(dt_post_inflation);
                evolve_fields(dt_post_inflation);

                numsteps++;
                report_step(numsteps);
            }
            break;
        }
        case INTEGRATOR_RK4: {
            const ScopedPostInflationRhsMode post_rhs_mode(true);
            while (a <= af_post_inflation) {
                rk4_step_inflation(dt_post_inflation, rk_scratch);
                numsteps++;
                report_step(numsteps);
            }
            break;
        }
        case INTEGRATOR_RK45:
        default: {
            const ScopedPostInflationRhsMode post_rhs_mode(true);
            double h = clamp_rk45_step(dt_post_inflation, rk45_max_dt);

            while (a <= af_post_inflation) {
                const double hmax = rk45_hmax_from_base_step(dt_post_inflation);
                h = clamp_rk45_step(h, hmax);

                h = rk45_accept_step(
                    h,
                    hmax,
                    [&](double h_trial, double h_limit, double& h_next) {
                        return rk45_step_inflation(h_trial, h_limit, h_next, rk_scratch);
                    },
                    [&](int attempts, double failed_h) {
                        std::fprintf(stderr,
                            "RK45 failed to converge in post-inflation at a=%e, t=%e (attempts=%d, h=%e)\n",
                            a, t, attempts, failed_h);
                        std::exit(1);
                    }
                );

                numsteps++;
                report_step(numsteps);
            }
            break;
        }
    }

    printf("Saving final inflaton data\n");
    fflush(output_);
    save_post_inflation(1);
    // Note: save_post_inflation() performs a temporary sync/desync for output.
    // Leave the integrator state unchanged here.
}

#endif
