// evolution_internal.h - Private interfaces shared by the evolution stages
//
// This header is internal to the C++ implementation. It keeps the public
// declarations in main.h small while allowing the evolution stages to share
// derivative helpers, RK coefficients, and integration interfaces.

#pragma once

#include "main.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

// Keep selected aggregate helpers outside already-large OpenMP lattice kernels.
// Inlining them increases code size and register pressure in the tensor loops.
#if defined(_MSC_VER)
#define INFLATIONEASY_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define INFLATIONEASY_NOINLINE __attribute__((noinline))
#else
#define INFLATIONEASY_NOINLINE
#endif

namespace evolution {

#if calculate_SIGW || post_inflation
// First derivative in code units. The explicit order-2 branch preserves the
// original default hot path; higher orders use the selected spatial stencil.
template <typename T>
inline T first_spatial_derivative(
    int axis, int i, int j, int k, const std::vector<T>& field)
{
    if constexpr (spatial::order == 2) {
        const double half_over_dx = 0.5 / dx;
        if (axis == 0) {
            if (i == 0 || i == N - 1) {
                return T((field[idx(spatial::increment(i), j, k)]
                        - field[idx(spatial::decrement(i), j, k)]) * half_over_dx);
            }
            return T((field[idx(i + 1, j, k)] - field[idx(i - 1, j, k)]) * half_over_dx);
        }
        if (axis == 1) {
            if (j == 0 || j == N - 1) {
                return T((field[idx(i, spatial::increment(j), k)]
                        - field[idx(i, spatial::decrement(j), k)]) * half_over_dx);
            }
            return T((field[idx(i, j + 1, k)] - field[idx(i, j - 1, k)]) * half_over_dx);
        }
        if (k == 0 || k == N - 1) {
            return T((field[idx(i, j, spatial::increment(k))]
                    - field[idx(i, j, spatial::decrement(k))]) * half_over_dx);
        }
        return T((field[idx(i, j, k + 1)] - field[idx(i, j, k - 1)]) * half_over_dx);
    } else {
        return spatial::first_derivative(axis, i, j, k, field, dx);
    }
}

#endif

// Simulation phase whose equations are selected by the shared field integrators.
enum class SimulationPhase {
    Inflation,
    PostInflation
};

// Reused stage/state buffers for the inflationary and post-inflationary RK
// methods. ensure_size() allocates only when the lattice size changes.
struct FieldRKWorkspace {
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

    void ensure_size(size_t site_count);
};

struct RK45StepResult {
    double accepted_step_size;
    double next_step_size;
};

constexpr int RK45_MAX_ATTEMPTS = 25;
constexpr double RK45_MIN_STEP_GUARD = 1.0 + 1e-12;

// Shared accept/reject driver. step_function leaves the state unchanged after
// rejection; failure_handler terminates or otherwise resolves a stalled step.
template <typename StepFunction, typename FailureHandler>
RK45StepResult rk45_accept_step(
    double h,
    double hmax,
    StepFunction&& step_function,
    FailureHandler&& failure_handler)
{
    bool accepted = false;
    int attempts = 0;
    double accepted_h = h;
    while (!accepted) {
        accepted_h = h;
        double h_suggested = h;
        accepted = step_function(h, hmax, h_suggested);
        h = h_suggested;
        ++attempts;

        if (!accepted &&
            (attempts > RK45_MAX_ATTEMPTS ||
             std::abs(h) <= rk45_min_dt * RK45_MIN_STEP_GUARD)) {
            failure_handler(attempts, h);
        }
    }
    return {accepted_h, h};
}

// Compile-time-visible Runge--Kutta coefficients. The field and deltaN
// integrators share these tables, while each retains its own state layout.
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

inline constexpr double RK4_A[4][4] = {
    {0.0, 0.0, 0.0, 0.0},
    {0.5, 0.0, 0.0, 0.0},
    {0.0, 0.5, 0.0, 0.0},
    {0.0, 0.0, 1.0, 0.0}
};

inline constexpr double RK4_B[4] = {
    1.0 / 6.0,
    1.0 / 3.0,
    1.0 / 3.0,
    1.0 / 6.0
};

inline constexpr double DP_A[7][7] = {
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {DormandPrince45Coefficients::a21, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {DormandPrince45Coefficients::a31,
     DormandPrince45Coefficients::a32,
     0.0, 0.0, 0.0, 0.0, 0.0},
    {DormandPrince45Coefficients::a41,
     DormandPrince45Coefficients::a42,
     DormandPrince45Coefficients::a43,
     0.0, 0.0, 0.0, 0.0},
    {DormandPrince45Coefficients::a51,
     DormandPrince45Coefficients::a52,
     DormandPrince45Coefficients::a53,
     DormandPrince45Coefficients::a54,
     0.0, 0.0, 0.0},
    {DormandPrince45Coefficients::a61,
     DormandPrince45Coefficients::a62,
     DormandPrince45Coefficients::a63,
     DormandPrince45Coefficients::a64,
     DormandPrince45Coefficients::a65,
     0.0, 0.0},
    {DormandPrince45Coefficients::a71,
     0.0,
     DormandPrince45Coefficients::a73,
     DormandPrince45Coefficients::a74,
     DormandPrince45Coefficients::a75,
     DormandPrince45Coefficients::a76,
     0.0}
};

inline constexpr double DP_B[7] = {
    DormandPrince45Coefficients::b1,
    0.0,
    DormandPrince45Coefficients::b3,
    DormandPrince45Coefficients::b4,
    DormandPrince45Coefficients::b5,
    DormandPrince45Coefficients::b6,
    0.0
};

inline constexpr double DP_E[7] = {
    DormandPrince45Coefficients::e1,
    0.0,
    DormandPrince45Coefficients::e3,
    DormandPrince45Coefficients::e4,
    DormandPrince45Coefficients::e5,
    DormandPrince45Coefficients::e6,
    DormandPrince45Coefficients::e7
};

double clamp_rk45_step(double h, double hmax);
double rk45_hmax_from_base_step(double base_step);

// Inflationary leapfrog kick and Runge--Kutta right-hand side.
void apply_inflation_leapfrog_kick(double step_size);

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
    double& daddt);

#if post_inflation
// Post-inflationary Runge--Kutta right-hand side.
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
    double& daddt);
#endif

// Shared fixed-step and adaptive field integrators.
void rk4_step_fields(
    double h,
    FieldRKWorkspace& workspace,
    SimulationPhase phase);

bool rk45_step_fields(
    double h,
    double hmax,
    double& h_next,
    FieldRKWorkspace& workspace,
    SimulationPhase phase);

} // namespace evolution
