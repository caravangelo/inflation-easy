// evolution_internal.h - Conventions and private interfaces of src/evolution/
//
// Start here when reading the evolution code. The run has three stages:
//   1. inflation       run_inflation_loop()       evolution.cpp, inflation.cpp
//   2. deltaN          run_deltaN_loop()          deltaN.cpp       (optional)
//   3. post-inflation  run_post_inflation_loop()  post_inflation.cpp (optional)
// integrators.cpp holds the leapfrog drift and the RK4/RK45 steppers shared by
// stages 1 and 3; linear_metric.cpp holds the optional metric correction of
// stage 1. initializeN() and initialize_post_inflation() (initialize.cpp)
// prepare the data for stages 2 and 3.
//
// Units and variables (notation of the paper). Reduced Planck units, M_Pl = 1,
// so the field is phi/M_Pl. As in LATTICEEASY, the evolution uses code units,
// defined by the rescaled variables
//     x~ = B x,    d tau~ = B a^s d tau,    V~ = V / B^2,
// where tau is conformal time, B = rescale_B and s = rescale_s (default 0).
// L and dx are in x~ units. In comments a prime denotes d/d tau~, and
// "Laplacian" means d^2/dx~^2 as approximated by the selected stencil.
// Outputs are converted back to reduced Planck units.
//
// Time steps. The run-time dt is a step in B*t (t: cosmic time). Because
// d tau~ = B a^(s-1) dt, fixed-step inflation uses dt*a^(s-1) in tau~,
// evaluating a at the start of each step. RK45 uses this as a reference for
// adaptive steps, not a fixed cosmic-time interval. dt_post_inflation is a
// step in tau~; the deltaN stage steps in local e-folds (dN).
//
// Global state (declared in main.h):
//   t          tau~ (in the deltaN stage it is only a running counter)
//   a, ad      scale factor and a';  ad2 = a'' (leapfrog bookkeeping)
//   astep      a at the end of the previous step; sets the next step size
//   f, fd      inflation:      phi and phi'
//              deltaN:         phi and dphi/dN_loc (N_loc: local e-folds)
//              post-inflation: Newtonian potential Phi and Phi'
//   hij, hijd  the six components of h_ij and h_ij', defined by
//              ds^2 = a^2 [-d tau^2 + (delta_ij + h_ij) dx^i dx^j]
//              and packed as [xx, yy, zz, xy, xz, yz] (sym_idx() in main.h).
//              All six are evolved without projection; output.cpp applies
//              the transverse-traceless (TT) projection to the spectra.
//   deltaN     accumulated local e-folds N(n) per site (zeta after saveN())
//   Ne, dN, Nend, phiref
//              deltaN loop: e-folds integrated so far, signed step, maximum
//              |Ne|, and field value of the final constant-phi slice
//   omega      post-inflationary equation of state w
//
// Time levels. RK4/RK45 keep fields and velocities at the same time. The
// leapfrog stores fields half a step ahead of velocities: a "kick" advances
// velocities using the fields, a "drift" advances fields at fixed velocity.
// save() and the end of run_inflation_loop() drift the fields back by half a
// step whenever synchronized values are needed.

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

// Midpoint kick for v' = force - friction*v, with force and friction
// evaluated halfway through the velocity step.
inline double damped_kick(double velocity, double force, double step, double half_drag) {
    return (velocity * (1.0 - half_drag) + step * force) / (1.0 + half_drag);
}

#if calculate_SIGW || post_inflation
// Periodic derivative along axis 0/1/2 = x/y/z, INCLUDING 1/dx. In contrast,
// spatial::laplacian() omits 1/dx^2, which its callers must supply.
// The order-2 branch preserves the default hot path; higher orders use the
// selected spatial stencil.
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

// Accept/reject loop shared by all RK45 drivers. step_function(h, hmax, h_next)
// attempts one step and returns true if it was accepted; after a rejection the
// evolved state is unchanged and h_next is a smaller trial step.
// failure_handler is called when the step cannot shrink further or too many
// attempts failed; the drivers abort the run there. Returns the accepted step
// and the proposed next step, both signed (backward deltaN runs use h < 0).
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

// Butcher tables of the Dormand-Prince 5(4) pair and of classical RK4, shared
// by the field and deltaN integrators. DP_B are the 5th-order weights used to
// advance the solution; DP_E = (5th - 4th order weights) gives the error
// estimate.
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

// RK45 step-size limits (see integrators.cpp).
double clamp_rk45_step(double h, double hmax);
double rk45_hmax_from_base_step(double base_step);

// Inflationary equations (inflation.cpp): one leapfrog kick, and the RHS
// evaluated on an RK stage state.
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
// Post-inflationary RHS evaluated on an RK stage state (post_inflation.cpp).
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

// Shared RK4 / RK45 steppers for stages 1 and 3 (integrators.cpp). h is a
// step in tau~; phase selects which right-hand side is used.
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
