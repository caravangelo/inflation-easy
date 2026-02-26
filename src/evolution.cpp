// evolution.cpp - Time evolution of fields and background
//
// This file implements the time integration of the scalar field,
// optional auxiliary systems (deltaN, tensor perturbations), and
// the background scale factor using the chosen update scheme.


#include "main.h"

// -------------------- Laplacians --------------------

// Helper for periodic indexing
inline int INCREMENT(int i) {
    return (i == N - 1) ? 0 : i + 1;
}

// Decrement index with periodic wrapping (i → i-1 mod N)
inline int DECREMENT(int i) {
    return (i == 0) ? N - 1 : i - 1;
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

bool g_use_post_inflation_rhs = false;

double clamp_rk45_step(double h, double hmax) {
    const double hmin = std::max(1e-16, rk45_min_dt);
    const double hhi = std::max(hmin, hmax);
    if (h < hmin) return hmin;
    if (h > hhi) return hhi;
    return h;
}

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

        auto d2_same_state = [&](int dim) {
            int ip = (dim==0? ((i==N-1)? (int)INCREMENT(i) : i+1) : i);
            int im = (dim==0? ((i==0)  ? (int)DECREMENT(i) : i-1) : i);
            int jp = (dim==1? ((j==N-1)? (int)INCREMENT(j) : j+1) : j);
            int jm = (dim==1? ((j==0)  ? (int)DECREMENT(j) : j-1) : j);
            int kp = (dim==2? ((k==N-1)? (int)INCREMENT(k) : k+1) : k);
            int km = (dim==2? ((k==0)  ? (int)DECREMENT(k) : k-1) : k);

            if (dim==0) return (f_state[idx(ip,j,k)] - 2.0*f_state[id] + f_state[idx(im,j,k)])/(dx*dx);
            if (dim==1) return (f_state[idx(i,jp,k)] - 2.0*f_state[id] + f_state[idx(i,jm,k)])/(dx*dx);
            return         (f_state[idx(i,j,kp)] - 2.0*f_state[id] + f_state[idx(i,j,km)])/(dx*dx);
        };

        auto d2_cross_state = [&](int d1, int d2) {
            int ip = ((d1==0||d2==0)? ((i==N-1)? (int)INCREMENT(i) : i+1) : i);
            int im = ((d1==0||d2==0)? ((i==0)  ? (int)DECREMENT(i) : i-1) : i);
            int jp = ((d1==1||d2==1)? ((j==N-1)? (int)INCREMENT(j) : j+1) : j);
            int jm = ((d1==1||d2==1)? ((j==0)  ? (int)DECREMENT(j) : j-1) : j);
            int kp = ((d1==2||d2==2)? ((k==N-1)? (int)INCREMENT(k) : k+1) : k);
            int km = ((d1==2||d2==2)? ((k==0)  ? (int)DECREMENT(k) : k-1) : k);

            const double fpp = f_state[idx(ip,jp,kp)];
            const double fpm = f_state[idx(ip,jm,km)];
            const double fmp = f_state[idx(im,jp,kp)];
            const double fmm = f_state[idx(im,jm,km)];
            return (fpp - fpm - fmp + fmm) / (4.0*dx*dx);
        };

        const double Hxx = d2_same_state(0);
        const double Hyy = d2_same_state(1);
        const double Hzz = d2_same_state(2);
        const double Hxy = d2_cross_state(0,1);
        const double Hxz = d2_cross_state(0,2);
        const double Hyz = d2_cross_state(1,2);

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

void rk4_step_inflation(double h, InflationRKScratch& scratch) {
    const size_t gs = f.size();
    scratch.ensure_size(gs);

    double ka[4], kad[4];

#if calculate_SIGW
    compute_inflation_rhs(f, fd, hij, hijd, scratch.kf[0], scratch.kfd[0], scratch.kh[0], scratch.khd[0], a, ad, ka[0], kad[0]);
#else
    compute_inflation_rhs(f, fd, scratch.kf[0], scratch.kfd[0], a, ad, ka[0], kad[0]);
#endif

    for (size_t id = 0; id < gs; ++id) {
        scratch.ftmp[id] = f[id] + 0.5 * h * scratch.kf[0][id];
        scratch.fdtmp[id] = fd[id] + 0.5 * h * scratch.kfd[0][id];
    }
#if calculate_SIGW
    for (int c = 0; c < 6; ++c) {
        for (size_t id = 0; id < gs; ++id) {
            scratch.htmp[c][id] = static_cast<float>(hij[c][id] + 0.5 * h * static_cast<double>(scratch.kh[0][c][id]));
            scratch.hdtmp[c][id] = static_cast<float>(hijd[c][id] + 0.5 * h * static_cast<double>(scratch.khd[0][c][id]));
        }
    }
    compute_inflation_rhs(scratch.ftmp, scratch.fdtmp, scratch.htmp, scratch.hdtmp, scratch.kf[1], scratch.kfd[1], scratch.kh[1], scratch.khd[1],
        a + 0.5 * h * ka[0], ad + 0.5 * h * kad[0], ka[1], kad[1]);
#else
    compute_inflation_rhs(scratch.ftmp, scratch.fdtmp, scratch.kf[1], scratch.kfd[1],
        a + 0.5 * h * ka[0], ad + 0.5 * h * kad[0], ka[1], kad[1]);
#endif

    for (size_t id = 0; id < gs; ++id) {
        scratch.ftmp[id] = f[id] + 0.5 * h * scratch.kf[1][id];
        scratch.fdtmp[id] = fd[id] + 0.5 * h * scratch.kfd[1][id];
    }
#if calculate_SIGW
    for (int c = 0; c < 6; ++c) {
        for (size_t id = 0; id < gs; ++id) {
            scratch.htmp[c][id] = static_cast<float>(hij[c][id] + 0.5 * h * static_cast<double>(scratch.kh[1][c][id]));
            scratch.hdtmp[c][id] = static_cast<float>(hijd[c][id] + 0.5 * h * static_cast<double>(scratch.khd[1][c][id]));
        }
    }
    compute_inflation_rhs(scratch.ftmp, scratch.fdtmp, scratch.htmp, scratch.hdtmp, scratch.kf[2], scratch.kfd[2], scratch.kh[2], scratch.khd[2],
        a + 0.5 * h * ka[1], ad + 0.5 * h * kad[1], ka[2], kad[2]);
#else
    compute_inflation_rhs(scratch.ftmp, scratch.fdtmp, scratch.kf[2], scratch.kfd[2],
        a + 0.5 * h * ka[1], ad + 0.5 * h * kad[1], ka[2], kad[2]);
#endif

    for (size_t id = 0; id < gs; ++id) {
        scratch.ftmp[id] = f[id] + h * scratch.kf[2][id];
        scratch.fdtmp[id] = fd[id] + h * scratch.kfd[2][id];
    }
#if calculate_SIGW
    for (int c = 0; c < 6; ++c) {
        for (size_t id = 0; id < gs; ++id) {
            scratch.htmp[c][id] = static_cast<float>(hij[c][id] + h * static_cast<double>(scratch.kh[2][c][id]));
            scratch.hdtmp[c][id] = static_cast<float>(hijd[c][id] + h * static_cast<double>(scratch.khd[2][c][id]));
        }
    }
    compute_inflation_rhs(scratch.ftmp, scratch.fdtmp, scratch.htmp, scratch.hdtmp, scratch.kf[3], scratch.kfd[3], scratch.kh[3], scratch.khd[3],
        a + h * ka[2], ad + h * kad[2], ka[3], kad[3]);
#else
    compute_inflation_rhs(scratch.ftmp, scratch.fdtmp, scratch.kf[3], scratch.kfd[3],
        a + h * ka[2], ad + h * kad[2], ka[3], kad[3]);
#endif

    const double sixth = h / 6.0;
    for (size_t id = 0; id < gs; ++id) {
        f[id] += sixth * (scratch.kf[0][id] + 2.0 * scratch.kf[1][id] + 2.0 * scratch.kf[2][id] + scratch.kf[3][id]);
        fd[id] += sixth * (scratch.kfd[0][id] + 2.0 * scratch.kfd[1][id] + 2.0 * scratch.kfd[2][id] + scratch.kfd[3][id]);
    }
#if calculate_SIGW
    for (int c = 0; c < 6; ++c) {
        for (size_t id = 0; id < gs; ++id) {
            hij[c][id] = static_cast<float>(hij[c][id] + sixth * (
                static_cast<double>(scratch.kh[0][c][id]) + 2.0 * static_cast<double>(scratch.kh[1][c][id]) +
                2.0 * static_cast<double>(scratch.kh[2][c][id]) + static_cast<double>(scratch.kh[3][c][id])));
            hijd[c][id] = static_cast<float>(hijd[c][id] + sixth * (
                static_cast<double>(scratch.khd[0][c][id]) + 2.0 * static_cast<double>(scratch.khd[1][c][id]) +
                2.0 * static_cast<double>(scratch.khd[2][c][id]) + static_cast<double>(scratch.khd[3][c][id])));
        }
    }
#endif

    a += sixth * (ka[0] + 2.0 * ka[1] + 2.0 * ka[2] + ka[3]);
    ad += sixth * (kad[0] + 2.0 * kad[1] + 2.0 * kad[2] + kad[3]);
    t += h;
}

bool rk45_step_inflation(double h, double hmax, double& h_next, InflationRKScratch& scratch) {
    const size_t gs = f.size();
    scratch.ensure_size(gs);

    static constexpr double a21 = 1.0 / 5.0;
    static constexpr double a31 = 3.0 / 40.0, a32 = 9.0 / 40.0;
    static constexpr double a41 = 44.0 / 45.0, a42 = -56.0 / 15.0, a43 = 32.0 / 9.0;
    static constexpr double a51 = 19372.0 / 6561.0, a52 = -25360.0 / 2187.0, a53 = 64448.0 / 6561.0, a54 = -212.0 / 729.0;
    static constexpr double a61 = 9017.0 / 3168.0, a62 = -355.0 / 33.0, a63 = 46732.0 / 5247.0, a64 = 49.0 / 176.0, a65 = -5103.0 / 18656.0;
    static constexpr double a71 = 35.0 / 384.0, a73 = 500.0 / 1113.0, a74 = 125.0 / 192.0, a75 = -2187.0 / 6784.0, a76 = 11.0 / 84.0;

    static constexpr double b1 = 35.0 / 384.0, b3 = 500.0 / 1113.0, b4 = 125.0 / 192.0, b5 = -2187.0 / 6784.0, b6 = 11.0 / 84.0;
    static constexpr double bs1 = 5179.0 / 57600.0, bs3 = 7571.0 / 16695.0, bs4 = 393.0 / 640.0, bs5 = -92097.0 / 339200.0, bs6 = 187.0 / 2100.0, bs7 = 1.0 / 40.0;

    static constexpr double e1 = b1 - bs1;
    static constexpr double e3 = b3 - bs3;
    static constexpr double e4 = b4 - bs4;
    static constexpr double e5 = b5 - bs5;
    static constexpr double e6 = b6 - bs6;
    static constexpr double e7 = -bs7;

    double ka[7], kad[7];

#if calculate_SIGW
    compute_inflation_rhs(f, fd, hij, hijd, scratch.kf[0], scratch.kfd[0], scratch.kh[0], scratch.khd[0], a, ad, ka[0], kad[0]);
#else
    compute_inflation_rhs(f, fd, scratch.kf[0], scratch.kfd[0], a, ad, ka[0], kad[0]);
#endif

    for (size_t id = 0; id < gs; ++id) {
        scratch.ftmp[id] = f[id] + h * (a21 * scratch.kf[0][id]);
        scratch.fdtmp[id] = fd[id] + h * (a21 * scratch.kfd[0][id]);
    }
#if calculate_SIGW
    for (int c = 0; c < 6; ++c) for (size_t id = 0; id < gs; ++id) {
        scratch.htmp[c][id] = static_cast<float>(hij[c][id] + h * (a21 * static_cast<double>(scratch.kh[0][c][id])));
        scratch.hdtmp[c][id] = static_cast<float>(hijd[c][id] + h * (a21 * static_cast<double>(scratch.khd[0][c][id])));
    }
    compute_inflation_rhs(scratch.ftmp, scratch.fdtmp, scratch.htmp, scratch.hdtmp, scratch.kf[1], scratch.kfd[1], scratch.kh[1], scratch.khd[1],
        a + h * (a21 * ka[0]), ad + h * (a21 * kad[0]), ka[1], kad[1]);
#else
    compute_inflation_rhs(scratch.ftmp, scratch.fdtmp, scratch.kf[1], scratch.kfd[1],
        a + h * (a21 * ka[0]), ad + h * (a21 * kad[0]), ka[1], kad[1]);
#endif

    for (size_t id = 0; id < gs; ++id) {
        scratch.ftmp[id] = f[id] + h * (a31 * scratch.kf[0][id] + a32 * scratch.kf[1][id]);
        scratch.fdtmp[id] = fd[id] + h * (a31 * scratch.kfd[0][id] + a32 * scratch.kfd[1][id]);
    }
#if calculate_SIGW
    for (int c = 0; c < 6; ++c) for (size_t id = 0; id < gs; ++id) {
        scratch.htmp[c][id] = static_cast<float>(hij[c][id] + h * (
            a31 * static_cast<double>(scratch.kh[0][c][id]) + a32 * static_cast<double>(scratch.kh[1][c][id])));
        scratch.hdtmp[c][id] = static_cast<float>(hijd[c][id] + h * (
            a31 * static_cast<double>(scratch.khd[0][c][id]) + a32 * static_cast<double>(scratch.khd[1][c][id])));
    }
    compute_inflation_rhs(scratch.ftmp, scratch.fdtmp, scratch.htmp, scratch.hdtmp, scratch.kf[2], scratch.kfd[2], scratch.kh[2], scratch.khd[2],
        a + h * (a31 * ka[0] + a32 * ka[1]), ad + h * (a31 * kad[0] + a32 * kad[1]), ka[2], kad[2]);
#else
    compute_inflation_rhs(scratch.ftmp, scratch.fdtmp, scratch.kf[2], scratch.kfd[2],
        a + h * (a31 * ka[0] + a32 * ka[1]), ad + h * (a31 * kad[0] + a32 * kad[1]), ka[2], kad[2]);
#endif

    for (size_t id = 0; id < gs; ++id) {
        scratch.ftmp[id] = f[id] + h * (a41 * scratch.kf[0][id] + a42 * scratch.kf[1][id] + a43 * scratch.kf[2][id]);
        scratch.fdtmp[id] = fd[id] + h * (a41 * scratch.kfd[0][id] + a42 * scratch.kfd[1][id] + a43 * scratch.kfd[2][id]);
    }
#if calculate_SIGW
    for (int c = 0; c < 6; ++c) for (size_t id = 0; id < gs; ++id) {
        scratch.htmp[c][id] = static_cast<float>(hij[c][id] + h * (
            a41 * static_cast<double>(scratch.kh[0][c][id]) +
            a42 * static_cast<double>(scratch.kh[1][c][id]) +
            a43 * static_cast<double>(scratch.kh[2][c][id])));
        scratch.hdtmp[c][id] = static_cast<float>(hijd[c][id] + h * (
            a41 * static_cast<double>(scratch.khd[0][c][id]) +
            a42 * static_cast<double>(scratch.khd[1][c][id]) +
            a43 * static_cast<double>(scratch.khd[2][c][id])));
    }
    compute_inflation_rhs(scratch.ftmp, scratch.fdtmp, scratch.htmp, scratch.hdtmp, scratch.kf[3], scratch.kfd[3], scratch.kh[3], scratch.khd[3],
        a + h * (a41 * ka[0] + a42 * ka[1] + a43 * ka[2]),
        ad + h * (a41 * kad[0] + a42 * kad[1] + a43 * kad[2]), ka[3], kad[3]);
#else
    compute_inflation_rhs(scratch.ftmp, scratch.fdtmp, scratch.kf[3], scratch.kfd[3],
        a + h * (a41 * ka[0] + a42 * ka[1] + a43 * ka[2]),
        ad + h * (a41 * kad[0] + a42 * kad[1] + a43 * kad[2]), ka[3], kad[3]);
#endif

    for (size_t id = 0; id < gs; ++id) {
        scratch.ftmp[id] = f[id] + h * (
            a51 * scratch.kf[0][id] + a52 * scratch.kf[1][id] + a53 * scratch.kf[2][id] + a54 * scratch.kf[3][id]);
        scratch.fdtmp[id] = fd[id] + h * (
            a51 * scratch.kfd[0][id] + a52 * scratch.kfd[1][id] + a53 * scratch.kfd[2][id] + a54 * scratch.kfd[3][id]);
    }
#if calculate_SIGW
    for (int c = 0; c < 6; ++c) for (size_t id = 0; id < gs; ++id) {
        scratch.htmp[c][id] = static_cast<float>(hij[c][id] + h * (
            a51 * static_cast<double>(scratch.kh[0][c][id]) +
            a52 * static_cast<double>(scratch.kh[1][c][id]) +
            a53 * static_cast<double>(scratch.kh[2][c][id]) +
            a54 * static_cast<double>(scratch.kh[3][c][id])));
        scratch.hdtmp[c][id] = static_cast<float>(hijd[c][id] + h * (
            a51 * static_cast<double>(scratch.khd[0][c][id]) +
            a52 * static_cast<double>(scratch.khd[1][c][id]) +
            a53 * static_cast<double>(scratch.khd[2][c][id]) +
            a54 * static_cast<double>(scratch.khd[3][c][id])));
    }
    compute_inflation_rhs(scratch.ftmp, scratch.fdtmp, scratch.htmp, scratch.hdtmp, scratch.kf[4], scratch.kfd[4], scratch.kh[4], scratch.khd[4],
        a + h * (a51 * ka[0] + a52 * ka[1] + a53 * ka[2] + a54 * ka[3]),
        ad + h * (a51 * kad[0] + a52 * kad[1] + a53 * kad[2] + a54 * kad[3]), ka[4], kad[4]);
#else
    compute_inflation_rhs(scratch.ftmp, scratch.fdtmp, scratch.kf[4], scratch.kfd[4],
        a + h * (a51 * ka[0] + a52 * ka[1] + a53 * ka[2] + a54 * ka[3]),
        ad + h * (a51 * kad[0] + a52 * kad[1] + a53 * kad[2] + a54 * kad[3]), ka[4], kad[4]);
#endif

    for (size_t id = 0; id < gs; ++id) {
        scratch.ftmp[id] = f[id] + h * (
            a61 * scratch.kf[0][id] + a62 * scratch.kf[1][id] + a63 * scratch.kf[2][id] +
            a64 * scratch.kf[3][id] + a65 * scratch.kf[4][id]);
        scratch.fdtmp[id] = fd[id] + h * (
            a61 * scratch.kfd[0][id] + a62 * scratch.kfd[1][id] + a63 * scratch.kfd[2][id] +
            a64 * scratch.kfd[3][id] + a65 * scratch.kfd[4][id]);
    }
#if calculate_SIGW
    for (int c = 0; c < 6; ++c) for (size_t id = 0; id < gs; ++id) {
        scratch.htmp[c][id] = static_cast<float>(hij[c][id] + h * (
            a61 * static_cast<double>(scratch.kh[0][c][id]) +
            a62 * static_cast<double>(scratch.kh[1][c][id]) +
            a63 * static_cast<double>(scratch.kh[2][c][id]) +
            a64 * static_cast<double>(scratch.kh[3][c][id]) +
            a65 * static_cast<double>(scratch.kh[4][c][id])));
        scratch.hdtmp[c][id] = static_cast<float>(hijd[c][id] + h * (
            a61 * static_cast<double>(scratch.khd[0][c][id]) +
            a62 * static_cast<double>(scratch.khd[1][c][id]) +
            a63 * static_cast<double>(scratch.khd[2][c][id]) +
            a64 * static_cast<double>(scratch.khd[3][c][id]) +
            a65 * static_cast<double>(scratch.khd[4][c][id])));
    }
    compute_inflation_rhs(scratch.ftmp, scratch.fdtmp, scratch.htmp, scratch.hdtmp, scratch.kf[5], scratch.kfd[5], scratch.kh[5], scratch.khd[5],
        a + h * (a61 * ka[0] + a62 * ka[1] + a63 * ka[2] + a64 * ka[3] + a65 * ka[4]),
        ad + h * (a61 * kad[0] + a62 * kad[1] + a63 * kad[2] + a64 * kad[3] + a65 * kad[4]), ka[5], kad[5]);
#else
    compute_inflation_rhs(scratch.ftmp, scratch.fdtmp, scratch.kf[5], scratch.kfd[5],
        a + h * (a61 * ka[0] + a62 * ka[1] + a63 * ka[2] + a64 * ka[3] + a65 * ka[4]),
        ad + h * (a61 * kad[0] + a62 * kad[1] + a63 * kad[2] + a64 * kad[3] + a65 * kad[4]), ka[5], kad[5]);
#endif

    for (size_t id = 0; id < gs; ++id) {
        scratch.ftmp[id] = f[id] + h * (
            a71 * scratch.kf[0][id] + a73 * scratch.kf[2][id] + a74 * scratch.kf[3][id] +
            a75 * scratch.kf[4][id] + a76 * scratch.kf[5][id]);
        scratch.fdtmp[id] = fd[id] + h * (
            a71 * scratch.kfd[0][id] + a73 * scratch.kfd[2][id] + a74 * scratch.kfd[3][id] +
            a75 * scratch.kfd[4][id] + a76 * scratch.kfd[5][id]);
    }
#if calculate_SIGW
    for (int c = 0; c < 6; ++c) for (size_t id = 0; id < gs; ++id) {
        scratch.htmp[c][id] = static_cast<float>(hij[c][id] + h * (
            a71 * static_cast<double>(scratch.kh[0][c][id]) +
            a73 * static_cast<double>(scratch.kh[2][c][id]) +
            a74 * static_cast<double>(scratch.kh[3][c][id]) +
            a75 * static_cast<double>(scratch.kh[4][c][id]) +
            a76 * static_cast<double>(scratch.kh[5][c][id])));
        scratch.hdtmp[c][id] = static_cast<float>(hijd[c][id] + h * (
            a71 * static_cast<double>(scratch.khd[0][c][id]) +
            a73 * static_cast<double>(scratch.khd[2][c][id]) +
            a74 * static_cast<double>(scratch.khd[3][c][id]) +
            a75 * static_cast<double>(scratch.khd[4][c][id]) +
            a76 * static_cast<double>(scratch.khd[5][c][id])));
    }
    compute_inflation_rhs(scratch.ftmp, scratch.fdtmp, scratch.htmp, scratch.hdtmp, scratch.kf[6], scratch.kfd[6], scratch.kh[6], scratch.khd[6],
        a + h * (a71 * ka[0] + a73 * ka[2] + a74 * ka[3] + a75 * ka[4] + a76 * ka[5]),
        ad + h * (a71 * kad[0] + a73 * kad[2] + a74 * kad[3] + a75 * kad[4] + a76 * kad[5]), ka[6], kad[6]);
#else
    compute_inflation_rhs(scratch.ftmp, scratch.fdtmp, scratch.kf[6], scratch.kfd[6],
        a + h * (a71 * ka[0] + a73 * ka[2] + a74 * ka[3] + a75 * ka[4] + a76 * ka[5]),
        ad + h * (a71 * kad[0] + a73 * kad[2] + a74 * kad[3] + a75 * kad[4] + a76 * kad[5]), ka[6], kad[6]);
#endif

    const double a5 = a + h * (b1 * ka[0] + b3 * ka[2] + b4 * ka[3] + b5 * ka[4] + b6 * ka[5]);
    const double ad5 = ad + h * (b1 * kad[0] + b3 * kad[2] + b4 * kad[3] + b5 * kad[4] + b6 * kad[5]);
    const double erra = h * (e1 * ka[0] + e3 * ka[2] + e4 * ka[3] + e5 * ka[4] + e6 * ka[5] + e7 * ka[6]);
    const double errad = h * (e1 * kad[0] + e3 * kad[2] + e4 * kad[3] + e5 * kad[4] + e6 * kad[5] + e7 * kad[6]);

    double err_acc = 0.0;
    std::size_t nvars = 2 * gs + 2;

    for (size_t id = 0; id < gs; ++id) {
        const double y5f = f[id] + h * (
            b1 * scratch.kf[0][id] + b3 * scratch.kf[2][id] + b4 * scratch.kf[3][id] +
            b5 * scratch.kf[4][id] + b6 * scratch.kf[5][id]);
        const double y5fd = fd[id] + h * (
            b1 * scratch.kfd[0][id] + b3 * scratch.kfd[2][id] + b4 * scratch.kfd[3][id] +
            b5 * scratch.kfd[4][id] + b6 * scratch.kfd[5][id]);
        const double errf = h * (
            e1 * scratch.kf[0][id] + e3 * scratch.kf[2][id] + e4 * scratch.kf[3][id] +
            e5 * scratch.kf[4][id] + e6 * scratch.kf[5][id] + e7 * scratch.kf[6][id]);
        const double errfd = h * (
            e1 * scratch.kfd[0][id] + e3 * scratch.kfd[2][id] + e4 * scratch.kfd[3][id] +
            e5 * scratch.kfd[4][id] + e6 * scratch.kfd[5][id] + e7 * scratch.kfd[6][id]);

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
            const double y5h = static_cast<double>(hij[c][id]) + h * (
                b1 * static_cast<double>(scratch.kh[0][c][id]) +
                b3 * static_cast<double>(scratch.kh[2][c][id]) +
                b4 * static_cast<double>(scratch.kh[3][c][id]) +
                b5 * static_cast<double>(scratch.kh[4][c][id]) +
                b6 * static_cast<double>(scratch.kh[5][c][id]));
            const double y5hd = static_cast<double>(hijd[c][id]) + h * (
                b1 * static_cast<double>(scratch.khd[0][c][id]) +
                b3 * static_cast<double>(scratch.khd[2][c][id]) +
                b4 * static_cast<double>(scratch.khd[3][c][id]) +
                b5 * static_cast<double>(scratch.khd[4][c][id]) +
                b6 * static_cast<double>(scratch.khd[5][c][id]));
            const double errh = h * (
                e1 * static_cast<double>(scratch.kh[0][c][id]) +
                e3 * static_cast<double>(scratch.kh[2][c][id]) +
                e4 * static_cast<double>(scratch.kh[3][c][id]) +
                e5 * static_cast<double>(scratch.kh[4][c][id]) +
                e6 * static_cast<double>(scratch.kh[5][c][id]) +
                e7 * static_cast<double>(scratch.kh[6][c][id]));
            const double errhd = h * (
                e1 * static_cast<double>(scratch.khd[0][c][id]) +
                e3 * static_cast<double>(scratch.khd[2][c][id]) +
                e4 * static_cast<double>(scratch.khd[3][c][id]) +
                e5 * static_cast<double>(scratch.khd[4][c][id]) +
                e6 * static_cast<double>(scratch.khd[5][c][id]) +
                e7 * static_cast<double>(scratch.khd[6][c][id]));

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
        fflush(output_);
    };

    if (integrator == INTEGRATOR_LEAPFROG) {
        evolve_fields(0.5 * dt); // First leapfrog step

        while (a <= af) {
            double dt_rescaled = dt * std::pow(astep, rescale_s - 1.0);
            evolve_derivs(dt_rescaled);
            evolve_fields(dt_rescaled);

            numsteps++;
            report_step(numsteps);
            astep = a;
        }
    } else if (integrator == INTEGRATOR_RK4) {
        while (a <= af) {
            const double dt_rescaled = dt * std::pow(astep, rescale_s - 1.0);
            rk4_step_inflation(dt_rescaled, rk_scratch);

            numsteps++;
            report_step(numsteps);
            astep = a;
        }
    } else {
        double h = clamp_rk45_step(dt * std::pow(astep, rescale_s - 1.0), rk45_max_dt);

        while (a <= af) {
            const double base_step = dt * std::pow(astep, rescale_s - 1.0);
            const double hmax = std::min(rk45_max_dt, std::max(rk45_min_dt, 2.0 * base_step));
            h = clamp_rk45_step(h, hmax);

            bool accepted = false;
            int attempts = 0;
            while (!accepted) {
                double h_suggested = h;
                accepted = rk45_step_inflation(h, hmax, h_suggested, rk_scratch);
                h = h_suggested;
                attempts++;

                if (!accepted && (attempts > 25 || h <= rk45_min_dt * (1.0 + 1e-12))) {
                    std::fprintf(stderr,
                        "RK45 failed to converge at a=%e, t=%e (attempts=%d, h=%e)\n",
                        a, t, attempts, h);
                    std::exit(1);
                }
            }

            numsteps++;
            report_step(numsteps);
            astep = a;
        }
    }

    printf("Saving final inflaton data\n");
    save(1);
    if (inflation_uses_staggered_derivatives()) {
        evolve_fields(-0.5 * dt * std::pow(astep, rescale_s - 1.0)); // sync field values with velocities
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

void rk4_step_deltaN(double h, DeltaNRKScratch& scratch) {
    const size_t gs = f.size();
    scratch.ensure_size(gs);

    double ka[4];
    compute_deltaN_rhs(f, fd, scratch.kf[0], scratch.kfd[0], scratch.kdn[0], ka[0]);

    for (size_t id = 0; id < gs; ++id) {
        scratch.ftmp[id] = f[id] + 0.5 * h * scratch.kf[0][id];
        scratch.fdtmp[id] = fd[id] + 0.5 * h * scratch.kfd[0][id];
        scratch.dntmp[id] = deltaN[id] + 0.5 * h * scratch.kdn[0][id];
    }
    compute_deltaN_rhs(scratch.ftmp, scratch.fdtmp, scratch.kf[1], scratch.kfd[1], scratch.kdn[1], ka[1]);

    for (size_t id = 0; id < gs; ++id) {
        scratch.ftmp[id] = f[id] + 0.5 * h * scratch.kf[1][id];
        scratch.fdtmp[id] = fd[id] + 0.5 * h * scratch.kfd[1][id];
        scratch.dntmp[id] = deltaN[id] + 0.5 * h * scratch.kdn[1][id];
    }
    compute_deltaN_rhs(scratch.ftmp, scratch.fdtmp, scratch.kf[2], scratch.kfd[2], scratch.kdn[2], ka[2]);

    for (size_t id = 0; id < gs; ++id) {
        scratch.ftmp[id] = f[id] + h * scratch.kf[2][id];
        scratch.fdtmp[id] = fd[id] + h * scratch.kfd[2][id];
        scratch.dntmp[id] = deltaN[id] + h * scratch.kdn[2][id];
    }
    compute_deltaN_rhs(scratch.ftmp, scratch.fdtmp, scratch.kf[3], scratch.kfd[3], scratch.kdn[3], ka[3]);

    const double sixth = h / 6.0;
    for (size_t id = 0; id < gs; ++id) {
        f[id] += sixth * (scratch.kf[0][id] + 2.0 * scratch.kf[1][id] + 2.0 * scratch.kf[2][id] + scratch.kf[3][id]);
        fd[id] += sixth * (scratch.kfd[0][id] + 2.0 * scratch.kfd[1][id] + 2.0 * scratch.kfd[2][id] + scratch.kfd[3][id]);
        deltaN[id] += sixth * (scratch.kdn[0][id] + 2.0 * scratch.kdn[1][id] + 2.0 * scratch.kdn[2][id] + scratch.kdn[3][id]);
    }

    a += sixth * (ka[0] + 2.0 * ka[1] + 2.0 * ka[2] + ka[3]);
    t += h;
}

bool rk45_step_deltaN(double h, double hmax, double& h_next, DeltaNRKScratch& scratch) {
    const size_t gs = f.size();
    scratch.ensure_size(gs);

    static constexpr double a21 = 1.0 / 5.0;
    static constexpr double a31 = 3.0 / 40.0, a32 = 9.0 / 40.0;
    static constexpr double a41 = 44.0 / 45.0, a42 = -56.0 / 15.0, a43 = 32.0 / 9.0;
    static constexpr double a51 = 19372.0 / 6561.0, a52 = -25360.0 / 2187.0, a53 = 64448.0 / 6561.0, a54 = -212.0 / 729.0;
    static constexpr double a61 = 9017.0 / 3168.0, a62 = -355.0 / 33.0, a63 = 46732.0 / 5247.0, a64 = 49.0 / 176.0, a65 = -5103.0 / 18656.0;
    static constexpr double a71 = 35.0 / 384.0, a73 = 500.0 / 1113.0, a74 = 125.0 / 192.0, a75 = -2187.0 / 6784.0, a76 = 11.0 / 84.0;

    static constexpr double b1 = 35.0 / 384.0, b3 = 500.0 / 1113.0, b4 = 125.0 / 192.0, b5 = -2187.0 / 6784.0, b6 = 11.0 / 84.0;
    static constexpr double bs1 = 5179.0 / 57600.0, bs3 = 7571.0 / 16695.0, bs4 = 393.0 / 640.0, bs5 = -92097.0 / 339200.0, bs6 = 187.0 / 2100.0, bs7 = 1.0 / 40.0;

    static constexpr double e1 = b1 - bs1;
    static constexpr double e3 = b3 - bs3;
    static constexpr double e4 = b4 - bs4;
    static constexpr double e5 = b5 - bs5;
    static constexpr double e6 = b6 - bs6;
    static constexpr double e7 = -bs7;

    double ka[7];
    compute_deltaN_rhs(f, fd, scratch.kf[0], scratch.kfd[0], scratch.kdn[0], ka[0]);

    for (size_t id = 0; id < gs; ++id) {
        scratch.ftmp[id] = f[id] + h * (a21 * scratch.kf[0][id]);
        scratch.fdtmp[id] = fd[id] + h * (a21 * scratch.kfd[0][id]);
        scratch.dntmp[id] = deltaN[id] + h * (a21 * scratch.kdn[0][id]);
    }
    compute_deltaN_rhs(scratch.ftmp, scratch.fdtmp, scratch.kf[1], scratch.kfd[1], scratch.kdn[1], ka[1]);

    for (size_t id = 0; id < gs; ++id) {
        scratch.ftmp[id] = f[id] + h * (a31 * scratch.kf[0][id] + a32 * scratch.kf[1][id]);
        scratch.fdtmp[id] = fd[id] + h * (a31 * scratch.kfd[0][id] + a32 * scratch.kfd[1][id]);
        scratch.dntmp[id] = deltaN[id] + h * (a31 * scratch.kdn[0][id] + a32 * scratch.kdn[1][id]);
    }
    compute_deltaN_rhs(scratch.ftmp, scratch.fdtmp, scratch.kf[2], scratch.kfd[2], scratch.kdn[2], ka[2]);

    for (size_t id = 0; id < gs; ++id) {
        scratch.ftmp[id] = f[id] + h * (a41 * scratch.kf[0][id] + a42 * scratch.kf[1][id] + a43 * scratch.kf[2][id]);
        scratch.fdtmp[id] = fd[id] + h * (a41 * scratch.kfd[0][id] + a42 * scratch.kfd[1][id] + a43 * scratch.kfd[2][id]);
        scratch.dntmp[id] = deltaN[id] + h * (a41 * scratch.kdn[0][id] + a42 * scratch.kdn[1][id] + a43 * scratch.kdn[2][id]);
    }
    compute_deltaN_rhs(scratch.ftmp, scratch.fdtmp, scratch.kf[3], scratch.kfd[3], scratch.kdn[3], ka[3]);

    for (size_t id = 0; id < gs; ++id) {
        scratch.ftmp[id] = f[id] + h * (
            a51 * scratch.kf[0][id] + a52 * scratch.kf[1][id] + a53 * scratch.kf[2][id] + a54 * scratch.kf[3][id]);
        scratch.fdtmp[id] = fd[id] + h * (
            a51 * scratch.kfd[0][id] + a52 * scratch.kfd[1][id] + a53 * scratch.kfd[2][id] + a54 * scratch.kfd[3][id]);
        scratch.dntmp[id] = deltaN[id] + h * (
            a51 * scratch.kdn[0][id] + a52 * scratch.kdn[1][id] + a53 * scratch.kdn[2][id] + a54 * scratch.kdn[3][id]);
    }
    compute_deltaN_rhs(scratch.ftmp, scratch.fdtmp, scratch.kf[4], scratch.kfd[4], scratch.kdn[4], ka[4]);

    for (size_t id = 0; id < gs; ++id) {
        scratch.ftmp[id] = f[id] + h * (
            a61 * scratch.kf[0][id] + a62 * scratch.kf[1][id] + a63 * scratch.kf[2][id] +
            a64 * scratch.kf[3][id] + a65 * scratch.kf[4][id]);
        scratch.fdtmp[id] = fd[id] + h * (
            a61 * scratch.kfd[0][id] + a62 * scratch.kfd[1][id] + a63 * scratch.kfd[2][id] +
            a64 * scratch.kfd[3][id] + a65 * scratch.kfd[4][id]);
        scratch.dntmp[id] = deltaN[id] + h * (
            a61 * scratch.kdn[0][id] + a62 * scratch.kdn[1][id] + a63 * scratch.kdn[2][id] +
            a64 * scratch.kdn[3][id] + a65 * scratch.kdn[4][id]);
    }
    compute_deltaN_rhs(scratch.ftmp, scratch.fdtmp, scratch.kf[5], scratch.kfd[5], scratch.kdn[5], ka[5]);

    for (size_t id = 0; id < gs; ++id) {
        scratch.ftmp[id] = f[id] + h * (
            a71 * scratch.kf[0][id] + a73 * scratch.kf[2][id] + a74 * scratch.kf[3][id] +
            a75 * scratch.kf[4][id] + a76 * scratch.kf[5][id]);
        scratch.fdtmp[id] = fd[id] + h * (
            a71 * scratch.kfd[0][id] + a73 * scratch.kfd[2][id] + a74 * scratch.kfd[3][id] +
            a75 * scratch.kfd[4][id] + a76 * scratch.kfd[5][id]);
        scratch.dntmp[id] = deltaN[id] + h * (
            a71 * scratch.kdn[0][id] + a73 * scratch.kdn[2][id] + a74 * scratch.kdn[3][id] +
            a75 * scratch.kdn[4][id] + a76 * scratch.kdn[5][id]);
    }
    compute_deltaN_rhs(scratch.ftmp, scratch.fdtmp, scratch.kf[6], scratch.kfd[6], scratch.kdn[6], ka[6]);

    const double a5 = a + h * (b1 * ka[0] + b3 * ka[2] + b4 * ka[3] + b5 * ka[4] + b6 * ka[5]);
    const double erra = h * (e1 * ka[0] + e3 * ka[2] + e4 * ka[3] + e5 * ka[4] + e6 * ka[5] + e7 * ka[6]);

    double err_acc = 0.0;
    std::size_t nvars = 3 * gs + 1;

    for (size_t id = 0; id < gs; ++id) {
        const double y5f = f[id] + h * (
            b1 * scratch.kf[0][id] + b3 * scratch.kf[2][id] + b4 * scratch.kf[3][id] +
            b5 * scratch.kf[4][id] + b6 * scratch.kf[5][id]);
        const double y5fd = fd[id] + h * (
            b1 * scratch.kfd[0][id] + b3 * scratch.kfd[2][id] + b4 * scratch.kfd[3][id] +
            b5 * scratch.kfd[4][id] + b6 * scratch.kfd[5][id]);
        const double y5dn = deltaN[id] + h * (
            b1 * scratch.kdn[0][id] + b3 * scratch.kdn[2][id] + b4 * scratch.kdn[3][id] +
            b5 * scratch.kdn[4][id] + b6 * scratch.kdn[5][id]);

        const double errf = h * (
            e1 * scratch.kf[0][id] + e3 * scratch.kf[2][id] + e4 * scratch.kf[3][id] +
            e5 * scratch.kf[4][id] + e6 * scratch.kf[5][id] + e7 * scratch.kf[6][id]);
        const double errfd = h * (
            e1 * scratch.kfd[0][id] + e3 * scratch.kfd[2][id] + e4 * scratch.kfd[3][id] +
            e5 * scratch.kfd[4][id] + e6 * scratch.kfd[5][id] + e7 * scratch.kfd[6][id]);
        const double errdn = h * (
            e1 * scratch.kdn[0][id] + e3 * scratch.kdn[2][id] + e4 * scratch.kdn[3][id] +
            e5 * scratch.kdn[4][id] + e6 * scratch.kdn[5][id] + e7 * scratch.kdn[6][id]);

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
        fflush(output_);
    };

    if (deltaN_integrator == INTEGRATOR_LEAPFROG) {
        while (Ne <= Nend) {
            evolve_derivsN(dN);
            evolve_fieldsN(dN);
            Ne += dN;
            numsteps++;
            report_step(numsteps);
        }
    } else if (deltaN_integrator == INTEGRATOR_RK4) {
        while (Ne <= Nend) {
            rk4_step_deltaN(dN, rk_scratch);
            Ne += dN;
            numsteps++;
            report_step(numsteps);
        }
    } else {
        double h = clamp_rk45_step(dN, rk45_max_dt);

        while (Ne <= Nend) {
            const double hmax = std::min(rk45_max_dt, std::max(rk45_min_dt, 2.0 * dN));
            h = clamp_rk45_step(h, hmax);

            bool accepted = false;
            int attempts = 0;
            while (!accepted) {
                double h_suggested = h;
                accepted = rk45_step_deltaN(h, hmax, h_suggested, rk_scratch);
                h = h_suggested;
                attempts++;

                if (!accepted && (attempts > 25 || h <= rk45_min_dt * (1.0 + 1e-12))) {
                    std::fprintf(stderr,
                        "RK45 failed to converge in deltaN loop at N=%e, t=%e (attempts=%d, h=%e)\n",
                        Ne, t, attempts, h);
                    std::exit(1);
                }
            }

            Ne += h;
            numsteps++;
            report_step(numsteps);
        }
    }

    saveN();
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
    int ip = (dim==0? ((i==N-1)? (int)INCREMENT(i) : i+1) : i);
    int im = (dim==0? ((i==0)  ? (int)DECREMENT(i) : i-1) : i);
    int jp = (dim==1? ((j==N-1)? (int)INCREMENT(j) : j+1) : j);
    int jm = (dim==1? ((j==0)  ? (int)DECREMENT(j) : j-1) : j);
    int kp = (dim==2? ((k==N-1)? (int)INCREMENT(k) : k+1) : k);
    int km = (dim==2? ((k==0)  ? (int)DECREMENT(k) : k-1) : k);

    if (dim==0) return (A[idx(ip,j,k)] - 2.0*A[idx(i,j,k)] + A[idx(im,j,k)])/(dx*dx);
    if (dim==1) return (A[idx(i,jp,k)] - 2.0*A[idx(i,j,k)] + A[idx(i,jm,k)])/(dx*dx);
    return         (A[idx(i,j,kp)] - 2.0*A[idx(i,j,k)] + A[idx(i,j,km)])/(dx*dx);
}

inline double d2_cross(int d1, int d2, int i, int j, int k, const std::vector<double>& A) {
    int ip = ((d1==0||d2==0)? ((i==N-1)? (int)INCREMENT(i) : i+1) : i);
    int im = ((d1==0||d2==0)? ((i==0)  ? (int)DECREMENT(i) : i-1) : i);
    int jp = ((d1==1||d2==1)? ((j==N-1)? (int)INCREMENT(j) : j+1) : j);
    int jm = ((d1==1||d2==1)? ((j==0)  ? (int)DECREMENT(j) : j-1) : j);
    int kp = ((d1==2||d2==2)? ((k==N-1)? (int)INCREMENT(k) : k+1) : k);
    int km = ((d1==2||d2==2)? ((k==0)  ? (int)DECREMENT(k) : k-1) : k);

    double fpp = A[idx(ip,jp,kp)];
    double fpm = A[idx(ip,jm,km)];
    double fmp = A[idx(im,jp,kp)];
    double fmm = A[idx(im,jm,km)];
    return (fpp - fpm - fmp + fmm) / (4.0*dx*dx);
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
        fflush(output_);
    };

    if (post_inflation_integrator == INTEGRATOR_LEAPFROG) {
        evolve_fields(0.5 * dt_post_inflation);

        while (a <= af_post_inflation) {
            const double dt_rescaled = dt_post_inflation;
            evolve_derivs_post_inflation(dt_rescaled);
            evolve_fields(dt_rescaled);

            numsteps++;
            report_step(numsteps);
        }
    } else if (post_inflation_integrator == INTEGRATOR_RK4) {
        g_use_post_inflation_rhs = true;
        while (a <= af_post_inflation) {
            rk4_step_inflation(dt_post_inflation, rk_scratch);
            numsteps++;
            report_step(numsteps);
        }
        g_use_post_inflation_rhs = false;
    } else {
        g_use_post_inflation_rhs = true;
        double h = clamp_rk45_step(dt_post_inflation, rk45_max_dt);

        while (a <= af_post_inflation) {
            const double hmax = std::min(rk45_max_dt, std::max(rk45_min_dt, 2.0 * dt_post_inflation));
            h = clamp_rk45_step(h, hmax);

            bool accepted = false;
            int attempts = 0;
            while (!accepted) {
                double h_suggested = h;
                accepted = rk45_step_inflation(h, hmax, h_suggested, rk_scratch);
                h = h_suggested;
                attempts++;

                if (!accepted && (attempts > 25 || h <= rk45_min_dt * (1.0 + 1e-12))) {
                    std::fprintf(stderr,
                        "RK45 failed to converge in post-inflation at a=%e, t=%e (attempts=%d, h=%e)\n",
                        a, t, attempts, h);
                    std::exit(1);
                }
            }

            numsteps++;
            report_step(numsteps);
        }
        g_use_post_inflation_rhs = false;
    }

    printf("Saving final inflaton data\n");
    save_post_inflation(1);
    // Note: save_post_inflation() performs a temporary sync/desync for output.
    // Leave the integrator state unchanged here.
}

#endif
