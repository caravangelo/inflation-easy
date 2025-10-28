// evolution.cpp - Core evolution algorithm for lattice fields

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

// Generic 7-point Laplacian for real-space arrays (works for double or float vectors)
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

#if calculate_SIGW
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
};

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
#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
    LOOP {
        fd[idx(i,j,k)] += d * (
            laplnorm * lapl<double>(i, j, k, f)
            - (2.0 + rescale_s) * ad * fd[idx(i,j,k)] / a
            - std::pow(a, 2.0 - 2.0 * rescale_s) * potential_derivative(i, j, k)
        );
    }

#if calculate_SIGW
    // Gravitational wave tensor evolution (store as float, compute RHS in double)
#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
    for (int i=0; i<N; ++i)
    for (int j=0; j<N; ++j)
    for (int k=0; k<N; ++k) {
        const int id = idx(i,j,k);

        // build gradients once
        GradPack G; build_grad_pack(i,j,k,G);

        // explicit six components (double)
        const double T_xx = stress_energy_fast(0,0,G);
        const double T_yy = stress_energy_fast(1,1,G);
        const double T_zz = stress_energy_fast(2,2,G);
        const double T_xy = stress_energy_fast(0,1,G);
        const double T_xz = stress_energy_fast(0,2,G);
        const double T_yz = stress_energy_fast(1,2,G);

        const double srcAmp = 2.0 * std::pow(a, -2.0 * rescale_s);

        // compute each RHS entirely in double; cast once on store
        {
          const double lap_h = static_cast<double>(lapl<float>(i,j,k,hij[0]));
          const double rhs = d * ( laplnorm * lap_h
                                 - (2.0 + rescale_s) * ad * static_cast<double>(hijd[0][id]) / a
                                 + srcAmp * T_xx );
          hijd[0][id] += static_cast<float>(rhs);
        }
        {
          const double lap_h = static_cast<double>(lapl<float>(i,j,k,hij[1]));
          const double rhs = d * ( laplnorm * lap_h
                                 - (2.0 + rescale_s) * ad * static_cast<double>(hijd[1][id]) / a
                                 + srcAmp * T_yy );
          hijd[1][id] += static_cast<float>(rhs);
        }
        {
          const double lap_h = static_cast<double>(lapl<float>(i,j,k,hij[2]));
          const double rhs = d * ( laplnorm * lap_h
                                 - (2.0 + rescale_s) * ad * static_cast<double>(hijd[2][id]) / a
                                 + srcAmp * T_zz );
          hijd[2][id] += static_cast<float>(rhs);
        }
        {
          const double lap_h = static_cast<double>(lapl<float>(i,j,k,hij[3]));
          const double rhs = d * ( laplnorm * lap_h
                                 - (2.0 + rescale_s) * ad * static_cast<double>(hijd[3][id]) / a
                                 + srcAmp * T_xy );
          hijd[3][id] += static_cast<float>(rhs);
        }
        {
          const double lap_h = static_cast<double>(lapl<float>(i,j,k,hij[4]));
          const double rhs = d * ( laplnorm * lap_h
                                 - (2.0 + rescale_s) * ad * static_cast<double>(hijd[4][id]) / a
                                 + srcAmp * T_xz );
          hijd[4][id] += static_cast<float>(rhs);
        }
        {
          const double lap_h = static_cast<double>(lapl<float>(i,j,k,hij[5]));
          const double rhs = d * ( laplnorm * lap_h
                                 - (2.0 + rescale_s) * ad * static_cast<double>(hijd[5][id]) / a
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

// -------------------- Main Evolution Loop --------------------

void run_evolution_loop(FILE* output_) {
  int numsteps = 0;

  evolve_fields(0.5 * dt); // First leapfrog step

  while (a <= af) {
    double dt_rescaled = dt * std::pow(astep, rescale_s - 1.0);
    evolve_derivs(dt_rescaled);
    evolve_fields(dt_rescaled);

    numsteps++;

    if (numsteps % output_freq == 0 && a < af) {
      save((numsteps % output_infrequent_freq == 0) ? 1 : 0);
    }

    if (screen_updates && numsteps % output_freq == 0) {
      printf("scale factor a = %f\n", a);
      printf("numsteps %i\n\n", numsteps);
    }

    fprintf(output_, "scale factor a = %f\n", a);
    fprintf(output_, "numsteps %i\n\n", numsteps);
    fflush(output_);

    astep = a;
  }

  printf("Saving final inflaton data\n");
  evolve_fields(-0.5 * dt * std::pow(astep, rescale_s - 1.0)); // sync field values wih velocities
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

// Determine reference φ value for ending deltaN integration
#ifndef phiref_manual
double get_phiref() {
  DECLARE_INDICES
  double fref = f[idx(0,0,0)];

#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
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
#endif

// Run main deltaN integration loop
void run_deltaN_loop(FILE* output_) {
  printf("Starting deltaN calculation\n");
  fprintf(output_, "Starting deltaN calculation\n");

  int numsteps = 0;
  Ne = 0.0;

  initializeN();
  evolve_fieldsN(0.5 * dN);

#ifdef phiref_manual
  phiref = phiref_manual;
#else
  phiref = get_phiref();
#endif

  while (Ne <= Nend) {
    evolve_derivsN(dN);
    evolve_fieldsN(dN);
    Ne += dN;
    numsteps++;

    if (screen_updates && numsteps % output_freq == 0) {
      printf("N = %f\n\n", Ne);
    }

    fprintf(output_, "N = %f\n\n", Ne);
    fflush(output_);
  }

  saveN();
  evolve_fieldsN(-0.5 * dN);
}
#endif

#if post_inflation
// -------------------- Post-inflation Evolution --------------------

// ===== derivative pack (uses existing dfdx, idx, INCREMENT/DECREMENT, N, dx) =====
struct DerivPack {
  double gf[3], gfd[3];  // ∂_i f, ∂_i fd
  double Hf[6];          // Hessian(f): [xx, yy, zz, xy, xz, yz] in your sym_idx order
  double f_here;
};

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


// ===== fast source using your sym_idx and comp_to_indices =====
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
#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
for (int i=0; i<N; ++i)
for (int j=0; j<N; ++j)
for (int k=0; k<N; ++k) {
    const int id = idx(i,j,k);

    // build derivatives once per site
    DerivPack P; 
    build_deriv_pack(i,j,k,P);

    // bare sources (not TT-projected), explicit pairs
    const double S_xx = stress_energy_post_inflation_fast(0,0,P);
    const double S_yy = stress_energy_post_inflation_fast(1,1,P);
    const double S_zz = stress_energy_post_inflation_fast(2,2,P);
    const double S_xy = stress_energy_post_inflation_fast(0,1,P);
    const double S_xz = stress_energy_post_inflation_fast(0,2,P);
    const double S_yz = stress_energy_post_inflation_fast(1,2,P);

    // RHS amplitude (C=2 convention)
    const double srcAmp = 2.0 * std::pow(a, -2.0 * rescale_s);

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
  evolve_derivs_post_inflation(0.5 * dt_post_inflation);

  while (a <= af_post_inflation) {
    double dt_rescaled = dt_post_inflation; // * std::pow(astep, rescale_s - 1);
    evolve_derivs_post_inflation(dt_rescaled);
    evolve_fields(dt_rescaled);

    numsteps++;

    if (numsteps % output_freq == 0) {
      save_post_inflation((numsteps % output_infrequent_freq == 0) ? 1 : 0);
    }

    if (screen_updates && numsteps % output_freq == 0) {
      printf("scale factor a = %f\n", a);
      printf("numsteps %i\n\n", numsteps);
    }

    fprintf(output_, "scale factor a = %f\n", a);
    fprintf(output_, "numsteps %i\n\n", numsteps);
    fflush(output_);

    astep = a;
  }

  printf("Saving final inflaton data\n");
  evolve_fields(-0.5 * dt * std::pow(astep, rescale_s - 1.0));
}

#endif