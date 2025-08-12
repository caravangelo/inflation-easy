// evolution.cpp - Core evolution algorithm for lattice fields

#include "main.h"

// -------------------- Index Wrapping for Periodic Boundary Conditions --------------------

// Increment index with periodic wrapping (i → i+1 mod N)
inline int INCREMENT(int i) {
  return (i == N - 1) ? 0 : i + 1;
}

// Decrement index with periodic wrapping (i → i-1 mod N)
inline int DECREMENT(int i) {
  return (i == 0) ? N - 1 : i - 1;
}

// -------------------- Laplacians --------------------

// Helper for periodic indexing
inline int PBC(int x) {
  if (x < 0) return x + N;
  if (x >= N) return x - N;
  return x;
}

inline float lapl(int i, int j, int k, const std::vector<float>& field) {
  // 3D 2nd order central finite difference Laplacian with periodic BC
  return (
    field[idx(PBC(i+1), j, k)] + field[idx(PBC(i-1), j, k)] +
    field[idx(i, PBC(j+1), k)] + field[idx(i, PBC(j-1), k)] +
    field[idx(i, j, PBC(k+1))] + field[idx(i, j, PBC(k-1))] -
    6.f * field[idx(i, j, k)]
  );
}

#if calculate_SIGW

// Central difference for spatial derivative of f in direction dim (0=x,1=y,2=z)
inline float dfdx(int dim, int i, int j, int k, const std::vector<float>& field) {
  if (dim == 0) {
    return (field[idx(PBC(i+1), j, k)] - field[idx(PBC(i-1), j, k)]) * 0.5f;
  } else if (dim == 1) {
    return (field[idx(i, PBC(j+1), k)] - field[idx(i, PBC(j-1), k)]) * 0.5f;
  } else { // dim == 2
    return (field[idx(i, j, PBC(k+1))] - field[idx(i, j, PBC(k-1))]) * 0.5f;
  }
}

// Stress-energy tensor component T_{lm} at site (i,j,k)
float stress_energy(int l, int m, int i, int j, int k) {
  // Derivatives of phi
  float dphi_l = dfdx(l, i, j, k, f);
  float dphi_m = dfdx(m, i, j, k, f);
  //float dt_phi = fd[idx(i, j, k)]; // conformal time derivative

  // Gradient magnitude squared: sum_n (dphi_n)^2
  //float grad_phi_sq = 0.f;
  //for (int n = 0; n < 3; ++n) {
  //  float dphi_n = dfdx(n, i, j, k, f);
  //  grad_phi_sq += dphi_n * dphi_n;
  //}

  // Prefactors from metric in tilde units
  // Recall T_{ij} = d_i phi d_j phi + delta_ij [ 1/2 a^{2s} (dt phi)^2 - 1/2 (grad phi)^2 ] - (a^2 / B^2) delta_ij V(phi)

  //float delta_lm = (l == m) ? 1.f : 0.f;

  float Tij = dphi_l * dphi_m;
              //+ delta_lm * (0.5f * powf(a, 2.f * rescale_s) * dt_phi * dt_phi - 0.5f * grad_phi_sq)
              //- delta_lm * pw2(a) * potential(f[idx(i,j,k)]);

  return Tij/pw2(rescale_B);
}

#endif
// -------------------- Energy Calculations --------------------

// Compute gradient energy density (averaged)
float gradient_energy() {
  DECLARE_INDICES
  float gradient = 0.;
  float norm = pw2(1. / (a * dx));
  LOOP gradient -= f[idx(i,j,k)] * lapl(i, j, k, f);
  return 0.5 * gradient * norm / (float)gridsize;
}

// Compute kinetic energy density (averaged)
float kin_energy() {
  DECLARE_INDICES
  float deriv_energy = 0.;
  LOOP deriv_energy += pw2(fd[idx(i,j,k)]);
  deriv_energy /= (float)gridsize;
  return 0.5 * pow(a, 2. * rescale_s - 2.) * deriv_energy;
}

// -------------------- Main Field Evolution --------------------

// Leapfrog update of field derivatives
void evolve_derivs(float d) {
    DECLARE_INDICES

    float laplnorm = pow(a, -2.0f * rescale_s) / pw2(dx);
    float sfev1    = rescale_s + 1.0f;
    float sfev2    = -2.0f * rescale_s + 2.0f;

    // Update second derivative of scale factor (ad2)
    ad2 = (-2.0f * ad - 2.0f * a / d / sfev1 * (
               1.0f - sqrt(1.0f + 2.0f * d * sfev1 * ad / a +
               pw2(d) * sfev1 * pow(a, sfev2) *
               (2.0f * gradient_energy() / 3.0f + potential_energy()))
           )) / d;

    ad += 0.5f * d * ad2;

    // Scalar field evolution
#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
    LOOP {
        fd[idx(i,j,k)] += d * (
            laplnorm * lapl(i, j, k, f)
            - (2.0f + rescale_s) * ad * fd[idx(i,j,k)] / a
            - pow(a, 2.0f - 2.0f * rescale_s) * potential_derivative(i, j, k)
        );
    }

#if calculate_SIGW
    // Gravitational wave tensor evolution
#if parallel_calculation
#pragma omp parallel for collapse(4)
#endif
    for (int comp = 0; comp < 6; ++comp)
    for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
    for (int k = 0; k < N; ++k) {
        auto [l, m] = comp_to_indices(comp);
        hijd[comp][idx(i,j,k)] += d * (
            laplnorm * lapl(i, j, k, hij[comp])
            - (2.0f + rescale_s) * ad * hijd[comp][idx(i,j,k)] / a
            + 2.0f * pow(a, 2.0f - 2.0f * rescale_s) * stress_energy(l, m, i, j, k)
        );
    }
#endif

    ad += 0.5f * d * ad2;
}

// Update fields using current derivatives
void evolve_fields(float d) {
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
    for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
    for (int k = 0; k < N; ++k) {
        hij[comp][idx(i,j,k)] += d * hijd[comp][idx(i,j,k)];
    }
#endif

    a += d * ad;
}

// -------------------- Main Evolution Loop --------------------

void run_evolution_loop(FILE* output_) {
  int numsteps = 0;

  while (a <= af) {
    float dt_rescaled = dt * pow(astep, rescale_s - 1);
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
  evolve_fields(-0.5 * dt * pow(astep, rescale_s - 1));
  save_last();
}

#if perform_deltaN

// -------------------- DeltaN Evolution --------------------

// Update fd in e-folding time coordinates
void evolve_derivsN(float d) {
  DECLARE_INDICES
#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
  LOOP {
    fd[idx(i,j,k)] += d * (-(3 - 0.5 * pw2(fd[idx(i,j,k)])) * (fd[idx(i,j,k)] + pot_ratio(i, j, k)));
  }
}

// Update f and deltaN grid
void evolve_fieldsN(float d) {
  DECLARE_INDICES
  t += d;

#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
  LOOP {
#if monotonic_potential
    if (abs(f[idx(i,j,k)]) > abs(phiref))
#elif antimonotonic_potential
    if (abs(f[idx(i,j,k)]) < abs(phiref))
#else
    if (potential(f[idx(i,j,k)]) > potential(phiref))
#endif
    deltaN[idx(i,j,k)] += d;

    f[idx(i,j,k)] += d * fd[idx(i,j,k)];
  }

  a += d * ad;
}

// Determine reference φ value for ending deltaN integration
#ifndef phiref_manual
float get_phiref() {
  DECLARE_INDICES
  float fref = f[idx(0,0,0)];

#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
  LOOP {
#if monotonic_potential
    if (abs(f[idx(i,j,k)]) < abs(fref))
#elif antimonotonic_potential
    if (abs(f[idx(i,j,k)]) > abs(fref))
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
  Ne = 0;

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
}
#endif