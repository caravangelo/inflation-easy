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

inline float lapl(int i, int j, int k, const std::vector<float>& field) {
 if (i == 0 || j == 0 || k == 0 || i == N - 1 || j == N - 1 || k == N - 1) {
    return (
      field[idx(i,j,INCREMENT(k))] + field[idx(i,j,DECREMENT(k))] +
      field[idx(i,INCREMENT(j),k)] + field[idx(i,DECREMENT(j),k)] +
      field[idx(INCREMENT(i),j,k)] + field[idx(DECREMENT(i),j,k)] -
      6. * field[idx(i,j,k)]
    );
  } else {
    return (
      field[idx(i,j,k+1)] + field[idx(i,j,k-1)] +
      field[idx(i,j+1,k)] + field[idx(i,j-1,k)] +
      field[idx(i+1,j,k)] + field[idx(i-1,j,k)] -
      6. * field[idx(i,j,k)]
    );
  }
}

#if calculate_SIGW

// Central difference for spatial derivative of field in direction dim (0=x,1=y,2=z)
inline float dfdx(int dim, int i, int j, int k, const std::vector<float>& field) {
  if (dim == 0) {
    if (i == 0 || i == N - 1) {
      return (field[idx(INCREMENT(i), j, k)] - field[idx(DECREMENT(i), j, k)]) * (0.5f / dx);
    } else {
      return (field[idx(i+1, j, k)] - field[idx(i-1, j, k)]) * (0.5f / dx);
    }
  } else if (dim == 1) {
    if (j == 0 || j == N - 1) {
      return (field[idx(i, INCREMENT(j), k)] - field[idx(i, DECREMENT(j), k)]) * (0.5f / dx);
    } else {
      return (field[idx(i, j+1, k)] - field[idx(i, j-1, k)]) * (0.5f / dx);
    }
  } else { // dim == 2
    if (k == 0 || k == N - 1) {
      return (field[idx(i, j, INCREMENT(k))] - field[idx(i, j, DECREMENT(k))]) * (0.5f / dx);
    } else {
      return (field[idx(i, j, k+1)] - field[idx(i, j, k-1)]) * (0.5f / dx);
    }
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

  return Tij; // Tij/pw2(rescale_B);
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

    // Scalar field evolution (unchanged)
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
    // Precompute RHS prefactor for the source term
    const float fac_source = 2.0f * pow(a, 2.0f - 2.0f * rescale_s);

    // Gravitational wave tensor evolution:
    // compute spatial derivatives once per cell, build trace-free Pi_{ij} and update all 6 components
#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
    LOOP {
        const int id = idx(i,j,k);

        // compute tilde-derivatives once
        float d0 = dfdx(0, i, j, k, f);
        float d1 = dfdx(1, i, j, k, f);
        float d2 = dfdx(2, i, j, k, f);

        float grad_sq = d0*d0 + d1*d1 + d2*d2;
        float one_third = (1.f/3.f) * grad_sq;

        // Laplacians of hij components (each uses the same laplnorm)
        // Note: lapl reads hij[comp] neighbors; kept as-is to match your stencil
        float lap_h0 = lapl(i, j, k, hij[0]); // xx
        float lap_h1 = lapl(i, j, k, hij[1]); // yy
        float lap_h2 = lapl(i, j, k, hij[2]); // zz
        float lap_h3 = lapl(i, j, k, hij[3]); // xy
        float lap_h4 = lapl(i, j, k, hij[4]); // xz
        float lap_h5 = lapl(i, j, k, hij[5]); // yz

        // friction factor for h' term (same as scalar)
        float fric_common = (2.0f + rescale_s) * ad / a;

        // compute source (trace-free Pi components in tilde indices)
        float Pi_xx = d0 * d0 - one_third;
        float Pi_yy = d1 * d1 - one_third;
        float Pi_zz = d2 * d2 - one_third;
        float Pi_xy = d0 * d1;
        float Pi_xz = d0 * d2;
        float Pi_yz = d1 * d2;

        // update hijd components (leapfrog style update of derivatives)
        hijd[0][id] += d * (
            laplnorm * lap_h0
            - fric_common * hijd[0][id] 
            + fac_source * Pi_xx
        );

        hijd[1][id] += d * (
            laplnorm * lap_h1
            - fric_common * hijd[1][id] 
            + fac_source * Pi_yy
        );

        hijd[2][id] += d * (
            laplnorm * lap_h2
            - fric_common * hijd[2][id] 
            + fac_source * Pi_zz
        );

        hijd[3][id] += d * (
            laplnorm * lap_h3
            - fric_common * hijd[3][id] 
            + fac_source * Pi_xy
        );

        hijd[4][id] += d * (
            laplnorm * lap_h4
            - fric_common * hijd[4][id] 
            + fac_source * Pi_xz
        );

        hijd[5][id] += d * (
            laplnorm * lap_h5
            - fric_common * hijd[5][id] 
            + fac_source * Pi_yz
        );
    }
#endif

    ad += 0.5f * d * ad2;
}

// Leapfrog update of field derivatives
void evolve_derivs_old(float d) {
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
    LOOP {
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
    LOOP {
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