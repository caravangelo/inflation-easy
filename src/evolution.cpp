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

// Laplacian of fd used in evolution equations (with edge handling)
inline float lapld(INDEXLIST) {
  if (i == 0 || j == 0 || k == 0 || i == N - 1 || j == N - 1 || k == N - 1) {
    return (
      fd[i][j][INCREMENT(k)] + fd[i][j][DECREMENT(k)] +
      fd[i][INCREMENT(j)][k] + fd[i][DECREMENT(j)][k] +
      fd[INCREMENT(i)][j][k] + fd[DECREMENT(i)][j][k] -
      6. * fd[i][j][k]
    );
  } else {
    return (
      fd[i][j][k + 1] + fd[i][j][k - 1] +
      fd[i][j + 1][k] + fd[i][j - 1][k] +
      fd[i + 1][j][k] + fd[i - 1][j][k] -
      6. * fd[i][j][k]
    );
  }
}

// Laplacian of f (same structure as above)
inline float lapl(INDEXLIST) {
  if (i == 0 || j == 0 || k == 0 || i == N - 1 || j == N - 1 || k == N - 1) {
    return (
      f[i][j][INCREMENT(k)] + f[i][j][DECREMENT(k)] +
      f[i][INCREMENT(j)][k] + f[i][DECREMENT(j)][k] +
      f[INCREMENT(i)][j][k] + f[DECREMENT(i)][j][k] -
      6. * f[i][j][k]
    );
  } else {
    return (
      f[i][j][k + 1] + f[i][j][k - 1] +
      f[i][j + 1][k] + f[i][j - 1][k] +
      f[i + 1][j][k] + f[i - 1][j][k] -
      6. * f[i][j][k]
    );
  }
}

// -------------------- Energy Calculations --------------------

// Compute gradient energy density (averaged)
float gradient_energy() {
  DECLARE_INDICES
  double gradient = 0.;
  float norm = pw2(1. / (a * dx));
  LOOP gradient -= f[i][j][k] * lapl(i, j, k);
  return 0.5 * gradient * norm / (float)gridsize;
}

// Compute kinetic energy density (averaged)
float kin_energy() {
  DECLARE_INDICES
  double deriv_energy = 0.;
  LOOP deriv_energy += pw2(fd[i][j][k]);
  deriv_energy /= (float)gridsize;
  return 0.5 * pow(a, 2. * rescale_s - 2.) * deriv_energy;
}

// -------------------- Main Field Evolution --------------------

// Leapfrog update of fd and ad
void evolve_derivs(float d) {
  DECLARE_INDICES

  float laplnorm = pow(a, -2. * rescale_s) / pw2(dx);
  float sfev1 = rescale_s + 1.;
  float sfev2 = -2. * rescale_s + 2.;

  // Update second derivative of scale factor (ad2)
  ad2 = (-2. * ad - 2. * a / d / sfev1 * (
            1. - sqrt(1. + 2. * d * sfev1 * ad / a +
            pw2(d) * sfev1 * pow(a, sfev2) *
            (2. * gradient_energy() / 3. + potential_energy()))
         )) / d;

  ad += 0.5 * d * ad2;

#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
  LOOP {
    fd[i][j][k] += d * (
      laplnorm * lapl(i, j, k)
      - (2. + rescale_s) * ad * fd[i][j][k] / a
      - pow(a, 2. - 2. * rescale_s) * potential_derivative(i, j, k)
    );
  }

  ad += 0.5 * d * ad2;
}

// Update f using current fd (standard leapfrog step)
void evolve_fields(float d) {
  DECLARE_INDICES
  t += d;

#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
  LOOP f[i][j][k] += d * fd[i][j][k];

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
    fd[i][j][k] += d * (-(3 - 0.5 * pw2(fd[i][j][k])) * (fd[i][j][k] + pot_ratio(i, j, k)));
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
    if (abs(f[i][j][k]) > abs(phiref))
#else
    if (potential(f[i][j][k]) > potential(phiref))
#endif
    deltaN[i][j][k] += d;

    f[i][j][k] += d * fd[i][j][k];
  }

  a += d * ad;
}

// Determine reference ϕ value for ending deltaN integration
#ifndef phiref_manual
float get_phiref() {
  DECLARE_INDICES
  float fref = f[0][0][0];

#if parallel_calculation
#pragma omp parallel for collapse(3)
#endif
  LOOP {
#if monotonic_potential
    if (abs(f[i][j][k]) < abs(fref))
#else
    if (potential(f[i][j][k]) < potential(fref))
#endif
    fref = f[i][j][k];
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
