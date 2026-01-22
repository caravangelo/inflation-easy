// initialize.cpp - Initialization routines for field values and simulation parameters

#include "main.h"
#include "ffteasy.hpp"

// -------------------- Random Number Generator --------------------
#define randa 16807
#define randm 2147483647
#define randq 127773
#define randr 2836

double rand_uniform(void) {
  if (seed < 1) return 0.33; // Fixed fallback for debugging
  static int i = 0;
  static int next = seed;
  if (!(next > 0)) {
    printf("Invalid seed used in random number function. Using seed=1\n");
    next = 1;
  }
  if (i == 0) for (i = 1; i < 100; i++) rand_uniform();

  next = randa * (next % randq) - randr * (next / randq);
  if (next < 0) next += randm;
  return (static_cast<double>(next) / static_cast<double>(randm));
}

#undef randa
#undef randm
#undef randq
#undef randr

// -------------------- High-Level Simulation Initialization --------------------

void initialize_simulation() {
  initialize();        // Basic checks and global param setup
  initializef();       // Set field configuration and apply FFT
  #if calculate_SIGW
  initializeGW();
  #endif
  output_parameters(); // Save run configuration
  save(1);             // First output
  t = t0;              // Set initial time
}

// -------------------- Field Mode Initialization --------------------
// Set amplitude and phase of vacuum mode using Rayleigh distribution
void set_mode(double p2, double p2lat, double m2, double *field, double *deriv, int real) {
    double phase, amplitude, rms_amplitude, omega;
    double re_f_left, im_f_left;
    static double norm = rescale_B * pow(L / pw2(dx), 1.5); // Normalization (see 2209.13616)
    double ii = L / (2.0 * pi) * sqrt(p2lat);
    static double hbterm = -hubble_init;
    static int tachyonic = 0;

    omega = (p2 + m2 > 0.0) ? sqrt(p2 + m2) : sqrt(p2);
    if (p2 + m2 <= 0.0 && tachyonic == 0) {
        printf("Warning: Tachyonic mode(s) may be initialized inaccurately\n");
        tachyonic = 1;
    }

    rms_amplitude = (omega > 0.0) ? norm / sqrt(2.0 * omega) : 0.0;

    // Apply high/low frequency cutoff
    if (high_cutoff_index > 0 && (ii > high_cutoff_index || ii < low_cutoff_index)) {
        rms_amplitude = 0.0;
    }

    amplitude = rms_amplitude * sqrt(log(1.0 / rand_uniform()));
    phase = 2.0 * pi * rand_uniform();

    re_f_left = amplitude * cos(phase);
    im_f_left = amplitude * sin(phase);

    field[0] = re_f_left;
    field[1] = im_f_left;
    deriv[0] = omega * im_f_left + hbterm * field[0];
    deriv[1] = -omega * re_f_left + hbterm * field[1];

    if (real == 1) {
        field[1] = 0.0;
        deriv[1] = 0.0;
    }
}

// -------------------- Basic Global Initialization --------------------

void initialize() {
  if (dt > dx / sqrt(3.0)) {
    printf("Time step too large. Courant condition violated: dt/dx = %f\n", dt / dx);
    exit(1);
  }

  printf("Generating initial conditions for new run at t = 0\n");

  t0 = 0.0;
  hubble_init = sqrt(potential(initial_field) / 3.0);

  if (hubble_init < 0.0) {
    printf("Error: Initial Hubble constant is negative or undefined\n");
    exit(1);
  }

  ad = hubble_init;

#if numerical_potential
  lstart.resize(N * N * N, 0);
#endif
}

// -------------------- Vacuum Fluctuation Initialization --------------------

void initializef() {
  f.resize(N * N * N, 0.0);
  fd.resize(N * N * N, 0.0);

  double p2, p2lat, mass_sq = 0.0;
  double dp2 = pw2(2.0 * pi / L);
  double initial_field_values = initial_field;
  double initial_field_derivs = initial_derivative;
  int i, j, k, iconj, jconj;
  double px, py, pz;
  int arraysize[] = {N, N, N};

  if (dt > dx / sqrt(3.0)) {
    printf("Time step too large. Courant condition violated: dt/dx = %f\n", dt / dx);
    exit(1);
  }

  hubble_init = sqrt((0.5 * pw2(initial_field_derivs) + potential(initial_field_values)) / 3.0);
  if (hubble_init < 0.0) {
    printf("Error: Initial Hubble constant is negative or undefined\n");
    exit(1);
  }
  ad = hubble_init;

  for (i = 0; i < N; i++) {
    px = (i <= N / 2 ? i : i - N);
    iconj = (i == 0 ? 0 : N - i);

    for (j = 0; j < N; j++) {
      py = (j <= N / 2 ? j : j - N);

      for (k = 1; k < N / 2; k++) {
        pz = k;
        p2lat = dp2 * (pw2(px) + pw2(py) + pw2(pz));
        p2 = 4.0 * pw2(N / L) * (pw2(sin(px * pi / N)) + pw2(sin(py * pi / N)) + pw2(sin(pz * pi / N)));
        set_mode(p2, p2lat, mass_sq, &f[idx(i,j,2*k)], &fd[idx(i,j,2*k)], 0);
      }

      if (j > N / 2 || (i > N / 2 && (j == 0 || j == N / 2))) {
        jconj = (j == 0 ? 0 : N - j);

        p2lat = dp2 * (pw2(px) + pw2(py));
        p2 = 4.0 * pw2(N / L) * (pw2(sin(px * pi / N)) + pw2(sin(py * pi / N)));
        set_mode(p2, p2lat, mass_sq, &f[idx(i,j,0)], &fd[idx(i,j,0)], 0);

        f[idx(iconj,jconj,0)] = f[idx(i,j,0)];
        f[idx(iconj,jconj,1)] = -f[idx(i,j,1)];
        fd[idx(iconj,jconj,0)] = fd[idx(i,j,0)];
        fd[idx(iconj,jconj,1)] = -fd[idx(i,j,1)];

        p2lat = dp2 * (pw2(px) + pw2(py) + pw2(N / 2));
        p2 = 4.0 * pw2(N / L) * (pw2(sin(px * pi / N)) + pw2(sin(py * pi / N)) + 1.0);
        set_mode(p2, p2lat, mass_sq, &fnyquist_p[i][2*j], &fdnyquist_p[i][2*j], 0);

        fnyquist_p[iconj][2*jconj]   = fnyquist_p[i][2*j];
        fnyquist_p[iconj][2*jconj+1] = -fnyquist_p[i][2*j+1];
        fdnyquist_p[iconj][2*jconj]   = fdnyquist_p[i][2*j];
        fdnyquist_p[iconj][2*jconj+1] = -fdnyquist_p[i][2*j+1];
      } else if ((i == 0 || i == N / 2) && (j == 0 || j == N / 2)) {
        p2lat = dp2 * (pw2(px) + pw2(py));
        p2 = 4.0 * pw2(N / L) * (pw2(sin(px * pi / N)) + pw2(sin(py * pi / N)));
        if (p2 > 0.0) set_mode(p2, p2lat, mass_sq, &f[idx(i,j,0)], &fd[idx(i,j,0)], 1);

        p2lat = dp2 * (pw2(px) + pw2(py) + pw2(N / 2));
        p2 = 4.0 * pw2(N / L) * (pw2(sin(px * pi / N)) + pw2(sin(py * pi / N)) + 1.0);
        set_mode(p2, p2lat, mass_sq, &fnyquist_p[i][2*j], &fdnyquist_p[i][2*j], 1);
      }
    }
  }

  f[idx(0,0,0)]  = 0.0;
  f[idx(0,0,1)]  = 0.0;
  fd[idx(0,0,0)] = 0.0;
  fd[idx(0,0,1)] = 0.0;

  // Inverse FFT to real space (double-precision buffers)
  fftrnd(f.data(),  &fnyquist_p[0][0],  3, arraysize, -1);
  fftrnd(fd.data(), &fdnyquist_p[0][0], 3, arraysize, -1);

  LOOP {
    f[idx(i,j,k)]  += initial_field_values;
    fd[idx(i,j,k)] += initial_field_derivs;
  }

  double vard = 0.0;
  LOOP vard += pw2(fd[idx(i,j,k)]);
  double deriv_energy_in = 0.5 * vard / static_cast<double>(gridsize);

  hubble_init = sqrt((deriv_energy_in + potential_energy() + gradient_energy()) / 3.0);
  if (hubble_init < 0.0) {
    printf("Error in calculating initial Hubble constant. Exiting.\n");
    exit(1);
  }

  ad = hubble_init;
  printf("Finished initial conditions\n");
}

// -------------------- DeltaN Initialization --------------------

#if perform_deltaN
void initializeN() {
  DECLARE_INDICES
  deltaN.resize(N * N * N, 0.0);

  double H_lat = 0.0;
  LOOP {
#if numerical_potential
    lstart[idx(i,j,k)] -= 100;
#endif
    H_lat = sqrt(potential(f[idx(i,j,k)]) / 3.0);
    fd[idx(i,j,k)] = fd[idx(i,j,k)] * pow(a, rescale_s - 1.0) / H_lat;
    deltaN[idx(i,j,k)] = 0.0;
  }
}
#endif

#if calculate_SIGW
void initializeGW() {
  const size_t gs = static_cast<size_t>(N) * N * N;
  for (int c = 0; c < 6; ++c) {
    hij[c].assign(gs, 0.0f);   // size to gs and zero all entries
    hijd[c].assign(gs, 0.0f);
  }
}
#endif

#if post_inflation
void initialize_post_inflation() {

  t0 = 0.0;

  const size_t gs = static_cast<size_t>(N) * N * N;
  for (int c = 0; c < 6; ++c) {
    hij[c].assign(gs, 0.0f);   // size to gs and zero all entries
    hijd[c].assign(gs, 0.0f);
  }

  DECLARE_INDICES
  
  #if perform_deltaN
  double Nmean = 0.;
  LOOP
  Nmean += deltaN[idx(i,j,k)];
  Nmean = Nmean / static_cast<double>(gridsize);
  LOOP f[idx(i,j,k)] = 2.0/3.0 * (deltaN[idx(i,j,k)] - Nmean); // Important: only valid for radiation!
  #else
  double fmean = 0.0;
  double fdmean = 0.0;
  LOOP
  { 
  fmean += f[idx(i,j,k)];
  fdmean += fd[idx(i,j,k)];

  }
  fmean = fmean / static_cast<double>(gridsize);
  fdmean = fdmean / static_cast<double>(gridsize);
  LOOP f[idx(i,j,k)] = 2.0/3.0 * (f[idx(i,j,k)] - fmean) / (fdmean * pow(a, rescale_s - 1.0) / (ad * pow(a, rescale_s - 2.0)));
  #endif

  fd.resize(N * N * N, 0.0);

  a = 1.0;
  ad = horizon_factor * N * 2.0 * pi / L; //setting ad is delicate, think 

  save_post_inflation(1);
}
#endif
