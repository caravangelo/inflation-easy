// initialize.cpp - Initialization routines for fields and background
//
// This file implements the initialization steps required before time
// evolution starts: consistency checks, background quantities, vacuum
// fluctuations for the scalar field, and optional auxiliary systems
// such as deltaN and tensor perturbations.


#include "main.h"
#include "ffteasy.hpp"

namespace {
// Construct H from an energy density only after checking the square-root input.
double checked_hubble(double energy_density, const char* context) {
    if (!std::isfinite(energy_density) || !(energy_density > 0.0)) {
        std::fprintf(stderr,
            "Cannot initialize the Hubble parameter from %s: energy density must be finite and positive (got %.17g).\n",
            context, energy_density);
        std::exit(EXIT_FAILURE);
    }

    const double value = std::sqrt(energy_density / 3.0);
    if (!std::isfinite(value) || !(value > 0.0)) {
        std::fprintf(stderr,
            "Cannot initialize the Hubble parameter from %s: result is not finite and positive.\n",
            context);
        std::exit(EXIT_FAILURE);
    }
    return value;
}
} // namespace

// -------------------- Random Number Generator --------------------
#define randa 16807
#define randm 2147483647
#define randq 127773
#define randr 2836

// Return a reproducible uniform deviate in (0,1) from the legacy minimal-standard RNG.
// The global seed fixes the sequence; seed < 1 selects a deterministic debug fallback.
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

// Execute initialization in dependency order: validate the background, create
// lattice fields, allocate optional sectors, and record the initial state.
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
// Draw one Gaussian vacuum mode in the packed complex FFT layout. `p2` is the
// lattice-dispersion momentum used in omega, while `p2lat` labels the DFT mode
// used by the optional spectral cutoffs. `real` enforces a self-conjugate mode.
void set_mode(
    double p2,
    double p2lat,
    const std::optional<double>& mass_squared,
    double *field,
    double *deriv,
    int real)
{
    double phase, amplitude, rms_amplitude, omega;
    double re_f_left, im_f_left;
    static double norm = rescale_B * pow(L / pw2(dx), 1.5); // Normalization (see 2209.13616)
    double ii = L / (2.0 * pi) * sqrt(p2lat);
    static double hbterm = -hubble_init;
    static int tachyonic = 0;

    double omega_squared = p2;
    if (mass_squared.has_value()) omega_squared += *mass_squared;

    omega = (omega_squared > 0.0) ? sqrt(omega_squared) : sqrt(p2);
    if (omega_squared <= 0.0 && tachyonic == 0) {
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

// Validate the stencil-dependent Courant condition and initialize the homogeneous expansion state.
void initialize() {
    const double courant_limit = spatial::courant_dt_over_dx_limit();
    const double courant_ratio = dt / dx;
    if (!std::isfinite(courant_ratio) || !(courant_ratio > 0.0)
        || courant_ratio > courant_limit) {
        printf("Time step too large for the order-%d stencil: dt/dx = %f, limit = %f\n",
               spatial::order, courant_ratio, courant_limit);
        exit(1);
    }

    printf("Generating initial conditions for new run at t = 0\n");

    t0 = 0.0;
    const double initial_background_energy =
        0.5 * pw2(initial_derivative) + potential(initial_field);
    hubble_init = checked_hubble(initial_background_energy, "the homogeneous initial energy");

    ad = hubble_init;

#if numerical_potential
    lstart.resize(N * N * N, 0);
#endif
}

// -------------------- Vacuum Fluctuation Initialization --------------------

// Generate Bunch-Davies scalar fluctuations in Fourier space, enforce Hermitian
// symmetry, transform to real space, and add the homogeneous initial values.
// The distinct DFT and effective lattice momenta must remain consistent with
// the spatial stencil and spectral-output conventions.
void initializef() {
    f.resize(N * N * N, 0.0);
    fd.resize(N * N * N, 0.0);

    double p2, p2lat;
    double dp2 = pw2(2.0 * pi / L);
    double initial_field_values = initial_field;
    double initial_field_derivs = initial_derivative;
    int i, j, k, iconj, jconj;
    int px, py, pz;
    int arraysize[] = {N, N, N};

    const double homogeneous_energy =
        0.5 * pw2(initial_field_derivs) + potential(initial_field_values);
    hubble_init = checked_hubble(homogeneous_energy, "the homogeneous initial energy");
    ad = hubble_init;

    for (i = 0; i < N; i++) {
        px = spatial::signed_mode(i);
        iconj = (i == 0 ? 0 : N - i);

        for (j = 0; j < N; j++) {
            py = spatial::signed_mode(j);

            for (k = 1; k < N / 2; k++) {
                pz = k;
                p2lat = dp2 * (pw2(px) + pw2(py) + pw2(pz));
                p2 = spatial::effective_momentum_squared(px, py, pz, dx);
                set_mode(p2, p2lat, initial_mass_squared, &f[idx(i,j,2*k)], &fd[idx(i,j,2*k)], 0);
            }

            if (j > N / 2 || (i > N / 2 && (j == 0 || j == N / 2))) {
                jconj = (j == 0 ? 0 : N - j);

                p2lat = dp2 * (pw2(px) + pw2(py));
                p2 = spatial::effective_momentum_squared(px, py, 0, dx);
                set_mode(p2, p2lat, initial_mass_squared, &f[idx(i,j,0)], &fd[idx(i,j,0)], 0);

                f[idx(iconj,jconj,0)] = f[idx(i,j,0)];
                f[idx(iconj,jconj,1)] = -f[idx(i,j,1)];
                fd[idx(iconj,jconj,0)] = fd[idx(i,j,0)];
                fd[idx(iconj,jconj,1)] = -fd[idx(i,j,1)];

                p2lat = dp2 * (pw2(px) + pw2(py) + pw2(N / 2));
                p2 = spatial::effective_momentum_squared(px, py, N / 2, dx);
                set_mode(p2, p2lat, initial_mass_squared, &fnyquist_p[i][2*j], &fdnyquist_p[i][2*j], 0);

                fnyquist_p[iconj][2*jconj]   = fnyquist_p[i][2*j];
                fnyquist_p[iconj][2*jconj+1] = -fnyquist_p[i][2*j+1];
                fdnyquist_p[iconj][2*jconj]   = fdnyquist_p[i][2*j];
                fdnyquist_p[iconj][2*jconj+1] = -fdnyquist_p[i][2*j+1];
            } else if ((i == 0 || i == N / 2) && (j == 0 || j == N / 2)) {
                p2lat = dp2 * (pw2(px) + pw2(py));
                p2 = spatial::effective_momentum_squared(px, py, 0, dx);
                if (p2 > 0.0) set_mode(p2, p2lat, initial_mass_squared, &f[idx(i,j,0)], &fd[idx(i,j,0)], 1);

                p2lat = dp2 * (pw2(px) + pw2(py) + pw2(N / 2));
                p2 = spatial::effective_momentum_squared(px, py, N / 2, dx);
                set_mode(p2, p2lat, initial_mass_squared, &fnyquist_p[i][2*j], &fdnyquist_p[i][2*j], 1);
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

    const double total_initial_energy =
        deriv_energy_in + potential_energy() + gradient_energy();
    hubble_init = checked_hubble(total_initial_energy, "the initialized lattice energy");

    ad = hubble_init;
    printf("Finished initial conditions\n");
}

// -------------------- DeltaN Initialization --------------------

#if perform_deltaN
// Prepare the separate-universe variables at every lattice site. The velocity
// is converted from code time to e-fold time before the deltaN loop begins.
void initializeN() {
    DECLARE_INDICES
    deltaN.resize(N * N * N, 0.0);

    LOOP {
#if numerical_potential
        lstart[idx(i,j,k)] -= 100;
#endif
        const size_t id = idx(i,j,k);
        const double local_potential = potential(f[id]);
        if (!std::isfinite(local_potential) || !(local_potential > 0.0)) {
            std::fprintf(stderr,
                "Cannot initialize deltaN at site (%d,%d,%d): local potential must be finite and positive (got %.17g).\n",
                i, j, k, local_potential);
            std::exit(EXIT_FAILURE);
        }
        const double velocity_in_cosmic_time =
            fd[id] * std::pow(a, rescale_s - 1.0);
        const double local_energy =
            0.5 * velocity_in_cosmic_time * velocity_in_cosmic_time
            + local_potential;
        if (!std::isfinite(velocity_in_cosmic_time)
            || !std::isfinite(local_energy) || !(local_energy > 0.0)) {
            std::fprintf(stderr,
                "Cannot initialize deltaN at site (%d,%d,%d): local homogeneous energy and velocity must be finite, with positive energy.\n",
                i, j, k);
            std::exit(EXIT_FAILURE);
        }
        const double H_lat = std::sqrt(local_energy / 3.0);
        const double velocity_in_efold_time = velocity_in_cosmic_time / H_lat;
        if (!std::isfinite(H_lat) || !(H_lat > 0.0)
            || !std::isfinite(velocity_in_efold_time)) {
            std::fprintf(stderr,
                "Cannot initialize deltaN at site (%d,%d,%d): non-finite local Hubble parameter or e-fold velocity.\n",
                i, j, k);
            std::exit(EXIT_FAILURE);
        }
        fd[id] = velocity_in_efold_time;
        deltaN[id] = 0.0;
    }
}
#endif

#if calculate_SIGW
// Allocate the six packed components of the symmetric tensor field and velocity.
void initializeGW() {
    const size_t gs = static_cast<size_t>(N) * N * N;
    for (int c = 0; c < 6; ++c) {
        hij[c].assign(gs, 0.0f);   // size to gs and zero all entries
        hijd[c].assign(gs, 0.0f);
    }
}
#endif

#if post_inflation
// Map the final inflationary scalar data to the Newtonian-potential initial
// condition used by the post-inflationary perfect-fluid evolution.
void initialize_post_inflation() {

    t0 = 0.0;
    t = 0.0;

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
    const double Phi_from_zeta = 3.0 * (1.0 + omega) / (5.0 + 3.0 * omega);
    LOOP f[idx(i,j,k)] = Phi_from_zeta * (deltaN[idx(i,j,k)] - Nmean);

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
    if (!std::isfinite(fmean) || !std::isfinite(fdmean)) {
        std::fprintf(stderr,
            "Cannot initialize the post-inflationary scalar from non-finite background field averages.\n");
        std::exit(EXIT_FAILURE);
    }
    const double Phi_from_zeta = 3.0 * (1.0 + omega) / (5.0 + 3.0 * omega);
    const double mapping_denominator =
        fdmean * std::pow(a, rescale_s - 1.0)
        / (ad * std::pow(a, rescale_s - 2.0));
    if (!std::isfinite(mapping_denominator) || mapping_denominator == 0.0) {
        std::fprintf(stderr,
            "Cannot initialize the post-inflationary scalar: the background field-velocity mapping denominator is zero or non-finite.\n");
        std::exit(EXIT_FAILURE);
    }
    LOOP f[idx(i,j,k)] = -Phi_from_zeta
        * (f[idx(i,j,k)] - fmean)
        / mapping_denominator;

#endif

    // The perfect-fluid evolution starts from a constant Newtonian potential.
    fd.assign(gs, 0.0);

    a = 1.0;
    ad = horizon_factor * N * 2.0 * pi / L;

    save_post_inflation(1);
}
#endif
