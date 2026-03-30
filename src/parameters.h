// parameters.h - Simulation configuration interface (hybrid model)
//
// This header documents two classes of parameters:
//
// 1) Compile-time parameters (edited here, then recompiled):
//    - feature toggles and lattice layout
//    - preprocessor-controlled branches
//
// 2) Run-time parameters (edited in params.txt, no recompilation):
//    - physical model values
//    - integration/output controls
//    - deltaN and post-inflation settings
//
// Defaults for run-time parameters are defined in runtime_parameters.cpp.
// If params.txt is present, those values override the defaults at startup.

#pragma once

// -------------------- Compile-time feature toggles --------------------

// Set to 1 to enable a numerical potential loaded from file, 0 for analytic expression.
#define numerical_potential 1

// 1: evolve the auxiliary deltaN system to compute zeta.
#define perform_deltaN 1

// 1: compute scalar-induced gravitational waves (SIGWs) during inflation.
#define calculate_SIGW 1

// 1: extend SIGW computation to a post-inflationary phase.
#define post_inflation 1

#if post_inflation && !calculate_SIGW
#error "post_inflation requires calculate_SIGW = 1"
#endif

// 1: enable OpenMP parallelism in selected loops (requires compiler support).
#define parallel_calculation 1

// -------------------- Compile-time lattice parameters --------------------

// Number of grid points per spatial dimension (total grid size N^3).
// Must be a power of 2 for the FFT routines used in the code.
const int N = 128;

// -------------------- Monotonicity hints (compile-time) --------------------

// These compile-time flags select which deltaN stopping-condition potential
// criterion is compiled in (field rolling direction toward/away from the end surface).
// They are NOT runtime params.txt keys.
//
// Default behavior keeps the historical auto-switch:
//   - numerical_potential = 1 -> monotonic=1, anti-monotonic=0
//   - numerical_potential = 0 -> monotonic=0, anti-monotonic=1
//
// For custom models, manually edit the two defines below if needed:
//   - 1/0 -> monotonic potential criterion: evolve while |phi| > |phi_ref|
//   - 0/1 -> anti-monotonic potential criterion: evolve while |phi| < |phi_ref|
//   - 0/0 -> generic potential criterion: evolve while V(phi) > V(phi_ref)
#if numerical_potential
#define monotonic_potential 1
#define antimonotonic_potential 0
#else
#define monotonic_potential 0
#define antimonotonic_potential 1
#endif

#if monotonic_potential && antimonotonic_potential
#error "Set at most one of monotonic_potential / antimonotonic_potential to 1."
#endif

// -------------------- Run-time parameters --------------------
// Defaults are defined in runtime_parameters.cpp and can be overridden by params.txt.
// Keep these as extern declarations only; do not assign values here.

// Random seed used to initialize vacuum fluctuations.
extern int seed;

// Field-rescaling exponent used by internal conventions.
// Keep at 0 unless a different scaling has been validated.
extern double rescale_s;

// Final scale factor for the main inflationary run (optional in params.txt; defaults to 2*N).
extern double af;

// Potential normalization.
extern double V0;

// Rescaling factor B (if not explicitly set in params.txt, defaults to sqrt(V0)).
extern double rescale_B;

#if !numerical_potential
// Spectral index parameter used by the default analytic hilltop branch.
extern double ns;
#endif

// Homogeneous initial field value and time derivative (code units).
extern double initial_field;
extern double initial_derivative;

// Comoving box size and main time step (code units).
extern double L;
extern double dt;

// Main inflation integrator selection.
// 0: leapfrog (default, staggered field derivatives)
// 1: RK4 (fixed-step)
// 2: RK45 (adaptive-step Dormand-Prince)
enum IntegratorType {
    INTEGRATOR_LEAPFROG = 0,
    INTEGRATOR_RK4 = 1,
    INTEGRATOR_RK45 = 2
};
extern int integrator;
#if post_inflation
// Post-inflation integrator selection (same enum values as `integrator`).
extern int post_inflation_integrator;
#endif
#if perform_deltaN
// deltaN-loop integrator selection (same enum values as `integrator`).
extern int deltaN_integrator;
#endif

// RK45 controls (used when either integrator selection is INTEGRATOR_RK45).
// If rk45_min_dt <= 0, runtime uses 1e-6 * dt.
// If rk45_max_dt <= 0, runtime uses dt.
extern double rk45_abs_tol;
extern double rk45_rel_tol;
extern double rk45_min_dt;
extern double rk45_max_dt;
extern double rk45_safety;

// Output cadences in number of integration steps.
extern int output_freq;
extern int output_infrequent_freq;

#if perform_deltaN
// deltaN integration controls.
extern double dN;
extern double Nend;

// Manual reference-field option for deltaN stopping hypersurface.
// If use_phiref_manual != 0, phiref_manual_value is used.
extern int use_phiref_manual;
extern double phiref_manual_value;
#endif

#if post_inflation
// Post-inflation controls (used only if post_inflation==1).
extern double horizon_factor;
extern double omega;
extern double dt_post_inflation;
// Optional in params.txt; defaults to 2*N if not provided.
extern double af_post_inflation;
#endif

// Optional Fourier-space cutoff controls.
extern double high_cutoff_index;
extern double low_cutoff_index;
extern int forcing_cutoff;

// Output switches (1 = on, 0 = off).
extern int output_spectra;
extern int output_histogram;
extern int output_energy;
extern int output_box3D;
extern int output_box2D;
extern int output_bispectrum;

#if perform_deltaN
// LOG-based zeta approximation controls.
extern int output_LOG;
extern double eta_log;
#endif

// Screen progress output and histogram resolution.
extern int screen_updates;
extern int nbins;

#if numerical_potential
// Interpolation look-back controls for tabulated numerical potential.
extern int int_err;
extern int int_errN;
#endif

// Derived run-time quantity: comoving lattice spacing dx = L/N.
extern double dx;

// Load run-time overrides from a key=value text file.
// Unknown keys are ignored with a warning. Missing file is allowed.
void load_runtime_parameters(const char* filename = "params.txt");
