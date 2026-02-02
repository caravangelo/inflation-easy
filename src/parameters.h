// parameters.h - Simulation configuration (compile-time / run-time constants)
//
// This header contains the knobs that define *what* the code runs and *how* it runs:
// - feature toggles (numerical potential, deltaN, SIGWs, post-inflation phase)
// - lattice size and stepping parameters
// - model parameters for either a numerical or analytic potential
// - output switches and numerical cutoffs
//
// Conventions:
// - Preprocessor flags (#define) enable/disable modules at compile time.
// - Constants (const ...) are read by the code as fixed parameters.
//

// -------------------- Feature toggles --------------------

// 1: load inflaton potential from file; 0: use analytic potential.
#define numerical_potential 1

// 1: evolve the auxiliary deltaN system to compute zeta.
#define perform_deltaN 1

// 1: compute scalar-induced gravitational waves (SIGWs) during inflation.
#define calculate_SIGW 1

// 1: extend SIGW computation to a post-inflationary phase.
#define post_inflation 1

// -------------------- Lattice and evolution parameters --------------------

// Number of grid points per spatial dimension (total grid size N^3).
// Must be a power of 2 for the FFT routines used in the code.
const int N = 128;

// Field-rescaling exponent used by the code's internal conventions.
// Keep at 0 unless you have validated an alternative scaling.
const double rescale_s = 0.0;

// End condition for the inflationary evolution (stop once a >= af).
const double af = 2*N;

// Seed for the random-number generator used to initialize vacuum fluctuations.
const int seed = 8;

// -------------------- Physical model parameters --------------------

#if numerical_potential
// Numerical-potential branch: the potential and its derivative are read from file.
// The parameters below set the overall normalization and initial conditions.

const double V0 = 3e-9; // Potential normalization (model-dependent)
const double rescale_B = std::sqrt(V0); // Field/energy rescaling factor (see documentation)

// Homogeneous initial conditions (field value and velocity in code units)
const double initial_field = 2.9181235049318586;
const double initial_derivative = -0.06727651095116181;

// Comoving box size (code units). Must be sub-horizon at initialization.
const double L = 10.0;

// Time step used during inflationary evolution.
const double dt = 0.001;

// Output cadences (in number of time steps).
// "Standard" quantities and "infrequent/heavier" quantities can be separated.
const int output_freq = 500;
const int output_infrequent_freq = 500;

#if perform_deltaN
// deltaN evolution parameters.
const double dN = 0.0001; // Step in e-folds for deltaN integration
const double Nend = 5.0;  // Final e-fold time for deltaN run

// Optional manual reference field value for terminating/defining deltaN.
//const double phiref_manual = 2.680654050913246;

// Monotonicity hints for the potential table.
// - monotonic_potential = 1 if V increases with |phi|
// - antimonotonic_potential = 1 if V decreases with |phi|
// Set both to 0 for a general non-monotonic potential (slower but robust).
#define monotonic_potential 1
#define antimonotonic_potential 0

#endif

#if post_inflation
// Post-inflation evolution controls (if enabled).
const double horizon_factor = 1.0;

// Equation of state parameter w (only tested for radiation: w = 1/3).
const double omega = 1.0/3.0;

const double dt_post_inflation = 0.001;
const double af_post_inflation = 2.*N*horizon_factor;
#endif

#else
// Analytic-potential branch: parameters below correspond to the chosen analytic V(phi).
// The defaults shown here are for a quadratic potential.

const double V0 = 3.338e-13; // Potential normalization
const double ns = 0.97;       // Spectral index (used by some analytic setups)
const double rescale_B = std::sqrt(V0);

const double initial_field = 0.0935;
const double initial_derivative = 0.000796;

const double L = 10.0;   // Comoving box size (must be sub-horizon initially)
const double dt = 0.0005;

// Output cadences (in number of time steps).
const int output_freq = 200;
const int output_infrequent_freq = 200;

#if perform_deltaN
// deltaN evolution parameters (analytic-potential branch).
const double dN = 0.0000001;
const double Nend = 0.001;

// Optional manual reference field value for terminating/defining deltaN.
// const double phiref_manual = value;

// Monotonicity hints for the analytic potential.
#define monotonic_potential 0
#define antimonotonic_potential 1

#endif
#endif

// -------------------- Momentum cutoff options --------------------

// If high_cutoff_index > 0: remove modes with k > (2*pi/L)*high_cutoff_index.
// If low_cutoff_index  > 0: remove modes with k < (2*pi/L)*low_cutoff_index.
const double high_cutoff_index = 0.0;
const double low_cutoff_index = 0.0;

// 1: enforce the cutoff continuously during evolution (not only at initialization).
const int forcing_cutoff = 0;

// -------------------- Output configuration --------------------

// Toggle individual outputs (1 = enabled, 0 = disabled).
const int output_spectra = 1;
const int output_histogram = 1;
const int output_energy = 1;
const int output_box3D = 0;
const int output_box2D = 0;
const int output_bispectrum = 0;

#if perform_deltaN
// 1: compute zeta using the log relation (see documentation for assumptions).
const int output_LOG = 0;

// Constant eta assumed by the log relation.
const double eta_log = -0.5;
#endif

// 1: print periodic progress updates to stdout.
const int screen_updates = 1;

// Number of bins used for field histograms.
const int nbins = 256;

// -------------------- Numerical-potential handling --------------------

#if numerical_potential
// Interpolation tolerances for the numerical potential tables.
// Increase if you encounter interpolation warnings/errors.
const int int_err = 5;
const int int_errN = 5;
#endif

// 1: enable OpenMP parallelism in selected loops (requires compiler support).
#define parallel_calculation 1
