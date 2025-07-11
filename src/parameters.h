// -------------------- InflationEasy Configuration Flags --------------------

// Set to 1 to enable a numerical potential loaded from file, 0 for analytic expression
#define numerical_potential 1

// Set to 1 to perform deltaN evolution to calculate zeta
#define perform_deltaN 1

/// Set to 1 to enable OpenMP parallelism in selected loops (if supported by the compiler)
#define parallel_calculation 1

// -------------------- Lattice and Evolution Parameters --------------------

// Number of points along each spatial dimension (total lattice points will be N^3)
const int N = 128; // Must be a power of 2

// Rescaling exponent (must be 0 unless you've tested otherwise)
const float rescale_s = 0;

// Final scale factor value (simulation stops when a >= af)
const float af = 1000.;

// Random seed for generating initial vacuum fluctuations
const int seed = 8;

// -------------------- Physical Model Parameters --------------------

#if numerical_potential
// Parameters used if the potential is loaded from file

const float V0 = 3e-9; // Potential parameter
const float rescale_B = sqrt(V0); // Rescaling factor (see doc)
const float initial_field = 2.9181235049318586; // Initial field value
const float initial_derivative = -0.06727651095116181; // Initial field velocity (in code units)
const float L = 10.; // Comoving box size in code units (must be sub-horizon at start)
const float dt = 0.001; // Time step
// Output frequency in time steps (standard and infrequent quantities)
const int output_freq = 100;
const int output_infrequent_freq = 100;

#if perform_deltaN
const float dN = 0.0001; // e-fold increment in deltaN evolution
const float Nend = 2.;   // Final e-folding time for deltaN run
// const float phiref_manual = value; // (Optional) manually set the reference ϕ value

// Set monotonic_potential = 1 if the potential increases with |ϕ|
// Set antimonotonic_potential = 1 if it decreases with |ϕ|
// Leave both set to 0 for general potentials, in which case deltaN evolution is slower. 
#define monotonic_potential 1
#define antimonotonic_potential 0

#endif

#else
// Parameters used for an analytic potential (defaults are for quadratic potential)

const float V0 = 3.338e-13; // Potential parameter
const float ns = 0.97; // Spectral index
const float rescale_B = sqrt(V0); // Rescaling factor (see doc)
const float initial_field = 0.0935; // Initial field value
const float initial_derivative = 0.000796; // Initial field velocity (in code units)
const float L = 10.; // Comoving box size in code units (must be sub-horizon at start)
const float dt = 0.0005; // Time step
// Output frequency in time steps (standard and infrequent quantities)
const int output_freq = 200;
const int output_infrequent_freq = 200;

#if perform_deltaN
const float dN = 0.0000001; // e-fold increment in deltaN evolution
const float Nend = 0.001; // Final e-folding time for deltaN run
// const float phiref_manual = value; // (Optional) manually set the reference ϕ value

// Set monotonic_potential = 1 if the potential increases with |ϕ|
// Set antimonotonic_potential = 1 if it decreases with |ϕ|
// Leave both set to 0 for general potentials, in which case deltaN evolution is slower.
#define monotonic_potential 0
#define antimonotonic_potential 1

#endif
#endif

// -------------------- Momentum Cutoff Options --------------------

// Set high_cutoff_index > 0 to cutoff modes with k > (2 pi/L)*high_index and k < (2 pi/L)*low_index
const float high_cutoff_index = 0;
const float low_cutoff_index = 0;

// Set to 1 to enforce the cutoff continuously, not just at initialization
const int forcing_cutoff = 0;

// -------------------- Output Configuration --------------------

// Toggle various outputs (1 = enabled)
const int output_spectra = 1;
const int output_histogram = 1;
const int output_energy = 1;
const int output_box3D = 0;
const int output_box2D = 1;
const int output_bispectrum = 0;

#if perform_deltaN
const int output_LOG = 0;     // Calculates zeta with the log relation (see documentation)
const float eta_log = -0.5;   // Value of constant eta assumed in the log formula
#endif

// Print time/step info to console
const int screen_updates = 1;

// Number of bins used in field histograms
const int nbins = 256;

// -------------------- Numerical Potential Handling --------------------

#if numerical_potential
const int int_err = 5;   // Increase if the code gives "interpolation error"
const int int_errN = 1;  // Same, but for deltaN interpolation
#endif
