// main.h - Global declarations and shared interfaces
//
// This header declares global variables, lattice fields, helper routines,
// and function interfaces that are shared across the simulation modules.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <stdexcept>
#include <utility>
#include <vector>
#include "parameters.h" // Simulation configuration parameters
#include "spatial_discretization.h"

// -------------------- Mathematical helpers --------------------

//#define float double // Uncomment to force double precision for float-based arrays

// Numerical value of pi.
const double pi = 2.0 * std::asin(1.0);

// Convenience helper returning x^2.
inline double pw2(double x) { return x * x; }

// -------------------- Indexing helpers --------------------

// Map 3D lattice indices (i,j,k) to a flat index for arrays of size N^3.
inline size_t idx(int i, int j, int k) {
    return spatial::index(i, j, k);
}

// Map a symmetric tensor index (i,j) with i,j in {0,1,2}
// to a packed component index in {0,...,5}.
inline int sym_idx(int i, int j) {
    if (i > j) std::swap(i, j);
    if (i == 0 && j == 0) return 0;
    if (i == 1 && j == 1) return 1;
    if (i == 2 && j == 2) return 2;
    if (i == 0 && j == 1) return 3;
    if (i == 0 && j == 2) return 4;
    if (i == 1 && j == 2) return 5;
    throw std::out_of_range("Invalid tensor indices");
}

// Inverse mapping of sym_idx: return tensor indices (l,m)
// associated with packed component comp = 0..5.
inline std::pair<int,int> comp_to_indices(int comp) {
    switch (comp) {
        case 0: return {0,0};
        case 1: return {1,1};
        case 2: return {2,2};
        case 3: return {0,1};
        case 4: return {0,2};
        case 5: return {1,2};
        default: throw std::out_of_range("Invalid component index");
    }
}

// -------------------- Derived lattice quantities --------------------

// Comoving lattice spacing (run-time, derived from L and N).
extern double dx;

// Total number of lattice sites.
const int gridsize = N * N * N;

// -------------------- Global runtime state --------------------

// Time and background evolution variables.
extern double t, t0;
extern double astep, a;
extern double ad, ad2;
extern double Ne;
extern double hubble_init;

// Reference field value used in deltaN evolution.
extern double phiref;

// Output naming and file handling.
extern char ext_[500];
extern char mode_[];

// -------------------- Lattice fields --------------------

// Scalar field and its time derivative.
extern std::vector<double> f;
extern std::vector<double> fd;

#if perform_deltaN
// Auxiliary field used for deltaN evolution.
extern std::vector<double> deltaN;
#endif

#if calculate_SIGW
// Tensor perturbations and their time derivatives, stored as the six
// components of a symmetric tensor.
extern std::vector<float> hij[6];
extern std::vector<float> hijd[6];
#endif

// Nyquist-frequency plane data required by the FFT layout.
extern double fnyquist_p[N][2 * N], fdnyquist_p[N][2 * N];
#if calculate_SIGW
extern float hijnyquist_p[6][N][2 * N], hijdnyquist_p[6][N][2 * N];
#endif

#if numerical_potential
// Tables and bookkeeping used to interpolate a numerical potential.
extern std::vector<int>    lstart;
extern std::vector<double> field_numerical;
extern std::vector<double> potential_numerical;
extern std::vector<double> potential_derivative_numerical;
#endif

// -------------------- Grid convenience macros --------------------

// Loop over all lattice sites.
#define LOOP for (i = 0; i < N; ++i) for (j = 0; j < N; ++j) for (k = 0; k < N; ++k)

// Standard argument list for functions operating on lattice indices.
#define INDEXLIST int i, int j, int k

// Local declaration of lattice indices.
#define DECLARE_INDICES int i, j, k;

// -------------------- Function declarations --------------------

// Initialization routines. These functions mutate the global simulation state.
/// Validate the initial configuration and initialize background quantities.
void initialize();
/// Allocate and initialize the scalar lattice in Fourier space, then transform it to real space.
void initializef();
#if calculate_SIGW
/// Allocate and zero the six tensor components used by the inflationary GW module.
void initializeGW();
#endif
/// Run the complete initialization sequence and write the initial output record.
void initialize_simulation();
/// Convert the final lattice state into initial data for the separate-universe deltaN evolution.
void initializeN();
/// Convert inflationary output into initial data for the optional post-inflationary stage.
void initialize_post_inflation();

// Field evolution routines. The step argument is expressed in the active code-time variable.
/// Advance the scalar, optional tensor, and scale-factor fields by one leapfrog drift.
void apply_leapfrog_drift(double step_size);

// Energy diagnostics.
/// Return the box-averaged scalar gradient energy density in code units.
double gradient_energy();
/// Return the box-averaged scalar kinetic energy density in code units.
double kin_energy();
/// Return the box-averaged scalar potential energy density in code units.
double potential_energy();

// Potential interface.
/// Evaluate the configured potential at a field value in code units.
double potential(double field_value);
/// Evaluate the configured potential derivative at a lattice site.
double potential_derivative(int i, int j, int k);
/// Evaluate the configured potential derivative without lattice-index bookkeeping.
double potential_derivative_from_value(double field_value);
/// Evaluate the potential and derivative together, optionally updating an interpolation hint.
void evaluate_potential_from_value(double field_value, int hint, int lookback, int* next_hint, double& pot, double& pot_deriv);
/// Return V'/V at a lattice site for the separate-universe equations.
double pot_ratio(int i, int j, int k);
/// Return V'/V at an arbitrary field value.
double pot_ratio_from_value(double field_value);

// Output routines.
/// Write the resolved run configuration to the metadata output.
void output_parameters();
/// Write inflationary outputs; expensive products are controlled by the argument and run-time flags.
void save(int force);
/// Write outputs that are required only at the end of inflation.
void save_last();
/// Write final deltaN products and record any incomplete-patch warning.
void saveN(FILE* output_log);
/// Write outputs for the post-inflationary stage.
void save_post_inflation(int force);

// Utility helpers.
/// Load a single-column numerical input file into a vector.
void load_vector(const std::string& filename, std::vector<double>& vec);
/// Create the main output directory if needed and reject conflicting non-directory paths.
bool ensure_results_directory();
/// Return the configured inflation integrator name for logs and metadata.
const char* integrator_name();
#if post_inflation
/// Return the configured post-inflation integrator name.
const char* post_inflation_integrator_name();
#endif
#if perform_deltaN
/// Return the configured deltaN integrator name.
const char* deltaN_integrator_name();
#endif
/// Report whether inflation uses staggered leapfrog derivatives.
inline bool inflation_uses_staggered_derivatives() {
    return integrator == INTEGRATOR_LEAPFROG;
}
#if perform_deltaN
/// Report whether the deltaN loop uses staggered leapfrog derivatives.
inline bool deltaN_uses_staggered_derivatives() {
    return deltaN_integrator == INTEGRATOR_LEAPFROG;
}
/// Return whether a separate-universe patch has not yet reached the selected hypersurface.
bool deltaN_patch_is_active(double field_value);
#endif
#if post_inflation
/// Report whether the post-inflationary loop uses staggered leapfrog derivatives.
inline bool post_inflation_uses_staggered_derivatives() {
    return post_inflation_integrator == INTEGRATOR_LEAPFROG;
}
#endif

// Main evolution drivers. Each owns the complete time loop for one simulation stage.
/// Evolve the nonlinear inflationary lattice to the requested final scale factor.
void run_inflation_loop(FILE* output_log);
#if perform_deltaN
/// Evolve each lattice site as an independent homogeneous patch and construct deltaN.
void run_deltaN_loop(FILE* output_log);
#endif
#if post_inflation
/// Evolve the post-inflationary scalar and tensor systems.
void run_post_inflation_loop(FILE* output_log);
#endif
