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

// -------------------- Mathematical helpers --------------------

//#define float double // Uncomment to force double precision for float-based arrays

// Numerical value of pi.
const double pi = 2.0 * std::asin(1.0);

// Convenience helper returning x^2.
inline double pw2(double x) { return x * x; }

// -------------------- Indexing helpers --------------------

// Map 3D lattice indices (i,j,k) to a flat index for arrays of size N^3.
inline size_t idx(int i, int j, int k) {
    return static_cast<size_t>(i) * N * N
        + static_cast<size_t>(j) * N
        + static_cast<size_t>(k);
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

// Comoving lattice spacing.
const double dx = L / (double)N;

// Total number of lattice sites.
const int gridsize = N * N * N;

// -------------------- Global runtime state --------------------

// Time and background evolution variables.
extern double t, t0;
extern double astep, a;
extern double ad, ad2;
extern double aterm;
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
// Tensor perturbations and their time derivatives
// stored as the six independent components of a symmetric tensor.
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

// Initialization routines.
void initialize();
void initializef();
void initializeGW();
void initialize_simulation();
void initializeN();
void initialize_post_inflation();

// Field evolution routines.
void evolve_fields(double d);
void evolve_derivs(double d);
void evolve_scale(double d);
void evolve_fieldsN(double d);
void evolve_derivsN(double d);

// Energy diagnostics.
double gradient_energy();
double kin_energy();
double potential_energy();

// Potential interface.
double potential(double field_value);
double potential_derivative(int i, int j, int k);
double pot_ratio(int i, int j, int k);

// Output routines.
void output_parameters();
void save(int force);
void save_last();
void saveN();
void save_post_inflation(int force);

// Utility helpers.
void load_vector(const std::string& filename, std::vector<double>& vec);
bool ensure_results_directory();

// Main evolution drivers.
void run_evolution_loop(FILE* output_);
#if perform_deltaN
void run_deltaN_loop(FILE* output_);
#endif
#if post_inflation
void run_post_inflation_loop(FILE* output_);
#endif
