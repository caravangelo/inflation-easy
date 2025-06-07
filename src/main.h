// main.h - Global declarations and shared definitions for InflationEasy
//
// This header declares global variables, constants, and functions
// used across the simulation. Global variable definitions are in main.cpp.

#pragma once

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <fstream>
#include <vector>
#include "parameters.h" // Adjustable simulation parameters

// -------------------- Math and Constants --------------------

#define float double // Work with double precision

const float pi = (float)(2. * asin(1.));

inline float pw2(float x) { return x * x; } // Square of a float

// ----------------- Function to access vector elements -----------------

inline size_t idx(int i, int j, int k) {
    return static_cast<size_t>(i) * N * N + static_cast<size_t>(j) * N + static_cast<size_t>(k);
}

// -------------------- Simulation Parameters --------------------

// Grid spacing (comoving distance between points)
const float dx = L / (float)N;

// Total number of points in the 3D grid
const int gridsize = N * N * N;

// -------------------- Global Dynamic Variables --------------------

// Time and scale factor evolution
extern float t, t0;             // Time and initial time
extern float astep, a;          // Scale factor and backup (for leapfrog)
extern float ad, ad2;           // First and second derivatives of scale factor
extern float aterm;             // Intermediate quantity used in evolution
extern float Ne;                // e-folding number used in deltaN evolution
extern float hubble_init;       // Initial Hubble parameter

// Reference field value for deltaN termination
extern float phiref;

// I/O formatting
extern char ext_[500];          // Filename extension for outputs
extern char mode_[];            // File open mode (write/append)

// -------------------- Lattice Fields --------------------

// Main field and its time derivative
extern std::vector<float> f;
extern std::vector<float> fd;

#if perform_deltaN
// deltaN field for separate evolution
extern std::vector<float> deltaN;
#endif

// Nyquist frequency components for FFT
extern float fnyquist_p[N][2 * N], fdnyquist_p[N][2 * N];

#if numerical_potential
// Interpolation for numerical potential
extern std::vector<int> lstart;                     // Interpolation index grid
extern std::vector<float> field_numerical;          // Grid of ϕ values
extern std::vector<float> potential_numerical;      // V(ϕ) values
extern std::vector<float> potential_derivative_numerical; // dV/dϕ values
#endif

// -------------------- Grid Macros --------------------
                             
#define LOOP for(i=0;i<N;i++) for(j=0;j<N;j++) for(k=0;k<N;k++) // Loop over full grid
#define INDEXLIST int i, int j, int k                   // Function parameter list for indexing
#define DECLARE_INDICES int i, j, k;                    // Local index variable declarations

// -------------------- Function Declarations --------------------

// Initialization
void initialize();               // Basic parameter checks and hubble_init
void initializef();              // Generate vacuum fluctuations
void initializeN();              // Setup for deltaN calculation
void initialize_simulation();    // Full initialization pipeline

// Field evolution
void evolve_fields(float d);     // Evolve f with step d
void evolve_derivs(float d);     // Evolve fd with step d
void evolve_scale(float d);      // (Unused but present)
void evolve_fieldsN(float d);    // Evolve f during deltaN run
void evolve_derivsN(float d);    // Evolve fd during deltaN run

// Energies
float gradient_energy();         // Compute spatial gradient energy
float kin_energy();              // Compute kinetic energy
float potential_energy();        // Compute total potential energy on grid

// Potential interface
float potential(float field_value);                 // V(ϕ) for single value
float potential_derivative(int i, int j, int k);    // ∂V/∂ϕ at a grid point
float pot_ratio(int i, int j, int k);               // ∂V/∂ϕ / V

// Output routines
void output_parameters();       // Write simulation parameters to file
void save(int force);           // Save observables (means, spectra, etc.)
void save_last();               // Final snapshot
void saveN();                   // Save for deltaN evolution

// Utilities
void fftrn(float f[], float fnyquist[], int ndims, int size[], int forward); // Real FFT
void load_vector(const std::string& filename, std::vector<float>& vec);      // Read floats from file
bool ensure_results_directory();                                             // Create output directory if needed

// Main evolution loops
void run_evolution_loop(FILE* output_);
#if perform_deltaN
void run_deltaN_loop(FILE* output_);
#endif
