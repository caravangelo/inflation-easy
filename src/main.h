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

// #define float double // Work with double precision   <-- REMOVED

const double pi = (double)(2. * asin(1.));

inline double pw2(double x) { return x * x; } // Square of a number

// ----------------- Function to access vector elements -----------------

inline size_t idx(int i, int j, int k) {
    return static_cast<size_t>(i) * N * N + static_cast<size_t>(j) * N + static_cast<size_t>(k);
}

// Map tensor index (i,j) to the unique component index (0 to 5)
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

// Returns the symmetric tensor indices (l,m) for component index comp = 0..5
inline std::pair<int,int> comp_to_indices(int comp) {
    switch(comp) {
        case 0: return {0,0};
        case 1: return {1,1};
        case 2: return {2,2};
        case 3: return {0,1};
        case 4: return {0,2};
        case 5: return {1,2};
        default: throw std::out_of_range("Invalid component index");
    }
}

// -------------------- Simulation Parameters --------------------

// Grid spacing (comoving distance between points)
const double dx = L / (double)N;

// Total number of points in the 3D grid
const int gridsize = N * N * N;

// -------------------- Global Dynamic Variables --------------------

// Time and scale factor evolution
extern double t, t0;             // Time and initial time
extern double astep, a;          // Scale factor and backup (for leapfrog)
extern double ad, ad2;           // First and second derivatives of scale factor
extern double aterm;             // Intermediate quantity used in evolution
extern double Ne;                // e-folding number used in deltaN evolution
extern double hubble_init;       // Initial Hubble parameter

// Reference field value for deltaN termination
extern double phiref;

// I/O formatting
extern char ext_[500];          // Filename extension for outputs
extern char mode_[];            // File open mode (write/append)

// -------------------- Lattice Fields --------------------

// Main field and its time derivative (double)
extern std::vector<double> f;
extern std::vector<double> fd;

#if perform_deltaN
// deltaN field for separate evolution (double)
extern std::vector<double> deltaN;
#endif

#if calculate_SIGW
// Gravitational waves (float)
extern std::vector<float> hij[6];
extern std::vector<float> hijd[6];
#endif

// Nyquist frequency components for FFT
extern double fnyquist_p[N][2 * N], fdnyquist_p[N][2 * N];
#if calculate_SIGW
extern float hijnyquist_p[6][N][2*N], hijdnyquist_p[6][N][2*N];
#endif

#if numerical_potential
// Interpolation for numerical potential (double)
extern std::vector<int>    lstart;                     
extern std::vector<double> field_numerical;            
extern std::vector<double> potential_numerical;        
extern std::vector<double> potential_derivative_numerical; 
#endif

// -------------------- Grid Macros --------------------
                             
#define LOOP for(i=0;i<N;i++) for(j=0;j<N;j++) for(k=0;k<N;k++) // Loop over full grid
#define INDEXLIST int i, int j, int k                   // Function parameter list for indexing
#define DECLARE_INDICES int i, j, k;                    // Local index variable declarations

// -------------------- Function Declarations --------------------

// Initialization
void initialize();               // Basic parameter checks and hubble_init
void initializef();              // Generate vacuum fluctuations
void initializeGW(); 
void initialize_simulation();    // Full initialization pipeline
void initializeN();              // Setup for deltaN calculation
void initialize_post_inflation();// Post inflation initialization pipeline

// Field evolution
void evolve_fields(double d);     // Evolve f with step d
void evolve_derivs(double d);     // Evolve fd with step d
void evolve_scale(double d);      // (Unused but present)
void evolve_fieldsN(double d);    // Evolve deltaN with step d
void evolve_derivsN(double d);    // Evolve deltaN derivative with step d

// Energies
double gradient_energy();         
double kin_energy();              
double potential_energy();        

// Potential interface
double potential(double field_value);                 
double potential_derivative(int i, int j, int k);     
double pot_ratio(int i, int j, int k);                

// Output routines
void output_parameters();               
void save(int force);                   
void save_last();                     
void saveN();                          
void save_post_inflation(int force);    

// Utilities
void load_vector(const std::string& filename, std::vector<double>& vec);      
bool ensure_results_directory();                                              

// Main evolution loops
void run_evolution_loop(FILE* output_);
#if perform_deltaN
void run_deltaN_loop(FILE* output_);
#endif
#if post_inflation
void run_post_inflation_loop(FILE* output_);
#endif
