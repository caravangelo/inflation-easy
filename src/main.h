// main.h - Global variable declarations and shared definitions
/*
This file contains the global variable declarations, function declarations,
and some definitions used in many of the routines. The global variables are
defined in the file latticeeasy.cpp.
*/

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <fstream>
#include <vector>

#define float double
const float pi = (float)(2.*asin(1.));
inline float pw2(float x) {return x*x;} // Useful macro for squaring floats

/////////////////////////////////INCLUDE ADJUSTABLE PARAMETERS///////////////////
#include "parameters.h"

/////////////////////////////////GLOBAL DYNAMIC VARIABLES////////////////////////
extern float t,t0; // Current time and initial time (t0=0 unless the run is a continuation of a previous one)
extern float astep,a,ad,ad2,aterm,ad0,Ne; // Scale factor and its derivatives (aterm is a combination of the others used in the equations of motion)
extern float hubble_init; // Initial value of the Hubble constant

extern char ext_[500]; // Extension for filenames - set once and used by all functions
extern char mode_[]; // Mode in which to open files, i.e. write ("w") or append ("a+"). Depends on the variable continue_run and on whether a previous grid image was found.

const float dx=L/(float)N; // Distance between adjacent gridpoints

extern float f[N][N][N],fd[N][N][N];
#if perform_deltaN
extern float deltaN[N][N][N];
#endif

extern float fnyquist_p[N][2*N],fdnyquist_p[N][2*N];
#if numerical_potential
extern int lstart[N][N][N];
extern std::vector<float> field_numerical;
extern std::vector<float> potential_numerical;
extern std::vector<float> potential_derivative_numerical;
#endif

const int gridsize=N*N*N; // Number of spatial points in the grid
#define FIELD f[i][j][k]
#define LOOP for(i=0;i<N;i++) for(j=0;j<N;j++) for(k=0;k<N;k++)
#define INDEXLIST int i, int j, int k
#define DECLARE_INDICES int i,j,k;

/////////////////////////////////FUNCTION DECLARATIONS///////////////////////////
// Initializes the field values and Hubble parameter at the start of the simulation
void initialize();
void initializef();
void initializeN();
float gradient_energy();
// Computes the total potential energy on the lattice
float potential_energy();
float kin_energy();

void evolve_fields(float d);
void evolve_scale(float d);
void evolve_derivs(float d);

void evolve_fieldsN(float d);
void evolve_derivsN(float d);

// Computes the homogeneous potential energy for a single field value (used during initialization)
float potential_energy_hom(float field_value);
float dvdf(int i,int j, int k);
float pot_ratio(int i,int j, int k);

// Outputs the simulation parameters to file for reference
void output_parameters(); // Output information about the run parameters
void save(int force, int last); // Calculate and save quantities (means, variances, spectra, etc.)
void saveN(); // Calculate and save quantities (means, variances, spectra, etc.)
void fftrn(float f[], float fnyquist[], int ndims, int size[], int forward); // Do a Fourier transform of an ndims dimensional array of real numbers
void load_vector(const std::string& filename, std::vector<float>& vec);