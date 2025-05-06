// main.cpp - Program entry point and top-level driver for InflationEasy
#include "main.h"
#include <iostream>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>

// Field values and derivatives
float f[N][N][N],fd[N][N][N];
#if perform_deltaN
float deltaN[N][N][N];
#endif
// Nyquist values needed for FFT
float fnyquist_p[N][2*N],fdnyquist_p[N][2*N];

float t,t0; // Current time and initial time
float astep=1.,a=1.,ad=0.,ad2=0.,ad0=0.,Ne=0.; // Scale factor and its derivatives (aterm is a combination of the others used in the equations of motion).
float hubble_init=0.; // Initial value of the Hubble constant
char mode_[10]="w"; // Mode in which to open files
char ext_[500]=".dat"; // Extension for filenames - set once and used by all output functions

// Quantities only needed if working with a numerical potential
#if numerical_potential
int lstart[N][N][N];
std::vector<float> field_numerical;
std::vector<float> potential_numerical;
std::vector<float> potential_derivative_numerical;
#endif

int main()
{
  
  int numsteps=0; // Quantities for counting how often to calculate and output derived quantities
  
  FILE *output_=fopen("results/output.txt","w"); // Outputs time. Used to remotely monitor progress
  
  if(seed<1) // The use of seed<1 turns off certain functions (random numbers, fourier transforms, gradients, and potential energy) and should only be used for debugging
  printf("Warning: The parameter seed has been set to %d, which will result in incorrect output. For correct output set seed to a positive integer.",seed);
  
  // Ensure "results" directory exists
  struct stat info;
  if (stat("results", &info) != 0) {
    // Directory does not exist, create it
    if (mkdir("results", 0777) != 0) {
      std::cerr << "Failed to create 'results' directory." << std::endl;
      return 1;
    }
  } else if (!(info.st_mode & S_IFDIR)) {
    std::cerr << "'results' exists but is not a directory." << std::endl;
    return 1;
  }
  
  #if numerical_potential
  load_vector("inputs/field_values.dat", field_numerical);
  load_vector("inputs/potential.dat", potential_numerical);
  load_vector("inputs/potential_derivative.dat", potential_derivative_numerical);
  
  std::ifstream inr("inputs/field_value.dat");
  if (inr.is_open()) {
    float element;
    while (inr >> element) {
      field_numerical.push_back(element);
    }
  }
  inr.close();
  
  std::ifstream inrd("inputs/potential.dat");
  if (inrd.is_open()) {
    float element;
    while (inrd >> element) {
      potential_numerical.push_back(element);
    }
  }
  inrd.close();
  
  std::ifstream ini("inputs/potential_derivative.dat");
  if (ini.is_open()) {
    float element;
    while (ini >> element) {
      potential_derivative_numerical.push_back(element);
    }
  }
  ini.close();
  #endif
  
  initialize(); // Initialize background
  initializef(); // Initialize flucuations
  
  output_parameters(); // Initialize some outputs
  save(1,0);
  t=t0;
  
  evolve_fields(.5*dt); // First time step of the leapfrog integrator
  
  while(a<=af) // Main time evolution loop
  {
    // evolve derivatives and fields, allow for adaptive time step
    evolve_derivs(dt*pow(astep,rescale_s-1));
    evolve_fields(dt*pow(astep,rescale_s-1));
    
    numsteps++;
    if(numsteps%output_freq == 0 && a<af)
    {
      if(numsteps%output_infrequent_freq == 0)
      save(1,0);
      else
      save(0,0);
    }
    
    if(screen_updates) // This option determines whether or not to update progress on the screen
    {
      if(numsteps%output_freq == 0)
      {
        
        printf("scale factor a = %f\n",a);
        printf("numsteps %i\n",numsteps);
        printf("\n");

      }
    }
    fprintf(output_,"scale factor a = %f\n",a);// Output progress to a file for monitoring progress
    fprintf(output_,"numsteps %i\n",numsteps);
    fprintf(output_,"\n");
    fflush(output_); // Make sure output file is always up to date

    astep = a;
  }
  
  printf("Saving final inflaton data\n");
  save(1,1);
  evolve_fields(-.5*dt*pow(astep,rescale_s-1));
  
  #if perform_deltaN
  
  // Starting the deltaN calculation to compute zeta. See arXiv:?????
  printf("Starting deltaN calculation\n");
  fprintf(output_,"Starting deltaN calculation\n");

  numsteps = 0;
  initializeN(); //convert the grid to e-folds time and initialize relevant quantities
  Ne=0;
  
  evolve_fieldsN(.5*dN);
  while(Ne<=Nend) // Main time evolution loop
  {
    evolve_derivsN(dN);
    evolve_fieldsN(dN);
    numsteps++;
    Ne+=dN;
    
    if(screen_updates) // This option determines whether or not to update progress on the screen
    {
      if(numsteps%output_freq == 0)
      {
        printf("N = %f\n",Ne);
        printf("\n");
       
      }
    }
    fprintf(output_,"N = %f\n",Ne);
    fprintf(output_,"\n");
    fflush(output_); // Make sure output file is always up to date
  }
  
  saveN(); // Calculate and save deltaN
  
  #endif
  
  output_parameters(); // Save run parameters and elapsed time
  fprintf(output_,"InflationEasy program finished\n");
  printf("InflationEasy program finished\n");
  printf("NUMSTEPS = %d\n", numsteps);
  
  return(0);
}

void load_vector(const std::string& filename, std::vector<float>& vec) {
  std::ifstream infile(filename);
  if (!infile) {
    std::cerr << "Could not open file: " << filename << std::endl;
    exit(1);
  }
  
  float value;
  while (infile >> value) {
    vec.push_back(value);
  }
  
  if (vec.empty()) {
    std::cerr << "File appears to be empty or malformed: " << filename << std::endl;
    exit(1);
  }
}