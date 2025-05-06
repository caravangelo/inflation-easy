// initialize.cpp - Field and parameter setup routines
//Done with initialize except potential
/*
This file contains the initialization function for setting field and parameter values.
*/

#include "main.h"

// Generate a uniform deviate between 0 and 1 using the Park-Miller minimum standard algorithm
#define randa 16807
#define randm 2147483647
#define randq 127773
#define randr 2836
float rand_uniform(void)
{
  if(seed<1) return(0.33); // *DEBUG* This is used to avoid randomness, for debugging purposes only.
  static int i=0;
  static int next=seed;
  if(!(next>0)) // Guard against 0, negative, or other invalid seeds
  {
    printf("Invalid seed used in random number function. Using seed=1\n");
    next=1;
  }
  if(i==0) // On the first call run through 100 calls. This allows small seeds without the first results necessarily being small.
  for(i=1;i<100;i++)
  rand_uniform();
  next = randa*(next%randq) - randr*(next/randq);
  if(next<0) next += randm;
  return ((float)next/(float)randm);
}
#undef randa
#undef randm
#undef randq
#undef randr

// Set the amplitude and phase of a mode of vacuum fluctuations
// Phase is set randomly (except for modes that must be real).
// Amplitude is set with a Rayleigh distribution about an rms value of 1/sqrt(omega) (times normalization terms).

void set_mode(float p2, float p2lat, float m2, float *field, float *deriv, int real)
{
  
  // IMORTANT: this function only works for sub-horizon box (Bunch-Davies initial conditions);
  float phase,amplitude,rms_amplitude,omega;
  float re_f_left,im_f_left;
  static float norm = rescale_B*pow(L/pw2(dx),1.5); // See 2209.13616 for normalization convention
  float ii;
  static float hbterm = hubble_init*(-1.);
  static int tachyonic=0; // Avoid printing the same error repeatedly
  
  // Momentum cutoff. If kcutoff!=0 then eliminate all initial modes with k>kcutoff.
  static float k2cutoff = (kcutoff<2.*pi*(float)N/L ? pw2(kcutoff) : 0.);
  if(k2cutoff>0. && p2>k2cutoff)
  {
    
    field[0]=0.;
    field[1]=0.;
    deriv[0]=0.;
    deriv[1]=0.;
    return;
  }
  ii = L/(2.*pi)*sqrt(p2lat);
  if(p2+m2>0.) // Check to avoid floating point errors
  omega=sqrt(p2+m2); // Omega = Sqrt(p^2 + m^2)
  else
  {
    if(tachyonic==0)
    printf("Warning: Tachyonic mode(s) may be initialized inaccurately\n");
    omega=sqrt(p2); // If p^2 + m^2 < 0 use m^2=0
    tachyonic=1;
  }
  if(omega>0. && ii<high_cutoff_index && ii>low_cutoff_index)
  rms_amplitude=norm/sqrt(2*omega);
  else
  rms_amplitude=0.;
  
  // Amplitude = RMS amplitude x Rayleigh distributed random number
  amplitude=rms_amplitude*sqrt(log(1./rand_uniform()));
  
  phase=2.*pi*rand_uniform(); // Set phase randomly
  re_f_left = amplitude*cos(phase);
  im_f_left = amplitude*sin(phase);
  
  field[0] = re_f_left; // Re(field)
  field[1] = im_f_left; // Im(field)
  deriv[0] = omega*(im_f_left) + hbterm*field[0]; // Field derivative
  deriv[1] = -omega*(re_f_left) + hbterm*field[1];
  if(real==1) // For real modes set the imaginary parts to zero
  {
    field[1] = 0.;
    deriv[1] = 0.;
  }
  
  return;
}

/////////////////////////////////////////////////////
// Externally called function(s)
/////////////////////////////////////////////////////

// Set initial parameters and field background values
// Initialize the field configuration and parameters
// Initializes the field values and Hubble parameter at the start of the simulation
void initialize()
{
 
  // Check to make sure time step is small enough to satisfy Courant condition, dt/dx < 1/Sqrt(ndims)
  if(dt>dx/sqrt(3.))
  {
    printf("Time step too large. The ratio dt/dx is currently %f but for stability should never exceed 1/sqrt(%d) (%f)\n",dt/dx,3,1./sqrt(3.));
    printf("Adjust dt to AT MOST %e, and preferably somewhat smaller than that.\n",dx/sqrt(3.));
    exit(1);
  }
  
  printf("Generating initial conditions for new run at t=0\n");
  
  t0=0;
  
  // Set initial value of Hubble constant - it will be updated below to account for gradients and kinetic terms
  hubble_init = sqrt((potential_energy_hom(initfield))/3.);
  
  if(!(hubble_init>=0.)) // Make sure Hubble isn't negative or undefined
  {
    printf("Error in calculating initial Hubble constant. Exiting.\n");
    exit(1);
  }
  ad=hubble_init;
  ad0=hubble_init;
  
  #if numerical_potential
  int i,j,k;
  LOOP
  lstart[i][j][k]=0;
  #endif
  
  return;
}

void initializef()
{
  float p2,p2lat; // Total squared momentum
  float dp2=pw2(2.*pi/L); // Square of grid spacing in momentum space
  float mass_sq; // Effective mass squared of fields
  float initial_field_values; // Initial value of fields (set to zero unless specified in parameters.h)
  float initial_field_derivs; // Initial value of field derivatives (set to zero unless specified in parameters.h)

  int i,j,k,iconj,jconj; // iconj and jconj are the indices where conjugate modes to i and j are stored
  float px,py,pz; // Components of momentum in units of grid spacing
  //float fnyquist[N][2*N],fdnyquist[N][2*N]; // Used by FFT routine to store modes with k=nyquist
  
  int arraysize[]={N,N,N}; // Array of grid size in all dimensions - used by FFT routine
  
  // Check to make sure time step is small enough to satisfy Courant condition, dt/dx < 1/Sqrt(ndims)
  if(dt>dx/sqrt(3.))
  {
    printf("Time step too large. The ratio dt/dx is currently %f but for stability should never exceed 1/sqrt(%d) (%f)\n",dt/dx,3,1./sqrt(3.));
    printf("Adjust dt to AT MOST %e, and preferably somewhat smaller than that.\n",dx/sqrt(3.));
    exit(1);
  }
  
  // Initiating scalar field.
  initial_field_values=initfield;
  initial_field_derivs= initderivs;
  
  // Set initial values of effective mass.
  mass_sq=0; // Assuming inflaton mass to be zero
  
  // Set initial value of Hubble constant
  
  hubble_init = sqrt((0.5*pw2(initial_field_derivs) + potential_energy_hom(initial_field_values))/3.);
  
  if(!(hubble_init>=0.)) // Make sure Hubble isn't negative or undefined
  {
    printf("Error in calculating initial Hubble constant. Exiting.\n");
    exit(1);
  }
  ad=hubble_init;
  //adGF=hubble_init;
  ad0=hubble_init;
  
  // Loop over gridpoints.
  // px and py go over all mode values in wrap-around order, rising from 0 to N/2 and then from -N/2+1 to -1
  for(i=0;i<N;i++)
  {
    px=(i<=N/2 ? i : i-N); // x-component of momentum of modes at x=i
    iconj=(i==0 ? 0 : N-i); // Index where complex conjugates of modes at x=i are stored (only used for k=0 or k=N/2)
    for(j=0;j<N;j++)
    {
      py=(j<=N/2 ? j : j-N); // y-component of momentum of modes at y=j
      
      // Set all modes with 0<k<N/2. The complex conjugates of these modes do not need to be set.
      for(k=1;k<N/2;k++)
      {
        pz=k; // z-component of momentum of modes at z=k
        p2lat=dp2*(pw2(px)+pw2(py)+pw2(pz)); // Total lattice momentum squared
        p2 = 4.*pw2(N/L)*(pw2(sin(px*pi/N))+pw2(sin(py*pi/N))+pw2(sin(pz*pi/N))); //physical momentum defined in 2209.13616
        set_mode(p2,p2lat,mass_sq,&f[i][j][2*k],&fd[i][j][2*k],0); // Set mode
        
      }
      
      // Set modes with k=0 or N/2.
      if(j>N/2 || (i>N/2 && (j==0 || j==N/2))) // The complex conjugates of these modes appear explicitly on the lattice and must be set to satisfy f(-p)=f(p)*
      {
        jconj=(j==0 ? 0 : N-j); // Index where complex conjugates of modes at y=j are stored
        // k=0
        p2lat=dp2*(pw2(px)+pw2(py)); // Total momentum squared
        p2 = 4.*pw2(N/L)*(pw2(sin(px*pi/N))+pw2(sin(py*pi/N)));
        set_mode(p2,p2lat,mass_sq,&f[i][j][0],&fd[i][j][0],0); // Set mode
        
        f[iconj][jconj][0]=f[i][j][0]; // Set complex conjugate mode
        f[iconj][jconj][1]=-f[i][j][1];
        fd[iconj][jconj][0]=fd[i][j][0];
        fd[iconj][jconj][1]=-fd[i][j][1];
        
        // k=N/2
        p2lat=dp2*(pw2(px)+pw2(py)+pw2(N/2)); // Total momentum squared
        p2 = 4.*pw2(N/L)*(pw2(sin(px*pi/N))+pw2(sin(py*pi/N))+1.);
        set_mode(p2,p2lat,mass_sq,&fnyquist_p[i][2*j],&fdnyquist_p[i][2*j],0); // Set mode
        
        fnyquist_p[iconj][2*jconj]=fnyquist_p[i][2*j]; // Set complex conjugate mode
        fnyquist_p[iconj][2*jconj+1]=-fnyquist_p[i][2*j+1];
        fdnyquist_p[iconj][2*jconj]=fdnyquist_p[i][2*j];
        fdnyquist_p[iconj][2*jconj+1]=-fdnyquist_p[i][2*j+1];
        
      }
      else if((i==0 || i==N/2) && (j==0 || j==N/2)) // The 8 "corners" of the lattice are set to real values
      {
        p2lat=dp2*(pw2(px)+pw2(py)); // Total momentum squared for k=0
        p2 = 4.*pw2(N/L)*(pw2(sin(px*pi/N))+pw2(sin(py*pi/N)));
        if(p2>0.) // Don't set the zeromode here (see below)
        set_mode(p2,p2lat,mass_sq,&f[i][j][0],&fd[i][j][0],1);
        
        p2lat=dp2*(pw2(px)+pw2(py)+pw2(N/2)); // Total momentum squared for k=N/2
        p2 = 4.*pw2(N/L)*(pw2(sin(px*pi/N))+pw2(sin(py*pi/N))+1.);
        set_mode(p2,p2lat,mass_sq,&fnyquist_p[i][2*j],&fdnyquist_p[i][2*j],1); // The last argument specifies that a real value should be set
        
      }
    } // End of loop over j (y-index on lattice)
  } // End of loop over i (x-index on lattice)
  f[0][0][0]=0.; // Set zeromode of field and derivative to zero (it gets added in position space)
  f[0][0][1]=0.;
  fd[0][0][0]=0.;
  fd[0][0][1]=0.;
  
  // *DEBUG* The option to not use the FFTs is for debugging purposes only. For actual simulations seed should always be positive
  if(seed>=0)
  {
    fftrn((float *)f,(float *)fnyquist_p,3,arraysize,-1); // Inverse Fourier transform of field
    fftrn((float *)fd,(float *)fdnyquist_p,3,arraysize,-1); // Inverse Fourier transform of field derivatives
  }
  
  // Add zeromode
  LOOP
  {
    f[i][j][k] += initial_field_values;
    fd[i][j][k] += initial_field_derivs;
  }
  
  float vard=0;
  float deriv_energy_in=0.;
  
  //Update the hubble constant to include all terms
  LOOP
  {
    vard += pw2(fd[i][j][k]); // Sum over squared field derivatives
  }
  // Convert sums to averages
  vard = vard/(float)gridsize;
  deriv_energy_in = .5*vard; // Extra terms account for model-dependent rescalings
  
  // Calculate and output gradient energy
  
  hubble_init = sqrt((deriv_energy_in + potential_energy() + gradient_energy())/3.);
  
  if(!(hubble_init>=0.)) // Make sure Hubble isn't negative or undefined
  {
    printf("Error in calculating initial Hubble constant. Exiting.\n");
    exit(1);
  }
  
  ad=hubble_init;
  ad0=hubble_init;
  
  printf("Finished initial conditions\n");
  
  return;
}

#if perform_deltaN

void initializeN()
{
  //Initialize the grid used in the delta N calculation, convert derivatives in e-folds time
  DECLARE_INDICES
  float H_lat=0;
  LOOP
  {
    #if numerical_potential
    lstart[i][j][k] = lstart[i][j][k]-100;
    #endif
    H_lat = sqrt((potential_energy_hom(f[i][j][k]))/3.);
    fd[i][j][k] = fd[i][j][k]*pow(a,rescale_s-1)/H_lat;
    deltaN[i][j][k] = 0.;
  }
  
  return;
}

#endif