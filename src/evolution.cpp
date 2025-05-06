// evolution.cpp - Core evolution algorithm for lattice fields

#include "main.h"

// Increments a grid location accounting for periodic wrapping
inline int INCREMENT(int i)
{
  return( (i==N-1) ? 0 : i+1 );
}

// Decrements a grid location accounting for periodic wrapping
inline int DECREMENT(int i)
{
  return( (i==0) ? N-1 : i-1 );
}

// Calculate the Laplacian of a field point in the bulk. (The result must be divided by dx^2 to give a Laplacian.)
inline float lapld(INDEXLIST)
{
  if(i==0 || j==0 || k==0 || i==N-1 || j==N-1 || k==N-1)
  return (fd[i][j][INCREMENT(k)] + fd[i][j][DECREMENT(k)]
  +fd[i][INCREMENT(j)][k] + fd[i][DECREMENT(j)][k]
  +fd[INCREMENT(i)][j][k] + fd[DECREMENT(i)][j][k]
  -6.*fd[i][j][k]);
  
  else return (fd[i][j][k+1] + fd[i][j][k-1]
  +fd[i][j+1][k] + fd[i][j-1][k]
  +fd[i+1][j][k] + fd[i-1][j][k]
  -6.*fd[i][j][k]);
}

inline float lapl(INDEXLIST)
{
  if(i==0 || j==0 || k==0 || i==N-1 || j==N-1 || k==N-1)
  return (f[i][j][INCREMENT(k)] + f[i][j][DECREMENT(k)]
  +f[i][INCREMENT(j)][k] + f[i][DECREMENT(j)][k]
  +f[INCREMENT(i)][j][k] + f[DECREMENT(i)][j][k]
  -6.*f[i][j][k]);
  
  else return (f[i][j][k+1] + f[i][j][k-1]
  +f[i][j+1][k] + f[i][j-1][k]
  +f[i+1][j][k] + f[i-1][j][k]
  -6.*f[i][j][k]);
}

/////////////////////////////////////////////////////
// Externally called function(s)
/////////////////////////////////////////////////////

float gradient_energy()
{
  int i=0,j=0,k=0;
  double gradient=0.;
  float norm=pw2(1./(a*dx));
  LOOP
  gradient -= f[i][j][k]*lapl(i,j,k);
  return(.5*gradient*norm/(float)gridsize);
}

float kin_energy()
{
  int i=0,j=0,k=0;
  double deriv_energy=0.;
  LOOP
  {
    deriv_energy += pw2(fd[i][j][k]);
  }
  deriv_energy = deriv_energy/(float)gridsize;
  return(.5*pow(a,2.*rescale_s-2.)*deriv_energy);
}

void evolve_derivs(float d)
{
  int i=0,j=0,k=0;
  float laplnorm = pow(a,-2.*rescale_s)/pw2(dx); // Norm of laplacian
  float sfev1=rescale_s+1.; // See documentation of LATTICEEASY for an explanation of these terms in the evolution equation for a
  float sfev2=-2.*rescale_s+2.;
  
  ad2 = (-2.*ad - 2.*a/d/sfev1*(1.-sqrt(1.+2.*d*sfev1*ad/a+pw2(d)*sfev1*pow(a,sfev2)*(2.*gradient_energy()/3.+potential_energy()))))/d;
  ad += .5*d*ad2; // Advance ad to time t for field evolution equations
  
  #if parallel_calculation
  #pragma omp parallel for collapse(3)
  #endif
  LOOP
  fd[i][j][k] += d*(laplnorm*lapl(i,j,k)-(2.+rescale_s)*ad*fd[i][j][k]/a-pow(a,2.-2.*rescale_s)*dvdf(i,j,k));
  
  ad += .5*d*ad2; // Advance ad to time t for field evolution equations
}

void evolve_fields(float d)
{
  DECLARE_INDICES
  
  t += d;
  
  #if parallel_calculation
  #pragma omp parallel for collapse(3)
  #endif
  LOOP
  f[i][j][k] += d*fd[i][j][k];
  
  a += d*ad;
  
}
#if perform_deltaN
//The following functions are only needed for the deltaN calculation
void evolve_derivsN(float d)
{
  int i=0,j=0,k=0;
  
  #if parallel_calculation
  #pragma omp parallel for collapse(3)
  #endif
  LOOP
  fd[i][j][k] += d*(-(3-0.5*pw2(fd[i][j][k]))*(fd[i][j][k] + pot_ratio(i,j,k)));
}

void evolve_fieldsN(float d)
{
  DECLARE_INDICES
  
  t += d;
  
  #if parallel_calculation
  #pragma omp parallel for collapse(3)
  #endif
  LOOP
  {
    if(f[i][j][k]>phiref)
    deltaN[i][j][k] += d;
    f[i][j][k] += d*fd[i][j][k];
  }
  
  a += d*ad;
  
}
#endif