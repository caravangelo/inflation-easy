// potential.cpp - Defines inflationary potential and its derivative
/*
This file contains the definition of inflationary potential. If the potential is a known analytical function, then change the last 4 function of this file acordingly (these functions are currently set for a quadratic potential), as well as the potential parameters in parameters.h
*/

#include "main.h"

#if numerical_potential

// Computes the total potential energy on the lattice
float potential_energy()
{
  DECLARE_INDICES
  int l;
  float potential=0.;
  // Loop over grid to calculate potential term
  LOOP
  {
    l=lstart[i][j][k];
    while(field_numerical[l] >= f[i][j][k])
    l++;
    
    if (field_numerical[l-1] < f[i][j][k])
    {
      printf("Interpolation Error \n");
      exit(1);
    }
    
    potential += (potential_numerical[l]+(f[i][j][k]-field_numerical[l])*(potential_numerical[l-1]-potential_numerical[l])/(field_numerical[l-1]-field_numerical[l]))/pw2(rescale_B);
    
    lstart[i][j][k] = l-int_err;
  }
  potential /= gridsize;
  
  return (potential);
}

float pot_ratio(int i, int j, int k)
{
  // This function is needed for the deltaN calculation
  int l;
  
  float pot;
  float pot_deriv;
  
  l=lstart[i][j][k];
  
  while(field_numerical[l] >= f[i][j][k])
  l++;
  
  if (field_numerical[l-1] < f[i][j][k])
  {
    printf("Interpolation Error \n");
    
    exit(1);
  }
  
  pot = potential_numerical[l]+(f[i][j][k]-field_numerical[l])*(potential_numerical[l-1]-potential_numerical[l])/(field_numerical[l-1]-field_numerical[l]);
  pot_deriv = potential_derivative_numerical[l]+(f[i][j][k]-field_numerical[l])*(potential_derivative_numerical[l-1]-potential_derivative_numerical[l])/(field_numerical[l-1]-field_numerical[l]);
  
  lstart[i][j][k] = l-int_errN;
  
  return pot_deriv/pot;
}

float dvdf(int i, int j, int k)
{
  int l = lstart[i][j][k];
  while(field_numerical[l] >= f[i][j][k])
  l++;
  return (potential_derivative_numerical[l]+(f[i][j][k]-field_numerical[l])*(potential_derivative_numerical[l-1]-potential_derivative_numerical[l])/(field_numerical[l-1]-field_numerical[l]))/pw2(rescale_B);
}

// Computes the homogeneous potential energy for a single field value (used during initialization)
float potential_energy_hom(float field_value)
{
  int l=0;
  while(field_numerical[l] >= field_value)
  l++;
  
  return (potential_numerical[l]+(field_value-field_numerical[l])*(potential_numerical[l-1]-potential_numerical[l])/(field_numerical[l-1]-field_numerical[l]))/pw2(rescale_B);
}

#else

// Computes the total potential energy on the lattice
float potential_energy()
{
  DECLARE_INDICES
  float potential=0.;
  // Loop over grid to calculate potential term (without parallelization)
  LOOP
  potential += .5*m2*pw2(f[i][j][k])/pw2(rescale_B);
  
  potential /= gridsize;
  
  return (potential);
}

float pot_ratio(int i, int j, int k)
{
  // This function is needed for the deltaN calculation
  
  float pot=.5*m2*pw2(f[i][j][k])/pw2(rescale_B);;
  float pot_deriv=dvdf(i,j,k);
  
  return pot_deriv/pot;
}

float dvdf(int i, int j, int k)
{
  return f[i][j][k]*m2/pw2(rescale_B);
}

// Computes the homogeneous potential energy for a single field value (used during initialization)
float potential_energy_hom(float field_value)
{
  return (.5*m2*pw2(field_value)/pw2(rescale_B));
}

#endif