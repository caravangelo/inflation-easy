#include <filesystem>
using namespace std::filesystem;
/*
This file contains functions for computing and writing simulation outputs,
such as field means, variances, spectra, and energy densities.
Output files are saved in the 'results/' directory with filenames based on the extension 'ext_'.
*/

#include "main.h"

char name_[550]; // Filenames - set differently by each function to open output files

// Zeroes out a mode and its derivative (used when applying a cutoff)
void kill_mode(float *field, float *deriv)
{
    field[0] = 0.;
    field[1] = 0.;
    deriv[0] = 0.;
    deriv[1] = 0.;

    return;
}

// Computes spatial averages and variances of the field and its derivative
void meansvars(int flush)
{
    static FILE *means_, *vars_, *velocity_;
    DECLARE_INDICES

    float av, var, vel;

    static int first = 1;
    if (first) // Open output files
    {
        snprintf(name_, sizeof(name_), "results/means%s", ext_);
        means_ = fopen(name_, mode_);
        snprintf(name_, sizeof(name_), "results/variance%s", ext_);
        vars_ = fopen(name_, mode_);
        snprintf(name_, sizeof(name_), "results/velocity%s", ext_);
        velocity_ = fopen(name_, mode_);
        first = 0;
    }

    fprintf(means_, "%f", t);
    fprintf(means_, " %e", a);
    fprintf(velocity_, "%f", t);
    fprintf(velocity_, " %e", a);
    fprintf(vars_, "%f", t);
    fprintf(vars_, " %e", a);

    av = 0.;
    vel = 0.;
    var = 0.;
    // Calculate field mean
    LOOP
    {
        av += f[idx(i,j,k)];
        vel += fd[idx(i,j,k)];
        var += pw2(f[idx(i,j,k)]);
    }
    av = av / (float)gridsize; // Convert sum to average
    vel = vel / (float)gridsize;

    vel = vel * pow(a, rescale_s - 1) * rescale_B;

    fprintf(means_, " %e", av);
    fprintf(velocity_, " %e", vel);
    fprintf(vars_, " %e", var - pw2(av));
    // Check for instability. See if the field has grown exponentially and become non-numerical at any point.
    if (av + FLT_MAX == av || (av != 0. && av / av != 1.))
    {
        printf("Unstable solution developed. Scalar field not numerical at t=%f\n", t);
        output_parameters();
        fflush(means_);
        fflush(vars_);
        exit(1);
    }

    fprintf(means_, "\n");
    fprintf(vars_, "\n");
    fprintf(velocity_, "\n");
    if (flush)
    {
        fflush(means_);
        fflush(vars_);
        fflush(velocity_);
    }
}

// Outputs the time and the physical quantities a, adot/a (i.e. Hubble), and adotdot
void scale(int flush)
{
    static FILE *sf_;

    static int first = 1;
    if (first) // Open output file
    {
        snprintf(name_, sizeof(name_), "results/sf%s", ext_);
        sf_ = fopen(name_, mode_);
        first = 0;
    }

    // Output a, H, and adotdot in physical units using rescalings
    fprintf(sf_, "%f %f %e %e\n",
            t, a,
            ad * rescale_B * pow(a, rescale_s - 2.),
            pw2(rescale_B) * pow(a, 2. * rescale_s - 2.) *
                (ad2 + (rescale_s - 1.) * pw2(ad) / a));

    if (flush)
        fflush(sf_);
}

// Outputs power spectrum of the field and applies a high-momentum cutoff if enabled
void spectraf()
{
    static FILE *spectra_, *spectratimes_; // Output files for power spectra and times at which spectra were taken
    const int maxnumbins = (int)(1.73205 * (N / 2)) + 1; // Number of bins (bin spacing=lattice spacing in Fourier space) = sqrt(NDIMS)*(N/2)+1. Set for 3D.
    int numpoints[maxnumbins]; // Number of points in each momentum bin
    float p[maxnumbins], f2[maxnumbins]; // Values for each bin: Momentum, |f_k|^2
    int numbins = (int)(sqrt(3.) * (N / 2)) + 1; // Actual number of bins for the number of dimensions
    float pmagnitude; // Total momentum (p) in units of lattice spacing, pmagnitude = Sqrt(px^2+py^2+pz^2)
    float pdisc;
    float dp = 2. * pi / L; // Size of grid spacing in momentum space
    float fp2; // Square magnitude of field (fp2) for a given mode
    int i, j, k, px, py, pz, iconj, jconj;
    float norm1 = pow(L / rescale_B, 3) / pow(N, 6); // Normalization to reduced Planck mass units

    int arraysize[] = {N, N, N};

    static int first = 1;
    if (first)
    {
        snprintf(name_, sizeof(name_), "results/spectra%s", ext_);
        spectra_ = fopen(name_, mode_);

        snprintf(name_, sizeof(name_), "results/spectratimes%s", ext_);
        spectratimes_ = fopen(name_, mode_);
        first = 0;
    }

    for (i = 0; i < numbins; i++)
        p[i] = dp * i;

    for (i = 0; i < numbins; i++)
    {
        numpoints[i] = 0;
        f2[i] = 0.;
    }

    fftrn(f.data(), (float *)fnyquist_p, 3, arraysize, 1);

    for (i = 0; i < N; i++)
    {
        px = (i <= N / 2 ? i : i - N);
        for (j = 0; j < N; j++)
        {
            py = (j <= N / 2 ? j : j - N);
            for (k = 1; k < N / 2; k++)
            {
                pz = k;
                pmagnitude = sqrt(pw2(px) + pw2(py) + pw2(pz));
                fp2 = pw2(f[idx(i,j,2 * k)]) + pw2(f[idx(i,j,2 * k + 1)]);
                numpoints[(int)pmagnitude] += 2;
                f2[(int)pmagnitude] += 2. * fp2;
            }
            for (k = 0; k <= N / 2; k += N / 2)
            {
                pz = k;
                pmagnitude = sqrt(pw2(px) + pw2(py) + pw2(pz));
                if (k == 0)
                {
                    fp2 = pw2(f[idx(i,j,0)]) + pw2(f[idx(i,j,1)]);
                }
                else
                {
                    fp2 = pw2(fnyquist_p[i][2 * j]) + pw2(fnyquist_p[i][2 * j + 1]);
                }
                numpoints[(int)pmagnitude]++;
                f2[(int)pmagnitude] += fp2;
            }
        }
    }

    for (i = 0; i < numbins; i++)
    {
        if (numpoints[i] > 0)
        {
            f2[i] = f2[i] / numpoints[i];
        }
        fprintf(spectra_, "%e %d %e\n", p[i], numpoints[i], norm1 * f2[i]);
    }

    if (high_cutoff_index > 0 && forcing_cutoff)
    {
        fftrn(fd.data(), (float *)fdnyquist_p, 3, arraysize, 1);
        for (i = 0; i < N; i++)
        {
            px = (i <= N / 2 ? i : i - N);
            iconj = (i == 0 ? 0 : N - i);
            for (j = 0; j < N; j++)
            {
                py = (j <= N / 2 ? j : j - N);
                for (k = 1; k < N / 2; k++)
                {
                    pz = k;
                    pdisc = sqrt(pw2(px) + pw2(py) + pw2(pz));
                    if (pdisc > high_cutoff_index || pdisc < low_cutoff_index)
                    {
                        kill_mode(&f[idx(i,j,2 * k)], &fd[idx(i,j,2 * k)]);
                    }
                }

                if (j > N / 2 || (i > N / 2 && (j == 0 || j == N / 2)))
                {
                    jconj = (j == 0 ? 0 : N - j);

                    pdisc = sqrt(pw2(px) + pw2(py));
                    if (pdisc > high_cutoff_index || pdisc < low_cutoff_index)
                    {
                        kill_mode(&f[idx(i,j,0)], &fd[idx(i,j,0)]);
                        f[idx(iconj,jconj,0)] = f[idx(i,j,0)];
                        f[idx(iconj,jconj,1)] = -f[idx(i,j,1)];
                        fd[idx(iconj,jconj,0)] = fd[idx(i,j,0)];
                        fd[idx(iconj,jconj,1)] = -fd[idx(i,j,1)];
                    }

                    pdisc = sqrt(pw2(px) + pw2(py) + pw2(N / 2.));
                    if (pdisc > high_cutoff_index || pdisc < low_cutoff_index)
                    {
                        kill_mode(&fnyquist_p[i][2 * j], &fdnyquist_p[i][2 * j]);
                        fnyquist_p[iconj][2 * jconj] = fnyquist_p[i][2 * j];
                        fnyquist_p[iconj][2 * jconj + 1] = -fnyquist_p[i][2 * j + 1];
                        fdnyquist_p[iconj][2 * jconj] = fdnyquist_p[i][2 * j];
                        fdnyquist_p[iconj][2 * jconj + 1] = -fdnyquist_p[i][2 * j + 1];
                    }
                }
                else if ((i == 0 || i == N / 2) && (j == 0 || j == N / 2))
                {
                    pdisc = sqrt(pw2(px) + pw2(py));
                    if (pdisc > high_cutoff_index || pdisc < low_cutoff_index)
                    {
                        kill_mode(&f[idx(i,j,0)], &fd[idx(i,j,0)]);
                    }

                    pdisc = sqrt(pw2(px) + pw2(py) + pw2(N / 2.));
                    if (pdisc > high_cutoff_index || pdisc < low_cutoff_index)
                    {
                        kill_mode(&fnyquist_p[i][2 * j], &fdnyquist_p[i][2 * j]);
                    }
                }
            }
        }
        fftrn(fd.data(), (float *)fdnyquist_p, 3, arraysize, -1);
    }

    fftrn(f.data(), (float *)fnyquist_p, 3, arraysize, -1);

    fprintf(spectra_, "\n");
    fflush(spectra_);
    fprintf(spectratimes_, "%f %e\n", t, a);
    fflush(spectratimes_);

    return;
}

//Outputs the 1D physical momentum, that takes into account the modified dispersion relation (see 2209.13616)
void get_modes()
{
  static FILE *modes_;
  const int maxnumbins=(int)(1.73205*(N/2))+1; // Number of bins (bin spacing=lattice spacing in Fourier space) = sqrt(NDIMS)*(N/2)+1. Set for 3D (i.e. biggest possible).
  int numpoints[maxnumbins]; // Number of points in each momentum bin
  float p_phys[maxnumbins]; // Values for each bin: Momentum, |f_k|^2, |f_k'|^2, n_k, and rho_k
  //float average=0.;
  //float zk[N][N][N]; // Values for each bin: Momentum, |f_k|^2, |f_k'|^2, n_k, and rho_k
  int numbins=(int)(sqrt(3.)*(N/2))+1; // Actual number of bins for the number of dimensions
  float pmagnitude,pphysical; // Total momentum (p) in units of lattice spacing, pmagnitude = Sqrt(px^2+py^2+pz^2). This also gives the bin index since bin spacing is set to equal lattice spacing.
  int i,j,k,px,py,pz; // px, py, and pz are components of momentum in units of grid spacing
  float norm1=rescale_B;

  static int first=1;
  if(first) // Open output files
  {

    snprintf(name_, sizeof(name_), "results/modes%s",ext_);
    modes_=fopen(name_,mode_);

    first=0;

  }

  for(i=0;i<numbins;i++) // Initialize all bins to 0
  {
    numpoints[i]=0; // Number of points in the bin
    p_phys[i]=0.; // |f_p|^2
  }
  // Loop runs over all gridpoints. All points with k<N/2 are in the array f, while points with k=N/2 are in the array fnyquist.
  // px and py go over all mode values in wrap-around order, rising from 0 to N/2 and then from -N/2+1 to -1
  for(i=0;i<N;i++)
  {
    px=(i<=N/2 ? i : i-N); // x-component of momentum of modes at x=i
    for(j=0;j<N;j++)
    {
      py=(j<=N/2 ? j : j-N); // y-component of momentum of modes at y=j
      // Modes with 0<k<N/2 are counted twice to account for the modes f(-p)=f(p)* that aren't explicitly included in the lattice
      for(k=1;k<N/2;k++)
      {
        pz=k; // z-component of momentum of modes at z=k
        pmagnitude=sqrt(pw2(px)+pw2(py)+pw2(pz)); // Magnitude of momentum of mode in units of momentum grid spacing
        pphysical=sqrt(4.*pw2(N/L)*(pw2(sin(px*pi/N))+pw2(sin(py*pi/N))+pw2(sin(pz*pi/N)))); //physical momentum defined in 2209.13616
        numpoints[(int)pmagnitude] += 2; // Iterate the count of points in this bin
        p_phys[(int)pmagnitude] += 2.*pphysical; // Add the power of this mode to the bin
      }
      // Modes with k=0 or k=N/2 are only counted once
      for(k=0;k<=N/2;k+=N/2) // "Loop" over the two values k=0 and k=N/2
      {
        pz=k;
        pmagnitude=sqrt(pw2(px)+pw2(py)+pw2(pz));
        pphysical=sqrt(4.*pw2(N/L)*(pw2(sin(px*pi/N))+pw2(sin(py*pi/N))+pw2(sin(pz*pi/N))));

        numpoints[(int)pmagnitude]++; // Iterate the count of points in this bin
        p_phys[(int)pmagnitude] += pphysical; // Add the power of this mode to the bin
      }
    } // End of loop over j
  } // End of loop over i
  for(i=0;i<numbins;i++)
  {
    if(numpoints[i]>0) // Convert sums to averages. (numpoints[i] should always be greater than zero.)
    {
      p_phys[i] = p_phys[i]/numpoints[i];
    }
    // Output the momentum, number of points, omega, and calculated spectra for each bin
    fprintf(modes_,"%e\n",norm1*p_phys[i]);
  }

  fprintf(modes_,"\n");
  fflush(modes_);

  return;
}

void bispectraf()
{
  //This function outputs the equilateral bispectrum only (see 2209.13616, Sec. 5.3.2).
  static FILE *bispectra_; // Output files for power spectra and times at which spectra were taken
  const int maxnumbins=(int)(1.73205*(N/2))+1; // Number of bins (bin spacing=lattice spacing in Fourier space) = sqrt(NDIMS)*(N/2)+1. Set for 3D (i.e. biggest possible).
  int numpoints[maxnumbins]; // Number of points in each momentum bin
  float p[maxnumbins],bisreal[maxnumbins],bisimag[maxnumbins]; // Values for each bin: Momentum, |f_k|^2, |f_k'|^2, n_k, and rho_k
  int numbins=(int)(sqrt(3.)*(N/2))+1; // Actual number of bins for the number of dimensions
  float dp=2.*pi/L; // Size of grid spacing in momentum space
  float mean=0.;
  int i1,j1,k1; // px, py, and pz are components of momentum in units of grid spacing
  int i2,j2,k2,i2n,j2n;
  int i3,j3,k3;
  float px1,py1,px2,py2,px3,py3;
  int i,j,k;
  float pzaus1,pzaus2,counts;
  float f1r,f1i,f2r,f2i,f3r,f3i;
  float norm1=pow(L/rescale_B,6)/pow(N,9);
  int arraysize[]={N,N,N}; // Array of grid size in all dimensions - used by FFT routine
  float kf;

  LOOP
  mean += f[idx(i,j,k)];

  mean = mean/(float)gridsize;

  LOOP
  f[idx(i,j,k)] -= mean;

  static int first=1;
  if(first) // Open output files
  {
    snprintf(name_, sizeof(name_), "results/bispectra%s",ext_);
    bispectra_=fopen(name_,mode_);
    first=0;
  }

  // Calculate magnitude of momentum in each bin
  for(k=0;k<numbins;k++)
    p[k]=dp*k;

  for(k=0;k<numbins;k++) // Initialize all bins to 0
  {
    numpoints[k]=0; // Number of points in the bin
    bisreal[k]=0.; // |f_p|^2
    bisimag[k]=0.; // |f_p|^2
  }

  fftrn(f.data(), (float *)fnyquist_p, 3, arraysize, 1); // Transform field values to Fourier space

  for(k=0;k<numbins;k++)
  {
    kf = k;
    for(i1=0;i1<N;i1++) for(j1=0;j1<N;j1++)
    {//for1
    px1=(i1<=N/2 ? i1 : i1-N);
    py1=(j1<=N/2 ? j1 : j1-N);
    pzaus1 = pw2(kf)-pw2(px1)-pw2(py1);
    if(pzaus1>=0)
    for(k1=(int)round(sqrt(pzaus1))-1;k1 < (int)round(sqrt(pzaus1))+2;k1++)
    if(abs(sqrt(pw2(px1)+pw2(py1)+pw2(k1))-kf) < 1.5 && k1 <= N/2 && k1 >= 0)
    for(i2=0;i2<N;i2++) for(j2=0;j2<N;j2++)
    {//for2
    px2=(i2<=N/2 ? i2 : i2-N);
    py2=(j2<=N/2 ? j2 : j2-N);
    pzaus2 = pw2(kf)-pw2(px2)-pw2(py2);
    if(pzaus2>=0)
    for(k2=(int)round(sqrt(pzaus2))-1;k2 < (int)round(sqrt(pzaus2))+2;k2++)
    if(abs(sqrt(pw2(px2)+pw2(py2)+pw2(k2))-kf) < 1.5 && k2 <= N/2 && k2 >= 0)
    {//if triangleapprox2
    px3 = px1 + px2;
    py3 = py1 + py2;
    k3 = k1 + k2;

    if(px3 <= N/2 and px3 > -N/2  and py3 <= N/2 and py3 > -N/2)
    {//if in lattice
    i3=(px3>=0 ? px3 : px3+N);
    j3=(py3>=0 ? py3 : py3+N);
    if(abs(sqrt(pw2(px3)+pw2(py3)+pw2(k3))-kf)< 1.5 && k3 <= N/2)
    {
      f1r = f[idx(i1,j1,2*k1)];
      f1i = f[idx(i1,j1,2*k1+1)];
      f2r = f[idx(i2,j2,2*k2)];
      f2i = f[idx(i2,j2,2*k2+1)];
      f3r = f[idx(i3,j3,2*k3)];
      f3i = f[idx(i3,j3,2*k3+1)];
      counts = 1.;
      if( k1 == int(N/2))
      {
        f1r = fnyquist_p[i1][2*j1];
        f1i = fnyquist_p[i1][2*j1+1];
      }
      if( k2 == int(N/2))
      {
        f2r = fnyquist_p[i2][2*j2];
        f2i = fnyquist_p[i2][2*j2+1];
      }
      if( k3 == int(N/2))
      {
        f3r = fnyquist_p[i3][2*j3];
        f3i = fnyquist_p[i3][2*j3+1];
      }
      if(k1 != (int)N/2 and k2 != (int)N/2 and (k1 != 0 or k2 != 0))
      counts = 2.;

      numpoints[k] += counts;
      bisreal[k] += counts*(f1r*f2r*f3r-f1i*f2i*f3r+f1i*f2r*f3i+f2i*f1r*f3i);
      bisimag[k] += counts*(-f1r*f2r*f3i+f1i*f2i*f3i+f1i*f2r*f3r+f2i*f1r*f3r);

    }//ifpzaus3 part 1
    if(k1 != 0 && k2!= 0)
    {
      k3 = k1-k2;
      if(k3>=0)
      if(abs(sqrt(pw2(px3)+pw2(py3)+pw2(k3))-kf) < 1.5 && k3 <= N/2)
      {
        i2n=(-px2 >=0 ? -px2 : -px2+N);
        j2n=(-py2 >=0 ? -py2 : -py2+N);

        f1r = f[idx(i1,j1,2*k1)];
        f1i = f[idx(i1,j1,2*k1+1)];
        f2r = f[idx(i2n,j2n,2*k2)];
        f2i = -f[idx(i2n,j2n,2*k2+1)];
        f3r = f[idx(i3,j3,2*k3)];
        f3i = f[idx(i3,j3,2*k3+1)];

        if( k1 == int(N/2))
        {
          f1r = fnyquist_p[i1][2*j1];
          f1i = fnyquist_p[i1][2*j1+1];
        }
        if( k2 == int(N/2))
        {
          f2r = fnyquist_p[i2n][2*j2n];
          f2i = -fnyquist_p[i2n][2*j2n+1];
        }
        if( k3 == int(N/2))
        {
          f3r = fnyquist_p[i3][2*j3];
          f3i = fnyquist_p[i3][2*j3+1];
        }

        numpoints[k] += 2.;
        bisreal[k] += 2.*(f1r*f2r*f3r-f1i*f2i*f3r+f1i*f2r*f3i+f2i*f1r*f3i);
        bisimag[k] += 2.*(-f1r*f2r*f3i+f1i*f2i*f3i+f1i*f2r*f3r+f2i*f1r*f3r);

      }//ifpzaus3 part 2

    }//if k1 and k2 neq0
  }//endif inda lattice
}//if triangleapprox2
}//for2
}//for1
//printf("%d, BISPECTRUM COUNTS = %d\n",k,numpoints[k]);
}//forkv

for(k=0;k<numbins;k++)
{
if(numpoints[k]>0) // Convert sums to averages.
{
bisreal[k] = bisreal[k]/numpoints[k];
bisimag[k] = bisimag[k]/numpoints[k];
}

fprintf(bispectra_,"%e %d %e %e\n",
p[k],numpoints[k],norm1*bisreal[k],norm1*bisimag[k]);
}

fftrn(f.data(), (float *)fnyquist_p, 3, arraysize, -1);

LOOP
f[idx(i,j,k)] += mean;

fprintf(bispectra_,"\n");
fflush(bispectra_);

return;
}

void box()
{
    static FILE *box_;
    int i,j,k;
    static int first=1;
    if(first) // Open output files
    {
        snprintf(name_, sizeof(name_), "results/box%s",ext_);
        box_ = fopen(name_,mode_);
        first=0;
    }

    for(i=0;i<N;i++) for(j=0;j<N;j++) for(k=0;k<N;k++)
    {
        fprintf(box_,"%.17g\n",f[idx(i,j,k)]);
    }
    fprintf(box_,"\n");
    fflush(box_);
}

void box2d()
{
    static FILE *snapshots_2d_phi_;
    int i, j, k;
    if(N < 64)
        i = 5;
    else
        i = 50;

    static int first=1;
    if(first) // Open output files
    {
        snprintf(name_, sizeof(name_), "results/snapshots_2d_phi%s",ext_);
        snapshots_2d_phi_ = fopen(name_,mode_);
        first=0;
    }

    for(j=0;j<N;j++) for(k=0;k<N;k++)
    {
        fprintf(snapshots_2d_phi_,"%.17g\n",f[idx(i,j,k)]);
    }
    fprintf(snapshots_2d_phi_,"\n");
    fflush(snapshots_2d_phi_);
}

void box2dot()
{
    static FILE *snapshots_2d_phidot_;
    int i, j, k;
    if(N < 64)
        i = 5;
    else
        i = 50;

    static int first=1;
    if(first) // Open output files
    {
        snprintf(name_, sizeof(name_), "results/snapshots_2d_phidot%s",ext_);
        snapshots_2d_phidot_ = fopen(name_,mode_);
        first=0;
    }

    for(j=0;j<N;j++) for(k=0;k<N;k++)
    {
        fprintf(snapshots_2d_phidot_,"%.17g\n",fd[idx(i,j,k)] * rescale_B);
    }
    fprintf(snapshots_2d_phidot_,"\n");
    fflush(snapshots_2d_phidot_);
}

void energy()
{
    static FILE *energy_, *conservation_;
    float deriv_energy, grad_energy, pot_energy;

    float totalE = 0.;
    static int first = 1;
    if (first) // Open output files
    {
        snprintf(name_, sizeof(name_), "results/energy%s", ext_);
        energy_ = fopen(name_, mode_);

        snprintf(name_, sizeof(name_), "results/conservation%s", ext_);
        conservation_ = fopen(name_, mode_);

        first = 0;
    }

    fprintf(energy_, "%f", t); // Output time
    fprintf(energy_, " %e", a);

    // Calculate and output kinetic (time derivative) energy
    deriv_energy = kin_energy();
    totalE += deriv_energy;
    fprintf(energy_, " %e", deriv_energy * rescale_B);

    // Calculate and output gradient energy
    grad_energy = gradient_energy();
    totalE += grad_energy;
    fprintf(energy_, " %e", grad_energy * rescale_B);

    // Calculate and output potential energy
    pot_energy = potential_energy();
    totalE += pot_energy;
    fprintf(energy_, " %e", pot_energy * rescale_B);

    fprintf(energy_, "\n");
    fflush(energy_);

    // Energy conservation
    fprintf(conservation_, "%e %e %e\n",
            t, a, 3. * pow(a, 2. * rescale_s - 4.) * pw2(ad) / (totalE));
    fflush(conservation_);
}

void readable_time(int t, FILE *info_)
{
    int tminutes = 60, thours = 60 * tminutes, tdays = 24 * thours;

    if (t == 0)
    {
        fprintf(info_, "less than 1 second\n");
        return;
    }

    // Days
    if (t > tdays)
    {
        fprintf(info_, "%d days", t / tdays);
        t = t % tdays;
        if (t > 0)
            fprintf(info_, ", ");
    }
    // Hours
    if (t > thours)
    {
        fprintf(info_, "%d hours", t / thours);
        t = t % thours;
        if (t > 0)
            fprintf(info_, ", ");
    }
    // Minutes
    if (t > tminutes)
    {
        fprintf(info_, "%d minutes", t / tminutes);
        t = t % tminutes;
        if (t > 0)
            fprintf(info_, ", ");
    }
    // Seconds
    if (t > 0)
        fprintf(info_, "%d seconds", t);

    fprintf(info_, "\n");
    return;
}

void histograms()
{
    static FILE *histogram_, *histogramtimes_;
    int i = 0, j = 0, k = 0;
    int binnum; // Index of bin for a given field value
    float binfreq[nbins]; // The frequency of field values occurring within each bin
    float bmin, bmax, df; // Minimum and maximum field values for each field and spacing (in field values) between bins
    int numpts; // Count the number of points in the histogram for each field. (Should be all lattice points unless explicit field limits are given.)

    static int first = 1;
    if (first) // Open output files
    {
        snprintf(name_, sizeof(name_), "results/histogram%s", ext_);
        histogram_ = fopen(name_, mode_);

        snprintf(name_, sizeof(name_), "results/histogramtimes%s", ext_);
        histogramtimes_ = fopen(name_, mode_);
        first = 0;
    }

    fprintf(histogramtimes_, "%f", t); // Output time at which histograms were recorded
    fprintf(histogramtimes_, " %e", a);

    i = 0;
    j = 0;
    k = 0;
    bmin = f[idx(i,j,k)];
    bmax = bmin;
    LOOP
    {
        bmin = (f[idx(i,j,k)] < bmin ? f[idx(i,j,k)] : bmin);
        bmax = (f[idx(i,j,k)] > bmax ? f[idx(i,j,k)] : bmax);
    }

    // Find the difference (in field value) between successive bins
    df = (bmax - bmin) / (float)(nbins); // bmin will be at the bottom of the first bin and bmax at the top of the last

    // Initialize all frequencies to zero
    for (i = 0; i < nbins; i++)
        binfreq[i] = 0.;

    // Iterate over grid to determine bin frequencies
    numpts = 0;
    LOOP
    {
        binnum = (int)((f[idx(i,j,k)] - bmin) / df); // Find index of bin for each value
        if (f[idx(i,j,k)] == bmax) // The maximal field value is at the top of the highest bin
            binnum = nbins - 1;
        if (binnum >= 0 && binnum < nbins) // Increment frequency in the appropriate bin
        {
            binfreq[binnum]++;
            numpts++;
        }
    } // End of loop over grid

    // Output results
    for (i = 0; i < nbins; i++)
        fprintf(histogram_, "%e\n", binfreq[i] / (float)numpts); // Output bin frequency normalized so the total equals 1
    fprintf(histogram_, "\n"); // Stick a blank line between times to make the file more readable
    fflush(histogram_);
    fprintf(histogramtimes_, " %e %e", bmin, df); // Output the starting point and stepsize for the bins at each time

    fprintf(histogramtimes_, "\n");
    fflush(histogramtimes_);
}

#if perform_deltaN

void spectraN()
{
    static FILE *spectraN_; // Output file for power spectrum
    const int maxnumbins = (int)(1.73205 * (N / 2)) + 1; // Number of bins (bin spacing=lattice spacing in Fourier space)
    int numpoints[maxnumbins]; // Number of points in each momentum bin
    float p[maxnumbins], f2[maxnumbins]; // Values for each bin: Momentum, |f_k|^2
    int numbins = (int)(sqrt(3.) * (N / 2)) + 1; // Actual number of bins for the number of dimensions
    float pmagnitude; // Total momentum (p) in units of lattice spacing
    float dp = 2. * pi / L; // Grid spacing in momentum space
    float fp2; // Square magnitude of field mode
    int i, j, k, px, py, pz;
    float norm1 = pow(L / rescale_B, 3) / pow(N, 6);
    int arraysize[] = {N, N, N}; // Array of grid size for FFT routine

    static int first = 1;
    if (first) // Open output file
    {
        snprintf(name_, sizeof(name_), "results/spectraN%s", ext_);
        spectraN_ = fopen(name_, mode_);
        first = 0;
    }

    // Initialize bin values
    for (i = 0; i < numbins; i++)
    {
        p[i] = dp * i;
        numpoints[i] = 0;
        f2[i] = 0.;
    }

    // Transform deltaN to Fourier space
    fftrn(deltaN.data(), (float *)fnyquist_p, 3, arraysize, 1);

    // Loop over lattice grid points
    for (i = 0; i < N; i++)
    {
        px = (i <= N / 2 ? i : i - N);
        for (j = 0; j < N; j++)
        {
            py = (j <= N / 2 ? j : j - N);

            // Modes with 0 < k < N/2 are counted twice
            for (k = 1; k < N / 2; k++)
            {
                pz = k;
                pmagnitude = sqrt(pw2(px) + pw2(py) + pw2(pz));
                fp2 = pw2(deltaN[idx(i,j,2 * k)]) + pw2(deltaN[idx(i,j,2 * k + 1)]);
                numpoints[(int)pmagnitude] += 2;
                f2[(int)pmagnitude] += 2. * fp2;
            }

            // Modes with k = 0 or k = N/2 are counted once
            for (k = 0; k <= N / 2; k += N / 2)
            {
                pz = k;
                pmagnitude = sqrt(pw2(px) + pw2(py) + pw2(pz));
                if (k == 0)
                {
                    fp2 = pw2(deltaN[idx(i,j,0)]) + pw2(deltaN[idx(i,j,1)]);
                }
                else
                {
                    fp2 = pw2(fnyquist_p[i][2 * j]) + pw2(fnyquist_p[i][2 * j + 1]);
                }
                numpoints[(int)pmagnitude]++;
                f2[(int)pmagnitude] += fp2;
            }
        }
    }

    // Output binned power spectrum
    for (i = 0; i < numbins; i++)
    {
        if (numpoints[i] > 0)
        {
            f2[i] = f2[i] / numpoints[i];
        }
        fprintf(spectraN_, "%e %d %e\n", p[i], numpoints[i], norm1 * f2[i]);
    }

    // Transform deltaN back to real space
    fftrn(deltaN.data(), (float *)fnyquist_p, 3, arraysize, -1);

    fprintf(spectraN_, "\n");
    fflush(spectraN_);

    return;
}

void boxN()
{
    static FILE *boxN_;
    int i, j, k;
    static int first = 1;
    if (first) // Open output file
    {
        snprintf(name_, sizeof(name_), "results/boxN%s", ext_);
        boxN_ = fopen(name_, mode_);
        first = 0;
    }

    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            for (k = 0; k < N; k++)
            {
                fprintf(boxN_, "%.17g\n", deltaN[idx(i,j,k)]);
            }
    fprintf(boxN_, "\n");
    fflush(boxN_);
}

void box2dN()
{
    static FILE *snapshots_2d_deltaN_;
    int i, j, k;
    if (N < 64)
        i = 5;
    else
        i = 50;
    static int first = 1;
    if (first) // Open output file
    {
        snprintf(name_, sizeof(name_), "results/snapshots_2d_deltaN%s", ext_);
        snapshots_2d_deltaN_ = fopen(name_, mode_);
        first = 0;
    }

    for (j = 0; j < N; j++)
        for (k = 0; k < N; k++)
        {
            fprintf(snapshots_2d_deltaN_, "%.17g\n", deltaN[idx(i,j,k)]);
        }
    fprintf(snapshots_2d_deltaN_, "\n");
    fflush(snapshots_2d_deltaN_);
}

void histogramsN()
{
    static FILE *histogramN_, *histogramtimesN_;
    int i=0, j=0, k=0;
    int binnum;
    float binfreq[nbins];
    float bmin, bmax, df;
    int numpts;

    static int first = 1;
    if (first)
    {
        snprintf(name_, sizeof(name_), "results/histogramN%s", ext_);
        histogramN_ = fopen(name_, mode_);

        snprintf(name_, sizeof(name_), "results/histogramtimesN%s", ext_);
        histogramtimesN_ = fopen(name_, mode_);
        first = 0;
    }

    fprintf(histogramtimesN_, "%f %e", t, a);

    bmin = deltaN[idx(0,0,0)];
    bmax = bmin;
    LOOP
    {
        bmin = (deltaN[idx(i,j,k)] < bmin ? deltaN[idx(i,j,k)] : bmin);
        bmax = (deltaN[idx(i,j,k)] > bmax ? deltaN[idx(i,j,k)] : bmax);
    }

    df = (bmax - bmin) / (float)(nbins);

    for (i = 0; i < nbins; i++)
        binfreq[i] = 0.;

    numpts = 0;
    LOOP
    {
        binnum = (int)((deltaN[idx(i,j,k)] - bmin) / df);
        if (deltaN[idx(i,j,k)] == bmax)
            binnum = nbins - 1;
        if (binnum >= 0 && binnum < nbins)
        {
            binfreq[binnum]++;
            numpts++;
        }
    }

    for (i = 0; i < nbins; i++)
        fprintf(histogramN_, "%e\n", binfreq[i] / (float)numpts);
    fprintf(histogramN_, "\n");
    fflush(histogramN_);

    fprintf(histogramtimesN_, " %e %e\n", bmin, df);
    fflush(histogramtimesN_);
}

void spectraLOG()
{
    static FILE *spectraLOG_; // Output files for power spectra and times at which spectra were taken
    const int maxnumbins = (int)(1.73205 * (N / 2)) + 1;
    int numpoints[maxnumbins];
    float p[maxnumbins], f2[maxnumbins];
    int numbins = (int)(sqrt(3.) * (N / 2)) + 1;
    float pmagnitude;
    float dp = 2. * pi / L;
    float fp2;
    int i, j, k, px, py, pz;
    float norm1 = pow(L / rescale_B, 3) / pow(N, 6);
    int arraysize[] = {N, N, N};

    static int first = 1;
    if (first)
    {
        snprintf(name_, sizeof(name_), "results/spectraLOG%s", ext_);
        spectraLOG_ = fopen(name_, mode_);
        first = 0;
    }

    for (i = 0; i < numbins; i++)
        p[i] = dp * i;

    for (i = 0; i < numbins; i++)
    {
        numpoints[i] = 0;
        f2[i] = 0.;
    }

    float factt = 0.;
    float fmean = 0.;

    LOOP
    {
        factt += fd[idx(i,j,k)] * pow(a, rescale_s - 1) / (ad * pow(a, rescale_s - 2.));
        fmean += f[idx(i,j,k)];
    }

    factt = factt / (float)gridsize;
    fmean = fmean / (float)gridsize;

    LOOP
    {
        deltaN[idx(i,j,k)] = 0;
        if ((1 - eta_log * (f[idx(i,j,k)] - fmean) / factt) > 0)
            deltaN[idx(i,j,k)] = 1./eta_log * log(1 - eta_log * (f[idx(i,j,k)] - fmean) / factt);
    }

    fftrn(deltaN.data(), (float *)fnyquist_p, 3, arraysize, 1);

    for (i = 0; i < N; i++)
    {
        px = (i <= N / 2 ? i : i - N);
        for (j = 0; j < N; j++)
        {
            py = (j <= N / 2 ? j : j - N);
            for (k = 1; k < N / 2; k++)
            {
                pz = k;
                pmagnitude = sqrt(pw2(px) + pw2(py) + pw2(pz));
                fp2 = pw2(deltaN[idx(i,j,2 * k)]) + pw2(deltaN[idx(i,j,2 * k + 1)]);
                numpoints[(int)pmagnitude] += 2;
                f2[(int)pmagnitude] += 2. * fp2;
            }
            for (k = 0; k <= N / 2; k += N / 2)
            {
                pz = k;
                pmagnitude = sqrt(pw2(px) + pw2(py) + pw2(pz));
                if (k == 0)
                {
                    fp2 = pw2(deltaN[idx(i,j,0)]) + pw2(deltaN[idx(i,j,1)]);
                }
                else
                {
                    fp2 = pw2(fnyquist_p[i][2 * j]) + pw2(fnyquist_p[i][2 * j + 1]);
                }
                numpoints[(int)pmagnitude]++;
                f2[(int)pmagnitude] += fp2;
            }
        }
    }

    for (i = 0; i < numbins; i++)
    {
        if (numpoints[i] > 0)
        {
            f2[i] = f2[i] / numpoints[i];
        }
        fprintf(spectraLOG_, "%e %d %e\n", p[i], numpoints[i], norm1 * f2[i]);
    }

    fftrn(deltaN.data(), (float *)fnyquist_p, 3, arraysize, -1);

    fprintf(spectraLOG_, "\n");
    fflush(spectraLOG_);

    return;
}

void histogramsLOG()
{
    static FILE *histogramLOG_, *histogramtimesLOG_;
    int i=0, j=0, k=0;
    int binnum;
    float binfreq[nbins];
    float bmin, bmax, df;
    int numpts;

    static int first = 1;
    if (first)
    {
        snprintf(name_, sizeof(name_), "results/histogramLOG%s", ext_);
        histogramLOG_ = fopen(name_, mode_);

        snprintf(name_, sizeof(name_), "results/histogramtimesLOG%s", ext_);
        histogramtimesLOG_ = fopen(name_, mode_);
        first = 0;
    }

    fprintf(histogramtimesLOG_, "%f %e", t, a);

    float factt = 0.;
    float fmean = 0.;

    LOOP
    {
        factt += fd[idx(i,j,k)] * pow(a, rescale_s - 1) / (ad * pow(a, rescale_s - 2.));
        fmean += f[idx(i,j,k)];
    }

    factt = factt / (float)gridsize;
    fmean = fmean / (float)gridsize;

    LOOP
    {
        deltaN[idx(i,j,k)] = 0;
        if ((1 - eta_log * (f[idx(i,j,k)] - fmean) / factt) > 0)
            deltaN[idx(i,j,k)] = 1./eta_log * log(1 - eta_log * (f[idx(i,j,k)] - fmean) / factt);
    }

    bmin = deltaN[idx(0,0,0)];
    bmax = bmin;
    LOOP
    {
        if (1 - eta_log * (f[idx(i,j,k)] - fmean) / factt > 0)
        {
            bmin = (deltaN[idx(i,j,k)] < bmin ? deltaN[idx(i,j,k)] : bmin);
            bmax = (deltaN[idx(i,j,k)] > bmax ? deltaN[idx(i,j,k)] : bmax);
        }
    }

    df = (bmax - bmin) / (float)(nbins);

    for (i = 0; i < nbins; i++)
        binfreq[i] = 0.;

    numpts = 0;
    LOOP
    {
        if (1 - eta_log * (f[idx(i,j,k)] - fmean) / factt > 0)
        {
            binnum = (int)((deltaN[idx(i,j,k)] - bmin) / df);
            if (deltaN[idx(i,j,k)] == bmax)
                binnum = nbins - 1;
            if (binnum >= 0 && binnum < nbins)
            {
                binfreq[binnum]++;
                numpts++;
            }
        }
    }

    for (i = 0; i < nbins; i++)
        fprintf(histogramLOG_, "%e\n", binfreq[i] / (float)numpts);
    fprintf(histogramLOG_, "\n");
    fflush(histogramLOG_);

    fprintf(histogramtimesLOG_, " %e %e\n", bmin, df);
    fflush(histogramtimesLOG_);
}

#endif

// Output information about the run parameters.
// This need only be called at the beginning and end of the run.
void output_parameters()
{
    static FILE *info_;
    static time_t tStart, tFinish; // Keep track of elapsed clock time

    static int first = 1;
    if (first) // At beginning of run output run parameters
    {
        snprintf(name_, sizeof(name_), "results/info%s", ext_);
        info_ = fopen(name_, mode_);

        fprintf(info_, "--------------------------\n");
        fprintf(info_, "General Program Information\n");
        fprintf(info_, "-----------------------------\n");
        fprintf(info_, "Grid size=%d^%d\n", N, 3);
        fprintf(info_, "L=%f\n", L);
        fprintf(info_, "f0=%f\n", initial_field);
        fprintf(info_, "fd0=%f\n", initial_derivative);
        fprintf(info_, "dt=%f, dt/dx=%f\n", dt, dt / dx);
        fprintf(info_, "rescale_s=%f\n", rescale_s);
        fprintf(info_, "rescale_B=%e\n", rescale_B);
        time(&tStart);
        fprintf(info_, "\nRun began at %s", ctime(&tStart)); // Output date in readable form
        first = 0;
    }
    else // If not at beginning record elapsed time for run
    {
        time(&tFinish);
        fprintf(info_, "Run ended at %s", ctime(&tFinish)); // Output ending date
        fprintf(info_, "\nRun from t=%f to t=%f took ", t0, t);
        readable_time((int)(tFinish - tStart), info_);
        fprintf(info_, "\n");
    }

    fflush(info_);
    return;
}

// Calculate and save quantities (means, variances, etc.). If force>0 all infrequent calculations will be performed
void save(int infrequent)
{
    if (t > 0.) // Synchronize field values and derivatives
        evolve_fields(-.5 * dt * pow(astep, rescale_s - 1));

    meansvars(infrequent);
    scale(infrequent);

    // Infrequent calculations
    if (infrequent)
    {
        if (output_box3D)
            box();
        if (output_box2D)
        {
            box2d();
            box2dot();
        }
        if (output_energy)
            energy();
        if (output_spectra)
            spectraf();
        if (output_histogram)
            histograms();
    }

    if (t > 0.) // Desynchronize field values and derivatives
        evolve_fields(.5 * dt * pow(astep, rescale_s - 1));
}

void save_last()
{
    get_modes();
    if (output_bispectrum)
        bispectraf();

#if perform_deltaN
    if (output_LOG)
    {
        //careful if you place these functions somewhere else, they use the vector deltaN as ausiliary variable
        spectraLOG();
        histogramsLOG();
    }
#endif
}

#if perform_deltaN
void saveN()
{
    if (t > 0.)
        evolve_fieldsN(-.5 * dN);

    float Nmean = 0.;
    DECLARE_INDICES

    LOOP
    {
        if (deltaN[idx(i,j,k)] <= (Nend - dN))
            Nmean += deltaN[idx(i,j,k)];
    }
    Nmean = Nmean / gridsize;

    LOOP
    {
        deltaN[idx(i,j,k)] -= Nmean;
    }

    histogramsN();

    if (output_box2D)
        box2dN();

    LOOP
    {
        if (deltaN[idx(i,j,k)] > (Nend - dN - Nmean))
            deltaN[idx(i,j,k)] = 0;
    }

    if (output_spectra)
        spectraN();

    if (t > 0.)
        evolve_fieldsN(.5 * dN);
}
#endif
