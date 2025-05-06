#define numerical_potential 1 // Set 1 to work with a numerical potential, 0 for analytical expression
#define perform_deltaN 1 // Set 1 to do the deltaN calculation to compute zeta
#define parallel_calculation 1 // Set 1 to do parallelize some of the loops with openmp

//Now we define the physical parameters. They differ in different models, so in the case of analytical potential we chose different values (chosen for a quadratic potential)
#if numerical_potential
const float V0 = 3e-9;//Potential parameter. Used also for rescaling.
const float initfield=2.9181235049318586;//initial field value
const float initderivs=-0.06727651095116181;//velocity in program units
const float rescale_B=sqrt(V0); //rescaling factor, see 2209.13616
const float L = 10.; //Comoving box size, to be chosen to make initial box sub-horizon
const float dt = .001; //time-step
#if perform_deltaN
const float phiref = 2.680654050913246;//final field hypersurface used for deltaN calculation. Model dependent.
const float dN = 0.0001;
const float Nend = 5.;
#endif
#else
//Values for quadratic potential
const float m2 = .263e-10;
const float initfield=14.5;
const float initderivs=-0.8157;//0.81522
const float rescale_B=sqrt(m2); //rescaling factor, see 2209.13616
const float L = 2.; //Comoving box size, to be chosen to make initial box sub-horizon
const float dt = .0001; //time-step
#if perform_deltaN
const float phiref = 13.516;//final field hypersurface used for deltaN calculation. Model dependent.
const float dN = 0.0000005;
const float Nend = .005;
#endif
#endif
const float rescale_s=0; //rescaling factor, see 2209.13616. IMPORTANT: the code has only been tested for rescale_s=0

const int N = 128; // Number of points along each edge of the cubical lattice. Should be a power of 2.

const float af = 1000.; // Final scale factor time
const int seed=8; // Random number seed. Should be a positive integer

#if numerical_potential
const int int_err=5; //interpolation errors only used for numerical potential. Leave like this unless the code gives "interpolation error", and if so try to increase this number a bit
const int int_errN=1; //interpolation errors for deltaN calculation only used for numerical potential. 
#endif

const float kcutoff=0;//Momentum for initial lowpass filter. Set to 0 to not filter
const int forcing_cutoff = 0;//set to 1 to apply a forced cutoff during the simulation with following values
const float high_cutoff_index=1000;
const float low_cutoff_index=0;

const int output_freq = 100;//output frequency
const int output_infrequent_freq =100;//output frequency of numerically expensive quantities
const int output_spectra = 1;
const int output_histogram = 1;
const int output_energy = 1;
const int output_box3D = 0;
const int output_box2D = 1;
const int output_bispectrum = 0;
#if perform_deltaN
const int output_LOG = 0;//output quantities for the log relation (see arXiv:?????)
#endif
const int screen_updates=1; // Set to 1 for time to be periodically output to screen (0 otherwise)
const int nbins=256; // Number of bins for histograms
const float histogram_min=0., histogram_max=0.; // Upper and lower limits of the histograms. To use all current field values set these two to be equal.
