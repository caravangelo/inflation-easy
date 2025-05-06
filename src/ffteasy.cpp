/*
FFTEASY consists of the four C functions fftc1, fftcn, fftr1, and
fftrn. FFTEASY is free. I am not in any way, shape, or form expecting
to make money off of these routines. I wrote them because I needed
them for some work I was doing and I'm putting them out on the
Internet in case other people might find them useful. Feel free to
download them, incorporate them into your code, modify them, translate
the comment lines into Swahili, or whatever else you want. What I do
want is the following:
1) Leave this notice (i.e. this entire paragraph beginning with
``FFTEASY consists of...'' and ending with my email address) in with
the code wherever you put it. Even if you're just using it in-house in
your department, business, or wherever else I would like these credits
to remain with it. This is partly so that people can...
2) Give me feedback. Did FFTEASY work great for you and help your
work?  Did you hate it? Did you find a way to improve it, or translate
it into another programming language? Whatever the case might be, I
would love to hear about it. Please let me know at the email address
below.
3) Finally, insofar as I have the legal right to do so I forbid you
to make money off of this code without my consent. In other words if
you want to publish these functions in a book or bundle them into
commercial software or anything like that contact me about it
first. I'll probably say yes, but I would like to reserve that right.

For any comments or questions you can reach me at
gfelder@email.smith.edu.
*/

/* These declarations are put here so you can easily cut and paste them into your program. */
void fftc1(float f[], int N, int skip, int forward);
void fftcn(float f[], int ndims, int size[], int forward);
void fftr1(float f[], int N, int forward);
void fftrn(float f[], float fnyquist[], int ndims, int size[], int forward);

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define float double
struct complex
{
  float real;
  float imag;
};

/*
Do a Fourier transform of an array of N complex numbers separated by steps of (complex) size skip.
The array f should be of length 2N*skip and N must be a power of 2.
Forward determines whether to do a forward transform (1) or an inverse one(-1)
*/
void fftc1(float f[], int N, int skip, int forward)
{
  int b,index1,index2,trans_size,trans;
  float pi2 = 4.*asin(1.);
  float pi2n,cospi2n,sinpi2n; /* Used in recursive formulas for Re(W^b) and Im(W^b) */
  struct complex wb; /* wk = W^k = e^(2 pi i b/N) in the Danielson-Lanczos formula for a transform of length N */
  struct complex temp1,temp2; /* Buffers for implementing recursive formulas */
  struct complex *c = (struct complex *)f; /* Treat f as an array of N complex numbers */
  
  /* Place the elements of the array c in bit-reversed order */
  for(index1=1,index2=0;index1<N;index1++) /* Loop through all elements of c */
  {
    for(b=N/2;index2>=b;b/=2) /* To find the next bit reversed array index subtract leading 1's from index2 */
    index2-=b;
    index2+=b; /* Next replace the first 0 in index2 with a 1 and this gives the correct next value */
    if(index2>index1) /* Swap each pair only the first time it is found */
    {
      temp1 = c[index2*skip];
      c[index2*skip] = c[index1*skip];
      c[index1*skip] = temp1;
    }
  }
  
  /* Next perform successive transforms of length 2,4,...,N using the Danielson-Lanczos formula */
  for(trans_size=2;trans_size<=N;trans_size*=2) /* trans_size = size of transform being computed */
  {
    pi2n = forward*pi2/(float)trans_size; /* +- 2 pi/trans_size */
    cospi2n = cos(pi2n); /* Used to calculate W^k in D-L formula */
    sinpi2n = sin(pi2n);
    wb.real = 1.; /* Initialize W^b for b=0 */
    wb.imag = 0.;
    for(b=0;b<trans_size/2;b++) /* Step over half of the elements in the transform */
    {
      for(trans=0;trans<N/trans_size;trans++) /* Iterate over all transforms of size trans_size to be computed */
      {
        index1 = (trans*trans_size+b)*skip; /* Index of element in first half of transform being computed */
        index2 = index1 + trans_size/2*skip; /* Index of element in second half of transform being computed */
        temp1 = c[index1];
        temp2 = c[index2];
        c[index1].real = temp1.real + wb.real*temp2.real - wb.imag*temp2.imag; /* Implement D-L formula */
        c[index1].imag = temp1.imag + wb.real*temp2.imag + wb.imag*temp2.real;
        c[index2].real = temp1.real - wb.real*temp2.real + wb.imag*temp2.imag;
        c[index2].imag = temp1.imag - wb.real*temp2.imag - wb.imag*temp2.real;
      }
      temp1 = wb;
      wb.real = cospi2n*temp1.real - sinpi2n*temp1.imag; /* Real part of e^(2 pi i b/trans_size) used in D-L formula */
      wb.imag = cospi2n*temp1.imag + sinpi2n*temp1.real; /* Imaginary part of e^(2 pi i b/trans_size) used in D-L formula */
    }
  }
  
  /* For an inverse transform divide by the number of grid points */
  if(forward<0.)
  for(index1=0;index1<skip*N;index1+=skip)
  {
    c[index1].real /= N;
    c[index1].imag /= N;
  }
}

/*
Do a Fourier transform of an ndims dimensional array of complex numbers
Array dimensions are given by size[0],...,size[ndims-1]. Note that these are sizes of complex arrays.
The array f should be of length 2*size[0]*...*size[ndims-1] and all sizes must be powers of 2.
Forward determines whether to do a forward transform (1) or an inverse one(-1)
*/
void fftcn(float f[], int ndims, int size[], int forward)
{
  int i,j,dim;
  int planesize=1,skip=1; /* These determine where to begin successive transforms and the skip between their elements (see below) */
  int totalsize=1; /* Total size of the ndims dimensional array */
  
  for(dim=0;dim<ndims;dim++) /* Determine total size of array */
  totalsize *= size[dim];
  
  for(dim=ndims-1;dim>=0;dim--) /* Loop over dimensions */
  {
    planesize *= size[dim]; /* Planesize = Product of all sizes up to and including size[dim] */
    for(i=0;i<totalsize;i+=planesize) /* Take big steps to begin loops of transforms */
    for(j=0;j<skip;j++) /* Skip sets the number of transforms in between big steps as well as the skip between elements */
    fftc1(f+2*(i+j),size[dim],skip,forward); /* 1-D Fourier transform. (Factor of two converts complex index to float index.) */
    skip *= size[dim]; /* Skip = Product of all sizes up to (but not including) size[dim] */
  }
}

/*
Do a Fourier transform of an array of N real numbers
N must be a power of 2
Forward determines whether to do a forward transform (>=0) or an inverse one(<0)
*/
void fftr1(float f[], int N, int forward)
{
  int b;
  float pi2n = 4.*asin(1.)/N,cospi2n=cos(pi2n),sinpi2n=sin(pi2n); /* pi2n = 2 Pi/N */
  struct complex wb; /* wb = W^b = e^(2 pi i b/N) in the Danielson-Lanczos formula for a transform of length N */
  struct complex temp1,temp2; /* Buffers for implementing recursive formulas */
  struct complex *c = (struct complex *)f; /* Treat f as an array of N/2 complex numbers */
  
  if(forward==1)
  fftc1(f,N/2,1,1); /* Do a transform of f as if it were N/2 complex points */
  
  wb.real = 1.; /* Initialize W^b for b=0 */
  wb.imag = 0.;
  for(b=1;b<N/4;b++) /* Loop over elements of transform. See documentation for these formulas */
  {
    temp1 = wb;
    wb.real = cospi2n*temp1.real - sinpi2n*temp1.imag; /* Real part of e^(2 pi i b/N) used in D-L formula */
    wb.imag = cospi2n*temp1.imag + sinpi2n*temp1.real; /* Imaginary part of e^(2 pi i b/N) used in D-L formula */
    temp1 = c[b];
    temp2 = c[N/2-b];
    c[b].real = .5*(temp1.real+temp2.real + forward*wb.real*(temp1.imag+temp2.imag) + wb.imag*(temp1.real-temp2.real));
    c[b].imag = .5*(temp1.imag-temp2.imag - forward*wb.real*(temp1.real-temp2.real) + wb.imag*(temp1.imag+temp2.imag));
    c[N/2-b].real = .5*(temp1.real+temp2.real - forward*wb.real*(temp1.imag+temp2.imag) - wb.imag*(temp1.real-temp2.real));
    c[N/2-b].imag = .5*(-temp1.imag+temp2.imag - forward*wb.real*(temp1.real-temp2.real) + wb.imag*(temp1.imag+temp2.imag));
  }
  temp1 = c[0];
  c[0].real = temp1.real+temp1.imag; /* Set b=0 term in transform */
  c[0].imag = temp1.real-temp1.imag; /* Put b=N/2 term in imaginary part of first term */
  
  if(forward==-1)
  {
    c[0].real *= .5;
    c[0].imag *= .5;
    fftc1(f,N/2,1,-1);
  }
}

/*
Do a Fourier transform of an ndims dimensional array of real numbers
Array dimensions are given by size[0],...,size[ndims-1]. All sizes must be powers of 2.
The (complex) nyquist frequency components are stored in fnyquist[size[0]][size[1]]...[2*size[ndims-2]]
Forward determines whether to do a forward transform (1) or an inverse one(-1)
*/
void fftrn(float f[], float fnyquist[], int ndims, int size[], int forward)
{
  int i,j,b;
  int index,indexneg=0; /* Positions in the 1-d arrays of points labeled by indices (i0,i1,...,i(ndims-1)); indexneg gives the position in the array of the corresponding negative frequency */
  int stepsize; // Used in calculating indexneg
  int N=size[ndims-1]; /* The size of the last dimension is used often enough to merit its own name. */
  double pi2n = 4.*asin(1.)/N,cospi2n=cos(pi2n),sinpi2n=sin(pi2n); /* pi2n = 2 Pi/N */
  struct complex wb; /* wb = W^b = e^(2 pi i b/N) in the Danielson-Lanczos formula for a transform of length N */
  struct complex temp1,temp2; /* Buffers for implementing recursive formulas */
  struct complex *c = (struct complex *)f, *cnyquist = (struct complex *)fnyquist; /* Treat f and fnyquist as arrays of complex numbers */
  int totalsize=1; /* Total number of complex points in array */
  int *indices= (int *) malloc(ndims*sizeof(int)); /* Indices for looping through array */
  if(!indices) /* Make sure memory was correctly allocated */
  {
    printf("Error allocating memory in fftrn routine. Exiting.\n");
    exit(1);
  }
  
  size[ndims-1] /= 2; /* Set size[] to be the sizes of f viewed as a complex array */
  for(i=0;i<ndims;i++)
  {
    totalsize *= size[i];
    indices[i] = 0;
  }
  
  if(forward==1) /* Forward transform */
  {
    fftcn(f,ndims,size,1); /* Do a transform of f as if it were N/2 complex points */
    for(i=0;i<totalsize/size[ndims-1];i++) /* Copy b=0 data into cnyquist so the recursion formulas below for b=0 and cnyquist don't overwrite data they later need */
    cnyquist[i] = c[i*size[ndims-1]]; /* Only copy points where last array index for c is 0 */
  }
  
  for(index=0;index<totalsize;index+=size[ndims-1]) /* Loop over all but last array index */
  {
    wb.real = 1.; /* Initialize W^b for b=0 */
    wb.imag = 0.;
    for(b=1;b<N/4;b++) /* Loop over elements of transform. See documentation for these formulas */
    {
      temp1 = wb;
      wb.real = cospi2n*temp1.real - sinpi2n*temp1.imag; /* Real part of e^(2 pi i b/N_real) used in D-L formula */
      wb.imag = cospi2n*temp1.imag + sinpi2n*temp1.real; /* Imaginary part of e^(2 pi i b/N_real) used in D-L formula */
      temp1 = c[index+b];
      temp2 = c[indexneg+N/2-b]; /* Note that N-b is NOT the negative frequency for b. Only nonnegative b momenta are stored. */
      c[index+b].real = .5*(temp1.real+temp2.real + forward*wb.real*(temp1.imag+temp2.imag) + wb.imag*(temp1.real-temp2.real));
      c[index+b].imag = .5*(temp1.imag-temp2.imag - forward*wb.real*(temp1.real-temp2.real) + wb.imag*(temp1.imag+temp2.imag));
      c[indexneg+N/2-b].real = .5*(temp1.real+temp2.real - forward*wb.real*(temp1.imag+temp2.imag) - wb.imag*(temp1.real-temp2.real));
      c[indexneg+N/2-b].imag = .5*(-temp1.imag+temp2.imag - forward*wb.real*(temp1.real-temp2.real) + wb.imag*(temp1.imag+temp2.imag));
    }
    temp1 = c[index];
    temp2 = cnyquist[indexneg/size[ndims-1]]; /* Index is smaller for cnyquist because it doesn't have the last dimension */
    /* Set b=0 term in transform */
    c[index].real = .5*(temp1.real+temp2.real + forward*(temp1.imag+temp2.imag));
    c[index].imag = .5*(temp1.imag-temp2.imag - forward*(temp1.real-temp2.real));
    /* Set b=N/2 transform. */
    cnyquist[indexneg/size[ndims-1]].real = .5*(temp1.real+temp2.real - forward*(temp1.imag+temp2.imag));
    cnyquist[indexneg/size[ndims-1]].imag = .5*(-temp1.imag+temp2.imag - forward*(temp1.real-temp2.real));
    
    /* Find indices for positive and single index for negative frequency. In each dimension indexneg[j]=0 if index[j]=0, indexneg[j]=size[j]-index[j] otherwise. */
    stepsize=size[ndims-1]; /* Amount to increment indexneg by as each individual index is incremented */
    for(j=ndims-2;indices[j]==size[j]-1 && j>=0;j--) /* If the rightmost indices are maximal reset them to 0. Indexneg goes from 1 to 0 in these dimensions. */
    {
      indices[j]=0;
      indexneg -= stepsize;
      stepsize *= size[j];
    }
    if(indices[j]==0) /* If index[j] goes from 0 to 1 indexneg[j] goes from 0 to size[j]-1 */
    indexneg += stepsize*(size[j]-1);
    else /* Otherwise increasing index[j] decreases indexneg by one unit. */
    indexneg -= stepsize;
    if(j>=0) /* This avoids writing outside the array bounds on the last pass through the array loop */
    indices[j]++;
  } /* End of i loop (over total array) */
  
  if(forward==-1) /* Inverse transform */
  fftcn(f,ndims,size,-1);
  
  size[ndims-1] *= 2; /* Give the user back the array size[] in its original condition */
  free(indices); /* Free up memory allocated for indices array */
}