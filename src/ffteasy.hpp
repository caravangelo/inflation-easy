/*
This code is an adaptation of FFTEASY by G. Felder
*/

#pragma once
#include <cstdlib>
#include <cstdio>
#include <cmath>

// Simple complex struct parametrized by scalar type T
template <typename T>
struct cpx {
  T real;
  T imag;
};

// ---------- Core 1D complex FFT on strided data ----------
template <typename T>
inline void fftc1(T f[], int N, int skip, int forward) {
  int b, index1, index2, trans_size, trans;
  const T pi2 = T(4) * std::asin(T(1));
  T pi2n, cospi2n, sinpi2n;
  cpx<T> wb, temp1, temp2;
  auto* c = reinterpret_cast<cpx<T>*>(f);

  // Bit-reversed permutation
  for (index1 = 1, index2 = 0; index1 < N; ++index1) {
    for (b = N/2; index2 >= b; b >>= 1) index2 -= b;
    index2 += b;
    if (index2 > index1) {
      temp1 = c[index2 * skip];
      c[index2 * skip] = c[index1 * skip];
      c[index1 * skip] = temp1;
    }
  }

  // Danielson–Lanczos iterations
  for (trans_size = 2; trans_size <= N; trans_size <<= 1) {
    pi2n = T(forward) * pi2 / T(trans_size);
    cospi2n = std::cos(pi2n);
    sinpi2n = std::sin(pi2n);
    wb.real = T(1);
    wb.imag = T(0);
    for (b = 0; b < trans_size/2; ++b) {
      for (trans = 0; trans < N/trans_size; ++trans) {
        index1 = (trans*trans_size + b) * skip;
        index2 = index1 + (trans_size/2) * skip;
        temp1 = c[index1];
        temp2 = c[index2];
        c[index1].real = temp1.real + wb.real*temp2.real - wb.imag*temp2.imag;
        c[index1].imag = temp1.imag + wb.real*temp2.imag + wb.imag*temp2.real;
        c[index2].real = temp1.real - wb.real*temp2.real + wb.imag*temp2.imag;
        c[index2].imag = temp1.imag - wb.real*temp2.imag - wb.imag*temp2.real;
      }
      temp1 = wb;
      wb.real = cospi2n*temp1.real - sinpi2n*temp1.imag;
      wb.imag = cospi2n*temp1.imag + sinpi2n*temp1.real;
    }
  }

  // Inverse normalization
  if (forward < 0) {
    for (index1 = 0; index1 < skip*N; index1 += skip) {
      c[index1].real /= T(N);
      c[index1].imag /= T(N);
    }
  }
}

// ---------- n-D complex FFT (sizes are complex sizes) ----------
template <typename T>
inline void fftcn(T f[], int ndims, int size[], int forward) {
  int i, j, dim;
  int planesize = 1, skip = 1;
  int totalsize = 1;

  for (dim = 0; dim < ndims; ++dim) totalsize *= size[dim];

  for (dim = ndims - 1; dim >= 0; --dim) {
    planesize *= size[dim];
    for (i = 0; i < totalsize; i += planesize)
      for (j = 0; j < skip; ++j)
        fftc1<T>(f + 2*(i + j), size[dim], skip, forward);
    skip *= size[dim];
  }
}

// ---------- 1D real FFT ----------
template <typename T>
inline void fftr1(T f[], int N, int forward) {
  int b;
  const T pi2n = (T(4)*std::asin(T(1))) / T(N);
  const T cospi2n = std::cos(pi2n);
  const T sinpi2n = std::sin(pi2n);
  cpx<T> wb, temp1, temp2;
  auto* c = reinterpret_cast<cpx<T>*>(f);

  if (forward == 1) fftc1<T>(f, N/2, 1, 1);

  wb.real = T(1); wb.imag = T(0);
  for (b = 1; b < N/4; ++b) {
    temp1 = wb;
    wb.real = cospi2n*temp1.real - sinpi2n*temp1.imag;
    wb.imag = cospi2n*temp1.imag + sinpi2n*temp1.real;
    temp1 = c[b];
    temp2 = c[N/2 - b];
    c[b].real     = T(0.5)*( temp1.real+temp2.real + T(forward)*wb.real*(temp1.imag+temp2.imag) + wb.imag*(temp1.real-temp2.real));
    c[b].imag     = T(0.5)*( temp1.imag-temp2.imag - T(forward)*wb.real*(temp1.real-temp2.real) + wb.imag*(temp1.imag+temp2.imag));
    c[N/2-b].real = T(0.5)*( temp1.real+temp2.real - T(forward)*wb.real*(temp1.imag+temp2.imag) - wb.imag*(temp1.real-temp2.real));
    c[N/2-b].imag = T(0.5)*( -temp1.imag+temp2.imag - T(forward)*wb.real*(temp1.real-temp2.real) + wb.imag*(temp1.imag+temp2.imag));
  }
  temp1 = c[0];
  c[0].real = temp1.real + temp1.imag;
  c[0].imag = temp1.real - temp1.imag;

  if (forward == -1) {
    c[0].real *= T(0.5);
    c[0].imag *= T(0.5);
    fftc1<T>(f, N/2, 1, -1);
  }
}

// ---------- n-D real FFT (sizes are real sizes; Nyquist plane split) ----------
template <typename T>
inline void fftrn(T f[], T fnyquist[], int ndims, int size[], int forward) {
  int i, j, b;
  int index, indexneg = 0;
  int stepsize;
  const int N = size[ndims - 1];
  const T pi2n = (T(4)*std::asin(T(1))) / T(N);
  const T cospi2n = std::cos(pi2n);
  const T sinpi2n = std::sin(pi2n);
  cpx<T> wb, temp1, temp2;
  auto* c        = reinterpret_cast<cpx<T>*>(f);
  auto* cnyquist = reinterpret_cast<cpx<T>*>(fnyquist);

  int totalsize = 1;
  int* indices = static_cast<int*>(std::malloc(ndims * sizeof(int)));
  if (!indices) { std::printf("Error allocating memory in fftrn routine.\n"); std::exit(1); }

  size[ndims - 1] /= 2;
  for (i = 0; i < ndims; ++i) { totalsize *= size[i]; indices[i] = 0; }

  if (forward == 1) {
    fftcn<T>(f, ndims, size, 1);
    for (i = 0; i < totalsize / size[ndims - 1]; ++i)
      cnyquist[i] = c[i * size[ndims - 1]];
  }

  for (index = 0; index < totalsize; index += size[ndims - 1]) {
    wb.real = T(1); wb.imag = T(0);
    for (b = 1; b < N/4; ++b) {
      temp1 = wb;
      wb.real = cospi2n*temp1.real - sinpi2n*temp1.imag;
      wb.imag = cospi2n*temp1.imag + sinpi2n*temp1.real;
      temp1 = c[index + b];
      temp2 = c[indexneg + N/2 - b];
      c[index + b].real          = T(0.5)*( temp1.real+temp2.real + T(forward)*wb.real*(temp1.imag+temp2.imag) + wb.imag*(temp1.real-temp2.real));
      c[index + b].imag          = T(0.5)*( temp1.imag-temp2.imag - T(forward)*wb.real*(temp1.real-temp2.real) + wb.imag*(temp1.imag+temp2.imag));
      c[indexneg + N/2 - b].real = T(0.5)*( temp1.real+temp2.real - T(forward)*wb.real*(temp1.imag+temp2.imag) - wb.imag*(temp1.real-temp2.real));
      c[indexneg + N/2 - b].imag = T(0.5)*( -temp1.imag+temp2.imag - T(forward)*wb.real*(temp1.real-temp2.real) + wb.imag*(temp1.imag+temp2.imag));
    }
    temp1 = c[index];
    temp2 = cnyquist[indexneg / size[ndims - 1]];
    c[index].real                           = T(0.5)*( temp1.real+temp2.real + T(forward)*(temp1.imag+temp2.imag));
    c[index].imag                           = T(0.5)*( temp1.imag-temp2.imag - T(forward)*(temp1.real-temp2.real));
    cnyquist[indexneg / size[ndims - 1]].real = T(0.5)*( temp1.real+temp2.real - T(forward)*(temp1.imag+temp2.imag));
    cnyquist[indexneg / size[ndims - 1]].imag = T(0.5)*( -temp1.imag+temp2.imag - T(forward)*(temp1.real-temp2.real));

    // advance multi-index and maintain indexneg
    stepsize = size[ndims - 1];
    for (j = ndims - 2; j >= 0 && indices[j] == size[j] - 1; --j) {
      indices[j] = 0;
      indexneg  -= stepsize;
      stepsize  *= size[j];
    }
    if (j >= 0) {
      if (indices[j] == 0) indexneg += stepsize * (size[j] - 1);
      else                 indexneg -= stepsize;
      ++indices[j];
    }
  }

  if (forward == -1) fftcn<T>(f, ndims, size, -1);

  size[ndims - 1] *= 2;
  std::free(indices);
}

// ---------- Public wrappers (C-style signatures) ----------

// FLOAT versions (for GW arrays)
inline void fftc1f(float* f, int N, int skip, int forward)                    { fftc1<float>(f,N,skip,forward); }
inline void fftcnf(float* f, int ndims, int size[], int forward)              { fftcn<float>(f,ndims,size,forward); }
inline void fftr1f(float* f, int N, int forward)                              { fftr1<float>(f,N,forward); }
inline void fftrnf(float* f, float* fny, int ndims, int size[], int forward)  { fftrn<float>(f,fny,ndims,size,forward); }

// DOUBLE versions (for scalar field, momentum, etc.)
inline void fftc1d(double* f, int N, int skip, int forward)                    { fftc1<double>(f,N,skip,forward); }
inline void fftcnd(double* f, int ndims, int size[], int forward)              { fftcn<double>(f,ndims,size,forward); }
inline void fftr1d(double* f, int N, int forward)                              { fftr1<double>(f,N,forward); }
inline void fftrnd(double* f, double* fny, int ndims, int size[], int forward) { fftrn<double>(f,fny,ndims,size,forward); }
