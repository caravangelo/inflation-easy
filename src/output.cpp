// output.cpp - Diagnostics and data output
//
// This file implements simulation I/O: parameter dumps, spectra, energy diagnostics,
// and field snapshots written to the results/ directory.
//
// Output routines intentionally preserve the plain-text schema documented in the
// manuscript and consumed by notebooks/plot.ipynb. Several FFT-based routines
// transform global arrays in place and restore them before returning; callers must
// not interleave those routines with time evolution.

#include <filesystem>
#include <cfloat>   // DBL_MAX, FLT_MAX
#include <limits>
using namespace std::filesystem;

#include "main.h"
#include "ffteasy.hpp"

// Shared filename buffer used only while opening persistent output streams.
char name_[550];

// Open an output file or terminate immediately with a path-specific error.
static FILE* open_output_or_die(const char* path, const char* mode) {
    FILE* fp = std::fopen(path, mode);
    if (!fp) {
        std::fprintf(stderr, "Failed to open output file: %s\n", path);
        std::exit(1);
    }
    return fp;
}

// Preserve the established slice for production grids while remaining valid
// for the very small lattices used in tests and exploratory runs.
static constexpr int snapshot_slice_index() {
    return N < 64 ? (N > 5 ? 5 : N - 1) : 50;
}

// Zero one packed complex mode and its derivative when applying a hard cutoff.
void kill_mode(double *field, double *deriv)
{
    field[0] = 0.;
    field[1] = 0.;
    deriv[0] = 0.;
    deriv[1] = 0.;
    return;
}

// Average the selected Laplacian's effective momentum over the same DFT shells
// and multiplicities used by all isotropic spectra.
static const std::vector<double>& effective_k_bin_centers()
{
    static std::vector<double> centers;
    if (!centers.empty()) return centers;

    const int numbins = static_cast<int>(std::sqrt(3.0) * (N / 2)) + 1;
    std::vector<int> counts(numbins, 0);
    std::vector<double> sums(numbins, 0.0);

    for (int i = 0; i < N; ++i) {
        const int px = spatial::signed_mode(i);
        for (int j = 0; j < N; ++j) {
            const int py = spatial::signed_mode(j);
            for (int k = 1; k < N / 2; ++k) {
                const int bin = spatial::shell_index(px, py, k);
                if (bin >= numbins) continue;
                const double keff = std::sqrt(spatial::effective_momentum_squared(px, py, k, dx));
                counts[bin] += 2;
                sums[bin] += 2.0 * keff;
            }
            for (int k = 0; k <= N / 2; k += N / 2) {
                const int bin = spatial::shell_index(px, py, k);
                if (bin >= numbins) continue;
                const double keff = std::sqrt(spatial::effective_momentum_squared(px, py, k, dx));
                counts[bin] += 1;
                sums[bin] += keff;
            }
        }
    }

    centers.resize(numbins, 0.0);
    for (int bin = 0; bin < numbins; ++bin) {
        if (counts[bin] > 0) centers[bin] = rescale_B * sums[bin] / counts[bin];
    }
    return centers;
}

// -----------------------------------------------------------------------------
// Power spectra helpers
// -----------------------------------------------------------------------------
//
// The routines below implement a single, shared pipeline for computing and
// writing a 3D isotropically-binned power spectrum from a real field stored in
// r2c FFT layout (with a separate Nyquist plane buffer). Individual output
// functions (e.g. for the inflaton, δN, or log-mapped δN) prepare the field in
// real space and then call the common backend.
//
// Conventions:
//   - Momentum bin index is the Euclidean norm of integer lattice wave numbers,
//     rounded to the nearest integer. Bins outside the allocated range are
//     ignored.
//   - The output for each bin is <|X_k|^2> averaged over the modes in that bin,
//     multiplied by the caller-provided normalization factor.
//
// -----------------------------------------------------------------------------
// Scalar power spectrum backend
// -----------------------------------------------------------------------------
//
// This routine expects `field_fft` to be in-place real-to-complex FFT layout
// (FFTW-style r2c, stored in the real array), with the Nyquist plane stored in
// `nyquist_plane` (same convention as the rest of the code base).
//
// The output is an isotropically binned power spectrum P(k) ∝ <|X_k|^2>, where
// bins are labeled by the rounded Euclidean norm of the integer lattice wave
// numbers. The caller supplies the overall normalization.

static void write_isotropic_spectrum_from_fft_r2c(
    FILE *out,
    const std::vector<double> &field_fft,
    const double (*nyquist_plane)[2 * N],
    const double norm1)
{
    const std::vector<double>& momenta = effective_k_bin_centers();
    const int numbins = static_cast<int>(momenta.size());
    std::vector<int> numpoints(numbins, 0);
    std::vector<double> f2(numbins, 0.0);

    // Match legacy binning + multiplicities + indexing exactly.
    for (int i = 0; i < N; i++)
    {
        const int px = spatial::signed_mode(i);

        for (int j = 0; j < N; j++)
        {
            const int py = spatial::signed_mode(j);

            // Interior modes: 1 <= k < N/2 carry a conjugate partner -> weight 2
            for (int k = 1; k < N / 2; k++)
            {
                const int pz = k;

                const int bin = spatial::shell_index(px, py, pz);
                if (bin < 0 || bin >= numbins) continue;

                const double re = field_fft[idx(i, j, 2 * k)];
                const double im = field_fft[idx(i, j, 2 * k + 1)];
                const double fp2 = re * re + im * im;

                numpoints[bin] += 2;
                f2[bin] += 2.0 * fp2;
            }

            // Special cases: k = 0 and k = N/2 (Nyquist plane stored separately for k=N/2)
            for (int k = 0; k <= N / 2; k += N / 2)
            {
                const int pz = k;

                const int bin = spatial::shell_index(px, py, pz);
                if (bin < 0 || bin >= numbins) continue;

                double fp2 = 0.0;
                if (k == 0)
                {
                    const double re = field_fft[idx(i, j, 0)];
                    const double im = field_fft[idx(i, j, 1)];
                    fp2 = re * re + im * im;
                }
                else
                {
                    const double re = nyquist_plane[i][2 * j];
                    const double im = nyquist_plane[i][2 * j + 1];
                    fp2 = re * re + im * im;
                }

                numpoints[bin] += 1;
                f2[bin] += fp2;
            }
        }
    }

    for (int i = 0; i < numbins; i++)
    {
        if (numpoints[i] > 0) f2[i] /= (double)numpoints[i];
        std::fprintf(out, "%e %d %e\n", momenta[i], numpoints[i], norm1 * f2[i]);
    }

    std::fprintf(out, "\n");
    std::fflush(out);
}


#if perform_deltaN
// Construct δN (linear mapping) in-place in deltaN.
static void build_deltaN_linear()
{
    DECLARE_INDICES

    double Nmean = 0.0;
    LOOP Nmean += deltaN[idx(i, j, k)];
    Nmean /= (double)gridsize;
    LOOP deltaN[idx(i, j, k)] -= Nmean;
}

// Construct δN using the log mapping used for the LOG spectrum output.
static void build_deltaN_log()
{
    DECLARE_INDICES

    double factt = 0.0;
    double fmean = 0.0;

    LOOP
    {
        factt += fd[idx(i, j, k)] * std::pow(a, rescale_s - 1.0) / (ad * std::pow(a, rescale_s - 2.0));
        fmean += f[idx(i, j, k)];
    }

    factt /= (double)gridsize;
    fmean /= (double)gridsize;

    if (!std::isfinite(factt) || std::abs(factt) < 1e-30) {
        LOOP deltaN[idx(i, j, k)] = 0.0;
        return;
    }

    LOOP
    {
        deltaN[idx(i, j, k)] = 0.0;

        const double arg = 1.0 - eta_log * (f[idx(i, j, k)] - fmean) / factt;
        if (arg > 0.0)
        {
            deltaN[idx(i, j, k)] = (1.0 / eta_log) * std::log(arg);
        }
    }
}
#endif



// Write the scalar mean, variance, and mean physical-time velocity.
void meansvars(int flush)
{
    static FILE *means_, *vars_, *velocity_;
    DECLARE_INDICES

    double av, var, vel;

    static int first = 1;
    if (first) // Open output files
    {
        snprintf(name_, sizeof(name_), "results/means%s", ext_);
        means_ = open_output_or_die(name_, mode_);
        snprintf(name_, sizeof(name_), "results/variance%s", ext_);
        vars_ = open_output_or_die(name_, mode_);
        snprintf(name_, sizeof(name_), "results/velocity%s", ext_);
        velocity_ = open_output_or_die(name_, mode_);
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
        av  += f[idx(i,j,k)];
        vel += fd[idx(i,j,k)];
    }
    const double inv_gridsize = 1.0 / static_cast<double>(gridsize);
    av  *= inv_gridsize;
    vel *= inv_gridsize;

    // Evaluate <(phi - <phi>)^2> directly to avoid cancellation when the
    // homogeneous field is much larger than its fluctuations.
    LOOP var += pw2(f[idx(i,j,k)] - av);
    var *= inv_gridsize;

    vel = vel * std::pow(a, rescale_s - 1.0) * rescale_B;

    fprintf(means_,    " %e", av);
    fprintf(velocity_, " %e", vel);
    fprintf(vars_,     " %e", var);

    // Check for instability. See if the field has grown exponentially and become non-numerical at any point.
    if (av + DBL_MAX == av || (av != 0. && av / av != 1.))
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

// Write the scale factor, Hubble rate, and scale-factor acceleration.
void scale(int flush)
{
    static FILE *sf_;

    static int first = 1;
    if (first) // Open output file
    {
        snprintf(name_, sizeof(name_), "results/sf%s", ext_);
        sf_ = open_output_or_die(name_, mode_);
        first = 0;
    }

    const double acceleration =
        std::pow(a, 3.0 - 2.0 * rescale_s)
        * (2.0 * gradient_energy() / 3.0 + potential_energy())
        - (rescale_s + 1.0) * pw2(ad) / a;

    // Output a, H, and adotdot in physical units using rescalings
    fprintf(sf_, "%f %f %e %e\n",
    t, a,
    ad * rescale_B * std::pow(a, rescale_s - 2.0),
    pw2(rescale_B) * std::pow(a, 2.0 * rescale_s - 2.0) *
    (acceleration + (rescale_s - 1.0) * pw2(ad) / a));

    if (flush)
    fflush(sf_);
}

// Write the scalar power spectrum. If forcing_cutoff is enabled, this routine
// also filters the live field and derivative in Fourier space before restoring them.
void spectraf()
{
    static FILE *spectra_, *spectratimes_;
    static int first = 1;

    if (first)
    {
        std::snprintf(name_, sizeof(name_), "results/spectra%s", ext_);
        spectra_ = open_output_or_die(name_, mode_);

        std::snprintf(name_, sizeof(name_), "results/spectratimes%s", ext_);
        spectratimes_ = open_output_or_die(name_, mode_);

        first = 0;
    }

    const double norm1 = std::pow(L / rescale_B, 3) / std::pow((double)N, 6);

    // Forward FFT to k-space for spectrum evaluation.


    int arraysize_spec[] = {N, N, N};


    fftrnd(f.data(), (double *)fnyquist_p, 3, arraysize_spec, 1);



    write_isotropic_spectrum_from_fft_r2c(spectra_, f, fnyquist_p, norm1);
    std::fprintf(spectratimes_, "%f %e\n", t, a);

    std::fflush(spectratimes_);
    // Preserve the original behaviour: optionally apply a hard k-space cutoff to (f, fd).
    // The optional forcing cutoff operates in Fourier space and, if enabled, permanently filters
    // both the field and its time derivative by inverse-transforming the filtered modes back to real space.
    // The implementation below is intentionally identical to the legacy mode-killing logic.
    int i, j, k, px, py, pz, iconj, jconj;
    double pdisc;
    int arraysize[] = {N, N, N};

    if (high_cutoff_index > 0 && forcing_cutoff)
    {
        // FFT of fd to k-space (double version)
        fftrnd(fd.data(), (double *)fdnyquist_p, 3, arraysize, 1);

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
                    pdisc = std::sqrt(pw2((double)px) + pw2((double)py) + pw2((double)pz));
                    if (pdisc > high_cutoff_index || pdisc < low_cutoff_index)
                    {
                        kill_mode(&f[idx(i,j,2 * k)], &fd[idx(i,j,2 * k)]);
                    }
                }

                if (j > N / 2 || (i > N / 2 && (j == 0 || j == N / 2)))
                {
                    jconj = (j == 0 ? 0 : N - j);

                    pdisc = std::sqrt(pw2((double)px) + pw2((double)py));
                    if (pdisc > high_cutoff_index || pdisc < low_cutoff_index)
                    {
                        kill_mode(&f[idx(i,j,0)], &fd[idx(i,j,0)]);
                        f[idx(iconj,jconj,0)] = f[idx(i,j,0)];
                        f[idx(iconj,jconj,1)] = -f[idx(i,j,1)];
                        fd[idx(iconj,jconj,0)] = fd[idx(i,j,0)];
                        fd[idx(iconj,jconj,1)] = -fd[idx(i,j,1)];
                    }

                    pdisc = std::sqrt(pw2((double)px) + pw2((double)py) + pw2((double)(N / 2.)));
                    if (pdisc > high_cutoff_index || pdisc < low_cutoff_index)
                    {
                        kill_mode(&fnyquist_p[i][2 * j], &fdnyquist_p[i][2 * j]);
                        fnyquist_p[iconj][2 * jconj]     = fnyquist_p[i][2 * j];
                        fnyquist_p[iconj][2 * jconj + 1] = -fnyquist_p[i][2 * j + 1];
                        fdnyquist_p[iconj][2 * jconj]     = fdnyquist_p[i][2 * j];
                        fdnyquist_p[iconj][2 * jconj + 1] = -fdnyquist_p[i][2 * j + 1];
                    }
                }
                else if ((i == 0 || i == N / 2) && (j == 0 || j == N / 2))
                {
                    pdisc = std::sqrt(pw2((double)px) + pw2((double)py));
                    if (pdisc > high_cutoff_index || pdisc < low_cutoff_index)
                    {
                        kill_mode(&f[idx(i,j,0)], &fd[idx(i,j,0)]);
                    }

                    pdisc = std::sqrt(pw2((double)px) + pw2((double)py) + pw2((double)(N / 2.)));
                    if (pdisc > high_cutoff_index || pdisc < low_cutoff_index)
                    {
                        kill_mode(&fnyquist_p[i][2 * j], &fdnyquist_p[i][2 * j]);
                    }
                }
            }
        }
        // Backward FFT of fd
        fftrnd(fd.data(), (double *)fdnyquist_p, 3, arraysize, -1);
    }

    // Backward FFT to restore the real-space field.
    fftrnd(f.data(), (double *)fnyquist_p, 3, arraysize, -1);

}



#if calculate_SIGW

//=============================================================================
// Gravitational-wave spectra (transverse-traceless projection)
//=============================================================================

namespace {

    // Open an output file the first time this routine is called.
    static inline void open_output_once(FILE *&fp, int &first, const char *relative_path)
    {
        if (!first) return;
        std::snprintf(name_, sizeof(name_), "%s%s", relative_path, ext_);
        fp = open_output_or_die(name_, mode_);
        first = 0;
    }

    // Conversion of tensor time-derivative modes to physical units before projection.
    enum class GwScaleMode {
        None,
        DoubleThenFloat,
        FloatMultiply
    };

    static inline void load_tensor_component(
        const std::vector<float> &buf,
        const float (*nyq)[2 * N],
        const bool on_nyquist_plane,
        const int i, const int j2,
        const size_t idx_mode,
        const GwScaleMode scale_mode,
        const double to_phys,
        float &re, float &im)
    {
        if (!on_nyquist_plane) {
            re = buf[idx_mode];
            im = buf[idx_mode + 1];
        } else {
            re = nyq[i][j2];
            im = nyq[i][j2 + 1];
        }

        if (scale_mode == GwScaleMode::None) return;

        if (scale_mode == GwScaleMode::DoubleThenFloat) {
            re = (float)((double)re * to_phys);
            im = (float)((double)im * to_phys);
        } else { // FloatMultiply
            const float tp = (float)to_phys;
            re *= tp;
            im *= tp;
        }
    }

    // Isotropically binned TT-projected spectrum over the shared output-shell range.
    // Output columns: comoving effective momentum (in reduced Planck units),
    // number of modes, and spectrum.
    static void write_gw_spectrum_impl(
        FILE *out,
        std::vector<float> (&h)[6],
        float (*hnyq)[N][2 * N],
        const GwScaleMode scale_mode,
        const double to_phys)
    {
        const int numbins = (int)(std::sqrt(3.0) * (N/2)) + 1;

        std::vector<int>   numpoints(numbins, 0);
        std::vector<float> p(numbins, 0.f), f2(numbins, 0.f);
        const std::vector<double>& bin_centers = effective_k_bin_centers();
        for (int i = 0; i < numbins; ++i) p[i] = static_cast<float>(bin_centers[i]);

        int arraysize[] = {N, N, N};
        for (int c = 0; c < 6; ++c) fftrnf(h[c].data(), (float*)hnyq[c], 3, arraysize, 1);

        // A self-conjugate Nyquist component has no unambiguous sign for the
        // real projector momentum. Modes containing one are retained only as
        // convention-dependent UV diagnostics; see DEVELOPER_GUIDE.md.
        for (int i = 0; i < N; ++i) {
            int px = spatial::signed_mode(i);
            for (int j = 0; j < N; ++j) {
                int py = spatial::signed_mode(j);

                // Interior r2c modes: 1 <= k < N/2 (count conjugate partner with weight 2).
                for (int k = 1; k < N/2; ++k) {
                    int pz = k;

                    float kx = static_cast<float>(spatial::effective_momentum_component(px, dx));
                    float ky = static_cast<float>(spatial::effective_momentum_component(py, dx));
                    float kz = static_cast<float>(spatial::effective_momentum_component(pz, dx));
                    float kt2 = kx*kx + ky*ky + kz*kz;
                    if (kt2 == 0.f) continue;

                    int bin = spatial::shell_index(px, py, pz);
                    if (bin >= numbins) continue;

                    size_t idx_mode = idx(i, j, 2*k);

                    float C_re[3][3] = {{0}}, C_im[3][3] = {{0}};
                    for (int l = 0; l < 3; ++l) for (int m = l; m < 3; ++m) {
                        int comp = sym_idx(l, m);
                        float re, im;
                        load_tensor_component(h[comp], hnyq[comp], false, 0, 0, idx_mode, scale_mode, to_phys, re, im);
                        C_re[l][m] = re; C_im[l][m] = im;
                        if (m != l) { C_re[m][l] = re; C_im[m][l] = im; }
                    }

                    float inv = 1.0f / std::sqrt(kt2);
                    float kh[3] = {kx*inv, ky*inv, kz*inv};

                    float P[3][3];
                    for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b)
                        P[a][b] = (a == b ? 1.f : 0.f) - kh[a]*kh[b];

                    float A_re[3][3] = {{0}}, A_im[3][3] = {{0}};
                    float T_re = 0.f, T_im = 0.f;

                    for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) {
                        float r = 0.f, im = 0.f;
                        for (int l = 0; l < 3; ++l) for (int m = 0; m < 3; ++m) {
                            float pal = P[a][l], pbm = P[b][m];
                            r  += pal * C_re[l][m] * pbm;
                            im += pal * C_im[l][m] * pbm;
                        }
                        A_re[a][b] = r; A_im[a][b] = im;
                    }

                    for (int l = 0; l < 3; ++l) for (int m = 0; m < 3; ++m) {
                        T_re += P[l][m] * C_re[l][m];
                        T_im += P[l][m] * C_im[l][m];
                    }

                    float fp2 = 0.f;
                    for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) {
                        float Xr = A_re[a][b] - 0.5f * P[a][b] * T_re;
                        float Xi = A_im[a][b] - 0.5f * P[a][b] * T_im;
                        fp2 += Xr*Xr + Xi*Xi;
                    }

                    numpoints[bin] += 2;
                    f2[bin] += 2.f * fp2;
                }

                // Special modes: k = 0 and k = N/2 (Nyquist plane). Use pz = +N/2 for k = N/2.
                for (int kk = 0; kk <= 1; ++kk) {
                    int k  = kk * (N/2);
                    int pz = k;

                    float kx = static_cast<float>(spatial::effective_momentum_component(px, dx));
                    float ky = static_cast<float>(spatial::effective_momentum_component(py, dx));
                    float kz = static_cast<float>(spatial::effective_momentum_component(pz, dx));
                    float kt2 = kx*kx + ky*ky + kz*kz;
                    if (kt2 == 0.f) continue;

                    int bin = spatial::shell_index(px, py, pz);
                    if (bin >= numbins) continue;

                    const int j2 = 2*j;

                    float C_re[3][3] = {{0}}, C_im[3][3] = {{0}};
                    for (int l = 0; l < 3; ++l) for (int m = l; m < 3; ++m) {
                        int comp = sym_idx(l, m);
                        float re, im;
                        if (k == 0) {
                            load_tensor_component(h[comp], hnyq[comp], false, 0, 0, idx(i, j, 0), scale_mode, to_phys, re, im);
                        } else {
                            load_tensor_component(h[comp], hnyq[comp], true, i, j2, 0, scale_mode, to_phys, re, im);
                        }
                        C_re[l][m] = re; C_im[l][m] = im;
                        if (m != l) { C_re[m][l] = re; C_im[m][l] = im; }
                    }

                    float inv = 1.0f / std::sqrt(kt2);
                    float kh[3] = {kx*inv, ky*inv, kz*inv};

                    float P[3][3];
                    for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b)
                        P[a][b] = (a == b ? 1.f : 0.f) - kh[a]*kh[b];

                    float A_re[3][3] = {{0}}, A_im[3][3] = {{0}};
                    float T_re = 0.f, T_im = 0.f;

                    for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) {
                        float r = 0.f, im = 0.f;
                        for (int l = 0; l < 3; ++l) for (int m = 0; m < 3; ++m) {
                            float pal = P[a][l], pbm = P[b][m];
                            r  += pal * C_re[l][m] * pbm;
                            im += pal * C_im[l][m] * pbm;
                        }
                        A_re[a][b] = r; A_im[a][b] = im;
                    }

                    for (int l = 0; l < 3; ++l) for (int m = 0; m < 3; ++m) {
                        T_re += P[l][m] * C_re[l][m];
                        T_im += P[l][m] * C_im[l][m];
                    }

                    float fp2 = 0.f;
                    for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) {
                        float Xr = A_re[a][b] - 0.5f * P[a][b] * T_re;
                        float Xi = A_im[a][b] - 0.5f * P[a][b] * T_im;
                        fp2 += Xr*Xr + Xi*Xi;
                    }

                    numpoints[bin] += 1;
                    f2[bin] += fp2;
                }
            }
        }

        for (int c = 0; c < 6; ++c) fftrnf(h[c].data(), (float*)hnyq[c], 3, arraysize, -1);

        const float norm1 = 0.5f * std::pow(L / rescale_B, 3) / std::pow((float)N, 6);
        for (int i = 0; i < numbins; ++i) {
            if (numpoints[i] > 0) f2[i] /= (float)numpoints[i];
            std::fprintf(out, "%e %d %e\n", p[i], numpoints[i], norm1 * f2[i]);
        }
        std::fprintf(out, "\n");
        std::fflush(out);
    }

} // anonymous namespace


// Write the TT-projected tensor power spectrum during inflation.
void spectraGW()
{
    static FILE *fp = nullptr;
    static int first = 1;
    open_output_once(fp, first, "results/spectraGW");
    write_gw_spectrum_impl(fp, hij, hijnyquist_p, GwScaleMode::None, 1.0);
}


// Write the TT-projected tensor-velocity spectrum in physical units during inflation.
void spectraGWdot()
{
    static FILE *fp = nullptr;
    static int first = 1;
    open_output_once(fp, first, "results/spectraGWdot");
    const double to_phys = rescale_B * std::pow(a, rescale_s - 1.0);
    write_gw_spectrum_impl(fp, hijd, hijdnyquist_p, GwScaleMode::DoubleThenFloat, to_phys);
}


// Write the TT-projected tensor power spectrum during post-inflationary evolution.
void spectraGW_post_inflation()
{
    static FILE *fp = nullptr;
    static int first = 1;
    open_output_once(fp, first, "results/post_inflation/spectraGW");
    write_gw_spectrum_impl(fp, hij, hijnyquist_p, GwScaleMode::None, 1.0);
}


// Write the post-inflationary tensor-velocity spectrum in physical units.
void spectraGWdot_post_inflation()
{
    static FILE *fp = nullptr;
    static int first = 1;
    open_output_once(fp, first, "results/post_inflation/spectraGWdot");
    const float to_phys = rescale_B * pow(a, rescale_s - 1.f);
    write_gw_spectrum_impl(fp, hijd, hijdnyquist_p, GwScaleMode::FloatMultiply, (double)to_phys);
}

#endif

// Write the comoving effective momentum (in reduced Planck units) associated
// with each isotropic DFT bin.
// The mapping uses the eigenvalue of the finite-difference Laplacian.
void get_modes()
{
    static FILE *modes_;
    const std::vector<double>& bin_centers = effective_k_bin_centers();

    static int first=1;
    if(first) // Open output files
    {
        snprintf(name_, sizeof(name_), "results/modes%s",ext_);
        modes_=open_output_or_die(name_, mode_);
        first=0;
    }

    for (double momentum : bin_centers) fprintf(modes_, "%e\n", momentum);

    fprintf(modes_,"\n");
    fflush(modes_);

    return;
}

// Write the equilateral scalar bispectrum using the code's established binning convention.
void bispectraf()
{
    static FILE *bispectra_; // Final-time equilateral scalar-bispectrum output
    const int numbins=(int)(std::sqrt(3.0)*(N/2))+1; // Actual number of bins for the number of dimensions
    std::vector<int> numpoints(numbins, 0); // Number of points in each momentum bin
    std::vector<double> bisreal(numbins, 0.0), bisimag(numbins, 0.0);
    const std::vector<double>& momenta = effective_k_bin_centers();
    double mean=0.;
    int i1,j1,k1;
    int i2,j2,k2,i2n,j2n;
    int i3,j3,k3;
    double px1,py1,px2,py2,px3,py3;
    int i,j,k;
    double pzaus1,pzaus2,counts;
    double f1r,f1i,f2r,f2i,f3r,f3i;
    const double norm1=std::pow(L/rescale_B,6)/std::pow((double)N,9);
    int arraysize[]={N,N,N}; // Array of grid size in all dimensions - used by FFT routine
    double kf;

    LOOP
    mean += f[idx(i,j,k)];

    mean = mean/(double)gridsize;

    LOOP
    f[idx(i,j,k)] -= mean;

    static int first=1;
    if(first) // Open output files
    {
        snprintf(name_, sizeof(name_), "results/bispectra%s",ext_);
        bispectra_=open_output_or_die(name_, mode_);
        first=0;
    }

    fftrnd(f.data(), (double *)fnyquist_p, 3, arraysize, 1); // Transform field values to Fourier space

    const auto load_mode = [&](int ii, int jj, int kk, bool conjugate,
                               double& re, double& im) {
        if (kk == N / 2) {
            re = fnyquist_p[ii][2 * jj];
            im = fnyquist_p[ii][2 * jj + 1];
        } else {
            re = f[idx(ii, jj, 2 * kk)];
            im = f[idx(ii, jj, 2 * kk + 1)];
        }
        if (conjugate) im = -im;
    };

    for(k=0;k<numbins;k++)
    {
        kf = (double)k;
        for(i1=0;i1<N;i1++) for(j1=0;j1<N;j1++)
        {//for1
            px1=(i1<=N/2 ? i1 : i1-N);
            py1=(j1<=N/2 ? j1 : j1-N);
            pzaus1 = pw2(kf)-pw2(px1)-pw2(py1);
            if(pzaus1>=0)
            for(k1=(int)std::round(std::sqrt(pzaus1))-1; k1 < (int)std::round(std::sqrt(pzaus1))+2; k1++)
            if(std::abs(std::sqrt(pw2(px1)+pw2(py1)+pw2((double)k1))-kf) < 1.5 && k1 <= N/2 && k1 >= 0)
            for(i2=0;i2<N;i2++) for(j2=0;j2<N;j2++)
            {//for2
                px2=(i2<=N/2 ? i2 : i2-N);
                py2=(j2<=N/2 ? j2 : j2-N);
                pzaus2 = pw2(kf)-pw2(px2)-pw2(py2);
                if(pzaus2>=0)
                for(k2=(int)std::round(std::sqrt(pzaus2))-1; k2 < (int)std::round(std::sqrt(pzaus2))+2; k2++)
                if(std::abs(std::sqrt(pw2(px2)+pw2(py2)+pw2((double)k2))-kf) < 1.5 && k2 <= N/2 && k2 >= 0)
                {//if triangleapprox2
                    px3 = px1 + px2;
                    py3 = py1 + py2;
                    k3  = k1 + k2;

                    if(px3 <= N/2 && px3 > -N/2  && py3 <= N/2 && py3 > -N/2)
                    {//if in lattice
                        i3=(px3>=0 ? px3 : px3+N);
                        j3=(py3>=0 ? py3 : py3+N);
                        if(std::abs(std::sqrt(pw2(px3)+pw2(py3)+pw2((double)k3))-kf)< 1.5 && k3 <= N/2)
                        {
                            load_mode(i1, j1, k1, false, f1r, f1i);
                            load_mode(i2, j2, k2, false, f2r, f2i);
                            load_mode(i3, j3, k3, false, f3r, f3i);
                            counts = 1.;
                            if(k1 != (int)N/2 && k2 != (int)N/2 && (k1 != 0 || k2 != 0))
                            counts = 2.;

                            numpoints[k] += (int)counts;
                            bisreal[k] += counts*(f1r*f2r*f3r - f1i*f2i*f3r + f1i*f2r*f3i + f2i*f1r*f3i);
                            bisimag[k] += counts*(-f1r*f2r*f3i + f1i*f2i*f3i + f1i*f2r*f3r + f2i*f1r*f3r);
                        }
                        if(k1 != 0 && k2!= 0)
                        {
                            k3 = k1-k2;
                            if(k3>=0)
                            if(std::abs(std::sqrt(pw2(px3)+pw2(py3)+pw2((double)k3))-kf) < 1.5 && k3 <= N/2)
                            {
                                i2n=(-px2 >=0 ? -px2 : -px2+N);
                                j2n=(-py2 >=0 ? -py2 : -py2+N);

                                load_mode(i1, j1, k1, false, f1r, f1i);
                                load_mode(i2n, j2n, k2, true, f2r, f2i);
                                load_mode(i3, j3, k3, false, f3r, f3i);

                                numpoints[k] += 2;
                                bisreal[k] += 2.0*(f1r*f2r*f3r - f1i*f2i*f3r + f1i*f2r*f3i + f2i*f1r*f3i);
                                bisimag[k] += 2.0*(-f1r*f2r*f3i + f1i*f2i*f3i + f1i*f2r*f3r + f2i*f1r*f3r);
                            }
                        }
                    } // if in lattice
                } // if triangleapprox2
            } // for2
        } // for1
    } // for k

    for(k=0;k<numbins;k++)
    {
        if(numpoints[k]>0) {
            bisreal[k] = bisreal[k]/numpoints[k];
            bisimag[k] = bisimag[k]/numpoints[k];
        }

        fprintf(bispectra_,"%e %d %e %e\n",
        momenta[k],numpoints[k],norm1*bisreal[k],norm1*bisimag[k]);
    }

    fftrnd(f.data(), (double *)fnyquist_p, 3, arraysize, -1);

    LOOP
    f[idx(i,j,k)] += mean;

    fprintf(bispectra_,"\n");
    fflush(bispectra_);

    return;
}

// Write a complete three-dimensional scalar-field snapshot.
void box()
{
    static FILE *box_;
    int i,j,k;
    static int first=1;
    if(first) // Open output files
    {
        snprintf(name_, sizeof(name_), "results/box%s",ext_);
        box_ = open_output_or_die(name_, mode_);
        first=0;
    }

    for(i=0;i<N;i++) for(j=0;j<N;j++) for(k=0;k<N;k++)
    {
        fprintf(box_,"%.17g\n",f[idx(i,j,k)]);
    }
    fprintf(box_,"\n");
    fflush(box_);
}

// Write a fixed two-dimensional slice of the scalar field.
void box2d()
{
    static FILE *snapshots_2d_phi_;
    const int i = snapshot_slice_index();
    int j, k;

    static int first=1;
    if(first) // Open output files
    {
        snprintf(name_, sizeof(name_), "results/snapshots_2d_phi%s",ext_);
        snapshots_2d_phi_ = open_output_or_die(name_, mode_);
        first=0;
    }

    for(j=0;j<N;j++) for(k=0;k<N;k++)
    {
        fprintf(snapshots_2d_phi_,"%.17g\n",f[idx(i,j,k)]);
    }
    fprintf(snapshots_2d_phi_,"\n");
    fflush(snapshots_2d_phi_);
}

// Write the matching two-dimensional slice of the scalar velocity.
void box2dot()
{
    static FILE *snapshots_2d_phidot_;
    const int i = snapshot_slice_index();
    int j, k;

    static int first=1;
    if(first) // Open output files
    {
        snprintf(name_, sizeof(name_), "results/snapshots_2d_phidot%s",ext_);
        snapshots_2d_phidot_ = open_output_or_die(name_, mode_);
        first=0;
    }

    for(j=0;j<N;j++) for(k=0;k<N;k++)
    {
        fprintf(snapshots_2d_phidot_, "%.17g\n",
                fd[idx(i,j,k)] * rescale_B * std::pow(a, rescale_s - 1.0));
    }
    fprintf(snapshots_2d_phidot_,"\n");
    fflush(snapshots_2d_phidot_);
}

// Write kinetic, gradient, and potential energies plus the Friedmann constraint ratio.
void energy()
{
    static FILE *energy_, *conservation_;
    double deriv_energy, grad_energy, pot_energy;

    double totalE = 0.;
    const double physical_energy_scale = pw2(rescale_B);
    static int first = 1;
    if (first) // Open output files
    {
        snprintf(name_, sizeof(name_), "results/energy%s", ext_);
        energy_ = open_output_or_die(name_, mode_);

        snprintf(name_, sizeof(name_), "results/conservation%s", ext_);
        conservation_ = open_output_or_die(name_, mode_);

        first = 0;
    }

    fprintf(energy_, "%f", t); // Output time
    fprintf(energy_, " %e", a);

    // Calculate and output kinetic (time derivative) energy
    deriv_energy = kin_energy();
    totalE += deriv_energy;
    fprintf(energy_, " %e", deriv_energy * physical_energy_scale);

    // Calculate and output gradient energy
    grad_energy = gradient_energy();
    totalE += grad_energy;
    fprintf(energy_, " %e", grad_energy * physical_energy_scale);

    // Calculate and output potential energy
    pot_energy = potential_energy();
    totalE += pot_energy;
    fprintf(energy_, " %e", pot_energy * physical_energy_scale);

    fprintf(energy_, "\n");
    fflush(energy_);

    // Energy conservation
    fprintf(conservation_, "%e %e %e\n",
    t, a, 3.0 * std::pow(a, 2.0 * rescale_s - 4.0) * pw2(ad) / (totalE));
    fflush(conservation_);
}

// Format an elapsed duration in days, hours, minutes, and seconds.
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

// Write a normalized one-point histogram of the scalar field and its bin metadata.
void histograms()
{
    static FILE *histogram_, *histogramtimes_;
    int i = 0, j = 0, k = 0;
    int binnum; // Index of bin for a given field value
    static std::vector<double> binfreq; // Reused histogram buffer to avoid per-call allocation
    double bmin, bmax, df; // Minimum and maximum field values for each field and bin spacing
    int numpts; // Count the number of points in the histogram for each field.
    static int first = 1;
    if (first) // Open output files
    {
        snprintf(name_, sizeof(name_), "results/histogram%s", ext_);
        histogram_ = open_output_or_die(name_, mode_);

        snprintf(name_, sizeof(name_), "results/histogramtimes%s", ext_);
        histogramtimes_ = open_output_or_die(name_, mode_);
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
    df = (bmax - bmin) / (double)(nbins); // bmin will be at the bottom of the first bin and bmax at the top of the last
    if (!std::isfinite(df) || df <= 0.0) df = 1.0;

    if ((int)binfreq.size() != nbins) binfreq.assign(nbins, 0.0);
    else std::fill(binfreq.begin(), binfreq.end(), 0.0);

    // Initialize all frequencies to zero
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
    if (numpts == 0) numpts = 1;
    for (i = 0; i < nbins; i++)
    fprintf(histogram_, "%e\n", binfreq[i] / (double)numpts); // Output bin frequency
    fprintf(histogram_, "\n"); // Stick a blank line between times to make the file more readable
    fflush(histogram_);
    fprintf(histogramtimes_, " %e %e", bmin, df); // Output the starting point and stepsize

    fprintf(histogramtimes_, "\n");
    fflush(histogramtimes_);
}

#if perform_deltaN


// Write the power spectrum of the mean-subtracted deltaN field.
void spectraN()
{
    static FILE *spectraN_;
    static int first = 1;

    if (first)
    {
        std::snprintf(name_, sizeof(name_), "results/spectraN%s", ext_);
        spectraN_ = open_output_or_die(name_, mode_);
        first = 0;
    }

    // δN is assumed to have been computed in real space by the evolution step.
    // The mean is removed before transforming to Fourier space.
    build_deltaN_linear();

    const double norm1 = std::pow(L / rescale_B, 3) / std::pow((double)N, 6);
    // Forward FFT to k-space for spectrum evaluation.

    int arraysize_spec[] = {N, N, N};

    fftrnd(deltaN.data(), (double *)fnyquist_p, 3, arraysize_spec, 1);


    write_isotropic_spectrum_from_fft_r2c(spectraN_, deltaN, fnyquist_p, norm1);

    // Backward FFT to restore the real-space field.
    fftrnd(deltaN.data(), (double *)fnyquist_p, 3, arraysize_spec, -1);

}



// Write the complete three-dimensional deltaN field.
void boxN()
{
    static FILE *boxN_;
    int i, j, k;
    static int first = 1;
    if (first) // Open output file
    {
        snprintf(name_, sizeof(name_), "results/boxN%s", ext_);
        boxN_ = open_output_or_die(name_, mode_);
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

// Write a fixed two-dimensional slice of the deltaN field.
void box2dN()
{
    static FILE *snapshots_2d_deltaN_;
    const int i = snapshot_slice_index();
    int j, k;
    static int first = 1;
    if (first) // Open output file
    {
        snprintf(name_, sizeof(name_), "results/snapshots_2d_deltaN%s", ext_);
        snapshots_2d_deltaN_ = open_output_or_die(name_, mode_);
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

// Write a normalized one-point histogram of completed deltaN patches.
void histogramsN(const std::vector<unsigned char>& completed)
{
    static FILE *histogramN_, *histogramtimesN_;
    int binnum;
    static std::vector<double> binfreq;
    double bmin, bmax, df;
    std::size_t numpts;

    static int first = 1;
    if (first)
    {
        snprintf(name_, sizeof(name_), "results/histogramN%s", ext_);
        histogramN_ = open_output_or_die(name_, mode_);

        snprintf(name_, sizeof(name_), "results/histogramtimesN%s", ext_);
        histogramtimesN_ = open_output_or_die(name_, mode_);
        first = 0;
    }

    fprintf(histogramtimesN_, "%f %e", t, a);

    bmin = std::numeric_limits<double>::infinity();
    bmax = -std::numeric_limits<double>::infinity();
    for (std::size_t id = 0; id < deltaN.size(); ++id) {
        if (!completed[id]) continue;
        bmin = std::min(bmin, deltaN[id]);
        bmax = std::max(bmax, deltaN[id]);
    }

    const bool single_value_histogram = (bmax == bmin);
    df = single_value_histogram ? 1.0 : (bmax - bmin) / (double)(nbins);
    if (!std::isfinite(df) || df <= 0.0) df = 1.0;

    if ((int)binfreq.size() != nbins) binfreq.assign(nbins, 0.0);
    else std::fill(binfreq.begin(), binfreq.end(), 0.0);

    numpts = 0;
    for (std::size_t id = 0; id < deltaN.size(); ++id) {
        if (!completed[id]) continue;
        binnum = single_value_histogram
            ? nbins - 1
            : (int)((deltaN[id] - bmin) / df);
        if (deltaN[id] == bmax) binnum = nbins - 1;
        if (binnum >= 0 && binnum < nbins)
        {
            binfreq[binnum]++;
            numpts++;
        }
    }

    if (numpts == 0) {
        std::fprintf(stderr, "No finite completed deltaN patches are available for the histogram.\n");
        std::exit(1);
    }

    for (int i = 0; i < nbins; i++)
    fprintf(histogramN_, "%e\n", binfreq[i] / (double)numpts);
    fprintf(histogramN_, "\n");
    fflush(histogramN_);

    fprintf(histogramtimesN_, " %e %e\n", bmin, df);
    fflush(histogramtimesN_);
}


// Write the spectrum of the optional logarithmic curvature mapping.
void spectraLOG()
{
    static FILE *spectraLOG_;
    static int first = 1;

    if (first)
    {
        std::snprintf(name_, sizeof(name_), "results/spectraLOG%s", ext_);
        spectraLOG_ = open_output_or_die(name_, mode_);
        first = 0;
    }

    // Build the log-mapped δN field in-place, then compute its spectrum.
    build_deltaN_log();

    const double norm1 = std::pow(L / rescale_B, 3) / std::pow((double)N, 6);
    // Forward FFT to k-space for spectrum evaluation.

    int arraysize_spec[] = {N, N, N};

    fftrnd(deltaN.data(), (double *)fnyquist_p, 3, arraysize_spec, 1);


    write_isotropic_spectrum_from_fft_r2c(spectraLOG_, deltaN, fnyquist_p, norm1);

    // Backward FFT to restore the real-space field.
    fftrnd(deltaN.data(), (double *)fnyquist_p, 3, arraysize_spec, -1);

}



// Write the one-point histogram of the optional logarithmic curvature mapping.
void histogramsLOG()
{
    static FILE *histogramLOG_, *histogramtimesLOG_;
    int i=0, j=0, k=0;
    int binnum;
    static std::vector<double> binfreq;
    double bmin, bmax, df;
    int numpts;

    static int first = 1;
    if (first)
    {
        snprintf(name_, sizeof(name_), "results/histogramLOG%s", ext_);
        histogramLOG_ = open_output_or_die(name_, mode_);

        snprintf(name_, sizeof(name_), "results/histogramtimesLOG%s", ext_);
        histogramtimesLOG_ = open_output_or_die(name_, mode_);
        first = 0;
    }

    fprintf(histogramtimesLOG_, "%f %e", t, a);

    double factt = 0.;
    double fmean = 0.;

    LOOP
    {
        factt += fd[idx(i,j,k)] * pow(a, rescale_s - 1) / (ad * pow(a, rescale_s - 2.));
        fmean += f[idx(i,j,k)];
    }

    factt = factt / (double)gridsize;
    fmean = fmean / (double)gridsize;

    if (!std::isfinite(factt) || std::abs(factt) < 1e-30) {
        for (i = 0; i < nbins; i++) fprintf(histogramLOG_, "%e\n", 0.0);
        fprintf(histogramLOG_, "\n");
        fflush(histogramLOG_);
        fprintf(histogramtimesLOG_, " %e %e\n", 0.0, 1.0);
        fflush(histogramtimesLOG_);
        return;
    }

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

    df = (bmax - bmin) / (double)(nbins);
    if (!std::isfinite(df) || df <= 0.0) df = 1.0;

    if ((int)binfreq.size() != nbins) binfreq.assign(nbins, 0.0);
    else std::fill(binfreq.begin(), binfreq.end(), 0.0);

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

    if (numpts == 0) numpts = 1;
    for (i = 0; i < nbins; i++)
    fprintf(histogramLOG_, "%e\n", binfreq[i] / (double)numpts);
    fprintf(histogramLOG_, "\n");
    fflush(histogramLOG_);

    fprintf(histogramtimesLOG_, " %e %e\n", bmin, df);
    fflush(histogramtimesLOG_);
}

#endif

// On its first call, record the resolved configuration and start time. On its
// second call, append the end time and elapsed wall-clock duration.
void output_parameters()
{
    static FILE *info_;
    static time_t tStart, tFinish; // Keep track of elapsed clock time

    static int first = 1;
    if (first) // At beginning of run output run parameters
    {
        snprintf(name_, sizeof(name_), "results/info%s", ext_);
        info_ = open_output_or_die(name_, mode_);

        fprintf(info_, "--------------------------\n");
        fprintf(info_, "General Program Information\n");
        fprintf(info_, "-----------------------------\n");
        fprintf(info_, "Grid size=%d^%d\n", N, 3);
        fprintf(info_, "L=%f\n", L);
        fprintf(info_, "f0=%f\n", initial_field);
        fprintf(info_, "fd0=%f\n", initial_derivative);
        fprintf(info_, "dt=%f, dt/dx=%f\n", dt, dt / dx);
        fprintf(info_, "spatial_stencil_order=%d\n", spatial::order);
        fprintf(info_, "inflation_integrator=%s\n", integrator_name());
#if perform_deltaN
        fprintf(info_, "deltaN_integrator=%s\n", deltaN_integrator_name());
#endif
#if post_inflation
        fprintf(info_, "post_inflation_integrator=%s\n", post_inflation_integrator_name());
#endif
        if (integrator == INTEGRATOR_RK45
#if perform_deltaN
            || deltaN_integrator == INTEGRATOR_RK45
#endif
#if post_inflation
            || post_inflation_integrator == INTEGRATOR_RK45
#endif
        ) {
            fprintf(info_, "rk45_abs_tol=%e\n", rk45_abs_tol);
            fprintf(info_, "rk45_rel_tol=%e\n", rk45_rel_tol);
            fprintf(info_, "rk45_min_dt=%e\n", rk45_min_dt);
            fprintf(info_, "rk45_max_dt=%e\n", rk45_max_dt);
            fprintf(info_, "rk45_safety=%e\n", rk45_safety);
        }
        fprintf(info_, "rescale_s=%f\n", rescale_s);
        fprintf(info_, "rescale_B=%e\n", rescale_B);
        if (initial_mass_squared.has_value()) {
            fprintf(info_, "initial_mass_squared=%.17g\n", *initial_mass_squared);
        } else {
            fprintf(info_, "initial_mass_squared=not_set\n");
        }
        fprintf(info_, "linear_metric_perturbations=%d\n", linear_metric_perturbations);
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

// Dispatch inflationary outputs. For leapfrog, temporarily synchronize fields
// and derivatives before diagnostics, then restore the staggered state.
void save(int infrequent)
{
    if (inflation_uses_staggered_derivatives() && t > 0.) // Synchronize field values and derivatives
    apply_leapfrog_drift(-.5 * dt * pow(astep, rescale_s - 1));

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
        {
            spectraf();
#if calculate_SIGW
            spectraGW();
            spectraGWdot();
#endif
        }
        if (output_histogram)
        histograms();
    }

    if (inflation_uses_staggered_derivatives() && t > 0.) // Desynchronize field values and derivatives
    apply_leapfrog_drift(.5 * dt * pow(astep, rescale_s - 1));
}

// Write products that are defined only for the final inflationary state.
void save_last()
{
    get_modes();
    if (output_bispectrum)
    bispectraf();

#if perform_deltaN
    if (output_LOG)
    {
        deltaN.resize(static_cast<std::size_t>(gridsize));
        // The logarithmic mapping outputs use deltaN as a temporary lattice field.
        spectraLOG();
        histogramsLOG();
    }
#endif
}

#if perform_deltaN
// Mean-subtract and mask incomplete patches before writing final deltaN products.
void saveN([[maybe_unused]] FILE* output_log)
{
    std::vector<unsigned char> completed(deltaN.size(), 0);
    std::size_t completed_count = 0;
    double Nmean = 0.0;

    for (std::size_t id = 0; id < deltaN.size(); ++id) {
        if (!std::isfinite(f[id]) || !std::isfinite(deltaN[id])) {
            std::fprintf(stderr, "Non-finite state encountered while finalizing deltaN outputs.\n");
            std::exit(1);
        }
        if (!deltaN_patch_is_active(f[id])) {
            completed[id] = 1;
            ++completed_count;
            Nmean += deltaN[id];
        }
    }

    if (completed_count == 0) {
        std::fprintf(stderr,
            "No deltaN patch reached the selected hypersurface within Nend=%g.\n",
            Nend);
        std::exit(1);
    }
    Nmean /= static_cast<double>(completed_count);

#if post_inflation
    if (completed_count != deltaN.size()) {
        std::fprintf(stderr,
            "The post-inflationary stage requires every deltaN patch to reach the selected hypersurface "
            "(%zu of %zu completed). Increase Nend.\n",
            completed_count, deltaN.size());
        std::exit(1);
    }
#else
    if (completed_count != deltaN.size()) {
        std::fprintf(stderr,
            "Warning: %zu of %zu deltaN patches did not reach the selected hypersurface; "
            "they are excluded from the histogram and masked in spatial outputs.\n",
            deltaN.size() - completed_count, deltaN.size());
        std::fprintf(output_log,
            "Warning: %zu of %zu deltaN patches did not reach the selected hypersurface; "
            "they are excluded from the histogram and masked in spatial outputs.\n",
            deltaN.size() - completed_count, deltaN.size());
        std::fflush(output_log);
    }
#endif

    for (std::size_t id = 0; id < deltaN.size(); ++id) {
        deltaN[id] = completed[id] ? deltaN[id] - Nmean : 0.0;
    }

    histogramsN(completed);

    if (output_box2D)
    box2dN();

    if (output_spectra)
    spectraN();

}
#else
// No-op implementation keeps the common output interface linkable when deltaN is disabled.
void saveN(FILE*) {}
#endif


// -------------------- Post-inflationary outputs --------------------

#if post_inflation

// Write post-inflationary scalar means, variances, and mean velocity.
void meansvars_post_inflation(int flush)
{
    static FILE *means_, *vars_, *velocity_;
    DECLARE_INDICES

    double av, var, vel;

    static int first = 1;
    if (first) // Open output files
    {
        snprintf(name_, sizeof(name_), "results/post_inflation/means%s", ext_);
        means_ = open_output_or_die(name_, mode_);
        snprintf(name_, sizeof(name_), "results/post_inflation/variance%s", ext_);
        vars_ = open_output_or_die(name_, mode_);
        snprintf(name_, sizeof(name_), "results/post_inflation/velocity%s", ext_);
        velocity_ = open_output_or_die(name_, mode_);
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
    }
    const double inv_gridsize = 1.0 / static_cast<double>(gridsize);
    av  *= inv_gridsize;
    vel *= inv_gridsize;

    // Use the centered form for an accurate variance when fluctuations are small.
    LOOP var += pw2(f[idx(i,j,k)] - av);
    var *= inv_gridsize;

    vel = vel * pow(a, rescale_s - 1) * rescale_B;

    fprintf(means_, " %e", av);
    fprintf(velocity_, " %e", vel);
    fprintf(vars_, " %e", var);
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

// Write post-inflationary background expansion quantities.
void scale_post_inflation(int flush)
{
    static FILE *sf_;

    static int first = 1;
    if (first) // Open output file
    {
        snprintf(name_, sizeof(name_), "results/post_inflation/sf%s", ext_);
        sf_ = open_output_or_die(name_, mode_);
        first = 0;
    }

    const double acceleration =
        -(rescale_s - 0.5 * (1.0 - 3.0 * omega)) * pw2(ad) / a;

    // Output a, H, and adotdot in physical units using rescalings
    fprintf(sf_, "%f %f %e %e\n",
    t, a,
    ad * rescale_B * pow(a, rescale_s - 2.),
    pw2(rescale_B) * pow(a, 2. * rescale_s - 2.) *
    (acceleration + (rescale_s - 1.) * pw2(ad) / a));

    if (flush)
    fflush(sf_);
}

// Write the post-inflationary scalar spectrum and apply the optional live cutoff.
void spectraf_post_inflation()
{
    static FILE *spectra_, *spectratimes_;
    static int first = 1;

    if (first)
    {
        std::snprintf(name_, sizeof(name_), "results/post_inflation/spectra%s", ext_);
        spectra_ = open_output_or_die(name_, mode_);

        std::snprintf(name_, sizeof(name_), "results/post_inflation/spectratimes%s", ext_);
        spectratimes_ = open_output_or_die(name_, mode_);

        first = 0;
    }

    const double norm1 = std::pow(L / rescale_B, 3) / std::pow((double)N, 6);

    // Forward FFT to k-space for spectrum evaluation.


    int arraysize_spec[] = {N, N, N};


    fftrnd(f.data(), (double *)fnyquist_p, 3, arraysize_spec, 1);



    write_isotropic_spectrum_from_fft_r2c(spectra_, f, fnyquist_p, norm1);
    std::fprintf(spectratimes_, "%f %e\n", t, a);

    std::fflush(spectratimes_);
    // Optional hard cutoff in Fourier space (identical mode filtering logic to the legacy implementation).
    int i, j, k, px, py, pz, iconj, jconj;
    double pdisc;
    int arraysize[] = {N, N, N};

    if (high_cutoff_index > 0 && forcing_cutoff)
    {
        fftrnd(fd.data(), (double *)fdnyquist_p, 3, arraysize, 1);
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
                    pdisc = sqrt(pw2((double)px) + pw2((double)py) + pw2((double)pz));
                    if (pdisc > high_cutoff_index || pdisc < low_cutoff_index)
                    {
                        kill_mode(&f[idx(i,j,2 * k)], &fd[idx(i,j,2 * k)]);
                    }
                }

                if (j > N / 2 || (i > N / 2 && (j == 0 || j == N / 2)))
                {
                    jconj = (j == 0 ? 0 : N - j);

                    pdisc = sqrt(pw2((double)px) + pw2((double)py));
                    if (pdisc > high_cutoff_index || pdisc < low_cutoff_index)
                    {
                        kill_mode(&f[idx(i,j,0)], &fd[idx(i,j,0)]);
                        f[idx(iconj,jconj,0)] = f[idx(i,j,0)];
                        f[idx(iconj,jconj,1)] = -f[idx(i,j,1)];
                        fd[idx(iconj,jconj,0)] = fd[idx(i,j,0)];
                        fd[idx(iconj,jconj,1)] = -fd[idx(i,j,1)];
                    }

                    pdisc = sqrt(pw2((double)px) + pw2((double)py) + pw2((double)N / 2.));
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
                    pdisc = sqrt(pw2((double)px) + pw2((double)py));
                    if (pdisc > high_cutoff_index || pdisc < low_cutoff_index)
                    {
                        kill_mode(&f[idx(i,j,0)], &fd[idx(i,j,0)]);
                    }

                    pdisc = sqrt(pw2((double)px) + pw2((double)py) + pw2((double)N / 2.));
                    if (pdisc > high_cutoff_index || pdisc < low_cutoff_index)
                    {
                        kill_mode(&fnyquist_p[i][2 * j], &fdnyquist_p[i][2 * j]);
                    }
                }
            }
        }
        fftrnd(fd.data(), (double *)fdnyquist_p, 3, arraysize, -1);
    }

    // Backward FFT to restore the real-space field.
    fftrnd(f.data(), (double *)fnyquist_p, 3, arraysize, -1);

}

// Write the post-inflationary scalar one-point histogram and bin metadata.
void histograms_post_inflation()
{
    static FILE *histogram_, *histogramtimes_;
    int i = 0, j = 0, k = 0;
    int binnum;
    static std::vector<double> binfreq;
    double bmin, bmax, df;
    int numpts;

    static int first = 1;
    if (first) // Open output files
    {
        snprintf(name_, sizeof(name_), "results/post_inflation/histogram%s", ext_);
        histogram_ = open_output_or_die(name_, mode_);

        snprintf(name_, sizeof(name_), "results/post_inflation/histogramtimes%s", ext_);
        histogramtimes_ = open_output_or_die(name_, mode_);
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
    df = (bmax - bmin) / (double)(nbins); // bmin will be at the bottom of the first bin and bmax at the top of the last
    if (!std::isfinite(df) || df <= 0.0) df = 1.0;

    if ((int)binfreq.size() != nbins) binfreq.assign(nbins, 0.0);
    else std::fill(binfreq.begin(), binfreq.end(), 0.0);

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
    if (numpts == 0) numpts = 1;
    for (i = 0; i < nbins; i++)
    fprintf(histogram_, "%e\n", binfreq[i] / (double)numpts); // Output bin frequency normalized so the total equals 1
    fprintf(histogram_, "\n"); // Stick a blank line between times to make the file more readable
    fflush(histogram_);
    fprintf(histogramtimes_, " %e %e", bmin, df); // Output the starting point and stepsize for the bins at each time

    fprintf(histogramtimes_, "\n");
    fflush(histogramtimes_);
}

// Dispatch post-inflationary outputs while respecting leapfrog staggering.
void save_post_inflation(int infrequent)
{
    if (post_inflation_uses_staggered_derivatives() && t > 0.) // Synchronize field values and derivatives
    apply_leapfrog_drift(-.5 * dt_post_inflation);

    meansvars_post_inflation(infrequent);
    scale_post_inflation(infrequent);

    // Infrequent calculations
    if (infrequent)
    {
        if (output_spectra)
        {
            spectraf_post_inflation();
#if calculate_SIGW
            spectraGW_post_inflation();
            spectraGWdot_post_inflation();
#endif
        }
        if (output_histogram)
        histograms_post_inflation();
    }

    if (post_inflation_uses_staggered_derivatives() && t > 0.) // Desynchronize field values and derivatives
    apply_leapfrog_drift(.5 * dt_post_inflation);
}

#endif
