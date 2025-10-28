void spectraGW()
{
    static FILE *spectraGW_ = nullptr;
    static int first = 1;

    const int   numbins   = (int)(std::sqrt(3.0) * (N/2)) + 1;
    const int   i_max_out = (int)std::floor(0.90 * numbins);   // keep up to 90% of Nyquist
    const float dp        = 2.f * (float)pi / (float)L;

    std::vector<int>   numpoints(numbins, 0);
    std::vector<float> p(numbins, 0.f), f2(numbins, 0.f);
    for (int i = 0; i < numbins; ++i) p[i] = dp * i;

    if (first) {
        snprintf(name_, sizeof(name_), "results/spectraGW%s", ext_);
        spectraGW_ = fopen(name_, mode_);
        first = 0;
    }

    // --- FFT to k-space (Nyquist planes in scratch) ---
    int arraysize[] = {N, N, N};
    for (int c = 0; c < 6; ++c) fftrnf(hij[c].data(), (float*)hijnyquist_p[c], 3, arraysize, 1);

    // --- loop over modes ---
    for (int i = 0; i < N; ++i) {
        int px = (i <= N/2 ? i : i - N);
        for (int j = 0; j < N; ++j) {
            int py = (j <= N/2 ? j : j - N);

            // interior (double-counted)
            for (int k = 1; k < N/2; ++k) {
                int pz = k;

                // lattice k and |k|
                float kx = (2.f/(float)dx) * std::sin((float)pi * px / (float)N);
                float ky = (2.f/(float)dx) * std::sin((float)pi * py / (float)N);
                float kz = (2.f/(float)dx) * std::sin((float)pi * pz / (float)N);
                float kt2 = kx*kx + ky*ky + kz*kz;
                if (kt2 == 0.f) continue;

                int bin = (int)lroundf(std::sqrt((float)(px*px + py*py + pz*pz)));
                if (bin >= i_max_out) continue;

                // read C_ij = h_ij(k) (complex) from hij (interior)
                int idx_mode = idx(i, j, 2*k);
                float C_re[3][3] = {{0}}, C_im[3][3] = {{0}};
                for (int l = 0; l < 3; ++l) for (int m = l; m < 3; ++m) {
                    int comp = sym_idx(l, m);
                    float re = hij[comp][idx_mode];
                    float im = hij[comp][idx_mode + 1];
                    C_re[l][m] = re; C_im[l][m] = im;
                    if (m != l) { C_re[m][l] = re; C_im[m][l] = im; }
                }

                // projector pieces
                float inv = 1.0f / std::sqrt(kt2);
                float kh[3] = {kx*inv, ky*inv, kz*inv};
                float P[3][3];
                for (int a=0;a<3;++a) for (int b=0;b<3;++b)
                    P[a][b] = (a==b?1.f:0.f) - kh[a]*kh[b];

                // A = P C P  and  T = Tr(P C)
                float A_re[3][3] = {{0}}, A_im[3][3] = {{0}};
                float T_re = 0.f, T_im = 0.f;
                for (int a=0;a<3;++a) for (int b=0;b<3;++b) {
                    float r=0.f, im=0.f;
                    for (int l=0;l<3;++l) for (int m=0;m<3;++m) {
                        float pal = P[a][l], pbm = P[b][m];
                        r  += pal * C_re[l][m] * pbm;
                        im += pal * C_im[l][m] * pbm;
                    }
                    A_re[a][b]=r; A_im[a][b]=im;
                }
                for (int l=0;l<3;++l) for (int m=0;m<3;++m) {
                    T_re += P[l][m]*C_re[l][m];
                    T_im += P[l][m]*C_im[l][m];
                }

                // Λ–contraction: sum |A - 0.5 P T|^2
                float fp2 = 0.f;
                for (int a=0;a<3;++a) for (int b=0;b<3;++b) {
                    float r  = A_re[a][b] - 0.5f*P[a][b]*T_re;
                    float im = A_im[a][b] - 0.5f*P[a][b]*T_im;
                    fp2 += r*r + im*im;
                }

                numpoints[bin] += 2;
                f2[bin]        += 2.f * fp2;
            }

            // 0 and Nyquist planes (single-counted)
            for (int k = 0; k <= N/2; k += N/2) {
                int pz = k;

                float kx = (2.f/(float)dx) * std::sin((float)pi * px / (float)N);
                float ky = (2.f/(float)dx) * std::sin((float)pi * py / (float)N);
                float kz = (2.f/(float)dx) * std::sin((float)pi * pz / (float)N);
                float kt2 = kx*kx + ky*ky + kz*kz;
                if (kt2 == 0.f) continue;

                int bin = (int)lroundf(std::sqrt((float)(px*px + py*py + pz*pz)));
                if (bin >= i_max_out) continue;

                float C_re[3][3] = {{0}}, C_im[3][3] = {{0}};
                if (k == 0) {
                    int idx_mode = idx(i, j, 0);
                    for (int l = 0; l < 3; ++l) for (int m = l; m < 3; ++m) {
                        int comp = sym_idx(l, m);
                        float re = hij[comp][idx_mode];
                        float im = hij[comp][idx_mode + 1];
                        C_re[l][m]=re; C_im[l][m]=im;
                        if (m!=l){ C_re[m][l]=re; C_im[m][l]=im; }
                    }
                } else { // Nyquist plane stored in scratch
                    int j2 = 2*j;
                    for (int l = 0; l < 3; ++l) for (int m = l; m < 3; ++m) {
                        int comp = sym_idx(l, m);
                        float re = hijnyquist_p[comp][i][j2];
                        float im = hijnyquist_p[comp][i][j2 + 1];
                        C_re[l][m]=re; C_im[l][m]=im;
                        if (m!=l){ C_re[m][l]=re; C_im[m][l]=im; }
                    }
                }

                float inv = 1.0f / std::sqrt(kt2);
                float kh[3] = {kx*inv, ky*inv, kz*inv};
                float P[3][3];
                for (int a=0;a<3;++a) for (int b=0;b<3;++b)
                    P[a][b] = (a==b?1.f:0.f) - kh[a]*kh[b];

                float A_re[3][3] = {{0}}, A_im[3][3] = {{0}};
                float T_re = 0.f, T_im = 0.f;
                for (int a=0;a<3;++a) for (int b=0;b<3;++b) {
                    float r=0.f, im=0.f;
                    for (int l=0;l<3;++l) for (int m=0;m<3;++m) {
                        float pal=P[a][l], pbm=P[b][m];
                        r  += pal*C_re[l][m]*pbm;
                        im += pal*C_im[l][m]*pbm;
                    }
                    A_re[a][b]=r; A_im[a][b]=im;
                }
                for (int l=0;l<3;++l) for (int m=0;m<3;++m) {
                    T_re += P[l][m]*C_re[l][m];
                    T_im += P[l][m]*C_im[l][m];
                }

                float fp2 = 0.f;
                for (int a=0;a<3;++a) for (int b=0;b<3;++b) {
                    float r  = A_re[a][b] - 0.5f*P[a][b]*T_re;
                    float im = A_im[a][b] - 0.5f*P[a][b]*T_im;
                    fp2 += r*r + im*im;
                }

                numpoints[bin] += 1;
                f2[bin]        += fp2;
            }
        }
    }

    // write (normalization: keep your convention if you had one)
    const float norm1 = std::pow(L / rescale_B, 3) / std::pow((float)N, 6); // replace with your previous norm if desired
    for (int i = 0; i < numbins; ++i) {
        if (numpoints[i] > 0) f2[i] /= numpoints[i];
        fprintf(spectraGW_, "%e %d %e\n", p[i], numpoints[i], norm1 * f2[i]);
    }
    fprintf(spectraGW_, "\n");
    fflush(spectraGW_);

    // restore to real space
    for (int c = 0; c < 6; ++c) fftrnf(hij[c].data(), (float*)hijnyquist_p[c], 3, arraysize, -1);
}

// --------------- a^4 P_{\dot h}(k) with Λ–contraction ---------------------
void spectraGW_dot()
{
    static FILE *spectraGWdot_ = nullptr;
    static int first = 1;

    const int   numbins   = (int)(std::sqrt(3.0) * (N/2)) + 1;
    const int   i_max_out = (int)std::floor(0.90 * numbins);   // keep up to 90% of Nyquist
    const float dp        = 2.f * (float)pi / (float)L;

    std::vector<int>   numpoints(numbins, 0);
    std::vector<float> p(numbins, 0.f), f2(numbins, 0.f);
    for (int i = 0; i < numbins; ++i) p[i] = dp * i;

    if (first) {
        snprintf(name_, sizeof(name_), "results/spectraGWdot%s", ext_);
        spectraGWdot_ = fopen(name_, mode_);
        first = 0;
    }

    // FFT to k-space for hdot_ij  (float path)
    int arraysize[] = {N, N, N};
    for (int c = 0; c < 6; ++c) fftrnf(hijd[c].data(), (float*)hijdnyquist_p[c], 3, arraysize, 1);

    // conformal -> physical time derivative factor (your convention)
    const double to_phys = rescale_B * std::pow(a, rescale_s - 1.0);

    for (int i = 0; i < N; ++i) {
        int px = (i <= N/2 ? i : i - N);
        for (int j = 0; j < N; ++j) {
            int py = (j <= N/2 ? j : j - N);

            for (int k = 1; k < N/2; ++k) {
                int pz = k;

                float kx = (2.f/(float)dx) * std::sin((float)pi * px / (float)N);
                float ky = (2.f/(float)dx) * std::sin((float)pi * py / (float)N);
                float kz = (2.f/(float)dx) * std::sin((float)pi * pz / (float)N);
                float kt2 = kx*kx + ky*ky + kz*kz;
                if (kt2 == 0.f) continue;

                int bin = (int)lroundf(std::sqrt((float)(px*px + py*py + pz*pz)));
                if (bin >= i_max_out) continue;

                int idx_mode = idx(i, j, 2*k);
                float C_re[3][3] = {{0}}, C_im[3][3] = {{0}};
                for (int l = 0; l < 3; ++l) for (int m = l; m < 3; ++m) {
                    int comp = sym_idx(l, m);
                    // multiply by to_phys (double) then cast to float once
                    float re = (float)( (double)hijd[comp][idx_mode]     * to_phys );
                    float im = (float)( (double)hijd[comp][idx_mode + 1] * to_phys );
                    C_re[l][m]=re; C_im[l][m]=im;
                    if (m!=l){ C_re[m][l]=re; C_im[m][l]=im; }
                }

                float inv = 1.0f / std::sqrt(kt2);
                float kh[3] = {kx*inv, ky*inv, kz*inv};
                float P[3][3];
                for (int a=0;a<3;++a) for (int b=0;b<3;++b)
                    P[a][b] = (a==b?1.f:0.f) - kh[a]*kh[b];

                float A_re[3][3] = {{0}}, A_im[3][3] = {{0}};
                float T_re = 0.f, T_im = 0.f;
                for (int a=0;a<3;++a) for (int b=0;b<3;++b) {
                    float r=0.f, im=0.f;
                    for (int l=0;l<3;++l) for (int m=0;m<3;++m) {
                        float pal=P[a][l], pbm=P[b][m];
                        r  += pal*C_re[l][m]*pbm;
                        im += pal*C_im[l][m]*pbm;
                    }
                    A_re[a][b]=r; A_im[a][b]=im;
                }
                for (int l=0;l<3;++l) for (int m=0;m<3;++m) {
                    T_re += P[l][m]*C_re[l][m];
                    T_im += P[l][m]*C_im[l][m];
                }

                float fp2 = 0.f;
                for (int a=0;a<3;++a) for (int b=0;b<3;++b) {
                    float r = A_re[a][b] - 0.5f*P[a][b]*T_re;
                    float im= A_im[a][b] - 0.5f*P[a][b]*T_im;
                    fp2 += r*r + im*im;
                }

                numpoints[bin] += 2;
                f2[bin]        += 2.f * fp2;
            }

            for (int k = 0; k <= N/2; k += N/2) {
                int pz = k;

                float kx = (2.f/(float)dx) * std::sin((float)pi * px / (float)N);
                float ky = (2.f/(float)dx) * std::sin((float)pi * py / (float)N);
                float kz = (2.f/(float)dx) * std::sin((float)pi * pz / (float)N);
                float kt2 = kx*kx + ky*ky + kz*kz;
                if (kt2 == 0.f) continue;

                int bin = (int)lroundf(std::sqrt((float)(px*px + py*py + pz*pz)));
                if (bin >= i_max_out) continue;

                float C_re[3][3] = {{0}}, C_im[3][3] = {{0}};
                if (k == 0) {
                    int idx_mode = idx(i, j, 0);
                    for (int l = 0; l < 3; ++l) for (int m = l; m < 3; ++m) {
                        int comp = sym_idx(l, m);
                        float re = (float)( (double)hijd[comp][idx_mode]     * to_phys );
                        float im = (float)( (double)hijd[comp][idx_mode + 1] * to_phys );
                        C_re[l][m]=re; C_im[l][m]=im;
                        if (m!=l){ C_re[m][l]=re; C_im[m][l]=im; }
                    }
                } else {
                    int j2 = 2*j;
                    for (int l = 0; l < 3; ++l) for (int m = l; m < 3; ++m) {
                        int comp = sym_idx(l, m);
                        float re = (float)( (double)hijdnyquist_p[comp][i][j2]     * to_phys );
                        float im = (float)( (double)hijdnyquist_p[comp][i][j2 + 1] * to_phys );
                        C_re[l][m]=re; C_im[l][m]=im;
                        if (m!=l){ C_re[m][l]=re; C_im[m][l]=im; }
                    }
                }

                float inv = 1.0f / std::sqrt(kt2);
                float kh[3] = {kx*inv, ky*inv, kz*inv};
                float P[3][3];
                for (int a=0;a<3;++a) for (int b=0;b<3;++b)
                    P[a][b] = (a==b?1.f:0.f) - kh[a]*kh[b];

                float A_re[3][3] = {{0}}, A_im[3][3] = {{0}};
                float T_re = 0.f, T_im = 0.f;
                for (int a=0;a<3;++a) for (int b=0;b<3;++b) {
                    float r=0.f, im=0.f;
                    for (int l=0;l<3;++l) for (int m=0;m<3;++m) {
                        float pal=P[a][l], pbm=P[b][m];
                        r  += pal*C_re[l][m]*pbm;
                        im += pal*C_im[l][m]*pbm;
                    }
                    A_re[a][b]=r; A_im[a][b]=im;
                }
                for (int l=0;l<3;++l) for (int m=0;m<3;++m) {
                    T_re += P[l][m]*C_re[l][m];
                    T_im += P[l][m]*C_im[l][m];
                }

                float fp2 = 0.f;
                for (int a=0;a<3;++a) for (int b=0;b<3;++b) {
                    float r = A_re[a][b] - 0.5f*P[a][b]*T_re;
                    float im= A_im[a][b] - 0.5f*P[a][b]*T_im;
                    fp2 += r*r + im*im;
                }

                numpoints[bin] += 1;
                f2[bin]        += fp2;
            }
        }
    }

    const float norm1 = std::pow(L / rescale_B, 3) / std::pow((float)N, 6); 
    for (int i = 0; i < numbins; ++i) {
        if (numpoints[i] > 0) f2[i] /= numpoints[i];
        fprintf(spectraGWdot_, "%e %d %e\n", p[i], numpoints[i], norm1 * f2[i]);
    }
    fprintf(spectraGWdot_, "\n");
    fflush(spectraGWdot_);

    // back to real space (float path)
    for (int c = 0; c < 6; ++c) fftrnf(hijd[c].data(), (float*)hijdnyquist_p[c], 3, arraysize, -1);
}
