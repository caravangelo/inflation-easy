// Linked with the production globals in a disposable 8^3 build.
void spectraf();
void spectraf_post_inflation();

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    const std::string phase = argv[1];
    forcing_cutoff = std::stoi(argv[2]);
    low_cutoff_index = 2;
    high_cutoff_index = 3;
    f.resize(gridsize);
    fd.resize(gridsize);
    // Give both fields a mean, a rejected mode, and a retained mode.
    for (int k = 0; k < N; ++k) {
        const double angle = 2 * pi * k / N;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                const auto id = idx(i,j,k);
                f[id] = 3 + std::cos(angle) + std::cos(2 * angle);
                fd[id] = -2 + std::sin(angle) + std::sin(2 * angle);
            }
        }
    }
    if (phase == "inflation") spectraf();
    else if (phase == "radiation") spectraf_post_inflation();
    else return 2;

    for (int k = 0; k < N; ++k) {
        const double angle = 2 * pi * k / N;
        const double expected_f = 3 + (!forcing_cutoff)*std::cos(angle) + std::cos(2 * angle);
        const double expected_fd = -2 + (!forcing_cutoff)*std::sin(angle) + std::sin(2 * angle);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                const auto id = idx(i,j,k);
                if (!(std::abs(f[id] - expected_f) < 1e-11 &&
                      std::abs(fd[id] - expected_fd) < 1e-11)) {
                    std::fprintf(stderr, "%s cutoff=%d (%d,%d,%d): f=%g fd=%g, expected %g %g\n",
                                 phase.c_str(), int(forcing_cutoff), i, j, k,
                                 f[id], fd[id], expected_f, expected_fd);
                    return 1;
                }
            }
        }
    }
    return 0;
}
