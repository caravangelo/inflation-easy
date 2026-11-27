// Appended to post_inflation.cpp in a disposable build to exercise its private
// kick directly. All drifts, kicks, RK stages and output use production code.
int main(int argc, char** argv) {
    if (argc != 4) return 2;
    const std::string test = argv[1];
    const bool use_rk = std::string(argv[2]) == "rk4";
    const int steps = std::stoi(argv[3]);
    const bool inflation = test == "inflation_tensor" || test == "metric";
    const double duration = inflation ? 0.4 : 1.0;
    const double step = duration / steps;
    t = 0.0;
    a = ad = rescale_B = 1.0;
    rescale_s = 0.0;
    omega = test == "background" ? 0.2 : 1.0 / 3.0;
    L = 2 * pi;
    dx = L / N;
    linear_metric_perturbations = test == "metric";
    field_numerical = {8, 4, 0, -4};
    potential_numerical.assign(4, test == "metric" ? 1.875 : 3.0);
    potential_derivative_numerical.assign(4, 0.0);
    lstart.assign(gridsize, 1);
    f.assign(gridsize, 2.0);
    fd.assign(gridsize, test == "metric" ? -1.5 : 0.0);
    for (int c = 0; c < 6; ++c) {
        hij[c].assign(gridsize, 0.0f);
        hijd[c].assign(gridsize, (test == "source" || test == "midpoint_source") ? 0.0f : c + 1.0f);
    }
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                const auto id = idx(i,j,k);
                const double mode = std::cos(2*pi*(i+j)/N);
                if (test == "metric") f[id] += 0.1 * mode;
                if (test == "radiation") fd[id] = 1.0;
                if (test == "midpoint_source") {
                    f[id] = 0.0;
                    fd[id] = std::sin(2*pi*(i+j+k)/N);
                }
                if (test == "source") {
                    f[id] = 0.1 * mode;
                    fd[id] = 0.15 * mode;
                }
            }
        }
    }

    if (test == "midpoint_source") {
        // Start at the kick, with zero fields and a nonzero velocity gradient.
        // Thus only the grad(Phi') product can generate a tensor velocity.
        apply_post_inflation_leapfrog_kick(step);
        for (int c = 0; c < 6; ++c) std::printf("%.17g ", hijd[c][0]);
        std::puts("");
        return 0;
    }
    if (test == "stop") {
        deltaN.assign(gridsize, 0.0);
        horizon_factor = 1.0 / N; // initialize_post_inflation sets a' = 1.
        dt_post_inflation = step;
        af_post_inflation = 1 + step / 4;
        output_histogram = output_spectra = 0;
        post_inflation_integrator = INTEGRATOR_LEAPFROG;
        run_post_inflation_loop(stdout);
        return 0;
    }

    evolution::FieldRKWorkspace workspace;
    const auto phase = inflation ? evolution::SimulationPhase::Inflation
                                 : evolution::SimulationPhase::PostInflation;
    if (!use_rk) apply_leapfrog_drift(step / 2);
    for (int n = 0; n < steps; ++n) {
        if (use_rk) {
            evolution::rk4_step_fields(step, workspace, phase);
        } else {
            if (inflation) evolution::apply_inflation_leapfrog_kick(step);
            else apply_post_inflation_leapfrog_kick(step);
            apply_leapfrog_drift(n + 1 == steps ? step / 2 : step);
        }
    }
    if (test == "metric") {
        // Difference at opposite phases cancels the homogeneous motion.
        std::printf("%.17g\n", 0.5 * (fd[0] - fd[idx(N/2,0,0)]));
    } else if (test == "source") {
        // A diagonal scalar wave drives the xy component. Its initial h and
        // h' vanish, so a measurable signal requires an active scalar source.
        std::printf("%.17g %.17g\n", hij[3][0], hijd[3][0]);
    } else if (test == "background") {
        std::printf("%.17g %.17g\n", a, ad);
    } else {
        for (int c = 0; c < 6; ++c) std::printf("%.17g ", hijd[c][0] / (c+1.0));
        if (test == "radiation") std::printf("%.17g", fd[0]);
        std::puts("");
    }
}
