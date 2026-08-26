// runtime_parameters.cpp - Run-time configuration defaults and parser
//
// Values in this file provide a complete fallback configuration. At program
// startup, load_runtime_parameters() applies recognized key-value overrides
// from params.txt, recomputes derived quantities, and sanitizes integrator
// controls. Compile-time switches remain exclusively in parameters.h.

#include "parameters.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <string>

// -------------------- Defaults --------------------

int seed = 8;

double rescale_s = 0.0;
// Default inflation stop scale factor.
double af = 2 * N;

#if numerical_potential
double V0 = 3e-9;
#else
double V0 = 3.338e-13;
#endif

double rescale_B = std::sqrt(V0);

int linear_metric_perturbations = 0;
std::optional<double> initial_mass_squared;

#if !numerical_potential
double ns = 0.97;
#endif

#if numerical_potential
double initial_field = 2.9181235049318586;
double initial_derivative = -0.06727651095116181;
double L = 10.0;
double dt = 0.001;
int output_freq = 500;
int output_infrequent_freq = 500;
#if perform_deltaN
double dN = 0.0001;
double Nend = 5.0;
int use_phiref_manual = 0;
double phiref_manual_value = 0.0;
#endif
#if post_inflation
double horizon_factor = 1.0;
double omega = 1.0 / 3.0;
double dt_post_inflation = 0.001;
double af_post_inflation = 2.0 * N;
#endif
#else
double initial_field = 0.0935;
double initial_derivative = 0.000796;
double L = 10.0;
double dt = 0.0005;
int output_freq = 200;
int output_infrequent_freq = 200;
#if perform_deltaN
double dN = 0.0000001;
double Nend = 0.001;
int use_phiref_manual = 0;
double phiref_manual_value = 0.0;
#endif
#if post_inflation
double horizon_factor = 1.0;
double omega = 1.0 / 3.0;
double dt_post_inflation = 0.001;
double af_post_inflation = 2.0 * N;
#endif
#endif

int integrator = INTEGRATOR_LEAPFROG;
#if perform_deltaN
int deltaN_integrator = INTEGRATOR_LEAPFROG;
#endif
#if post_inflation
int post_inflation_integrator = INTEGRATOR_LEAPFROG;
#endif
double rk45_abs_tol = 1e-8;
double rk45_rel_tol = 1e-6;
double rk45_min_dt = -1.0;
double rk45_max_dt = -1.0;
double rk45_safety = 0.9;

double high_cutoff_index = 0.0;
double low_cutoff_index = 0.0;
int forcing_cutoff = 0;

int output_spectra = 1;
int output_histogram = 1;
int output_energy = 1;
int output_box3D = 0;
int output_box2D = 0;
int output_bispectrum = 0;

#if perform_deltaN
int output_LOG = 0;
double eta_log = -0.5;
#endif

int screen_updates = 1;
int nbins = 256;

#if numerical_potential
int int_err = 5;
int int_errN = 5;
#endif

double dx = L / static_cast<double>(N);

namespace {
// Remove leading and trailing ASCII whitespace from a parser token.
std::string trim(const std::string& s) {
    size_t b = 0;
    while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
    size_t e = s.size();
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
    return s.substr(b, e - b);
}

// Parse a complete integer token; partial conversions are rejected.
bool parse_int(const std::string& s, int& out) {
    char* end = nullptr;
    errno = 0;
    long v = std::strtol(s.c_str(), &end, 10);
    if (errno != 0 || !end || end == s.c_str() || *end != '\0'
        || v < std::numeric_limits<int>::min()
        || v > std::numeric_limits<int>::max()) return false;
    out = static_cast<int>(v);
    return true;
}

// Parse a complete floating-point token; partial conversions are rejected.
bool parse_double(const std::string& s, double& out) {
    char* end = nullptr;
    errno = 0;
    double v = std::strtod(s.c_str(), &end);
    if (errno != 0 || !end || end == s.c_str() || *end != '\0') return false;
    out = v;
    return true;
}

// Accept documented integrator names, common aliases, or the corresponding enum value.
bool parse_integrator(const std::string& raw, int& out) {
    std::string s = raw;
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });

    if (s == "leapfrog" || s == "lf") {
        out = INTEGRATOR_LEAPFROG;
        return true;
    }
    if (s == "rk4") {
        out = INTEGRATOR_RK4;
        return true;
    }
    if (s == "rk45" || s == "rkf45" || s == "dopri5") {
        out = INTEGRATOR_RK45;
        return true;
    }

    int iv = 0;
    if (parse_int(s, iv) && iv >= INTEGRATOR_LEAPFROG && iv <= INTEGRATOR_RK45) {
        out = iv;
        return true;
    }
    return false;
}

// Fall back to leapfrog if an integrator value is outside the supported enum.
void sanitize_integrator_choice(int& value) {
    if (value < INTEGRATOR_LEAPFROG || value > INTEGRATOR_RK45) {
        value = INTEGRATOR_LEAPFROG;
    }
}

// Apply the same validation to every simulation stage compiled into the executable.
void sanitize_all_integrator_choices() {
    sanitize_integrator_choice(integrator);
#if perform_deltaN
    sanitize_integrator_choice(deltaN_integrator);
#endif
#if post_inflation
    sanitize_integrator_choice(post_inflation_integrator);
#endif
}

// Enforce positive tolerances and a consistent allowed RK45 step interval.
void sanitize_rk45_controls() {
    if (!(rk45_min_dt > 0.0)) rk45_min_dt = std::abs(dt) * 1e-6;
    if (!(rk45_max_dt > 0.0)) rk45_max_dt = std::abs(dt);
    if (rk45_min_dt <= 0.0) rk45_min_dt = 1e-16;
    if (rk45_max_dt < rk45_min_dt) rk45_max_dt = rk45_min_dt;
    if (!(rk45_abs_tol > 0.0)) rk45_abs_tol = 1e-8;
    if (!(rk45_rel_tol > 0.0)) rk45_rel_tol = 1e-6;
    if (!(rk45_safety > 0.0 && rk45_safety < 1.0)) rk45_safety = 0.9;
}

[[noreturn]] void invalid_runtime_configuration(const char* message) {
    std::fprintf(stderr, "Invalid run-time configuration: %s\n", message);
    std::exit(EXIT_FAILURE);
}

void require_finite(const char* name, double value) {
    if (!std::isfinite(value)) {
        std::fprintf(stderr, "Invalid run-time configuration: %s must be finite.\n", name);
        std::exit(EXIT_FAILURE);
    }
}

void require_positive(const char* name, double value) {
    require_finite(name, value);
    if (!(value > 0.0)) {
        std::fprintf(stderr, "Invalid run-time configuration: %s must be positive.\n", name);
        std::exit(EXIT_FAILURE);
    }
}

// Validate the resolved configuration once, before any lattice allocation or evolution.
// Bounds here are limited to values required for well-defined arithmetic and loops.
void validate_runtime_parameters() {
    require_finite("rescale_s", rescale_s);
    const double singularity_tolerance = 16.0 * std::numeric_limits<double>::epsilon();
    if (integrator == INTEGRATOR_LEAPFROG
        && std::abs(rescale_s + 1.0) <= singularity_tolerance) {
        invalid_runtime_configuration(
            "rescale_s is too close to -1 for the leapfrog scale-factor update; use RK4/RK45 or another rescaling.");
    }
    require_positive("V0", V0);
    require_positive("rescale_B", rescale_B);
    require_finite("initial_field", initial_field);
    require_finite("initial_derivative", initial_derivative);
    if (initial_mass_squared.has_value()) {
        require_finite("initial_mass_squared", *initial_mass_squared);
    }
#if !numerical_potential
    require_finite("ns", ns);
#endif

    require_positive("L", L);
    require_positive("dt", dt);
    require_positive("dx", dx);
    require_finite("af", af);
    if (af < 1.0) {
        invalid_runtime_configuration("af must be at least 1, the initial scale factor.");
    }

    if (output_freq <= 0) {
        invalid_runtime_configuration("output_freq must be a positive number of steps.");
    }
    if (output_infrequent_freq <= 0) {
        invalid_runtime_configuration("output_infrequent_freq must be a positive number of steps.");
    }
    if (nbins <= 0) {
        invalid_runtime_configuration("nbins must be positive.");
    }
    const std::size_t lattice_sites =
        static_cast<std::size_t>(N) * static_cast<std::size_t>(N) * static_cast<std::size_t>(N);
    if (static_cast<std::size_t>(nbins) > lattice_sites) {
        invalid_runtime_configuration("nbins must not exceed the number of lattice sites N^3.");
    }

    require_positive("rk45_abs_tol", rk45_abs_tol);
    require_positive("rk45_rel_tol", rk45_rel_tol);
    require_positive("rk45_min_dt", rk45_min_dt);
    require_positive("rk45_max_dt", rk45_max_dt);
    require_finite("rk45_safety", rk45_safety);
    if (!(rk45_safety > 0.0 && rk45_safety < 1.0)) {
        invalid_runtime_configuration("rk45_safety must lie strictly between 0 and 1.");
    }
    if (rk45_min_dt > rk45_max_dt) {
        invalid_runtime_configuration("rk45_min_dt must not exceed rk45_max_dt.");
    }

    require_finite("low_cutoff_index", low_cutoff_index);
    require_finite("high_cutoff_index", high_cutoff_index);
    if (low_cutoff_index < 0.0 || high_cutoff_index < 0.0) {
        invalid_runtime_configuration("cutoff indices must be nonnegative.");
    }
    if (high_cutoff_index > 0.0 && low_cutoff_index > high_cutoff_index) {
        invalid_runtime_configuration("low_cutoff_index must not exceed high_cutoff_index.");
    }

#if perform_deltaN
    require_finite("dN", dN);
    if (dN == 0.0) {
        invalid_runtime_configuration("dN must be nonzero; use its sign to select the integration direction.");
    }
    require_finite("Nend", Nend);
    if (Nend < 0.0) {
        invalid_runtime_configuration("Nend is a nonnegative integration magnitude.");
    }
    if (dN < 0.0 && !use_phiref_manual) {
        invalid_runtime_configuration("backward deltaN evolution (dN < 0) requires use_phiref_manual = 1.");
    }
    if (use_phiref_manual) {
        require_finite("phiref_manual_value", phiref_manual_value);
    }
    require_finite("eta_log", eta_log);
    if (output_LOG && eta_log == 0.0) {
        invalid_runtime_configuration("eta_log must be nonzero when output_LOG is enabled.");
    }
#endif

#if post_inflation
    require_positive("horizon_factor", horizon_factor);
    require_finite("omega", omega);
    const double one_plus_omega = 1.0 + omega;
    const double five_plus_three_omega = 5.0 + 3.0 * omega;
    const double denominator_tolerance = 16.0 * std::numeric_limits<double>::epsilon();
    if (std::abs(one_plus_omega) <= denominator_tolerance) {
        invalid_runtime_configuration("omega is too close to -1, where the post-inflation equations are singular.");
    }
    if (std::abs(five_plus_three_omega) <= denominator_tolerance) {
        invalid_runtime_configuration("omega is too close to -5/3, where the zeta-to-Phi map is singular.");
    }
    const double Phi_from_zeta = 3.0 * one_plus_omega / five_plus_three_omega;
    if (!std::isfinite(Phi_from_zeta)) {
        invalid_runtime_configuration("omega produces a non-finite zeta-to-Phi conversion factor.");
    }
    require_positive("dt_post_inflation", dt_post_inflation);
    require_finite("af_post_inflation", af_post_inflation);
    if (af_post_inflation < 1.0) {
        invalid_runtime_configuration("af_post_inflation must be at least 1, the post-inflation initial scale factor.");
    }
#endif
}
} // namespace

// Load optional key-value overrides. Unknown, malformed, or compile-time-only
// entries are reported and ignored; defaults remain active for omitted keys.
void load_runtime_parameters(const char* filename) {
    std::ifstream in(filename);
    if (!in.good()) {
        // Optional file: keep defaults if not present.
        af = 2.0 * N;
        dx = L / static_cast<double>(N);
#if post_inflation
        af_post_inflation = 2.0 * N;
#endif
        sanitize_rk45_controls();
        sanitize_all_integrator_choices();
        validate_runtime_parameters();
        return;
    }

    bool rescale_B_overridden = false;
    bool af_overridden = false;
    bool af_post_overridden = false;
    bool rk45_min_overridden = false;
    bool rk45_max_overridden = false;
#if !perform_deltaN
    bool ignored_deltaN_runtime_keys = false;
#endif
#if !post_inflation
    bool ignored_post_inflation_runtime_keys = false;
#endif

    std::string line;
    int lineno = 0;
    while (std::getline(in, line)) {
        ++lineno;
        const auto hash = line.find('#');
        if (hash != std::string::npos) line.erase(hash);
        line = trim(line);
        if (line.empty()) continue;

        const auto eq = line.find('=');
        if (eq == std::string::npos) {
            std::fprintf(stderr, "Ignoring malformed line %d in %s: %s\n", lineno, filename, line.c_str());
            continue;
        }

        std::string key = trim(line.substr(0, eq));
        std::string val = trim(line.substr(eq + 1));

        int ival = 0;
        double dval = 0.0;

        if (key == "seed" && parse_int(val, ival)) seed = ival;
        else if (key == "rescale_s" && parse_double(val, dval)) rescale_s = dval;
        else if (key == "af" && parse_double(val, dval)) { af = dval; af_overridden = true; }
        else if (key == "V0" && parse_double(val, dval)) V0 = dval;
        else if (key == "rescale_B" && parse_double(val, dval)) { rescale_B = dval; rescale_B_overridden = true; }
#if !numerical_potential
        else if (key == "ns" && parse_double(val, dval)) ns = dval;
#endif
        else if (key == "initial_field" && parse_double(val, dval)) initial_field = dval;
        else if (key == "initial_derivative" && parse_double(val, dval)) initial_derivative = dval;
        else if (key == "initial_mass_squared" && parse_double(val, dval)
              && std::isfinite(dval)) {
            initial_mass_squared = dval;
        }
        else if (key == "linear_metric_perturbations" && parse_int(val, ival)) {
            linear_metric_perturbations = (ival != 0);
        }
        else if (key == "L" && parse_double(val, dval)) L = dval;
        else if (key == "dt" && parse_double(val, dval)) dt = dval;
        else if (key == "inflation_integrator") {
            int parsed_integrator = INTEGRATOR_LEAPFROG;
            if (parse_integrator(val, parsed_integrator)) {
                integrator = parsed_integrator;
            } else {
                std::fprintf(stderr, "Ignoring invalid inflation_integrator '%s' on line %d.\n", val.c_str(), lineno);
            }
        }
#if perform_deltaN
        else if (key == "deltaN_integrator") {
            int parsed_integrator = INTEGRATOR_LEAPFROG;
            if (parse_integrator(val, parsed_integrator)) {
                deltaN_integrator = parsed_integrator;
            } else {
                std::fprintf(stderr, "Ignoring invalid deltaN_integrator '%s' on line %d.\n", val.c_str(), lineno);
            }
        }
#endif
#if post_inflation
        else if (key == "post_inflation_integrator") {
            int parsed_integrator = INTEGRATOR_LEAPFROG;
            if (parse_integrator(val, parsed_integrator)) {
                post_inflation_integrator = parsed_integrator;
            } else {
                std::fprintf(stderr, "Ignoring invalid post_inflation_integrator '%s' on line %d.\n", val.c_str(), lineno);
            }
        }
#endif
        else if (key == "rk45_abs_tol" && parse_double(val, dval)) rk45_abs_tol = dval;
        else if (key == "rk45_rel_tol" && parse_double(val, dval)) rk45_rel_tol = dval;
        else if (key == "rk45_min_dt" && parse_double(val, dval)) { rk45_min_dt = dval; rk45_min_overridden = true; }
        else if (key == "rk45_max_dt" && parse_double(val, dval)) { rk45_max_dt = dval; rk45_max_overridden = true; }
        else if (key == "rk45_safety" && parse_double(val, dval)) rk45_safety = dval;
        else if (key == "output_freq" && parse_int(val, ival)) output_freq = ival;
        else if (key == "output_infrequent_freq" && parse_int(val, ival)) output_infrequent_freq = ival;
#if perform_deltaN
        else if (key == "dN" && parse_double(val, dval)) dN = dval;
        else if (key == "Nend" && parse_double(val, dval)) Nend = dval;
        else if (key == "use_phiref_manual" && parse_int(val, ival)) use_phiref_manual = (ival != 0);
        else if (key == "phiref_manual_value" && parse_double(val, dval)) phiref_manual_value = dval;
#endif
#if post_inflation
        else if (key == "horizon_factor" && parse_double(val, dval)) horizon_factor = dval;
        else if (key == "omega" && parse_double(val, dval)) omega = dval;
        else if (key == "dt_post_inflation" && parse_double(val, dval)) dt_post_inflation = dval;
        else if (key == "af_post_inflation" && parse_double(val, dval)) { af_post_inflation = dval; af_post_overridden = true; }
#endif
        else if (key == "high_cutoff_index" && parse_double(val, dval)) high_cutoff_index = dval;
        else if (key == "low_cutoff_index" && parse_double(val, dval)) low_cutoff_index = dval;
        else if (key == "forcing_cutoff" && parse_int(val, ival)) forcing_cutoff = (ival != 0);
        else if (key == "output_spectra" && parse_int(val, ival)) output_spectra = (ival != 0);
        else if (key == "output_histogram" && parse_int(val, ival)) output_histogram = (ival != 0);
        else if (key == "output_energy" && parse_int(val, ival)) output_energy = (ival != 0);
        else if (key == "output_box3D" && parse_int(val, ival)) output_box3D = (ival != 0);
        else if (key == "output_box2D" && parse_int(val, ival)) output_box2D = (ival != 0);
        else if (key == "output_bispectrum" && parse_int(val, ival)) output_bispectrum = (ival != 0);
#if perform_deltaN
        else if (key == "output_LOG" && parse_int(val, ival)) output_LOG = (ival != 0);
        else if (key == "eta_log" && parse_double(val, dval)) eta_log = dval;
#endif
        else if (key == "screen_updates" && parse_int(val, ival)) screen_updates = (ival != 0);
        else if (key == "nbins" && parse_int(val, ival)) nbins = ival;
#if numerical_potential
        else if (key == "int_err" && parse_int(val, ival)) int_err = std::max(1, ival);
        else if (key == "int_errN" && parse_int(val, ival)) int_errN = std::max(1, ival);
#endif
#if !perform_deltaN
        else if (key == "dN" || key == "Nend" || key == "use_phiref_manual"
              || key == "phiref_manual_value" || key == "output_LOG"
              || key == "eta_log" || key == "deltaN_integrator") {
            ignored_deltaN_runtime_keys = true;
        }
#endif
#if !post_inflation
        else if (key == "horizon_factor" || key == "omega" || key == "dt_post_inflation"
              || key == "af_post_inflation" || key == "post_inflation_integrator") {
            ignored_post_inflation_runtime_keys = true;
        }
#endif
        else if (key == "N" || key == "numerical_potential" || key == "perform_deltaN"
              || key == "calculate_SIGW" || key == "post_inflation" || key == "parallel_calculation"
              || key == "monotonic_potential" || key == "antimonotonic_potential") {
            std::fprintf(stderr, "Ignoring compile-time parameter '%s' in %s (requires recompilation).\n", key.c_str(), filename);
        } else {
            std::fprintf(stderr, "Ignoring unknown or invalid parameter '%s' on line %d.\n", key.c_str(), lineno);
        }
    }

    if (!rescale_B_overridden) rescale_B = std::sqrt(V0);
    if (!af_overridden) af = 2.0 * N;
    if (!af_post_overridden) {
#if post_inflation
        af_post_inflation = 2.0 * N;
#endif
    }
    dx = L / static_cast<double>(N);

    if (!rk45_min_overridden || rk45_min_dt <= 0.0) rk45_min_dt = std::abs(dt) * 1e-6;
    if (!rk45_max_overridden || rk45_max_dt <= 0.0) rk45_max_dt = std::abs(dt);
    sanitize_rk45_controls();
    sanitize_all_integrator_choices();
    validate_runtime_parameters();
#if !perform_deltaN
    if (ignored_deltaN_runtime_keys) {
        std::fprintf(stderr, "Ignoring deltaN runtime parameters in %s because perform_deltaN=0.\n", filename);
    }
#endif
#if !post_inflation
    if (ignored_post_inflation_runtime_keys) {
        std::fprintf(stderr, "Ignoring post-inflation runtime parameters in %s because post_inflation=0.\n", filename);
    }
#endif
}

// Stable text representation used in logs and reproducibility metadata.
const char* integrator_name() {
    switch (integrator) {
        case INTEGRATOR_LEAPFROG: return "leapfrog";
        case INTEGRATOR_RK4: return "rk4";
        case INTEGRATOR_RK45: return "rk45";
        default: return "unknown";
    }
}

#if perform_deltaN
// Stable text representation of the deltaN integrator.
const char* deltaN_integrator_name() {
    switch (deltaN_integrator) {
        case INTEGRATOR_LEAPFROG: return "leapfrog";
        case INTEGRATOR_RK4: return "rk4";
        case INTEGRATOR_RK45: return "rk45";
        default: return "unknown";
    }
}
#endif

#if post_inflation
// Stable text representation of the post-inflation integrator.
const char* post_inflation_integrator_name() {
    switch (post_inflation_integrator) {
        case INTEGRATOR_LEAPFROG: return "leapfrog";
        case INTEGRATOR_RK4: return "rk4";
        case INTEGRATOR_RK45: return "rk45";
        default: return "unknown";
    }
}
#endif
