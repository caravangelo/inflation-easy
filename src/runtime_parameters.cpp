// runtime_parameters.cpp - Run-time configuration defaults and parser

#include "parameters.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
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
std::string trim(const std::string& s) {
    size_t b = 0;
    while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
    size_t e = s.size();
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
    return s.substr(b, e - b);
}

bool parse_int(const std::string& s, int& out) {
    char* end = nullptr;
    errno = 0;
    long v = std::strtol(s.c_str(), &end, 10);
    if (errno != 0 || !end || *end != '\0') return false;
    out = static_cast<int>(v);
    return true;
}

bool parse_double(const std::string& s, double& out) {
    char* end = nullptr;
    errno = 0;
    double v = std::strtod(s.c_str(), &end);
    if (errno != 0 || !end || *end != '\0') return false;
    out = v;
    return true;
}

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
} // namespace

void load_runtime_parameters(const char* filename) {
    std::ifstream in(filename);
    if (!in.good()) {
        // Optional file: keep defaults if not present.
        af = 2.0 * N;
        dx = L / static_cast<double>(N);
#if post_inflation
        af_post_inflation = 2.0 * N;
#endif
        if (!(rk45_min_dt > 0.0)) rk45_min_dt = std::abs(dt) * 1e-6;
        if (!(rk45_max_dt > 0.0)) rk45_max_dt = std::abs(dt);
        if (rk45_max_dt < rk45_min_dt) rk45_max_dt = rk45_min_dt;
        if (!(rk45_abs_tol > 0.0)) rk45_abs_tol = 1e-8;
        if (!(rk45_rel_tol > 0.0)) rk45_rel_tol = 1e-6;
        if (!(rk45_safety > 0.0 && rk45_safety < 1.0)) rk45_safety = 0.9;
        if (integrator < INTEGRATOR_LEAPFROG || integrator > INTEGRATOR_RK45) {
            integrator = INTEGRATOR_LEAPFROG;
        }
#if perform_deltaN
        if (deltaN_integrator < INTEGRATOR_LEAPFROG || deltaN_integrator > INTEGRATOR_RK45) {
            deltaN_integrator = INTEGRATOR_LEAPFROG;
        }
#endif
#if post_inflation
        if (post_inflation_integrator < INTEGRATOR_LEAPFROG || post_inflation_integrator > INTEGRATOR_RK45) {
            post_inflation_integrator = INTEGRATOR_LEAPFROG;
        }
#endif
        return;
    }

    bool rescale_B_overridden = false;
    bool af_overridden = false;
    bool af_post_overridden = false;
    bool rk45_min_overridden = false;
    bool rk45_max_overridden = false;

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
        else if (key == "L" && parse_double(val, dval)) L = dval;
        else if (key == "dt" && parse_double(val, dval)) dt = dval;
        else if (key == "integrator") {
            int parsed_integrator = INTEGRATOR_LEAPFROG;
            if (parse_integrator(val, parsed_integrator)) {
                integrator = parsed_integrator;
            } else {
                std::fprintf(stderr, "Ignoring invalid integrator '%s' on line %d.\n", val.c_str(), lineno);
            }
        }
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
    if (rk45_min_dt <= 0.0) rk45_min_dt = 1e-16;
    if (rk45_max_dt < rk45_min_dt) rk45_max_dt = rk45_min_dt;
    if (!(rk45_abs_tol > 0.0)) rk45_abs_tol = 1e-8;
    if (!(rk45_rel_tol > 0.0)) rk45_rel_tol = 1e-6;
    if (!(rk45_safety > 0.0 && rk45_safety < 1.0)) rk45_safety = 0.9;
    if (integrator < INTEGRATOR_LEAPFROG || integrator > INTEGRATOR_RK45) {
        integrator = INTEGRATOR_LEAPFROG;
    }
#if perform_deltaN
    if (deltaN_integrator < INTEGRATOR_LEAPFROG || deltaN_integrator > INTEGRATOR_RK45) {
        deltaN_integrator = INTEGRATOR_LEAPFROG;
    }
#endif
#if post_inflation
    if (post_inflation_integrator < INTEGRATOR_LEAPFROG || post_inflation_integrator > INTEGRATOR_RK45) {
        post_inflation_integrator = INTEGRATOR_LEAPFROG;
    }
#endif
}

const char* integrator_name() {
    switch (integrator) {
        case INTEGRATOR_LEAPFROG: return "leapfrog";
        case INTEGRATOR_RK4: return "rk4";
        case INTEGRATOR_RK45: return "rk45";
        default: return "unknown";
    }
}

#if perform_deltaN
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
const char* post_inflation_integrator_name() {
    switch (post_inflation_integrator) {
        case INTEGRATOR_LEAPFROG: return "leapfrog";
        case INTEGRATOR_RK4: return "rk4";
        case INTEGRATOR_RK45: return "rk45";
        default: return "unknown";
    }
}
#endif
