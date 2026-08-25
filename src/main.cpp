// main.cpp - Program entry point and top-level driver
//
// This file defines the global variables declared in main.h and implements
// the top-level control flow of the simulation: initialization, directory
// setup, output management, and dispatch to the appropriate evolution loops.


#include "main.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include <fstream>

// Field values and derivatives
std::vector<double> f;
std::vector<double> fd;
#if perform_deltaN
std::vector<double> deltaN;
#endif

#if calculate_SIGW
std::vector<float> hij[6];
std::vector<float> hijd[6];
#endif

double fnyquist_p[N][2*N], fdnyquist_p[N][2*N];
#if calculate_SIGW
float hijnyquist_p[6][N][2*N], hijdnyquist_p[6][N][2*N];
#endif

double t, t0;
double astep = 1.0, a = 1.0, ad = 0.0, ad2 = 0.0, ad0 = 0.0, Ne = 0.0;
double phiref = 0.0;       // Reference field value for deltaN calculation
double hubble_init = 0.0;
char mode_[10] = "w";
char ext_[500] = ".dat";

#if numerical_potential
std::vector<int>    lstart;
std::vector<double> field_numerical;
std::vector<double> potential_numerical;
std::vector<double> potential_derivative_numerical;

// Load the three tabulated inputs required by the numerical-potential branch.
void load_numerical_inputs() {
    load_vector("inputs/field_values.dat", field_numerical);
    load_vector("inputs/potential.dat", potential_numerical);
    load_vector("inputs/potential_derivative.dat", potential_derivative_numerical);

    const size_t table_size = field_numerical.size();
    if (table_size < 2) {
        std::cerr << "Numerical potential tables must contain at least two entries.\n";
        exit(EXIT_FAILURE);
    }
    if (potential_numerical.size() != table_size
        || potential_derivative_numerical.size() != table_size) {
        std::cerr << "Numerical potential tables must contain the same number of entries.\n";
        exit(EXIT_FAILURE);
    }
    for (size_t i = 1; i < table_size; ++i) {
        if (!(field_numerical[i] < field_numerical[i - 1])) {
            std::cerr << "inputs/field_values.dat must be strictly descending; entries "
                      << i << " and " << i + 1 << " are out of order or repeated.\n";
            exit(EXIT_FAILURE);
        }
    }
}
#endif

// Ensure the main output path exists and is a directory.
bool ensure_results_directory() {
    struct stat info;
    if (stat("results", &info) != 0) {
        if (mkdir("results", 0777) != 0) {
            std::cerr << "Failed to create 'results' directory.\n";
            return false;
        }
    } else if (!(info.st_mode & S_IFDIR)) {
        std::cerr << "'results' exists but is not a directory.\n";
        return false;
    }
    return true;
}

// Ensure the optional post-inflation output path exists and is a directory.
bool ensure_post_inflation_directory() {
    const char* path = "results/post_inflation";
    struct stat info;
    if (stat(path, &info) != 0) {
        if (mkdir(path, 0777) != 0) {
            std::cerr << "Failed to create '" << path << "' directory.\n";
            return false;
        }
    } else if (!(info.st_mode & S_IFDIR)) {
        std::cerr << "'" << path << "' exists but is not a directory.\n";
        return false;
    }
    return true;
}

// Load a nonempty whitespace-delimited column of finite floating-point values.
void load_vector(const std::string& filename, std::vector<double>& vec) {
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Could not open file: " << filename << std::endl;
        exit(1);
    }

    vec.clear();
    double value = 0.0;
    while (infile >> value) {
        if (!std::isfinite(value)) {
            std::cerr << "Non-finite value at entry " << vec.size() + 1
                      << " in " << filename << ".\n";
            exit(EXIT_FAILURE);
        }
        vec.push_back(value);
    }

    if (!infile.eof()) {
        std::cerr << "Malformed numerical value near entry " << vec.size() + 1
                  << " in " << filename << ".\n";
        exit(EXIT_FAILURE);
    }

    if (vec.empty()) {
        std::cerr << "File appears to be empty or malformed: " << filename << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Top-level stage driver. Compile-time switches determine which optional stages
// are present; params.txt supplies run-time settings before any allocation.
int main() {
    load_runtime_parameters("params.txt");

#if numerical_potential
    load_numerical_inputs();
#endif

    if (seed < 1) {
        printf("Warning: The parameter seed has been set to %d, which will result in incorrect output.\n", seed);
    }

    if (!ensure_results_directory()) {
        return 1;
    }
#if post_inflation
    if (!ensure_post_inflation_directory()) {
        return 1;
    }
#endif

    FILE* output_ = fopen("results/output.txt", "w");
    if (!output_) {
        std::cerr << "Failed to open output file.\n";
        return 1;
    }

    initialize_simulation();     // Initialize the simulation
    run_evolution_loop(output_); // Run the main evolution loop

#if perform_deltaN
    run_deltaN_loop(output_);    // Run the deltaN evolution loop
#endif

#if post_inflation
    run_post_inflation_loop(output_);    // Run the post inflation evolution loop
#endif

    output_parameters();
    fprintf(output_, "InflationEasy program finished\n");
    printf("InflationEasy program finished\n");

    fclose(output_);
    return 0;
}
