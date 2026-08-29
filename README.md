
# InflationEasy

![](https://github.com/user-attachments/assets/7af3e20c-ec15-4f93-8764-85e422bbe8d7)

**InflationEasy** is a lattice code specifically developed for cosmological inflation. It simulates the nonlinear dynamics of a scalar field on a three-dimensional lattice in an expanding FLRW universe using finite-difference spatial derivatives. Building on the well-known [LATTICEEASY](http://www.felderbooks.com/latticeeasy/) by Gary Felder and Igor Tkachev, the code incorporates several features tailored to inflationary applications, including a nonperturbative $\delta N$ calculation of the curvature perturbation $\zeta$, optional linear scalar metric corrections, and the calculation of scalar-induced gravitational waves generated during inflation and at subsequent horizon re-entry.

More information is available in the associated publication: [arXiv:2506.11797](https://arxiv.org/abs/2506.11797). **Note:** The associated paper is currently under review. Version 1.1.1 includes features introduced after the current arXiv version, such as selectable higher-order spatial stencils and optional linear metric corrections. In case of discrepancies with the arXiv version, this documentation takes precedence.

## Key Features

- **Lattice-Based Simulation:** Evolves a scalar field on a discrete lattice in an expanding universe.
- **Flexible Potential Handling:** Supports both analytical and numerical representations of the inflationary potential.  
  **Note:** The code typically runs faster when using an analytical potential.
- **Gravitational-Waves:** Can compute the scalar-induced gravitational-wave background sourced during inflation and after inflation at horizon re-entry.
- **Linear Metric Corrections:** Optionally includes the leading linear scalar metric correction to the inflaton equation of motion.
- **Selectable Spatial Accuracy:** Supports spatial stencils up to sixth order, with consistent discretization corrections in the initialization and outputs.
- **OpenMP Parallelization:** Optionally leverages OpenMP for accelerated computation on multi-core systems.
- **Comprehensive Output:** Produces detailed outputs including field statistics, background quantities, and power spectra.

## Quickstart Notebook (Recommended For Non-Coders)

Use `notebooks/quickstart.ipynb` for a guided, self-running workflow that:

- sets up a small analytical run (default `N=64`, editable in the notebook),
- compiles and runs the code,
- stores outputs in `pedagogical_runs/...`,
- and produces key summary plots.

This is the recommended first entry point for new users.

The notebooks are optional. The numerical evolution is performed entirely by the C++ executable, which can be configured, compiled, and run directly from the command line.

## Prerequisites

- A C++17-compliant compiler (e.g., GCC or Clang).
- (Optional) OpenMP for parallel execution.
- (Optional) Python 3 for the automated test targets and notebook post-processing; it is not required to build or run the C++ simulation itself.
- (Optional) Jupyter for running the notebooks interactively.

## Building the Code

To compile the program, simply run:

```bash
make
```

### Optional Performance Flags

The default build is portable across machines. Two optional flags can be enabled for local speedups:

```bash
# CPU-specific code generation (non-portable binary)
make ENABLE_NATIVE=1

# Link-time optimization
make ENABLE_LTO=1

# Combined
make ENABLE_NATIVE=1 ENABLE_LTO=1
```

`ENABLE_NATIVE=1` may produce faster binaries on the build machine, but the executable may not run on different CPU architectures.

## Running the Simulation

### Input Setup

Two default examples are supported:

1. **Example A (current repository default): Numerical potential**  
   Implements the ultra-slow-roll (USR) potential described as Case I in [arXiv:2410.23942](https://arxiv.org/abs/2410.23942), with $\mathcal{P}_{\zeta,\rm{tree}}^{\rm{max}} = 10^{-2}$.
   This potential is built by prescribing a phenomenological SR-USR-SR profile for $\eta(N)$ and numerically reconstructing $V(\phi)$ (reverse-engineering approach, see also [arXiv:2207.10056](https://arxiv.org/abs/2207.10056)). For this reason, the default USR case is provided as tabulated input files in `inputs/` rather than as a closed analytic formula.

   **Note:** This option is slower than the analytical one. For faster runs (e.g., on a laptop), prefer an analytical potential.

2. **Example B: Analytical hilltop potential**  
   The hilltop potential:
   
 $$V(\phi) = V_0 \left(1 - \frac{1-n_s}{2}\frac{\phi^2}{2 M_{\rm Pl}^2}\right).$$  
   
Switch between these via the compile-time `numerical_potential` flag in `parameters.h` (requires recompilation).
To keep run-time values consistent with the selected potential mode, use the matching preset:

- Example A (Numerical): `params.numerical.txt`
- Example B (Analytical hilltop): `params.analytic.txt`

### Configuration Model

InflationEasy supports two classes of parameters:

- Compile-time parameters (`src/parameters.h`): potential representation (numerical vs. analytical), optional-module switches, the number of lattice points per spatial direction (`N`), and spatial stencil order.
- Run-time parameters (`params.txt`): potential parameters, time steps, output options, and most scan parameters.

A ready-to-use `params.txt` is included at the repository root.
In the current repository state, `params.txt` matches **Example A (Numerical)**.
You can edit it directly; values there override the defaults compiled into the executable.
Invalid values that would make an integration or output operation undefined are reported before the lattice is initialized.
Two preset profiles are also included for convenience:

```bash
# Switch active profile to Example A (Numerical)
cp params.numerical.txt params.txt

# Switch active profile to Example B (Analytical hilltop)
cp params.analytic.txt params.txt
```

The spatial discretization is selected in `src/parameters.h`:

```cpp
#define SPATIAL_STENCIL_ORDER 2
```

Choose `2` (default), `4`, or `6`, then recompile. The selected order is used consistently for the scalar and tensor Laplacians, the directional derivatives entering the GW sector, vacuum initialization, effective-momentum output and binning, and the Fourier-space TT projector. Higher orders use wider stencils and a correspondingly stricter stability limit; initialization checks `dt/dx` against the selected limit.

### Essential Runtime Parameters (Quick Guide)

The most commonly adjusted run-time keys in `params.txt` are summarized below.

#### Inflationary Initialization and Evolution

- `initial_field`: homogeneous initial value of the inflaton, in code units.
- `initial_derivative`: homogeneous initial inflaton velocity, in code units.
- `dt`: inflationary time step.
- `af`: final scale factor for the inflationary evolution. If omitted, it defaults to `2*N`.
- `initial_mass_squared`: optional mass-squared contribution to the initial inflaton mode frequency, in the same internal code units as $k_{\rm eff}^2$. If omitted, the initializer uses $\omega_k^2=k_{\rm eff}^2$; if supplied, it uses $\omega_k^2=k_{\rm eff}^2+m_{\rm init}^2$, with `initial_mass_squared` providing $m_{\rm init}^2$.
- `linear_metric_perturbations`: set to `1` to add the leading linear scalar metric correction during the inflationary evolution. It defaults to `0` and, as a run-time option, does not require recompilation.
- `output_bispectrum`: set to `1` to compute the equilateral scalar-field bispectrum at the final inflationary time. The resulting `bispectra.dat` file contains only equilateral configurations and does not contain general triangle configurations. This calculation is computationally expensive for large lattices and is disabled by default.

#### $\delta N$ Stage

These parameters are used when `perform_deltaN=1`:

- `dN`, `Nend`: integration step and maximum integration magnitude. `dN` must be nonzero; a negative value integrates backward and requires `use_phiref_manual=1`, while `Nend` remains nonnegative.

Important: `monotonic_potential` / `antimonotonic_potential` select the compile-time deltaN stopping potential criterion in `src/parameters.h`; they are not `params.txt` keys. For forward integration, the implemented criteria are: monotonic -> evolve while `|phi| > |phi_ref|`, anti-monotonic -> evolve while `|phi| < |phi_ref|`, and if both are `0`, generic potential fallback -> evolve while `V(phi) > V(phi_ref)`.

If `Nend` is reached before every patch crosses the selected surface, the code records a warning in both the terminal and `results/output.txt`. The deltaN histogram then uses only completed patches, while incomplete sites are set to zero in spatial and spectral products. The post-inflationary stage requires every patch to complete and stops with an error otherwise.

#### Post-inflationary Stage

These parameters are used when `post_inflation=1`:

- `dt_post_inflation`, `af_post_inflation`: time step and final scale factor (`af_post_inflation` defaults to `2*N` if omitted).
- `omega`: constant equation-of-state parameter $w$.

#### Integrator Selection

- `inflation_integrator`, `deltaN_integrator`, `post_inflation_integrator`: choose `leapfrog`, `rk4`, or `rk45` for the corresponding stage (each defaults to `leapfrog` if omitted).
- `rk45_abs_tol`, `rk45_rel_tol`, `rk45_min_dt`, `rk45_max_dt`, `rk45_safety`: used only by stages whose integrator is set to `rk45`.

### Custom Potentials

To define a custom potential:

- For an analytical potential, modify the relevant functions in `potential.cpp`.
- For a numerical potential, place `field_values.dat`, `potential.dat`, and `potential_derivative.dat` in the `inputs/` directory. These must be finite one-value-per-line tables of equal length, with a strictly descending field grid.
- Adjust physical and numerical run parameters in `params.txt` (or in defaults inside `runtime_parameters.cpp`).

### Running the Code

After compilation, run the simulation via:

```bash
./inflation_easy
```

Output will appear in the `results/` directory. A runtime log is saved in `results/output.txt`, along with energy densities, field values, spectra, and more.

Histogram files store normalized probabilities per bin. The example notebooks use the bin metadata in the corresponding `histogramtimes*.dat` files to construct bin centers and divide by the bin width when plotting probability densities.

The documentation and outputs use natural units, $\hbar=c=1$, and dimensionful quantities are expressed in units of the reduced Planck mass $M_{\rm Pl}\equiv(8\pi G)^{-1/2}$. Thus $M_{\rm Pl}=1$ and $8\pi G=1$.

## Code Structure

### Source files (`src/`)

The top-level source files provide program orchestration and numerical
infrastructure:

- `main.cpp`: Entry point of the program; orchestrates the simulation workflow.
- `main.h`: Declares the shared simulation state and core functions.
- `initialize.cpp`: Prepares the scalar, expansion, and optional tensor initial conditions.
- `spatial_discretization.h`: Defines the selectable spatial operators and their Fourier-space effective momenta.
- `ffteasy.hpp`: Provides the FFT routines inherited from LATTICEEASY.
- `output.cpp`: Handles writing results to disk, including observables and diagnostics.
- `potential.cpp`: Defines the inflationary potential, either analytically or via input files.
- `parameters.h`: Compile-time configuration for the potential representation (numerical vs. analytical), optional modules, number of lattice points per spatial direction, and spatial stencil order.
  **Important:** Edit this file only for settings that require recompilation.
- `runtime_parameters.cpp`: Run-time defaults and parser for `params.txt` overrides.

The evolution subsystem is grouped under `src/evolution/` by responsibility:

- `evolution/evolution.cpp`: Runs the inflationary time-evolution loop and selects its integrator.
- `evolution/inflation.cpp`: Implements the inflationary equations, energy diagnostics, and leapfrog updates.
- `evolution/integrators.cpp`: Provides the shared leapfrog drift and RK4/RK45 stepping machinery.
- `evolution/deltaN.cpp`: Implements the separate-universe $\delta N$ evolution and stopping conditions.
- `evolution/post_inflation.cpp`: Implements the post-inflationary scalar and tensor evolution.
- `evolution/linear_metric.cpp` / `evolution/linear_metric.h`: Implement the optional linear scalar metric correction.
- `evolution/evolution_internal.h`: Declares interfaces shared only within the evolution subsystem.

### Input files (`inputs/`)
These files are only required when using a **numerical potential** (`numerical_potential = 1` in `parameters.h`):

- `field_values.dat`: Field values at which the potential is defined.
- `potential.dat`: Corresponding potential values.
- `potential_derivative.dat`: First derivative of the potential.

Each file should contain a single column of finite values, one per line. The three tables must have the same length, contain at least two entries, and use a strictly descending field grid. The executable reports malformed or inconsistent tables before starting the simulation. Analytical potentials do not require any input files.

### Output files (`results/`)
- Simulation results, logs, spectra, and other diagnostics are written here.

### Notebooks (`notebooks/`)
- `quickstart.ipynb`: Guided end-to-end run notebook (recommended first notebook).
- `plot.ipynb`: Post-processing and visualization notebook for simulation outputs.

### Tests (`tests/`)
- `spatial_discretization_test.cpp`: Checks each spatial stencil against its Fourier eigenvalue.
- `release_smoke.py`: Runs the clean-tree smoke, sanitizer, and pre-release configuration matrices.
- `refactor_equivalence.py`: Compares deterministic outputs from a source-only refactor against a selected Git revision.

Users extending the C++ implementation should read [`DEVELOPER_GUIDE.md`](DEVELOPER_GUIDE.md), which documents the responsibilities of these modules and the relationships that must remain consistent among potentials, spatial stencils, integrators, optional modules, and outputs.

## Reproducibility Notes

- The random seed is controlled by `seed` in `params.txt` (or by defaults in `src/runtime_parameters.cpp`).
- For reproducible results, record:
  - commit hash (`git rev-parse HEAD`)
  - compiler and version (`c++ --version`)
  - full `src/parameters.h` (compile-time)
  - full `params.txt` used for the run (run-time)
  - whether OpenMP was enabled
- Main run metadata is written by the code to `results/info.dat`.

## Developer Notes

The complete module map and modification checklist are provided in [`DEVELOPER_GUIDE.md`](DEVELOPER_GUIDE.md).

If you modify the evolution modules listed above, keep these invariants unchanged unless you intentionally redesign the algorithm:

- Periodic finite-difference stencils for all lattice derivatives.
- Leapfrog staggering semantics (half-step synchronization only at output boundaries).
- RK45 acceptance logic based on weighted RMS error with `rk45_abs_tol`/`rk45_rel_tol`.
- Existing output schema (`results/*.dat` and `results/post_inflation/*.dat`) used by analysis scripts.

### Automated Tests

Run the fast spatial-operator and simulation smoke tests with:

```bash
make test
```

Before preparing a release, run the focused memory checks and the broader
integrator, stencil, and optional-feature matrix:

```bash
make test-sanitizers
make test-release
```

For a source-only architectural change, compare every deterministic output in
the dedicated equivalence matrix against the intended baseline, for example:

```bash
python3 tests/refactor_equivalence.py --reference v1.1.0 --current-worktree
```

`make test-spatial` remains available when only the second-, fourth-, and
sixth-order spatial operators and their Fourier eigenvalues need to be checked.

## Jupyter Notebooks

Two notebooks are included:

- `quickstart.ipynb`: recommended first notebook for a guided run + summary plots.
- `plot.ipynb`: dedicated post-processing notebook for detailed visualization of outputs.

To use it:

1. Ensure Python 3 and Jupyter are installed.
2. Run:

```bash
jupyter notebook notebooks/quickstart.ipynb
```

For post-processing existing outputs, you can also run:

```bash
jupyter notebook notebooks/plot.ipynb
```

## Citing This Work

If you use *InflationEasy* in your research, please cite the associated code paper:  
[arXiv:2506.11797](https://arxiv.org/abs/2506.11797)

ASCL entry: [ascl:2603.004](https://ascl.net/2603.004)

A machine-readable citation file is provided in `CITATION.cff`.

Please cite also these additional references where *InflationEasy* was developed and applied:

- [arXiv:2102.06378](https://arxiv.org/abs/2102.06378)
- [arXiv:2209.13616](https://arxiv.org/abs/2209.13616)
- [arXiv:2403.12811](https://arxiv.org/abs/2403.12811)
- [arXiv:2410.23942](https://arxiv.org/abs/2410.23942)
- [arXiv:2506.11795](https://arxiv.org/abs/2506.11795)
- [arXiv:2604.03628](https://arxiv.org/abs/2604.03628)

## License

This project is released under the MIT License. Portions of the code are adapted from LATTICEEASY.
