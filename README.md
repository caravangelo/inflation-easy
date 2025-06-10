
# InflationEasy

**InflationEasy** is a C++ program designed to simulate single-field inflation on a 3D lattice in an expanding FLRW universe. It is inspired by and partially adapted from [LATTICEEASY](http://www.felderbooks.com/latticeeasy/) by Gary Felder and Igor Tkachev.

More information is available in the associated publication: [arXiv:?????](https://arxiv.org/abs/????).

## Key Features

- **Lattice-Based Simulation:** Evolves a scalar field on a discrete lattice in an expanding universe.
- **Flexible Potential Handling:** Supports both analytical and numerical representations of the inflationary potential.  
  **Note:** The code typically runs faster when using an analytical potential.
- **OpenMP Parallelization:** Optionally leverages OpenMP for accelerated computation on multi-core systems.
- **Comprehensive Output:** Produces detailed outputs including field statistics, background quantities, and power spectra.

## Code Structure

### Source files (`src/`)
- `main.cpp`: Entry point of the program; orchestrates the simulation workflow.
- `initialize.cpp`: Sets initial field conditions and computes the initial Hubble parameter.
- `evolution.cpp`: Implements the core algorithm for time evolution of the scalar field.
- `output.cpp`: Handles writing results to disk, including observables and diagnostics.
- `potential.cpp`: Defines the inflationary potential, either analytically or via input files.
- `parameters.h`: Central header file for configuring physical and numerical parameters.  
  **Important:** Edit this file to configure your simulation setup.

### Input files (`inputs/`)
These files are only required when using a **numerical potential** (`numerical_potential = 1` in `parameters.h`):

- `field_values.dat`: Field values at which the potential is defined.
- `potential.dat`: Corresponding potential values.
- `potential_derivative.dat`: First derivative of the potential.

Each file should contain a single column of values, one per line. Analytical potentials do not require any input files.

### Output files (`results/`)
- Simulation results, logs, spectra, and other diagnostics are written here.

### Notebook (`notebooks/`)
- `plot.ipynb`: A Jupyter notebook for post-processing and visualizing simulation outputs.

## Prerequisites

- A C++17-compliant compiler (e.g., GCC or Clang).
- (Optional) OpenMP for parallel execution.
- (Optional) Python 3 and Jupyter for running the notebook.

## Building the Code

To compile the program, simply run:

```bash
make
```

## Running the Simulation

### Input Setup

Two default potentials are supported:

1. **Numerical potential (default):**  
   Implements the ultra-slow-roll (USR) potential described as Case I in [arXiv:2410.23942](https://arxiv.org/abs/2410.23942), with $\mathcal{P}_{\zeta,	ext{tree}}^{	ext{max}} = 10^{-2}$.

   **Note:** This option is slower than the analytical one. For faster runs (e.g., on a laptop), prefer an analytical potential.

2. **Analytical potential:**  
   The quadratic potential $V(\phi) = \frac{1}{2}m^2\phi^2$.  
   In this case, the (optional) $\delta N$ calculation for nonlinear $\zeta$ (see arXiv:????.?????) is typically unnecessary due to the small amplitude of perturbations.

Switch between these via the `numerical_potential` flag in `parameters.h`.

### Custom Potentials

To define a custom potential:

- For an analytical potential, modify the relevant functions in `potential.cpp`.
- For a numerical potential, place `field_values.dat`, `potential.dat`, and `potential_derivative.dat` in the `inputs/` directory. These must be one-value-per-line.
- Adjust physical and numerical parameters in `parameters.h`.

### Running the Code

After compilation, run the simulation via:

```bash
./inflation_easy
```

Output will appear in the `results/` directory. A runtime log is saved in `output.txt`, along with energy densities, field values, spectra, and more.

All quantities are given in **reduced Planck units**, where $M_{\mathrm{Pl}}^{\text{red}} = \frac{M_{\mathrm{Pl}}}{\sqrt{8\pi}} = 1$. This sets $\hbar = c = 1$ and $8\pi G = 1$, simplifying the equations.

## Jupyter Notebook

An example notebook, `plot.ipynb`, is included to help visualize output data.

To use it:

1. Ensure Python 3 and Jupyter are installed.
2. Run:

```bash
jupyter notebook plot.ipynb
```

The notebook reads data from `results/` and plots quantities like the field evolution and power spectra.

## Citing This Work

If you use *InflationEasy* in your research, please cite the associated code paper:  
[arXiv:?????](https://arxiv.org/abs/????)

Please cite also these additional references where *InflationEasy* was developed and applied:

- [arXiv:2102.06378](https://arxiv.org/abs/2102.06378)
- [arXiv:2209.13616](https://arxiv.org/abs/2209.13616)
- [arXiv:2403.12811](https://arxiv.org/abs/2403.12811)
- [arXiv:2410.23942](https://arxiv.org/abs/2410.23942)
- [arXiv:????](https://arxiv.org/abs/????)

## License

This project is released under the MIT License. Portions of the code are adapted from LATTICEEASY.
