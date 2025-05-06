# InflationEasy

**InflationEasy** is a C++ program designed to simulate single-field inflation on a lattice. It is inspired by and partially adapted from [LATTICEEASY](http://www.felderbooks.com/latticeeasy/) by Gary Felder and Igor Tkachev.

## Key Features

- **Lattice-Based Simulation:** Evolves a scalar field on a discrete lattice in an expanding universe.
- **Flexible Potential Handling:** Supports both analytical and numerical representations of the inflationary potential.  
  **Note:** The code runs significantly faster when using an analytical potential.
- **OpenMP Parallelization:** Optionally leverages OpenMP for accelerated computation on multi-core systems.
- **Comprehensive Output:** Generates detailed output including field statistics, background quantities, and power spectra.

## Code Structure

- `main.cpp`: Program entry point; contains the `main` function and global variable declarations.
- `initialize.cpp`: Sets initial field conditions and computes the initial Hubble parameter.
- `evolution.cpp`: Contains the core algorithms for time evolution of the scalar field.
- `output.cpp`: Handles writing simulation results and statistics to disk.
- `potential.cpp`: Defines the inflationary potential, either analytically or from input files.
- `parameters.h`: Central header for simulation and model parameters.  
  **Important:** Users should carefully edit this file to configure simulations.
- `inputs/`: Directory for input files (required for numerical potentials):
  - `field_values.dat`: Field values.
  - `potential.dat`: Corresponding potential values.
  - `potential_derivative.dat`: Derivatives of the potential.
- `results/`: Output directory for simulation data and logs.
- `plot.ipynb`: A Jupyter notebook that provides a basic example of how to process and visualize the output data.

## Prerequisites

- A C++17-compliant compiler (e.g., GCC or Clang).
- (Optional) OpenMP library for parallel execution.
- (Optional) Python 3 and Jupyter for running the example notebook `plot.ipynb`.

## Building the Code

To compile the program, run:

```bash
make
```

The `Makefile` automatically enables OpenMP if the `parallel_calculation` flag is set to `1` in `parameters.h`. Ensure OpenMP is available on your system if you enable this option.

## Running the Simulation

### Input Setup

Two potential options are currently supported:

1. **Numerical potential (default):**  
   The default setup implements the USR potential described as Case I in [arXiv:2410.23942](https://arxiv.org/abs/2410.23942), with $\mathcal{P}_{\zeta,\text{tree}}^{\text{max}} = 10^{-2}.$

   **Note:** The code is slower when using a numerical potential instead of an analytical potential. If you need a faster code (for example if you are running it on a laptop), use an analytical potential instead.

2. **Analytical potential:**  
   A simple quadratic potential $V(\phi) = \frac{1}{2}m^2\phi^2$.  
   In this case, the (optional) $\delta N$ calculation to obtain the fully nonlinear $\zeta$ (see arXiv:????.?????) is slower and typically unnecessary due to the small amplitude of perturbations.

To switch between the two, modify the `numerical_potential` flag in `parameters.h`.

### Custom Potentials

To define your own potential:

- If using an analytical potential, modify the functions in `potential.cpp`.
- If using a numerical potential, ensure the files `field_values.dat`, `potential.dat`, and `potential_derivative.dat` are placed in the `inputs/` directory. These files must contain one value per line.
- Adjust relevant parameters in `parameters.h`.

### Executing the Program

After compilation, run the program from the terminal:

```bash
./inflation_easy
```

Simulation output will be written to the `results/` directory. A log of progress is saved to `output.txt`, and additional files contain energy densities, field values, power spectra and other useful statistics.

All quantities in the simulation are expressed in **reduced Planck units**, where the reduced Planck mass $M_{\mathrm{Pl}}^{\text{red}} = \frac{M_{\mathrm Pl}}{\sqrt{8\pi}} =1$. This means that fundamental constants are set such that $\hbar = c = 1$ and $8\pi G = 1$, simplifying the equations.


## Processing example with Jupyter

An example notebook `plot.ipynb` is provided to demonstrate how to read and visualize the simulation output using Python. This is optional, but can be useful for quick data inspection.

To run the notebook:

1. Make sure you have Python 3 and Jupyter installed.
2. Launch the notebook:

```bash
jupyter notebook plot.ipynb
```

The notebook loads output files from the `results/` directory and generates basic plots such as background evolution and field statistics.

## Citing This Work

If you use InflationEasy in your research, please cite the following papers, where this code was developed and tested:

- [arXiv:2102.06378](https://arxiv.org/abs/2102.06378)
- [arXiv:2209.13616](https://arxiv.org/abs/2209.13616)
- [arXiv:2403.12811](https://arxiv.org/abs/2403.12811)
- [arXiv:2410.23942](https://arxiv.org/abs/2410.23942)

## License

This project is released under the MIT License. Portions of the code are adapted from LATTICEEASY.
