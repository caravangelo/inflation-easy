# InflationEasy Developer Guide

This guide describes the relationships that must remain consistent among the C++ modules. It is intended for users who want to implement a model or modify the algorithm, rather than only change run parameters.

## Execution stages

The executable reads `params.txt` once at startup and then runs the enabled stages in order:

1. initialize the inflationary background and lattice fields;
2. evolve the nonlinear inflationary lattice;
3. optionally evolve each site with the separate-universe deltaN module;
4. optionally evolve the post-inflationary scalar and tensor systems;
5. write plain-text outputs for independent post-processing.

Compile-time switches in `src/parameters.h` determine which stages and data structures are present. Run-time settings in `params.txt` select physical values, integrators, stopping times, and outputs.

## Module responsibilities

- `src/main.cpp` owns startup, input loading, output-directory creation, and stage ordering.
- `src/parameters.h` declares compile-time switches and the run-time configuration interface.
- `src/runtime_parameters.cpp` defines defaults, parses `params.txt`, and validates derived controls.
- `src/initialize.cpp` constructs Fourier-space vacuum modes, enforces Hermitian symmetry, transforms them to real space, and prepares optional stages.
- `src/potential.cpp` is the only potential interface used by the evolution code. It supports either analytic functions or descending tabulated data.
- `src/spatial_discretization.h` defines the second-, fourth-, and sixth-order centered spatial operators, their Fourier eigenvalues, momentum-shell convention, and stencil-dependent stability bound.
- `src/output.cpp` computes diagnostics and writes the stable text-output schema.
- `src/ffteasy.hpp` provides the inherited in-place FFT implementation and separate Nyquist-plane convention.

The time-evolution implementation is grouped under `src/evolution/`:

- `src/evolution/evolution.cpp` owns the inflationary control loop, including integrator selection and output cadence.
- `src/evolution/inflation.cpp` contains the inflationary right-hand sides, leapfrog updates, and energy diagnostics.
- `src/evolution/integrators.cpp` contains the leapfrog drift and the RK4/RK45 machinery shared by the inflationary and post-inflationary phases.
- `src/evolution/deltaN.cpp` contains the separate-universe equations, stopping logic, and $\delta N$ control loop.
- `src/evolution/post_inflation.cpp` contains the post-inflationary scalar and tensor equations and control loop.
- `src/evolution/evolution_internal.h` documents the units, variables and time-level conventions shared by the evolution code, and declares the derivative wrappers and integration interfaces used only within this subsystem.
- `src/evolution/linear_metric.cpp` and `src/evolution/linear_metric.h` implement the optional linear scalar metric correction from spatially averaged lattice quantities.

## Reading the evolution code

Start with the comment block at the top of `src/evolution/evolution_internal.h`. It lists the three stages and their drivers, the code units (`x~ = B x`, `d tau~ = B a^s d tau`, `V~ = V/B^2`), the meaning of the global state in each stage, and the leapfrog/RK time-level conventions. Each file under `src/evolution/` then opens with the equations it implements, written in the notation of the paper.

A run proceeds as follows:

1. `initialize_simulation()` sets up the background and the vacuum fluctuations.
2. `run_inflation_loop()` evolves the lattice to the target scale factor `af`, stopping at the first step beyond it. It calls the kernels in `inflation.cpp` through the leapfrog or the RK steppers in `integrators.cpp`; `save()` writes periodic outputs, and the loop synchronizes the leapfrog state before `save_last()` writes the end-of-inflation outputs.
3. `run_deltaN_loop()` (optional) converts each site to e-fold time in `initializeN()`, integrates every site to `phiref`, and `saveN()` subtracts the mean local expansion to form zeta.
4. `run_post_inflation_loop()` (optional) maps zeta to the Newtonian potential in `initialize_post_inflation()` and evolves the potential and the tensors.

The table below links the equations of the manuscript to their implementation. Equation numbers refer to the resubmitted manuscript of 19 October 2026 (Overleaf commit `ef83d27`) and may change in later versions; the source comments do not depend on them.

| Manuscript equation | Implementation |
| --- | --- |
| Code units, Eq. (13); equations of motion in code units, Eq. (14) | `inflation.cpp`: `apply_inflation_leapfrog_kick()`, `compute_inflation_rhs()` |
| Linear metric mass shift, Eq. (6), added as in Eq. (1) | `linear_metric.cpp`: `compute_linear_metric_correction()` |
| Effective momentum, Eq. (12), and its higher-order versions | `spatial_discretization.h`: `effective_component_squared()` |
| Constant-step rescaling `dt a^(s-1)` (Sec. 3.3) | `evolution.cpp`: `inflation_code_time_step()` |
| Local e-fold equation, Eq. (21); zeta from N, Eq. (20) | `deltaN.cpp`: `compute_deltaN_rhs()`, `apply_deltaN_leapfrog_kick()`; `output.cpp`: `saveN()` |
| Tensor equation with inflationary source, Eqs. (31), (39) | `inflation.cpp`: tensor part of both kernels |
| Bardeen equation, Eq. (32); initial data, Eq. (33) | `post_inflation.cpp`; `initialize.cpp`: `initialize_post_inflation()` |
| Post-inflationary source, Eq. (34) | `post_inflation.cpp` (the code accumulates `-S_ij` and multiplies by `-2`) |
| TT projector, Eq. (41); spectra normalization, Eq. (35) | `output.cpp`: `write_gw_spectrum_impl()` |

## Numerical consistency

### Spatial operators and Fourier modes

`SPATIAL_STENCIL_ORDER` in `parameters.h` selects a centered stencil of order 2, 4, or 6 at compile time. `spatial_discretization.h` defines the corresponding real-space Laplacian and directional derivatives, together with the Laplacian-derived effective momentum. `initialize.cpp` uses the selected Laplacian eigenvalue for the vacuum frequencies and checks its Courant bound; `output.cpp` uses the same effective momentum for spectral labels, shell averages, and the tensor projector. Keep these uses synchronized when changing a stencil.

GW spectra are written over the same output shells as the scalar spectra. A self-conjugate Nyquist component has no unambiguous sign for the real projector momentum, so projected values for modes containing such a component are convention-dependent UV diagnostics rather than physical predictions. Near the lattice UV boundary, the first-derivative and Laplacian symbols differ and finite differences do not obey the continuum product rule exactly, so cancellations in the post-inflationary source require explicit resolution and stencil-order convergence checks.

The order-2 kernels in the evolution modules intentionally retain the original hot path. This preserves default-run performance, while orders 4 and 6 delegate to the shared spatial-discretization helper. The tensor sector requires both first derivatives and same-axis or mixed second derivatives, so a spatial-order change must cover all three operator families rather than only the scalar Laplacian.

Real-to-complex transforms operate in place and store the final Nyquist plane in a separate buffer. Any routine that transforms a live field for diagnostics must restore the real-space field before returning.

### Time integration

Leapfrog stores fields and derivatives at half-step offsets. Output routines temporarily synchronize them and restore the staggered state afterward. RK4 and RK45 store fields and derivatives at the same time. An RK45 rejection must not modify the accepted global state.

The inflation, deltaN, and post-inflation stages have independent integrator selections. The first and third stages share the field RK machinery in `src/evolution/integrators.cpp`; the $\delta N$ stage uses its own state layout and reuses the same RK coefficients and acceptance policy. A change to shared RK staging must therefore be tested in every compiled stage.

### Potentials and rescaling

For an analytic model, update both `analytic_potential()` and `analytic_potential_derivative()` in `potential.cpp`. For a numerical model, provide consistent `field_values.dat`, `potential.dat`, and `potential_derivative.dat` tables. The executable requires equal lengths of at least two finite entries, rejects malformed tokens, and requires the field grid to be strictly descending. The grid must cover the complete simulated trajectory.

Potential values and derivatives returned by the interface are in code units. A model change must be consistent with `V0`, `rescale_B`, the homogeneous initial field and velocity, and the chosen box and time step. Verify that the numerical trajectory never leaves a tabulated potential domain.

The optional `initial_mass_squared` run-time value contributes directly to the initial frequency as $\omega_k^2=k_{\rm eff}^2+m_{\rm init}^2$. Absence of this key is distinct from an explicitly supplied value: the default initialization contains no mass-squared contribution. Negative values are accepted; modes with non-positive total frequency squared retain the existing tachyonic-mode warning and massless fallback.

### deltaN stopping surface

The deltaN loop treats lattice sites as independent homogeneous patches. The compile-time monotonicity switches select the stopping test; `phiref` defines the final constant-field surface. `dN` is signed: positive values integrate forward and negative values integrate backward, while `Nend` is always a nonnegative integration magnitude. Backward integration requires a manually supplied `phiref`, and its stopping comparison is reversed consistently for every integrator. A new stopping prescription must be implemented consistently in both the leapfrog and RK right-hand sides.

### Tensor storage and projection

Tensor fields use six packed symmetric components in the order `[xx, yy, zz, xy, xz, yz]`. They are evolved without a transverse-traceless projection in real space; the projection is applied when spectra are constructed. Changes to component packing must update evolution, FFT buffers, and output contractions together.

At the beginning of the post-inflationary evolution, the code resets the scale factor to $a=1$ and sets $\partial_{\tilde\tau}a=f_{\rm hor}\,2\pi N/L$, where $f_{\rm hor}$ is the run-time parameter `horizon_factor`. This parameter controls the initial Hubble scale relative to the lattice resolution. The default places the full resolved range outside the Hubble scale, but the most ultraviolet modes are not necessarily deeply super-Hubble. Applications that require the complete range to satisfy $k/(aH)\ll1$ should increase `horizon_factor` and test sensitivity to this choice.

### Optional linear metric correction

`linear_metric_perturbations` is a run-time switch. When it is off, the main lattice loop performs no metric average or metric arithmetic. When it is on, the correction is computed from the current inflationary state and added only to the scalar acceleration. It does not alter the deltaN or post-inflationary equations.

### Outputs and notebooks

The C++ executable is the numerical program. It writes documented plain-text files under `results/`; the notebooks do not participate in the evolution. `notebooks/quickstart.ipynb` optionally automates configuration, compilation, and a small example, while `notebooks/plot.ipynb` reads existing outputs for visualization.

Downstream scripts rely on filenames, column order, normalization, time records, and blank-line block separators. Treat changes to this schema as interface changes and update the manuscript, README, and notebooks together.

Histogram outputs contain normalized probabilities per bin, not probability densities. Their companion `histogramtimes*.dat` files record the bin minimum and width needed to construct bin centers and densities.

When `output_bispectrum=1`, `bispectra.dat` contains only the equilateral scalar-field bispectrum at the final inflationary state; it is not a general three-momentum bispectrum output. Its filename and column layout are part of the existing output schema.

## Validation checklist

For any model or algorithm change:

1. rebuild from a clean tree with the intended `parameters.h` switches;
2. record the complete `params.txt`, compiler, thread count, and commit hash;
3. check time-step convergence and the Friedmann-constraint output;
4. vary lattice size and box size while retaining the physical scales of interest;
5. inspect the resolved momentum interval throughout the active dynamics;
6. run all affected integrators and optional stages;
7. run `make test` for the spatial-operator checks and clean-tree smoke matrix;
8. before a release, run `make test-sanitizers` and `make test-release` for the broader integrator, stencil, and optional-feature coverage.

For a source-only refactor, also run `tests/refactor_equivalence.py` against the
pre-refactor revision and require byte-for-byte agreement of deterministic outputs.

New physics modules should include a focused validation case and a brief statement of their regime of applicability.
