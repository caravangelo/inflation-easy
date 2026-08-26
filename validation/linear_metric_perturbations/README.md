# Linear metric correction validation

This directory records a paired validation of the optional linear scalar metric correction included in InflationEasy v1.1.0.

## Run configuration

- Lattice: `128^3`
- Random seed: `8`
- Time integration: default leapfrog settings
- Physical profiles: unmodified `params.numerical.txt` and `params.analytic.txt`
- Pairing: the off/on runs differ only in `linear_metric_perturbations`
- Scope: `perform_deltaN=1`; the unrelated `calculate_SIGW` and `post_inflation` modules were compiled out
- Output-only changes: histograms, energy files, and 2D snapshots were disabled

The numerical off/on runs took 3:35 and 3:38, respectively. The analytic off/on runs each took 1:24.

## Result

For the analytic slow-roll example, both estimators agree between the off/on runs at the few-times-`10^-5` relative level or better across the retained spectrum.

For the numerical USR example, the broad peak is also stable: its relative change is `3.14e-5` for linear zeta and `3.00e-5` for deltaN zeta. After excluding the first six nonzero infrared bins and the final ten ultraviolet bins, the maximum relative change is below `4.7e-4` for both estimators.

The infrared edge is not unchanged. The first nonzero bin changes by about `-16.3%`, followed by `-4.45%` and `-1.68%` in the next two bins. These modes start only marginally inside the horizon in the default setup, so this feature should be checked against an earlier start/larger box before interpreting it physically.

The CSV files contain the full off/on spectra and signed residuals. Regenerate the tables and plots with:

```bash
python3 analyze_spectra.py --runs-root /path/to/paired/runs
```
