# Linear metric perturbation validation

Paired default-example runs at $128^3$ with seed 8. The GW and post-inflation modules were compiled out because this validation concerns only the scalar and deltaN spectra.

| Example | Spectrum | Peak residual | Median | Central-band maximum | First nonzero bin |
|---|---|---:|---:|---:|---:|
| numerical | linear zeta | 3.141416e-05 | 7.422225e-07 | 4.686272e-04 | -1.629254e-01 |
| numerical | deltaN zeta | 2.998637e-05 | 7.802947e-06 | 4.676849e-04 | -1.636952e-01 |
| analytic | linear zeta | 2.233170e-05 | 1.074782e-05 | 1.840192e-05 | 2.264264e-05 |
| analytic | deltaN zeta | 2.139441e-05 | 1.158931e-05 | 3.451729e-05 | 2.263550e-05 |

Residuals are `(P_on/P_off) - 1`. The median excludes the final 10 UV bins; the central-band maximum additionally excludes the first six nonzero IR bins.
