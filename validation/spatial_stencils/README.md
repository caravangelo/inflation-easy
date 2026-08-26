# Higher-order spatial-stencil validation

The selectable spatial operators were validated against Appendix C of
[Caravano et al., arXiv:2102.06378](https://arxiv.org/abs/2102.06378), excluding
the isotropic stencils as intended for InflationEasy.

## Figure 9: modified dispersion relation

For `N=128` and `L=1.4/m`, the high-momentum endpoints obtained with the same
rounded shells and real-FFT multiplicities used by the code are:

| stencil order | final `k_eff/m` |
| ---: | ---: |
| 2 | 316.646340 |
| 4 | 365.611078 |
| 6 | 389.204649 |

The second- and fourth-order curves reproduce Fig. 9; the sixth-order curve is
the corresponding extension. `make test-spatial` checks these values together
with the real-space plane-wave eigenvalues of every derivative family.

## Figure 10: initialized power spectrum

The full `128^3` quadratic-potential test used
`m=0.51e-5`, `phi_in=14.5`, `phi'_in=-0.8152 m`, `L=1.4/m`, and stopped at
`a=1000`, following Appendix C. Binning against the naive lattice momentum
separates the spectra for different stencils, whereas using each stencil's
`k_eff` collapses them as in Fig. 10. Over the central resolved interval, the
corrected median spectra and relative RMS variations were:

| stencil order | median spectrum | relative RMS |
| ---: | ---: | ---: |
| 2 | `1.94782254e-11` | 1.537% |
| 4 | `1.93492585e-11` | 1.456% |
| 6 | `1.93445195e-11` | 1.499% |

The fourth- and sixth-order medians agree to approximately 0.025%.

## Default-path performance

The order-2 evolution wrappers retain the original kernels. Interleaved runs of
the original revision and the new default on the same `N=64` configuration
differed by 0.35% with eight OpenMP threads; single-thread orderings differed by
-1.03% and +0.10%. These variations are consistent with timing noise rather
than an added default-path cost.
