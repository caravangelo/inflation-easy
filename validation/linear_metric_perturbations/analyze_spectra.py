#!/usr/bin/env python3
"""Compare final zeta spectra with the linear metric correction off and on."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


@dataclass
class SpectrumPair:
    k: np.ndarray
    counts: np.ndarray
    off: np.ndarray
    on: np.ndarray

    @property
    def residual(self) -> np.ndarray:
        out = np.full_like(self.off, np.nan)
        np.divide(self.on, self.off, out=out, where=self.off != 0.0)
        return out - 1.0


def load_final_spectra(run_dir: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    results = run_dir / "results"
    k = np.atleast_1d(np.loadtxt(results / "modes.dat"))
    inflaton = np.atleast_2d(np.loadtxt(results / "spectra.dat"))
    delta_n = np.atleast_2d(np.loadtxt(results / "spectraN.dat"))
    sf = np.atleast_2d(np.loadtxt(results / "sf.dat"))
    velocity = np.atleast_2d(np.loadtxt(results / "velocity.dat"))

    n_bins = len(k)
    final_inflaton = inflaton[-n_bins:]
    phase_space = k**3 / (2.0 * np.pi**2)
    linear_conversion = (sf[-1, 2] / velocity[-1, 2]) ** 2

    linear_zeta = final_inflaton[:, 2] * phase_space * linear_conversion
    delta_n_zeta = delta_n[:, 2] * phase_space
    return k, final_inflaton[:, 1], linear_zeta, delta_n_zeta


def load_pair(runs_root: Path, example: str, spectrum_index: int) -> SpectrumPair:
    off = load_final_spectra(runs_root / f"{example}_off")
    on = load_final_spectra(runs_root / f"{example}_on")
    if not np.array_equal(off[0], on[0]):
        raise ValueError(f"Mode grids differ for {example}")
    return SpectrumPair(k=off[0], counts=off[1], off=off[spectrum_index], on=on[spectrum_index])


def valid_mask(pair: SpectrumPair) -> np.ndarray:
    return (pair.k > 0.0) & (pair.counts > 0.0) & np.isfinite(pair.residual)


def metrics(pair: SpectrumPair) -> dict[str, float]:
    mask = valid_mask(pair)
    no_uv = mask.copy()
    no_uv[-10:] = False
    central = no_uv.copy()
    central[:7] = False
    peak_index = np.flatnonzero(no_uv)[np.argmax(pair.off[no_uv])]
    residual = pair.residual
    return {
        "peak_k": float(pair.k[peak_index]),
        "peak_off": float(pair.off[peak_index]),
        "peak_on": float(pair.on[peak_index]),
        "peak_relative_residual": float(residual[peak_index]),
        "median_abs_relative_residual_without_last_10_bins": float(np.median(np.abs(residual[no_uv]))),
        "max_abs_relative_residual_central_band": float(np.max(np.abs(residual[central]))),
        "first_nonzero_bin_relative_residual": float(residual[np.flatnonzero(mask)[0]]),
    }


def write_csv(path: Path, pair: SpectrumPair) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["k", "mode_count", "P_zeta_off", "P_zeta_on", "relative_residual_on_over_off_minus_1"])
        for row in zip(pair.k, pair.counts, pair.off, pair.on, pair.residual):
            writer.writerow([f"{float(value):.17e}" for value in row])


def write_figure(path: Path, example: str, linear: SpectrumPair, delta_n: SpectrumPair) -> None:
    fig, (ax_spectrum, ax_residual) = plt.subplots(
        2, 1, figsize=(7.2, 6.6), sharex=True, gridspec_kw={"height_ratios": [2.0, 1.0]}
    )
    mask = valid_mask(linear)

    ax_spectrum.loglog(linear.k[mask], linear.off[mask], color="C0", label=r"linear $\zeta$, correction off")
    ax_spectrum.loglog(linear.k[mask], linear.on[mask], "--", color="C0", label=r"linear $\zeta$, correction on")
    ax_spectrum.loglog(delta_n.k[mask], delta_n.off[mask], color="C2", label=r"$\delta N$ $\zeta$, correction off")
    ax_spectrum.loglog(delta_n.k[mask], delta_n.on[mask], "--", color="C2", label=r"$\delta N$ $\zeta$, correction on")
    ax_spectrum.set_ylabel(r"$\mathcal{P}_{\zeta}(k)$")
    ax_spectrum.set_title(f"{example.capitalize()} default example ($128^3$, seed 8)")
    ax_spectrum.grid(True, which="both", alpha=0.2)
    ax_spectrum.legend(fontsize=8, ncol=2)

    ax_residual.plot(linear.k[mask], 100.0 * linear.residual[mask], color="C0", label=r"linear $\zeta$")
    ax_residual.plot(delta_n.k[mask], 100.0 * delta_n.residual[mask], color="C2", label=r"$\delta N$ $\zeta$")
    ax_residual.axhline(0.0, color="black", lw=0.8)
    ax_residual.set_xscale("log")
    ax_residual.set_yscale("symlog", linthresh=1e-3)
    ax_residual.set_xlabel(r"$k$ [$M_{\mathrm{Pl}}$]")
    ax_residual.set_ylabel("residual [%]")
    ax_residual.grid(True, which="both", alpha=0.2)
    ax_residual.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rows: list[tuple[str, str, dict[str, float]]] = []
    for example in ("numerical", "analytic"):
        linear = load_pair(args.runs_root, example, 2)
        delta_n = load_pair(args.runs_root, example, 3)
        write_csv(args.output_dir / f"{example}_linear_zeta.csv", linear)
        write_csv(args.output_dir / f"{example}_deltaN_zeta.csv", delta_n)
        write_figure(args.output_dir / f"{example}_spectra_residuals.png", example, linear, delta_n)
        rows.append((example, "linear zeta", metrics(linear)))
        rows.append((example, "deltaN zeta", metrics(delta_n)))

    with (args.output_dir / "summary.md").open("w", encoding="utf-8") as stream:
        stream.write("# Linear metric perturbation validation\n\n")
        stream.write("Paired default-example runs at $128^3$ with seed 8. The GW and post-inflation modules were compiled out because this validation concerns only the scalar and deltaN spectra.\n\n")
        stream.write("| Example | Spectrum | Peak residual | Median | Central-band maximum | First nonzero bin |\n")
        stream.write("|---|---|---:|---:|---:|---:|\n")
        for example, spectrum, values in rows:
            stream.write(
                f"| {example} | {spectrum} | {values['peak_relative_residual']:.6e} | "
                f"{values['median_abs_relative_residual_without_last_10_bins']:.6e} | "
                f"{values['max_abs_relative_residual_central_band']:.6e} | "
                f"{values['first_nonzero_bin_relative_residual']:.6e} |\n"
            )
        stream.write("\nResiduals are `(P_on/P_off) - 1`. The median excludes the final 10 UV bins; the central-band maximum additionally excludes the first six nonzero IR bins.\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
