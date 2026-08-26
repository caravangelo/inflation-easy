#!/usr/bin/env python3
"""Strict clean-tree smoke and release matrix for InflationEasy.

The suite copies only release inputs from the current working tree into
temporary directories.  It never builds in, or writes results to, the source
checkout.  The default ``ci`` tier is deliberately compact; ``release`` adds a
broader configuration matrix intended for a tagged release candidate.
"""

from __future__ import annotations

import argparse
import itertools
import math
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Optional, Union


BUILD_TIMEOUT = 120
RUN_TIMEOUT = 45


@dataclass(frozen=True)
class BuildConfig:
    name: str
    numerical: int = 1
    delta_n: int = 1
    gw: int = 1
    post: int = 1
    stencil: int = 2
    parallel: int = 1
    n: int = 8


@dataclass(frozen=True)
class RunConfig:
    name: str
    inflation_integrator: str = "leapfrog"
    delta_n_integrator: str = "leapfrog"
    post_integrator: str = "leapfrog"
    linear_metric: int = 0
    initial_mass: Optional[float] = None
    output_log: int = 0
    output_bispectrum: int = 0
    output_spectra: int = 1
    output_histogram: int = 1
    output_energy: int = 1
    output_box2d: int = 1
    output_box3d: int = 0
    delta_n_step: Optional[float] = None
    delta_n_end: Optional[float] = None
    phiref_manual: Optional[float] = None
    low_cutoff_index: float = 0.0
    high_cutoff_index: float = 0.0
    expect_incomplete_delta_n: bool = False
    expect_degenerate_delta_n_histogram: bool = False
    expect_backward_crossing: bool = False
    expect_no_delta_n_steps: bool = False
    homogeneous_diagnostics: bool = False
    omega: float = 1.0 / 3.0


def command(
    args: list[str],
    *,
    cwd: Path,
    timeout: int,
    env: Optional[dict[str, str]] = None,
    expect_success: bool = True,
) -> subprocess.CompletedProcess[str]:
    try:
        proc = subprocess.run(
            args,
            cwd=cwd,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as exc:
        output = exc.stdout or ""
        raise AssertionError(
            f"Command timed out after {timeout}s: {' '.join(args)}\n{output}"
        ) from exc

    if expect_success and proc.returncode != 0:
        raise AssertionError(
            f"Command failed ({proc.returncode}): {' '.join(args)}\n{proc.stdout}"
        )
    if not expect_success and proc.returncode == 0:
        raise AssertionError(
            f"Command unexpectedly succeeded: {' '.join(args)}\n{proc.stdout}"
        )
    return proc


def release_files(repo: Path) -> list[Path]:
    """Return the deterministic build/run allowlist, with no Git dependency."""
    root_files = (
        Path("Makefile"),
        Path("params.txt"),
        Path("params.numerical.txt"),
        Path("params.analytic.txt"),
    )
    build_globs = (
        "src/*.cpp",
        "src/*.h",
        "src/*.hpp",
        "inputs/*.dat",
    )

    paths = list(root_files)
    for pattern in build_globs:
        paths.extend(path.relative_to(repo) for path in sorted(repo.glob(pattern)))

    missing = [relative for relative in root_files if not (repo / relative).is_file()]
    if missing:
        raise AssertionError(f"Required release inputs are missing: {missing}")
    if not any(relative.parts[0] == "src" for relative in paths):
        raise AssertionError("No source build inputs were found")
    if not any(relative.parts[0] == "inputs" for relative in paths):
        raise AssertionError("No numerical-potential input tables were found")
    if len(paths) != len(set(paths)):
        raise AssertionError("The release input allowlist contains duplicate paths")
    return sorted(paths)


def copy_clean_tree(repo: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=False)
    for relative in release_files(repo):
        source = repo / relative
        if not source.is_file():
            raise AssertionError(f"Allow-listed release input is missing: {relative}")
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)


def replace_checked(path: Path, pattern: str, replacement: str) -> None:
    text = path.read_text(encoding="utf-8")
    updated, count = re.subn(pattern, replacement, text, flags=re.MULTILINE)
    if count != 1:
        raise AssertionError(
            f"Expected exactly one configuration substitution in {path}: {pattern!r}; got {count}"
        )
    path.write_text(updated, encoding="utf-8")


def configure_build(work: Path, cfg: BuildConfig) -> None:
    header = work / "src" / "parameters.h"
    values = {
        "numerical_potential": cfg.numerical,
        "perform_deltaN": cfg.delta_n,
        "calculate_SIGW": cfg.gw,
        "post_inflation": cfg.post,
        "parallel_calculation": cfg.parallel,
    }
    for macro, value in values.items():
        replace_checked(
            header,
            rf"^(\s*#define\s+{re.escape(macro)}\s+)[01]\s*$",
            rf"\g<1>{value}",
        )
    replace_checked(
        header,
        r"^(\s*const\s+int\s+N\s*=\s*)\d+(\s*;)\s*$",
        rf"\g<1>{cfg.n}\g<2>",
    )
    replace_checked(
        header,
        r"^(\s*#define\s+SPATIAL_STENCIL_ORDER\s+)[246]\s*$",
        rf"\g<1>{cfg.stencil}",
    )


def make_args(*, sanitizers: bool) -> list[str]:
    args = ["make", "-j2"]
    if sanitizers:
        flags = (
            "-std=c++17 -O1 -g -Wall -Wextra "
            "-fsanitize=address,undefined -fno-omit-frame-pointer"
        )
        args.extend([f"CXXFLAGS={flags}", "LDFLAGS=-fsanitize=address,undefined", "LIBS=-lm"])
        if platform.system() == "Darwin":
            args.extend(["LLVM_PREFIX=", "CXX=/usr/bin/clang++"])
    return args


def build(work: Path, *, sanitizers: bool = False, expect_success: bool = True) -> subprocess.CompletedProcess[str]:
    proc = command(
        make_args(sanitizers=sanitizers),
        cwd=work,
        timeout=BUILD_TIMEOUT,
        expect_success=expect_success,
    )
    if expect_success and re.search(r"(?mi)\bwarning:", proc.stdout):
        raise AssertionError(f"Compiler warning in strict build:\n{proc.stdout}")
    return proc


def runtime_values(build_cfg: BuildConfig, run_cfg: RunConfig) -> dict[str, str]:
    if build_cfg.numerical:
        values: dict[str, str] = {
            "V0": "3e-9",
            "initial_field": "2.9181235049318586",
            "initial_derivative": "-0.06727651095116181",
        }
        d_n = "0.0005"
        n_end = "0.01"
    else:
        values = {
            "V0": "3.338e-13",
            "ns": "0.97",
            "initial_field": "0.0935",
            "initial_derivative": "0.000796",
        }
        d_n = "0.0001"
        n_end = "0.01"

    values.update(
        {
            "seed": "8",
            "rescale_s": "0.0",
            "af": "1.003",
            "dt": "0.0005",
            "linear_metric_perturbations": str(run_cfg.linear_metric),
            "L": "10.0",
            "output_freq": "1",
            "output_infrequent_freq": "1",
            "output_spectra": str(run_cfg.output_spectra),
            "output_histogram": str(run_cfg.output_histogram),
            "output_energy": str(run_cfg.output_energy),
            "output_box3D": str(run_cfg.output_box3d),
            "output_box2D": str(run_cfg.output_box2d),
            "output_bispectrum": str(run_cfg.output_bispectrum),
            "screen_updates": "0",
            "nbins": "16",
            "high_cutoff_index": f"{run_cfg.high_cutoff_index:.17g}",
            "low_cutoff_index": f"{run_cfg.low_cutoff_index:.17g}",
            "forcing_cutoff": "0",
            "inflation_integrator": run_cfg.inflation_integrator,
            "rk45_abs_tol": "1e-7",
            "rk45_rel_tol": "1e-5",
            "rk45_min_dt": "1e-12",
            "rk45_max_dt": "0.0005",
            "rk45_safety": "0.9",
        }
    )
    if run_cfg.initial_mass is not None:
        values["initial_mass_squared"] = f"{run_cfg.initial_mass:.17g}"

    if build_cfg.delta_n:
        values.update(
            {
                "dN": (
                    f"{run_cfg.delta_n_step:.17g}"
                    if run_cfg.delta_n_step is not None
                    else d_n
                ),
                "Nend": (
                    f"{run_cfg.delta_n_end:.17g}"
                    if run_cfg.delta_n_end is not None
                    else n_end
                ),
                "use_phiref_manual": "1" if run_cfg.phiref_manual is not None else "0",
                "output_LOG": str(run_cfg.output_log),
                "eta_log": "-0.5",
                "deltaN_integrator": run_cfg.delta_n_integrator,
            }
        )
        if run_cfg.phiref_manual is not None:
            values["phiref_manual_value"] = f"{run_cfg.phiref_manual:.17g}"
    if build_cfg.post:
        values.update(
            {
                "horizon_factor": "1.0",
                "omega": f"{run_cfg.omega:.17g}",
                "dt_post_inflation": "0.0005",
                "af_post_inflation": "1.003",
                "post_inflation_integrator": run_cfg.post_integrator,
            }
        )
    return values


def write_params(path: Path, values: dict[str, str]) -> None:
    body = "# Generated by tests/release_smoke.py\n"
    body += "".join(f"{key} = {value}\n" for key, value in values.items())
    path.write_text(body, encoding="utf-8")


def expected_outputs(build_cfg: BuildConfig, run_cfg: RunConfig) -> set[str]:
    expected = {
        "info.dat",
        "output.txt",
        "means.dat",
        "variance.dat",
        "velocity.dat",
        "sf.dat",
        "modes.dat",
    }
    if run_cfg.output_energy:
        expected.update({"energy.dat", "conservation.dat"})
    if run_cfg.output_spectra:
        expected.update({"spectra.dat", "spectratimes.dat"})
        if build_cfg.gw:
            expected.update({"spectraGW.dat", "spectraGWdot.dat"})
    if run_cfg.output_histogram:
        expected.update({"histogram.dat", "histogramtimes.dat"})
    if run_cfg.output_box3d:
        expected.add("box.dat")
    if run_cfg.output_box2d:
        expected.update({"snapshots_2d_phi.dat", "snapshots_2d_phidot.dat"})
    if run_cfg.output_bispectrum:
        expected.add("bispectra.dat")

    if build_cfg.delta_n:
        # deltaN histograms are part of the final deltaN product in the current interface.
        expected.update({"histogramN.dat", "histogramtimesN.dat"})
        if run_cfg.output_spectra:
            expected.add("spectraN.dat")
        if run_cfg.output_box2d:
            expected.add("snapshots_2d_deltaN.dat")
        if run_cfg.output_log:
            expected.update({"spectraLOG.dat", "histogramLOG.dat", "histogramtimesLOG.dat"})

    if build_cfg.post:
        expected.update(
            {
                "post_inflation/means.dat",
                "post_inflation/variance.dat",
                "post_inflation/velocity.dat",
                "post_inflation/sf.dat",
            }
        )
        if run_cfg.output_spectra:
            expected.update(
                {
                    "post_inflation/spectra.dat",
                    "post_inflation/spectratimes.dat",
                    "post_inflation/spectraGW.dat",
                    "post_inflation/spectraGWdot.dat",
                }
            )
        if run_cfg.output_histogram:
            expected.update(
                {"post_inflation/histogram.dat", "post_inflation/histogramtimes.dat"}
            )
    return expected


def numeric_rows(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    for lineno, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        raw = raw.strip()
        if not raw:
            continue
        try:
            row = [float(token) for token in raw.split()]
        except ValueError as exc:
            raise AssertionError(f"Non-numeric token in {path}:{lineno}: {raw}") from exc
        if not all(math.isfinite(value) for value in row):
            raise AssertionError(f"Non-finite value in {path}:{lineno}: {raw}")
        rows.append(row)
    if not rows:
        raise AssertionError(f"Numeric output is empty: {path}")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise AssertionError(f"Inconsistent numeric column count: {path}")
    return rows


def assert_monotonic_time(path: Path) -> None:
    rows = numeric_rows(path)
    times = [row[0] for row in rows]
    if any(right < left for left, right in zip(times, times[1:])):
        raise AssertionError(f"First column is not monotonic: {path}")
    if len(times) > 1 and not any(right > left for left, right in zip(times, times[1:])):
        raise AssertionError(f"First column never advances: {path}")


def validate_histogram(results: Path, data_rel: str, times_rel: str, nbins: int) -> None:
    data = numeric_rows(results / data_rel)
    times = numeric_rows(results / times_rel)
    if len(data) != nbins * len(times):
        raise AssertionError(
            f"Histogram shape mismatch for {data_rel}: {len(data)} rows, "
            f"expected {nbins} * {len(times)}"
        )
    values = [row[0] for row in data]
    for index in range(len(times)):
        total = sum(values[index * nbins : (index + 1) * nbins])
        if not math.isclose(total, 1.0, rel_tol=2e-5, abs_tol=2e-5):
            raise AssertionError(f"Histogram block {index} in {data_rel} sums to {total}")


def validate_degenerate_delta_n_histogram(results: Path, nbins: int) -> None:
    """Check the zero-width case where every completed deltaN value is identical."""
    probabilities = [row[0] for row in numeric_rows(results / "histogramN.dat")]
    if len(probabilities) != nbins:
        raise AssertionError("Degenerate deltaN histogram must contain one complete bin block")
    populated = [value for value in probabilities if value != 0.0]
    if len(populated) != 1 or not math.isclose(populated[0], 1.0, abs_tol=2e-6):
        raise AssertionError(
            "Degenerate deltaN histogram must place unit probability in exactly one bin"
        )

    metadata = numeric_rows(results / "histogramtimesN.dat")
    if len(metadata) != 1 or len(metadata[0]) != 4:
        raise AssertionError("Degenerate deltaN histogram metadata has an unexpected shape")
    bmin, width = metadata[0][2], metadata[0][3]
    if not math.isclose(bmin, 0.0, abs_tol=1e-15) or not (width > 0.0):
        raise AssertionError(
            f"Degenerate deltaN histogram requires bmin=0 and a finite positive width; "
            f"got bmin={bmin}, width={width}"
        )

    snapshot = [row[0] for row in numeric_rows(results / "snapshots_2d_deltaN.dat")]
    if any(value != 0.0 for value in snapshot):
        raise AssertionError("Mean-subtracted degenerate deltaN output is not identically zero")

    spectrum = numeric_rows(results / "spectraN.dat")
    if any(len(row) < 3 or row[2] != 0.0 for row in spectrum):
        raise AssertionError("Degenerate deltaN field must have a zero power spectrum")


def validate_spectrum_blocks(
    results: Path,
    relative: str,
    expected_modes: list[float],
    block_count: int,
    require_full_range: bool = False,
) -> None:
    rows = numeric_rows(results / relative)
    rows_per_block = len(expected_modes)
    expected_rows = rows_per_block * block_count
    if len(rows) != expected_rows:
        raise AssertionError(
            f"Spectrum shape mismatch for {relative}: {len(rows)} rows, "
            f"expected {rows_per_block} * {block_count}"
        )

    for block in range(block_count):
        offset = block * rows_per_block
        labels = [row[0] for row in rows[offset : offset + rows_per_block]]
        for index, (actual, expected) in enumerate(zip(labels, expected_modes)):
            if not math.isclose(actual, expected, rel_tol=2e-6, abs_tol=2e-12):
                raise AssertionError(
                    f"Spectrum label mismatch in {relative}, block {block}, row {index}: "
                    f"got {actual}, expected modes.dat value {expected}"
                )
            if require_full_range and index > 0 and rows[offset + index][1] <= 0:
                raise AssertionError(
                    f"GW spectrum does not populate the full lattice range in "
                    f"{relative}, block {block}, row {index}"
                )


def validate_energy_constraint(results: Path) -> None:
    energy = numeric_rows(results / "energy.dat")
    conservation = numeric_rows(results / "conservation.dat")
    scale_factor = numeric_rows(results / "sf.dat")
    if not (len(energy) == len(conservation) == len(scale_factor)):
        raise AssertionError(
            "energy.dat, conservation.dat, and sf.dat must describe the same output times"
        )

    for index, (energy_row, conservation_row, sf_row) in enumerate(
        zip(energy, conservation, scale_factor)
    ):
        if len(energy_row) != 5 or len(conservation_row) != 3 or len(sf_row) != 4:
            raise AssertionError(f"Unexpected energy/constraint output width at row {index}")
        if not (
            math.isclose(energy_row[0], conservation_row[0], abs_tol=5e-7)
            and math.isclose(energy_row[0], sf_row[0], abs_tol=5e-7)
            and math.isclose(energy_row[1], conservation_row[1], rel_tol=2e-6)
            and math.isclose(energy_row[1], sf_row[1], rel_tol=2e-6)
        ):
            raise AssertionError(f"Energy/constraint output times do not align at row {index}")

        density = sum(energy_row[2:5])
        friedmann_density = 3.0 * sf_row[2] * sf_row[2]
        if density <= 0.0:
            raise AssertionError(f"Non-positive total energy density at row {index}: {density}")
        ratio = friedmann_density / density
        reported_ratio = conservation_row[2]
        if not math.isclose(ratio, reported_ratio, rel_tol=2e-5, abs_tol=2e-5):
            raise AssertionError(
                f"conservation.dat is inconsistent with energy.dat and sf.dat at row {index}: "
                f"computed {ratio}, reported {reported_ratio}"
            )
        if not math.isclose(reported_ratio, 1.0, rel_tol=1e-3, abs_tol=1e-3):
            raise AssertionError(
                f"Friedmann constraint differs from unity at row {index}: {reported_ratio}"
            )


def validate_acceleration_outputs(
    results: Path,
    run_cfg: RunConfig,
    lattice_size: int,
) -> None:
    energy = numeric_rows(results / "energy.dat")
    scale_factor = numeric_rows(results / "sf.dat")
    for index, (energy_row, sf_row) in enumerate(zip(energy, scale_factor)):
        a = sf_row[1]
        hubble = sf_row[2]
        gradient = energy_row[3]
        potential = energy_row[4]
        # Cosmic-time form of the scale-factor RHS written by scale().
        expected = a * (potential + 2.0 * gradient / 3.0 - 2.0 * hubble * hubble)
        if not math.isclose(sf_row[3], expected, rel_tol=3e-5, abs_tol=2e-18):
            raise AssertionError(
                f"Inflationary acceleration output disagrees with its instantaneous RHS "
                f"at row {index}: got {sf_row[3]}, expected {expected}"
            )

    velocity = numeric_rows(results / "velocity.dat")
    snapshots = numeric_rows(results / "snapshots_2d_phidot.dat")
    values_per_snapshot = lattice_size * lattice_size
    if len(snapshots) != values_per_snapshot * len(velocity):
        raise AssertionError("Homogeneous velocity snapshots do not align with velocity.dat")
    if not any(abs(row[1] - 1.0) > 1e-6 for row in velocity):
        raise AssertionError("Homogeneous velocity conversion was not tested at a != 1")
    for block, velocity_row in enumerate(velocity):
        expected_velocity = velocity_row[2]
        begin = block * values_per_snapshot
        for snapshot_row in snapshots[begin : begin + values_per_snapshot]:
            if not math.isclose(
                snapshot_row[0], expected_velocity, rel_tol=3e-6, abs_tol=2e-15
            ):
                raise AssertionError(
                    f"Physical velocity conversion mismatch in snapshot block {block}: "
                    f"got {snapshot_row[0]}, expected {expected_velocity}"
                )

    post_sf_path = results / "post_inflation/sf.dat"
    if post_sf_path.exists():
        for index, row in enumerate(numeric_rows(post_sf_path)):
            a = row[1]
            hubble = row[2]
            expected = -0.5 * (1.0 + 3.0 * run_cfg.omega) * a * hubble * hubble
            if not math.isclose(row[3], expected, rel_tol=3e-5, abs_tol=2e-18):
                raise AssertionError(
                    f"Post-inflation acceleration output disagrees with the constant-w RHS "
                    f"at row {index}: got {row[3]}, expected {expected}"
                )


def delta_n_reports(results: Path) -> list[float]:
    text = (results / "output.txt").read_text(encoding="utf-8")
    return [
        float(match.group(1))
        for match in re.finditer(r"(?m)^N\s*=\s*([-+0-9.eE]+)\s*$", text)
    ]


def validate_delta_n_reporting(
    results: Path,
    run_cfg: RunConfig,
    process_output: str,
    total_sites: int,
) -> None:
    reports = delta_n_reports(results)
    if run_cfg.expect_no_delta_n_steps:
        if reports:
            raise AssertionError("The deltaN loop advanced despite a zero Nend budget")
    elif not reports:
        raise AssertionError("The deltaN loop did not report any accepted steps")

    configured_step = run_cfg.delta_n_step
    expected_sign = -1.0 if configured_step is not None and configured_step < 0.0 else 1.0
    if reports and not any(expected_sign * value > 0.0 for value in reports):
        raise AssertionError("The deltaN progress report has the wrong integration direction")
    if any(expected_sign * value < -1e-12 for value in reports):
        raise AssertionError("The deltaN progress report changes integration direction")
    if run_cfg.delta_n_end is not None:
        budget_tolerance = 2e-6
        if any(abs(value) > run_cfg.delta_n_end + budget_tolerance for value in reports):
            raise AssertionError("The deltaN loop exceeded the configured Nend budget")

    if run_cfg.expect_backward_crossing:
        if configured_step is None or configured_step >= 0.0 or run_cfg.phiref_manual is None:
            raise AssertionError("The backward-crossing test is not configured consistently")
        box_values = [row[0] for row in numeric_rows(results / "box.dat")]
        initial_slice = box_values[-total_sites:]
        if len(initial_slice) != total_sites:
            raise AssertionError("The backward-crossing test lacks a complete final inflationary box")
        if not all(abs(value) > abs(run_cfg.phiref_manual) for value in initial_slice):
            raise AssertionError(
                "The backward-crossing test did not start with every patch active"
            )

    # In the adaptive loop, the accumulated e-fold value must advance by the
    # accepted step, which is also the increment applied to the evolved state.
    # A proposed next step can differ substantially and must not be reported as
    # already accepted.
    if reports and run_cfg.delta_n_integrator == "rk45":
        main_time = numeric_rows(results / "sf.dat")[-1][0]
        delta_n_time = numeric_rows(results / "histogramtimesN.dat")[-1][0]
        state_increment = delta_n_time - main_time
        if not math.isclose(reports[-1], state_increment, rel_tol=2e-5, abs_tol=2e-6):
            raise AssertionError(
                "RK45 deltaN bookkeeping does not match the accepted state increment: "
                f"reported N={reports[-1]}, state increment={state_increment}"
            )

    warning = re.search(
        r"Warning:\s+(\d+) of (\d+) deltaN patches did not reach the selected hypersurface",
        process_output,
    )
    if run_cfg.expect_incomplete_delta_n:
        if warning is None:
            raise AssertionError("The deliberately incomplete deltaN case did not report masking")
        incomplete = int(warning.group(1))
        reported_total = int(warning.group(2))
        if reported_total != total_sites or not (0 < incomplete < reported_total):
            raise AssertionError(
                f"Unexpected incomplete-patch report: {incomplete} of {reported_total}"
            )
        if warning.group(0) not in (results / "output.txt").read_text(encoding="utf-8"):
            raise AssertionError("The incomplete-patch warning was not persisted in output.txt")
        completed = reported_total - incomplete

        histogram = [row[0] for row in numeric_rows(results / "histogramN.dat")]
        populated_bins = 0
        for probability in histogram:
            count = probability * completed
            if not math.isclose(count, round(count), rel_tol=2e-5, abs_tol=2e-5):
                raise AssertionError(
                    "deltaN histogram probabilities are not normalized by the completed-patch count"
                )
            if probability > 0.0:
                populated_bins += 1
        if populated_bins < 2:
            raise AssertionError(
                "The incomplete deltaN case needs at least two populated bins to test normalization"
            )

        snapshot = [row[0] for row in numeric_rows(results / "snapshots_2d_deltaN.dat")]
        if not any(value == 0.0 for value in snapshot):
            raise AssertionError("Incomplete deltaN patches were not visibly masked in the output slice")
        if not any(value != 0.0 for value in snapshot):
            raise AssertionError("The incomplete deltaN output slice contains no completed signal")
    elif warning is not None:
        raise AssertionError(f"Unexpected incomplete deltaN evolution:\n{process_output}")


def validate_outputs(work: Path, build_cfg: BuildConfig, run_cfg: RunConfig) -> None:
    results = work / "results"
    if not results.is_dir():
        raise AssertionError("Simulation did not create results/")

    actual = {
        path.relative_to(results).as_posix()
        for path in results.rglob("*")
        if path.is_file()
    }
    expected = expected_outputs(build_cfg, run_cfg)
    if actual != expected:
        missing = sorted(expected - actual)
        unexpected = sorted(actual - expected)
        raise AssertionError(f"Output inventory mismatch; missing={missing}, unexpected={unexpected}")

    for relative in sorted(actual):
        if relative in {"info.dat", "output.txt"}:
            if (results / relative).stat().st_size == 0:
                raise AssertionError(f"Empty text output: {relative}")
            continue
        numeric_rows(results / relative)

    time_files = ["sf.dat"]
    if run_cfg.output_spectra:
        time_files.append("spectratimes.dat")
    if run_cfg.output_histogram:
        time_files.append("histogramtimes.dat")
    if build_cfg.post:
        time_files.append("post_inflation/sf.dat")
        if run_cfg.output_spectra:
            time_files.append("post_inflation/spectratimes.dat")
        if run_cfg.output_histogram:
            time_files.append("post_inflation/histogramtimes.dat")
    for relative in time_files:
        assert_monotonic_time(results / relative)

    modes = [row[0] for row in numeric_rows(results / "modes.dat")]
    expected_mode_count = math.floor(math.sqrt(3.0) * (build_cfg.n // 2)) + 1
    if len(modes) != expected_mode_count:
        raise AssertionError(
            f"modes.dat has {len(modes)} rows, expected {expected_mode_count}"
        )
    if not math.isclose(modes[0], 0.0, abs_tol=1e-15) or any(
        right <= left for left, right in zip(modes, modes[1:])
    ):
        raise AssertionError("modes.dat must start at zero and be strictly increasing")
    if run_cfg.output_spectra:
        main_times = len(numeric_rows(results / "spectratimes.dat"))
        validate_spectrum_blocks(results, "spectra.dat", modes, main_times)
        if build_cfg.gw:
            validate_spectrum_blocks(
                results, "spectraGW.dat", modes, main_times, require_full_range=True
            )
            validate_spectrum_blocks(
                results, "spectraGWdot.dat", modes, main_times, require_full_range=True
            )
        if build_cfg.delta_n:
            validate_spectrum_blocks(results, "spectraN.dat", modes, 1)
        if build_cfg.post:
            post_times = len(numeric_rows(results / "post_inflation/spectratimes.dat"))
            validate_spectrum_blocks(
                results, "post_inflation/spectra.dat", modes, post_times
            )
            validate_spectrum_blocks(
                results,
                "post_inflation/spectraGW.dat",
                modes,
                post_times,
                require_full_range=True,
            )
            validate_spectrum_blocks(
                results,
                "post_inflation/spectraGWdot.dat",
                modes,
                post_times,
                require_full_range=True,
            )
    if build_cfg.delta_n and run_cfg.output_log:
        validate_spectrum_blocks(results, "spectraLOG.dat", modes, 1)

    if run_cfg.output_energy:
        validate_energy_constraint(results)
    if run_cfg.homogeneous_diagnostics:
        if not run_cfg.output_energy or not run_cfg.output_box2d:
            raise AssertionError("Homogeneous diagnostics require energy and 2D outputs")
        validate_acceleration_outputs(results, run_cfg, build_cfg.n)

    if run_cfg.output_histogram:
        validate_histogram(results, "histogram.dat", "histogramtimes.dat", 16)
        if build_cfg.post:
            validate_histogram(
                results,
                "post_inflation/histogram.dat",
                "post_inflation/histogramtimes.dat",
                16,
            )
    if build_cfg.delta_n:
        validate_histogram(results, "histogramN.dat", "histogramtimesN.dat", 16)
        if run_cfg.expect_degenerate_delta_n_histogram:
            if not run_cfg.output_spectra or not run_cfg.output_box2d:
                raise AssertionError(
                    "Degenerate deltaN regression requires spectrum and 2D outputs"
                )
            validate_degenerate_delta_n_histogram(results, 16)
        if run_cfg.output_log:
            validate_histogram(results, "histogramLOG.dat", "histogramtimesLOG.dat", 16)

    info = (results / "info.dat").read_text(encoding="utf-8")
    required_info = [
        f"Grid size={build_cfg.n}^3",
        f"spatial_stencil_order={build_cfg.stencil}",
        f"inflation_integrator={run_cfg.inflation_integrator}",
        f"linear_metric_perturbations={run_cfg.linear_metric}",
    ]
    if build_cfg.delta_n:
        required_info.append(f"deltaN_integrator={run_cfg.delta_n_integrator}")
    if build_cfg.post:
        required_info.append(f"post_inflation_integrator={run_cfg.post_integrator}")
    for marker in required_info:
        if marker not in info:
            raise AssertionError(f"Missing configuration marker in info.dat: {marker}")

    mass_match = re.search(r"(?m)^initial_mass_squared=(.+)$", info)
    if mass_match is None:
        raise AssertionError("Missing initial_mass_squared marker in info.dat")
    mass_value = mass_match.group(1).strip()
    if run_cfg.initial_mass is None:
        if mass_value != "not_set":
            raise AssertionError("An omitted initial mass must remain unset")
    else:
        try:
            reported_mass = float(mass_value)
        except ValueError as exc:
            raise AssertionError(f"Invalid initial mass marker: {mass_value}") from exc
        if not math.isclose(reported_mass, run_cfg.initial_mass, rel_tol=1e-14, abs_tol=1e-16):
            raise AssertionError(
                f"Initial mass marker mismatch: got {reported_mass}, expected {run_cfg.initial_mass}"
            )


def run_simulation(
    work: Path,
    build_cfg: BuildConfig,
    run_cfg: RunConfig,
    *,
    sanitizer_env: bool = False,
) -> None:
    results = work / "results"
    if results.exists():
        shutil.rmtree(results)
    write_params(work / "params.txt", runtime_values(build_cfg, run_cfg))
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "2" if build_cfg.parallel else "1"
    if sanitizer_env:
        detect_leaks = "0" if platform.system() == "Darwin" else "1"
        env["ASAN_OPTIONS"] = f"detect_leaks={detect_leaks}:halt_on_error=1"
        env["UBSAN_OPTIONS"] = "halt_on_error=1:print_stacktrace=1"
    proc = command(["./inflation_easy"], cwd=work, timeout=RUN_TIMEOUT, env=env)
    validate_outputs(work, build_cfg, run_cfg)
    if build_cfg.delta_n:
        validate_delta_n_reporting(
            results,
            run_cfg,
            proc.stdout,
            build_cfg.n * build_cfg.n * build_cfg.n,
        )


def invalid_runtime_case(
    work: Path,
    build_cfg: BuildConfig,
    run_cfg: RunConfig,
    key: str,
    value: str,
) -> None:
    values = runtime_values(build_cfg, run_cfg)
    values[key] = value
    write_params(work / "params.txt", values)
    results = work / "results"
    if results.exists():
        shutil.rmtree(results)
    proc = command(
        ["./inflation_easy"],
        cwd=work,
        timeout=10,
        env={**os.environ, "OMP_NUM_THREADS": "1"},
        expect_success=False,
    )
    if key.lower() not in proc.stdout.lower():
        raise AssertionError(
            f"Invalid {key} failed without a parameter-specific diagnostic:\n{proc.stdout}"
        )


def ignored_runtime_token_case(
    work: Path,
    build_cfg: BuildConfig,
    run_cfg: RunConfig,
    key: str,
    value: str,
) -> None:
    """Malformed tokens must be ignored rather than parsed as zero."""
    values = runtime_values(build_cfg, run_cfg)
    values[key] = value
    write_params(work / "params.txt", values)
    results = work / "results"
    if results.exists():
        shutil.rmtree(results)
    proc = command(
        ["./inflation_easy"],
        cwd=work,
        timeout=RUN_TIMEOUT,
        env={**os.environ, "OMP_NUM_THREADS": "1"},
    )
    expected = f"Ignoring unknown or invalid parameter '{key}'"
    if expected not in proc.stdout:
        raise AssertionError(
            f"Malformed {key} token was not reported as ignored:\n{proc.stdout}"
        )
    validate_outputs(work, build_cfg, run_cfg)


def invalid_numerical_table_case(
    work: Path,
    build_cfg: BuildConfig,
    run_cfg: RunConfig,
    replacements: dict[str, str],
    diagnostic: Union[str, tuple[str, ...]],
) -> None:
    """A malformed numerical-potential table must fail before initialization."""
    paths = {name: work / "inputs" / name for name in replacements}
    originals = {name: path.read_bytes() for name, path in paths.items()}
    try:
        for name, contents in replacements.items():
            paths[name].write_text(contents, encoding="utf-8")
        write_params(work / "params.txt", runtime_values(build_cfg, run_cfg))
        results = work / "results"
        if results.exists():
            shutil.rmtree(results)
        proc = command(
            ["./inflation_easy"],
            cwd=work,
            timeout=10,
            env={**os.environ, "OMP_NUM_THREADS": "1"},
            expect_success=False,
        )
        accepted_diagnostics = (diagnostic,) if isinstance(diagnostic, str) else diagnostic
        if not any(message.lower() in proc.stdout.lower() for message in accepted_diagnostics):
            raise AssertionError(
                "Malformed numerical tables lacked any accepted diagnostic "
                f"{accepted_diagnostics!r}:\n{proc.stdout}"
            )
    finally:
        for name, contents in originals.items():
            paths[name].write_bytes(contents)


class Matrix:
    def __init__(self, repo: Path, root: Path, *, sanitizers: bool = False):
        self.repo = repo
        self.root = root
        self.sanitizers = sanitizers
        self.built: dict[BuildConfig, Path] = {}
        self.build_count = 0
        self.run_count = 0

    def ensure(self, cfg: BuildConfig) -> Path:
        if cfg in self.built:
            return self.built[cfg]
        work = self.root / f"build-{len(self.built):02d}-{cfg.name}"
        copy_clean_tree(self.repo, work)
        configure_build(work, cfg)
        build(work, sanitizers=self.sanitizers)
        self.built[cfg] = work
        self.build_count += 1
        return work

    def run(self, cfg: BuildConfig, variant: RunConfig) -> None:
        work = self.ensure(cfg)
        run_simulation(work, cfg, variant, sanitizer_env=self.sanitizers)
        self.run_count += 1
        print(f"PASS run: {cfg.name}/{variant.name}")


def compile_failure(repo: Path, root: Path, cfg: BuildConfig, diagnostic: str) -> None:
    work = root / f"compile-failure-{cfg.name}"
    copy_clean_tree(repo, work)
    configure_build(work, cfg)
    proc = build(work, expect_success=False)
    if diagnostic.lower() not in proc.stdout.lower():
        raise AssertionError(
            f"Expected compile diagnostic containing {diagnostic!r}:\n{proc.stdout}"
        )
    print(f"PASS compile rejection: {cfg.name}")


def delta_n_local_energy_initialization_regression(repo: Path, root: Path) -> None:
    """Verify that initializeN uses the full local homogeneous energy in H.

    The probe is injected only into a clean temporary copy.  It reports the
    velocity, potential, and converted dphi/dN for one site immediately before
    the auxiliary deltaN evolution, allowing the test to check the Friedmann
    conversion directly without adding a diagnostic interface to production
    builds.
    """
    work = root / "deltaN-local-energy-initialization"
    copy_clean_tree(repo, work)
    cfg = BuildConfig(
        "deltaN-local-energy-initialization",
        delta_n=1,
        gw=0,
        post=0,
        parallel=0,
    )
    configure_build(work, cfg)

    replace_checked(
        work / "src" / "initialize.cpp",
        r"^(\s*)fd\[id\] = velocity_in_efold_time;\s*$",
        r'''\g<1>if (id == 0) {
\g<1>    const double probe_velocity =
\g<1>        fd[id] * std::pow(a, rescale_s - 1.0);
\g<1>    std::fprintf(stderr,
\g<1>        "DELTA_N_LOCAL_ENERGY_PROBE %.17g %.17g %.17g\\n",
\g<1>        probe_velocity, local_potential, velocity_in_efold_time);
\g<1>}
\g<1>fd[id] = velocity_in_efold_time;''',
    )
    build(work)

    run_cfg = RunConfig(
        "deltaN-local-energy-initialization",
        output_spectra=0,
        output_histogram=0,
        output_energy=0,
        output_box2d=0,
    )
    values = runtime_values(cfg, run_cfg)
    values["initial_derivative"] = "-0.1"
    write_params(work / "params.txt", values)
    proc = command(
        ["./inflation_easy"],
        cwd=work,
        timeout=RUN_TIMEOUT,
        env={**os.environ, "OMP_NUM_THREADS": "1"},
    )

    matches = re.findall(
        r"DELTA_N_LOCAL_ENERGY_PROBE\s+([-+0-9.eE]+)\s+"
        r"([-+0-9.eE]+)\s+([-+0-9.eE]+)",
        proc.stdout,
    )
    if len(matches) != 1:
        raise AssertionError(
            "Expected exactly one deltaN local-energy probe record; "
            f"got {len(matches)}:\n{proc.stdout}"
        )
    velocity, potential_value, converted = map(float, matches[0])
    local_energy = 0.5 * velocity * velocity + potential_value
    expected = velocity / math.sqrt(local_energy / 3.0)
    potential_only = velocity / math.sqrt(potential_value / 3.0)
    if not math.isclose(converted, expected, rel_tol=5e-13, abs_tol=5e-15):
        raise AssertionError(
            "deltaN velocity conversion does not use the full local homogeneous energy: "
            f"got {converted}, expected {expected}"
        )
    relative_gap = abs(expected - potential_only) / max(abs(expected), 1e-300)
    if relative_gap < 1e-4:
        raise AssertionError("The deltaN local-energy regression is numerically vacuous")
    print("PASS regression: deltaN initialization uses kinetic plus potential energy")


def run_ci_tier(repo: Path, root: Path) -> tuple[int, int]:
    matrix = Matrix(repo, root)
    full = BuildConfig("full-numerical")
    analytic = BuildConfig("full-analytic", numerical=0)
    scalar = BuildConfig("scalar-only", delta_n=0, gw=0, post=0, parallel=0)
    tiny_scalar = BuildConfig(
        "tiny-scalar-slices", delta_n=0, gw=0, post=0, parallel=0, n=4
    )
    post_no_dn = BuildConfig("post-without-deltaN", delta_n=0)
    order4 = BuildConfig("full-order4", stencil=4)
    order6 = BuildConfig("full-order6", stencil=6)
    delta_n_no_post = BuildConfig("deltaN-without-post", gw=0, post=0)

    for integrator in ("leapfrog", "rk4", "rk45"):
        matrix.run(
            full,
            RunConfig(
                f"all-{integrator}",
                inflation_integrator=integrator,
                delta_n_integrator=integrator,
                post_integrator=integrator,
            ),
        )
    matrix.run(
        full,
        RunConfig(
            "mixed-integrators",
            inflation_integrator="rk45",
            delta_n_integrator="leapfrog",
            post_integrator="rk4",
        ),
    )
    matrix.run(full, RunConfig("metric-and-mass", linear_metric=1, initial_mass=0.01022))
    matrix.run(full, RunConfig("log-output", output_log=1))
    matrix.run(full, RunConfig("bispectrum-output", output_bispectrum=1))
    matrix.run(analytic, RunConfig("analytic-negative-mass", initial_mass=-0.015))
    matrix.run(
        analytic,
        RunConfig(
            "zero-deltaN-budget",
            delta_n_end=0.0,
            low_cutoff_index=100.0,
            high_cutoff_index=100.0,
            expect_no_delta_n_steps=True,
            expect_degenerate_delta_n_histogram=True,
            homogeneous_diagnostics=True,
        ),
    )
    for integrator in ("leapfrog", "rk4", "rk45"):
        matrix.run(
            analytic,
            RunConfig(
                f"backward-deltaN-{integrator}",
                inflation_integrator=integrator,
                delta_n_integrator=integrator,
                post_integrator=integrator,
                delta_n_step=-0.0002,
                delta_n_end=0.01,
                phiref_manual=0.0935,
                low_cutoff_index=100.0,
                high_cutoff_index=100.0,
                output_box3d=1,
                expect_backward_crossing=True,
                homogeneous_diagnostics=True,
            ),
        )
    matrix.run(
        delta_n_no_post,
        RunConfig(
            "incomplete-deltaN-masking",
            delta_n_step=0.00001,
            delta_n_end=0.0001,
            expect_incomplete_delta_n=True,
        ),
    )
    matrix.run(scalar, RunConfig("scalar-only"))
    matrix.run(tiny_scalar, RunConfig("smallest-supported-slice-output"))
    matrix.run(
        scalar,
        RunConfig(
            "minimal-output-toggles",
            output_spectra=0,
            output_histogram=0,
            output_energy=0,
            output_box2d=0,
        ),
    )
    matrix.run(post_no_dn, RunConfig("post-without-deltaN"))
    matrix.run(order4, RunConfig("order4"))
    matrix.run(order6, RunConfig("order6"))

    scalar_work = matrix.ensure(scalar)
    for key in ("output_freq", "output_infrequent_freq", "nbins"):
        invalid_runtime_case(scalar_work, scalar, RunConfig("invalid"), key, "0")
        print(f"PASS runtime rejection: {key}=0")
    for key in ("output_spectra", "dt"):
        ignored_runtime_token_case(
            scalar_work, scalar, RunConfig(f"empty-{key}"), key, ""
        )
        matrix.run_count += 1
        print(f"PASS malformed-token rejection: {key}=<empty>")

    compile_failure(
        repo,
        root,
        BuildConfig("non-power-of-two-N", delta_n=0, gw=0, post=0, n=12),
        "power of two",
    )
    compile_failure(
        repo,
        root,
        BuildConfig("one-point-N", delta_n=0, gw=0, post=0, n=1),
        "at least 2",
    )
    compile_failure(
        repo,
        root,
        BuildConfig("post-without-GW", delta_n=0, gw=0, post=1),
        "requires calculate_SIGW",
    )
    delta_n_local_energy_initialization_regression(repo, root)
    return matrix.build_count + 1, matrix.run_count + 1


def run_release_tier(repo: Path, root: Path) -> tuple[int, int]:
    matrix = Matrix(repo, root)

    # Every valid compile-time feature combination at the default stencil.
    for numerical, delta_n, gw_state in itertools.product((0, 1), (0, 1), ("off", "inflation", "post")):
        gw = int(gw_state != "off")
        post = int(gw_state == "post")
        cfg = BuildConfig(
            f"features-p{numerical}-d{delta_n}-{gw_state}",
            numerical=numerical,
            delta_n=delta_n,
            gw=gw,
            post=post,
        )
        matrix.run(cfg, RunConfig("leapfrog"))

    full = BuildConfig("release-full-numerical")
    for inflation, delta_n, post in itertools.product(
        ("leapfrog", "rk4", "rk45"), repeat=3
    ):
        matrix.run(
            full,
            RunConfig(
                f"integrators-{inflation}-{delta_n}-{post}",
                inflation_integrator=inflation,
                delta_n_integrator=delta_n,
                post_integrator=post,
            ),
        )

    for numerical, order in itertools.product((0, 1), (4, 6)):
        cfg = BuildConfig(f"p{numerical}-order{order}", numerical=numerical, stencil=order)
        for integrator in ("leapfrog", "rk4", "rk45"):
            matrix.run(
                cfg,
                RunConfig(
                    f"all-{integrator}",
                    inflation_integrator=integrator,
                    delta_n_integrator=integrator,
                    post_integrator=integrator,
                ),
            )

    matrix.run(replace(full, name="release-full-serial", parallel=0), RunConfig("serial"))
    matrix.run(full, RunConfig("metric-on", linear_metric=1))
    matrix.run(full, RunConfig("mass-zero", initial_mass=0.0))
    matrix.run(full, RunConfig("mass-positive", initial_mass=0.01022))
    matrix.run(full, RunConfig("log-output", output_log=1))
    matrix.run(full, RunConfig("bispectrum-output", output_bispectrum=1))
    matrix.run(full, RunConfig("box3d-output", output_box2d=0, output_box3d=1))

    analytic = BuildConfig("release-full-analytic", numerical=0)
    matrix.run(analytic, RunConfig("mass-negative", initial_mass=-0.015))
    matrix.run(analytic, RunConfig("metric-on", linear_metric=1))
    matrix.run(
        analytic,
        RunConfig(
            "zero-deltaN-budget",
            delta_n_end=0.0,
            low_cutoff_index=100.0,
            high_cutoff_index=100.0,
            expect_no_delta_n_steps=True,
            expect_degenerate_delta_n_histogram=True,
            homogeneous_diagnostics=True,
        ),
    )
    for integrator in ("leapfrog", "rk4", "rk45"):
        matrix.run(
            analytic,
            RunConfig(
                f"backward-deltaN-{integrator}",
                inflation_integrator=integrator,
                delta_n_integrator=integrator,
                post_integrator=integrator,
                delta_n_step=-0.0002,
                delta_n_end=0.01,
                phiref_manual=0.0935,
                low_cutoff_index=100.0,
                high_cutoff_index=100.0,
                output_box3d=1,
                expect_backward_crossing=True,
                homogeneous_diagnostics=True,
            ),
        )

    matrix.run(
        BuildConfig("release-deltaN-without-post", gw=0, post=0),
        RunConfig(
            "incomplete-deltaN-masking",
            delta_n_step=0.00001,
            delta_n_end=0.0001,
            expect_incomplete_delta_n=True,
        ),
    )

    scalar = BuildConfig("release-scalar-invalid-tests", delta_n=0, gw=0, post=0, parallel=0)
    scalar_work = matrix.ensure(scalar)
    for key in ("output_freq", "output_infrequent_freq", "nbins", "dt"):
        invalid_runtime_case(scalar_work, scalar, RunConfig("invalid"), key, "0")

    table_files = {
        name: (scalar_work / "inputs" / name).read_text(encoding="utf-8").splitlines()
        for name in ("field_values.dat", "potential.dat", "potential_derivative.dat")
    }
    invalid_numerical_table_case(
        scalar_work,
        scalar,
        RunConfig("invalid-table-length"),
        {"potential_derivative.dat": "\n".join(table_files["potential_derivative.dat"][:-1]) + "\n"},
        "same number of entries",
    )
    invalid_numerical_table_case(
        scalar_work,
        scalar,
        RunConfig("invalid-table-minimum"),
        {name: lines[0] + "\n" for name, lines in table_files.items()},
        "at least two entries",
    )
    invalid_numerical_table_case(
        scalar_work,
        scalar,
        RunConfig("invalid-table-nonfinite"),
        {"potential.dat": "nan\n" + "\n".join(table_files["potential.dat"][1:]) + "\n"},
        ("non-finite value", "malformed numerical value"),
    )
    repeated_field = list(table_files["field_values.dat"])
    repeated_field[1] = repeated_field[0]
    invalid_numerical_table_case(
        scalar_work,
        scalar,
        RunConfig("invalid-table-order"),
        {"field_values.dat": "\n".join(repeated_field) + "\n"},
        "strictly descending",
    )

    compile_failure(
        repo,
        root,
        BuildConfig("release-non-power-of-two-N", delta_n=0, gw=0, post=0, n=12),
        "power of two",
    )
    compile_failure(
        repo,
        root,
        BuildConfig("release-one-point-N", delta_n=0, gw=0, post=0, n=1),
        "at least 2",
    )
    return matrix.build_count, matrix.run_count


def run_sanitizer_tier(repo: Path, root: Path) -> tuple[int, int]:
    matrix = Matrix(repo, root, sanitizers=True)
    cfg = BuildConfig("sanitizer-full", parallel=0)
    matrix.run(cfg, RunConfig("log-output", output_log=1))
    matrix.run(cfg, RunConfig("bispectrum-output", output_bispectrum=1))
    matrix.run(
        cfg,
        RunConfig(
            "degenerate-deltaN-histogram",
            delta_n_end=0.0,
            low_cutoff_index=100.0,
            high_cutoff_index=100.0,
            expect_no_delta_n_steps=True,
            expect_degenerate_delta_n_histogram=True,
        ),
    )
    matrix.run(
        cfg,
        RunConfig(
            "adaptive-integrators",
            inflation_integrator="rk45",
            delta_n_integrator="rk45",
            post_integrator="rk45",
        ),
    )
    return matrix.build_count, matrix.run_count


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--tier", choices=("ci", "release", "sanitizers"), default="ci")
    parser.add_argument("--keep-temp", action="store_true")
    args = parser.parse_args(argv)
    repo = args.repo.resolve()

    root = Path(tempfile.mkdtemp(prefix=f"inflationeasy-{args.tier}-", dir="/tmp"))
    succeeded = False

    try:
        if args.tier == "ci":
            builds, runs = run_ci_tier(repo, root)
        elif args.tier == "release":
            builds, runs = run_release_tier(repo, root)
        else:
            builds, runs = run_sanitizer_tier(repo, root)
        print(f"PASS {args.tier}: {builds} successful builds, {runs} successful simulations")
        succeeded = True
        return 0
    except Exception:
        print(f"FAILED {args.tier}; temporary artifacts: {root}", file=sys.stderr)
        raise
    finally:
        if succeeded and not args.keep_temp:
            shutil.rmtree(root)


if __name__ == "__main__":
    raise SystemExit(main())
