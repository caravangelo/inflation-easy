#!/usr/bin/env python3
"""Strict output-equivalence audit for source-only evolution refactors.

The audit exports a reference revision and the current revision into separate
temporary trees, builds representative compile-time configurations, runs both
with one OpenMP thread, and compares every emitted result file byte for byte.
Only the wall-clock lines in ``info.dat`` are removed before comparison.

This is deliberately stricter than the release smoke suite.  The smoke suite
checks physical and file-format invariants within one revision; this utility
checks that an architectural refactor has not changed any deterministic output
relative to the selected baseline.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


TESTS_DIR = Path(__file__).resolve().parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

import release_smoke as smoke  # noqa: E402


@dataclass(frozen=True)
class EquivalenceCase:
    """One compile-time build and run-time configuration to compare."""

    name: str
    build: smoke.BuildConfig
    run: smoke.RunConfig


FULL = smoke.BuildConfig("equiv-full", parallel=1, n=8)

CASES = (
    EquivalenceCase(
        "full-leapfrog-all-outputs",
        FULL,
        smoke.RunConfig(
            "full-leapfrog-all-outputs",
            output_log=1,
            output_bispectrum=1,
            output_box3d=1,
        ),
    ),
    EquivalenceCase(
        "full-mixed-integrators-metric-mass",
        FULL,
        smoke.RunConfig(
            "full-mixed-integrators-metric-mass",
            inflation_integrator="rk45",
            delta_n_integrator="leapfrog",
            post_integrator="rk4",
            linear_metric=1,
            initial_mass=0.01022,
        ),
    ),
    EquivalenceCase(
        "full-adaptive-integrators",
        FULL,
        smoke.RunConfig(
            "full-adaptive-integrators",
            inflation_integrator="rk45",
            delta_n_integrator="rk45",
            post_integrator="rk45",
        ),
    ),
    EquivalenceCase(
        "analytic-scalar-rk4",
        smoke.BuildConfig(
            "equiv-analytic-scalar",
            numerical=0,
            delta_n=0,
            gw=0,
            post=0,
            parallel=0,
            n=8,
        ),
        smoke.RunConfig(
            "analytic-scalar-rk4",
            inflation_integrator="rk4",
            linear_metric=1,
            initial_mass=-0.015,
        ),
    ),
    EquivalenceCase(
        "analytic-backward-deltaN-rk45",
        smoke.BuildConfig(
            "equiv-analytic-deltaN",
            numerical=0,
            delta_n=1,
            gw=0,
            post=0,
            parallel=0,
            n=8,
        ),
        smoke.RunConfig(
            "analytic-backward-deltaN-rk45",
            inflation_integrator="rk4",
            delta_n_integrator="rk45",
            delta_n_step=-0.0002,
            delta_n_end=0.01,
            phiref_manual=0.0935,
            low_cutoff_index=100.0,
            high_cutoff_index=100.0,
            output_box3d=1,
            expect_backward_crossing=True,
        ),
    ),
    EquivalenceCase(
        "post-inflation-without-deltaN",
        smoke.BuildConfig(
            "equiv-post-without-deltaN",
            delta_n=0,
            gw=1,
            post=1,
            parallel=0,
            n=8,
        ),
        smoke.RunConfig(
            "post-inflation-without-deltaN",
            inflation_integrator="rk4",
            post_integrator="rk45",
        ),
    ),
    EquivalenceCase(
        "inflationary-GW-order4",
        smoke.BuildConfig(
            "equiv-inflationary-GW-order4",
            delta_n=0,
            gw=1,
            post=0,
            stencil=4,
            parallel=0,
            n=8,
        ),
        smoke.RunConfig("inflationary-GW-order4"),
    ),
    EquivalenceCase(
        "full-order6-rk4",
        smoke.BuildConfig(
            "equiv-full-order6",
            stencil=6,
            parallel=0,
            n=8,
        ),
        smoke.RunConfig(
            "full-order6-rk4",
            inflation_integrator="rk4",
            delta_n_integrator="rk4",
            post_integrator="rk4",
        ),
    ),
)

CASE_BY_NAME = {case.name: case for case in CASES}


@dataclass(frozen=True)
class RunResult:
    stdout: bytes
    stderr: bytes
    results: Path


def resolve_ref(repo: Path, ref: str) -> str:
    proc = subprocess.run(
        ["git", "rev-parse", "--verify", f"{ref}^{{commit}}"],
        cwd=repo,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return proc.stdout.strip()


def export_git_ref(repo: Path, ref: str, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=False)
    proc = subprocess.run(
        ["git", "archive", "--format=tar", ref],
        cwd=repo,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    with tarfile.open(fileobj=io.BytesIO(proc.stdout), mode="r:") as archive:
        archive.extractall(destination)


def export_worktree(repo: Path, destination: Path) -> None:
    """Copy only release inputs, allowing an uncommitted refactor to be audited."""
    smoke.copy_clean_tree(repo, destination)


def build_tree(
    source_tree: Path,
    destination: Path,
    config: smoke.BuildConfig,
    jobs: int,
) -> Path:
    smoke.copy_clean_tree(source_tree, destination)
    smoke.configure_build(destination, config)
    proc = subprocess.run(
        ["make", f"-j{jobs}"],
        cwd=destination,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=smoke.BUILD_TIMEOUT,
    )
    if proc.returncode != 0:
        raise AssertionError(
            f"Build failed for {config.name} in {destination}:\n{proc.stdout}"
        )
    if re.search(r"(?mi)\bwarning:", proc.stdout):
        raise AssertionError(
            f"Compiler warning for {config.name} in {destination}:\n{proc.stdout}"
        )
    return destination


def run_case(work: Path, case: EquivalenceCase) -> RunResult:
    results = work / "results"
    if results.exists():
        shutil.rmtree(results)

    smoke.write_params(
        work / "params.txt",
        smoke.runtime_values(case.build, case.run),
    )
    env = os.environ.copy()
    env.update(
        {
            "LC_ALL": "C",
            "OMP_DYNAMIC": "FALSE",
            "OMP_NUM_THREADS": "1",
            "TZ": "UTC",
        }
    )
    try:
        proc = subprocess.run(
            ["./inflation_easy"],
            cwd=work,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=smoke.RUN_TIMEOUT,
        )
    except subprocess.TimeoutExpired as exc:
        raise AssertionError(
            f"Simulation timed out for {case.name} in {work}"
        ) from exc

    if proc.returncode != 0:
        output = (proc.stdout + proc.stderr).decode("utf-8", errors="replace")
        raise AssertionError(
            f"Simulation failed for {case.name} in {work}:\n{output}"
        )

    smoke.validate_outputs(work, case.build, case.run)
    if case.build.delta_n:
        smoke.validate_delta_n_reporting(
            results,
            case.run,
            (proc.stdout + proc.stderr).decode("utf-8", errors="replace"),
            case.build.n ** 3,
        )
    return RunResult(proc.stdout, proc.stderr, results)


def normalize_info(contents: bytes) -> bytes:
    """Remove only the three wall-clock-dependent records from info.dat."""
    kept: list[bytes] = []
    prefixes = (b"Run began at ", b"Run ended at ", b"Run from t=")
    for line in contents.splitlines(keepends=True):
        if line.startswith(prefixes):
            continue
        kept.append(line)
    return b"".join(kept)


def result_inventory(results: Path) -> dict[str, Path]:
    return {
        path.relative_to(results).as_posix(): path
        for path in sorted(results.rglob("*"))
        if path.is_file()
    }


def digest(contents: bytes) -> str:
    return hashlib.sha256(contents).hexdigest()[:16]


def first_differing_line(left: bytes, right: bytes) -> str:
    left_lines = left.decode("utf-8", errors="replace").splitlines()
    right_lines = right.decode("utf-8", errors="replace").splitlines()
    limit = max(len(left_lines), len(right_lines))
    for index in range(limit):
        left_line = left_lines[index] if index < len(left_lines) else "<missing>"
        right_line = right_lines[index] if index < len(right_lines) else "<missing>"
        if left_line != right_line:
            return (
                f"first difference at line {index + 1}: "
                f"current={left_line!r}, reference={right_line!r}"
            )
    return "contents differ at the byte level"


def compare_bytes(label: str, current: bytes, reference: bytes) -> list[str]:
    if current == reference:
        return []
    return [
        f"{label}: {first_differing_line(current, reference)}; "
        f"sha256 current={digest(current)}, reference={digest(reference)}"
    ]


def compare_case(
    case: EquivalenceCase,
    current: RunResult,
    reference: RunResult,
) -> list[str]:
    failures: list[str] = []
    failures.extend(compare_bytes("stdout", current.stdout, reference.stdout))
    failures.extend(compare_bytes("stderr", current.stderr, reference.stderr))

    current_files = result_inventory(current.results)
    reference_files = result_inventory(reference.results)
    current_names = set(current_files)
    reference_names = set(reference_files)
    if current_names != reference_names:
        failures.append(
            "output inventory differs: "
            f"missing={sorted(reference_names - current_names)}, "
            f"unexpected={sorted(current_names - reference_names)}"
        )

    for relative in sorted(current_names & reference_names):
        current_contents = current_files[relative].read_bytes()
        reference_contents = reference_files[relative].read_bytes()
        if relative == "info.dat":
            current_contents = normalize_info(current_contents)
            reference_contents = normalize_info(reference_contents)
        failures.extend(
            compare_bytes(relative, current_contents, reference_contents)
        )

    return [f"[{case.name}] {failure}" for failure in failures]


def selected_cases(names: Optional[list[str]]) -> list[EquivalenceCase]:
    if not names:
        return list(CASES)
    ordered: list[EquivalenceCase] = []
    seen: set[str] = set()
    for name in names:
        if name not in seen:
            ordered.append(CASE_BY_NAME[name])
            seen.add(name)
    return ordered


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument(
        "--reference",
        default="v1.1.0",
        help="Immutable pre-refactor Git revision (default: v1.1.0)",
    )
    parser.add_argument(
        "--current",
        default="HEAD",
        help="Current Git revision to compare (default: HEAD)",
    )
    parser.add_argument(
        "--current-worktree",
        action="store_true",
        help="Compare the current release-input working tree instead of --current",
    )
    parser.add_argument(
        "--case",
        action="append",
        choices=sorted(CASE_BY_NAME),
        help="Run one named case; repeat for multiple cases (default: all)",
    )
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument(
        "--allow-identical",
        action="store_true",
        help="Permit identical current/reference commits for a harness self-check",
    )
    parser.add_argument("--keep-temp", action="store_true")
    args = parser.parse_args(argv)

    repo = args.repo.resolve()
    if args.jobs < 1:
        parser.error("--jobs must be positive")

    try:
        reference_sha = resolve_ref(repo, args.reference)
        current_sha = None if args.current_worktree else resolve_ref(repo, args.current)
    except subprocess.CalledProcessError as exc:
        diagnostic = exc.stderr.strip() if exc.stderr else str(exc)
        print(f"Unable to resolve comparison revision: {diagnostic}", file=sys.stderr)
        return 2

    if current_sha == reference_sha and not args.allow_identical:
        print(
            "Current and reference revisions are identical; use --allow-identical "
            "only for a deliberate harness self-check.",
            file=sys.stderr,
        )
        return 2

    cases = selected_cases(args.case)
    root = Path(tempfile.mkdtemp(prefix="inflationeasy-refactor-equivalence-", dir="/tmp"))
    succeeded = False
    try:
        reference_source = root / "reference-source"
        current_source = root / "current-source"
        export_git_ref(repo, reference_sha, reference_source)
        if args.current_worktree:
            export_worktree(repo, current_source)
            current_label = "working tree"
        else:
            assert current_sha is not None
            export_git_ref(repo, current_sha, current_source)
            current_label = current_sha[:12]

        reference_builds: dict[smoke.BuildConfig, Path] = {}
        current_builds: dict[smoke.BuildConfig, Path] = {}
        failures: list[str] = []

        for case in cases:
            if case.build not in reference_builds:
                index = len(reference_builds)
                reference_builds[case.build] = build_tree(
                    reference_source,
                    root / f"reference-build-{index:02d}-{case.build.name}",
                    case.build,
                    args.jobs,
                )
                current_builds[case.build] = build_tree(
                    current_source,
                    root / f"current-build-{index:02d}-{case.build.name}",
                    case.build,
                    args.jobs,
                )

            reference_result = run_case(reference_builds[case.build], case)
            current_result = run_case(current_builds[case.build], case)
            case_failures = compare_case(case, current_result, reference_result)
            if case_failures:
                failures.extend(case_failures)
                print(f"FAIL {case.name}")
            else:
                print(f"PASS {case.name}")

        if failures:
            print("Strict refactor equivalence failed:", file=sys.stderr)
            for failure in failures:
                print(f" - {failure}", file=sys.stderr)
            print(f"Temporary artifacts retained at: {root}", file=sys.stderr)
            return 1

        print(
            "PASS strict refactor equivalence: "
            f"{len(reference_builds)} builds per revision, {len(cases)} paired runs; "
            f"current={current_label}, reference={reference_sha[:12]}"
        )
        succeeded = True
        return 0
    except Exception:
        print(f"FAILED; temporary artifacts retained at: {root}", file=sys.stderr)
        raise
    finally:
        if succeeded and not args.keep_temp:
            shutil.rmtree(root)


if __name__ == "__main__":
    raise SystemExit(main())
