#!/usr/bin/env python3
"""
Regression test: current branch vs main at N=16 (leapfrog).

Runs both revisions in isolated temp directories under /tmp, using identical
runtime parameters, then compares numeric output files.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

import numpy as np


NUMERIC_FILES = [
    "sf.dat",
    "means.dat",
    "velocity.dat",
    "variance.dat",
    "energy.dat",
    "conservation.dat",
    "spectra.dat",
    "spectratimes.dat",
    "histogram.dat",
    "histogramtimes.dat",
    "spectraN.dat",
    "histogramN.dat",
    "histogramtimesN.dat",
    "post_inflation/sf.dat",
    "post_inflation/means.dat",
    "post_inflation/velocity.dat",
    "post_inflation/variance.dat",
    "post_inflation/spectra.dat",
    "post_inflation/spectratimes.dat",
    "post_inflation/histogram.dat",
    "post_inflation/histogramtimes.dat",
]


def run(cmd: list[str], cwd: Path | None = None, env: dict[str, str] | None = None) -> None:
    subprocess.run(cmd, cwd=cwd, env=env, check=True)


def export_git_ref(repo: Path, ref: str, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        ["git", "archive", ref],
        cwd=repo,
        check=True,
        stdout=subprocess.PIPE,
    )
    with tempfile.NamedTemporaryFile(delete=False) as tmp:
        tmp.write(proc.stdout)
        tmp_path = Path(tmp.name)
    try:
        with tarfile.open(tmp_path) as tf:
            tf.extractall(path=out_dir)
    finally:
        tmp_path.unlink(missing_ok=True)


def patch_n_to_16(parameters_h: Path) -> None:
    txt = parameters_h.read_text()
    txt_new = re.sub(r"const int N = \d+;", "const int N = 16;", txt)
    if txt == txt_new:
        raise RuntimeError("Failed to patch N in parameters.h")
    parameters_h.write_text(txt_new)


def write_test_params(work: Path, source_params: Path) -> None:
    txt = source_params.read_text()
    # Strip any explicit per-loop integrator keys so we can enforce below.
    txt = re.sub(
        r"^\s*(inflation_integrator|deltaN_integrator|post_inflation_integrator)\s*=.*\n",
        "",
        txt,
        flags=re.MULTILINE,
    )
    txt += (
        "\n"
        "inflation_integrator = leapfrog\n"
        "deltaN_integrator = leapfrog\n"
        "post_inflation_integrator = leapfrog\n"
    )
    (work / "params.txt").write_text(txt)


def build_and_run(work: Path, inputs_dir: Path) -> None:
    if not inputs_dir.exists():
        raise RuntimeError(f"inputs directory not found: {inputs_dir}")
    work_inputs = work / "inputs"
    if not work_inputs.exists():
        os.symlink(inputs_dir, work_inputs)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    run(["make", "-j4"], cwd=work, env=env)
    run(["./inflation_easy"], cwd=work, env=env)


def load_numeric(path: Path) -> np.ndarray:
    return np.loadtxt(path)


def compare_results(a: Path, b: Path, rtol: float, atol: float) -> list[str]:
    failures: list[str] = []
    for rel in NUMERIC_FILES:
        pa = a / "results" / rel
        pb = b / "results" / rel
        if not pa.exists() or not pb.exists():
            continue
        da = load_numeric(pa)
        db = load_numeric(pb)
        if da.shape != db.shape:
            failures.append(f"{rel}: shape mismatch {da.shape} vs {db.shape}")
            continue
        if da.size == 0 and db.size == 0:
            continue
        if not np.allclose(da, db, rtol=rtol, atol=atol, equal_nan=True):
            diff = np.abs(da - db)
            max_abs = float(np.nanmax(diff))
            denom = np.maximum(np.abs(db), atol)
            max_rel = float(np.nanmax(diff / denom))
            failures.append(f"{rel}: max_abs={max_abs:.3e}, max_rel={max_rel:.3e}")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=".", help="Path to git repo")
    ap.add_argument("--main-ref", default="main", help="Reference branch/tag to compare against")
    ap.add_argument("--params", default="params.numerical.txt", help="Runtime params template from repo root")
    ap.add_argument("--rtol", type=float, default=1e-10)
    ap.add_argument("--atol", type=float, default=1e-12)
    args = ap.parse_args()

    repo = Path(args.repo).resolve()
    params_src = (repo / args.params).resolve()
    inputs_dir = (repo / "inputs").resolve()

    if not params_src.exists():
        print(f"Missing params template: {params_src}", file=sys.stderr)
        return 2

    with tempfile.TemporaryDirectory(prefix="ie-reg-main-n16-", dir="/tmp") as td:
        root = Path(td)
        work_head = root / "head"
        work_main = root / "main"

        export_git_ref(repo, "HEAD", work_head)
        export_git_ref(repo, args.main_ref, work_main)

        patch_n_to_16(work_head / "src" / "parameters.h")
        patch_n_to_16(work_main / "src" / "parameters.h")

        write_test_params(work_head, params_src)
        write_test_params(work_main, params_src)

        build_and_run(work_head, inputs_dir)
        build_and_run(work_main, inputs_dir)

        failures = compare_results(work_head, work_main, rtol=args.rtol, atol=args.atol)
        if failures:
            print("Regression mismatch vs main at N=16:")
            for f in failures:
                print(f" - {f}")
            print(f"Temp artifacts kept at: {root}")
            return 1

        print("PASS: current branch reproduces main numerically at N=16 (leapfrog, same params).")
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
