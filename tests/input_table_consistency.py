#!/usr/bin/env python3
"""Check that the default USR force table differentiates its potential table."""

from __future__ import annotations

import argparse
import math
from pathlib import Path


def values(path: Path) -> list[float]:
    result = [float(line) for line in path.read_text().splitlines() if line.strip()]
    if not result or not all(math.isfinite(value) for value in result):
        raise AssertionError(f"Empty or non-finite input table: {path}")
    return result


def check(inputs: Path) -> tuple[int, float]:
    field = values(inputs / "field_values.dat")
    potential = values(inputs / "potential.dat")
    force = values(inputs / "potential_derivative.dat")
    if len(field) < 3 or len(field) != len(potential) or len(field) != len(force):
        raise AssertionError("The three USR tables must have the same length of at least three")
    if any(left <= right for left, right in zip(field, field[1:])):
        raise AssertionError("The field grid must be strictly descending")

    force_scale = max(abs(value) for value in force)
    if force_scale == 0.0:
        raise AssertionError("The force table is identically zero")
    largest_relative_error = 0.0
    for index in range(1, len(field) - 1):
        left = field[index] - field[index - 1]
        right = field[index + 1] - field[index]
        left_secant = (potential[index] - potential[index - 1]) / left
        right_secant = (potential[index + 1] - potential[index]) / right
        numerical_derivative = (right * left_secant + left * right_secant) / (left + right)
        scale = max(abs(force[index]), force_scale * 1e-3)
        largest_relative_error = max(
            largest_relative_error,
            abs(numerical_derivative - force[index]) / scale,
        )

    # The dense table is sampled from a smooth spline; this leaves a generous
    # margin above its finite-difference truncation error while rejecting the
    # original derivative table constructed on an uneven field grid.
    if largest_relative_error >= 1e-3:
        raise AssertionError(
            f"V' is inconsistent with V(phi): maximum scaled error "
            f"{largest_relative_error:.3g} exceeds 0.001"
        )
    return len(field), largest_relative_error


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", type=Path, default=Path(__file__).resolve().parents[1] / "inputs")
    args = parser.parse_args()
    count, error = check(args.inputs)
    print(f"PASS USR input tables: {count} rows, maximum scaled derivative error {error:.3g}")


if __name__ == "__main__":
    main()
