"""Check Fourier cutoffs at initialization and during scalar spectrum output."""
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

import release_smoke as smoke


def main():
    repo = Path(__file__).resolve().parents[1]
    with TemporaryDirectory() as temp:
        work = Path(temp) / "probe"
        smoke.copy_clean_tree(repo, work)
        smoke.configure_build(work, smoke.BuildConfig("zero-mode", n=8, parallel=0))
        driver = work / "src/main.cpp"
        smoke.replace_checked(driver, r"^int main\(", "int simulation_main(")
        with driver.open("a") as output:
            output.write("\n" + (repo / "tests/spectrum_zero_mode.cpp").read_text())
        smoke.build(work)
        (work / "results/post_inflation").mkdir(parents=True)
        failures = []
        def initial_modes(low, high):
            # Each process starts the production RNG at seed 8.
            output = subprocess.check_output(
                [work / "inflation_easy", "init", "0", str(low), str(high)], cwd=work, text=True)
            return [tuple(map(float, line.split())) for line in output.splitlines()]

        baseline = initial_modes(0, 0)
        assert len(baseline) == 3 and all(any(mode) for mode in baseline)
        # With bounds at 2 the probe has modes below, on, and above each bound.
        for low, high in ((0, 0), (2, 0), (0, 2), (2, 2)):
            modes = initial_modes(low, high)
            assert len(modes) == len(baseline)
            for n, mode in enumerate(modes, 1):
                expected = baseline[n-1] if n >= low and (not high or n <= high) else (0.0,)*4
                if mode != expected:
                    failures.append(("init", low, high, n))
            for phase in ("inflation", "radiation"):
                for force in ("0", "1"):
                    result = subprocess.run(
                        [work / "inflation_easy", phase, force, str(low), str(high)], cwd=work)
                    if result.returncode:
                        failures.append((phase, force, low, high))
        if failures:
            raise AssertionError(f"Fourier cutoff regression failed: {failures}")
        print("PASS Fourier initialization and live filtering in both phases")


if __name__ == "__main__":
    main()
