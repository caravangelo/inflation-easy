"""Exercise both scalar spectrum writers with live filtering enabled and disabled."""
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
        for phase in ("inflation", "radiation"):
            for force in ("0", "1"):
                result = subprocess.run([work / "inflation_easy", phase, force], cwd=work)
                if result.returncode:
                    failures.append((phase, force, result.returncode))
        if failures:
            raise AssertionError(f"Zero-mode regression failed: {failures}")
        print("PASS scalar spectrum filtering at every site in both phases")


if __name__ == "__main__":
    main()
