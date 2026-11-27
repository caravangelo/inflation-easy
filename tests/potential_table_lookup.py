"""Check tabulated potential values and range errors for cached and uncached lookup."""
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

import release_smoke as smoke


def main():
    repo = Path(__file__).resolve().parents[1]
    with TemporaryDirectory() as temp:
        work = Path(temp) / "probe"
        smoke.copy_clean_tree(repo, work)
        smoke.configure_build(work, smoke.BuildConfig("table", parallel=0))
        driver = work / "src/main.cpp"
        smoke.replace_checked(driver, r"^int main\(", "int simulation_main(")
        with driver.open("a") as output:
            output.write("\n" + (repo / "tests/potential_table_lookup.cpp").read_text())
        smoke.build(work)
        for case in ("below", "above"):
            result = subprocess.run([work / "inflation_easy", case], cwd=work,
                                    capture_output=True, text=True)
            if result.returncode != 1 or "outside the numerical potential range" not in result.stderr:
                raise AssertionError((case, result.returncode, result.stderr))
        subprocess.run([work / "inflation_easy", "knots"], cwd=work, check=True)
        print("PASS potential table knots, lookup hints, and both range errors")


if __name__ == "__main__":
    main()
