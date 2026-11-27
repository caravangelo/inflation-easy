"""Measure temporal order for issue #12 on a fixed 8³ grid at common final times."""
import math
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

import release_smoke as smoke


def main():
    repo = Path(__file__).resolve().parents[1]
    failures = []
    with TemporaryDirectory() as temp:
        work = Path(temp) / "probe"
        smoke.copy_clean_tree(repo, work)
        smoke.configure_build(work, smoke.BuildConfig("probe", parallel=0))
        smoke.replace_checked(work / "src/main.cpp", r"^int main\(", "int simulation_main(")
        with (work / "src/evolution/post_inflation.cpp").open("a") as source:
            source.write("\n" + (repo / "tests/leapfrog_convergence.cpp").read_text())
        smoke.build(work)
        (work / "results/post_inflation").mkdir(parents=True)

        def run(case, method, steps):
            return subprocess.check_output(
                [work / "inflation_easy", case, method, str(steps)], cwd=work, text=True)

        def values(case, method, steps):
            return [float(token) for token in run(case, method, steps).split()]

        # Flat potential and zero scalar velocity give de Sitter a = 1/(1-t).
        # In radiation a = 1+t. For w=0.2, a''=q*a'^2/a and a'=a^q.
        q = 0.2
        final_a = (2-q)**(1/(1-q))
        exact = {"inflation_tensor": [0.6**2]*6, "radiation": [0.25]*6 + [1/16],
                 "background": [final_a, final_a**q]}
        for case in ("inflation_tensor", "metric", "radiation", "source", "background"):
            target = exact.get(case)
            coarse_ref = None
            if target is None:
                coarse_ref = values(case, "rk4", 128)
                target = values(case, "rk4", 256)
            samples = [values(case, "lf", n) for n in (8, 16, 32)]
            assert all(len(sample) == len(target) for sample in samples)
            for component, (reference, *observed) in enumerate(zip(target, *samples)):
                errors = [abs(value-reference) for value in observed]
                orders = [math.log2(e1/e2) for e1, e2 in zip(errors, errors[1:])]
                if coarse_ref is not None:
                    # Check RK4 resolution against the finest tested LF error.
                    if abs(coarse_ref[component]-reference) > errors[-1] / 10:
                        failures.append((case, component, "reference resolution"))
                if component == 0 or case in ("source", "background") or component == 6:
                    print(f"{case}[{component}] errors {errors}, orders {orders}")
                if not all(1.7 < order < 2.4 for order in orders):
                    failures.append((case, component, orders))

        # Isolate the velocity-gradient source at the origin in one kick.
        # Radiation gives scalar friction 4 and tensor friction 2 at a=a'=1.
        step = 1/8
        dx = 2*math.pi/8
        grad_mid = math.sin(dx) / dx / (1 + 2*step)
        expected = 2*step*grad_mid**2 / (1 + step)
        source = values("midpoint_source", "lf", 8)
        if len(source) != 6 or not all(math.isclose(v, expected, rel_tol=1e-6) for v in source):
            failures.append(("midpoint source", source, expected))

        run("stop", "lf", 8)
        rows = (work / "results/post_inflation/sf.dat").read_text().splitlines()
        time, saved_a = map(float, rows[-1].split()[:2])
        print(f"stopping target {1+1/32}, saved scale factor {saved_a} at time {time}")
        if (time, saved_a) != (1/8, 1+1/8):
            failures.append(("stop", time, saved_a))
    if failures:
        raise AssertionError(failures)


if __name__ == "__main__":
    main()
