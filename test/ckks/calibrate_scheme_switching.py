#!/usr/bin/env python3
"""Reproduce scheme-switching security screens using Parameter-Selection.

Run with Sage Python. This screens LWE instances, not complete protocol
security or correctness. In particular, a zero-padded ring secret needs
additional structural analysis before any passing estimate is meaningful.
"""

import argparse
import concurrent.futures
import datetime
import json
import math
from pathlib import Path
import subprocess
import sys


CASES = {
    "ckks-integer4": (16384, 256, 3.2, 3),
    "ckks-prototype": (32768, 1000, 3.2, 3),
    "ckks-larger-ring": (65536, 1000, 3.2, 3),
    "ckks-shorter-chain": (32768, 832, 3.2, 3),
    "ckks-880-chain": (32768, 880, 3.2, 3),
    "clpx-reverse-half832-low-noise": (832, 32, 2**13, 2),
    "clpx-reverse-half832": (832, 32, 2**14, 2),
    "tfhe1024-default": (1024, 32, 2**7, 3),
    "tfhe760-alpha17p125": (760, 32, 2**(32-17.125), 2),
    "tfhe760-alpha17p25": (760, 32, 2**(32-17.25), 2),
    "tfhe760-alpha18": (760, 32, 2**14, 2),
    "tfhe760-fresh": (760, 32, 2**15, 2),
    "tfhe2048-default": (2048, 64, 2**13, 3),
    # The highest active radix-256 row has message 2^40*s_CKKS.
    # Reducing that public row modulo 2^40 removes its unknown message.
    # The padded secret has only 760 unknown binary coefficients.
    "compact760-message-eliminated": (760, 40, 3.2, 2),
    "compact2048-full-modulus-screen": (2048, 46, 3.2, 2),
    "tlwe2048-small-noise": (2048, 32, 3.2, 2),
    "tlwe4096-small-noise": (4096, 32, 3.2, 2),
}
ATTACKS = ("primal_usvp", "primal_bdd", "dual", "dual_hybrid")


def revision(path):
    try:
        return subprocess.check_output(
            ["git", "-C", str(path), "rev-parse", "HEAD"], text=True
        ).strip()
    except FileNotFoundError:
        # Parameter-Selection's Sage image does not include the git binary.
        gitdir = path / ".git"
        if gitdir.is_file():
            gitdir = (path / gitdir.read_text().strip().removeprefix("gitdir: ")).resolve()
        head = (gitdir / "HEAD").read_text().strip()
        if not head.startswith("ref: "):
            return head
        ref = head.removeprefix("ref: ")
        if (gitdir / ref).exists():
            return (gitdir / ref).read_text().strip()
        for line in (gitdir / "packed-refs").read_text().splitlines():
            if line.endswith(" " + ref):
                return line.split()[0]
        raise RuntimeError(f"cannot resolve revision for {path}")


def estimate(task):
    root, name, attack = task
    sys.path.insert(0, str(Path(root) / "python" / "lattice-estimator"))
    from estimator import LWE, ND, RC
    from sage.all import log

    n, log_q, sigma, alphabet = CASES[name]
    params = LWE.Parameters(
        n=n, q=2**log_q, Xs=ND.UniformMod(alphabet),
        Xe=ND.DiscreteGaussian(sigma), tag=name,
    )
    try:
        result = getattr(LWE, attack)(params, red_cost_model=RC.BDGL16)
        bits = float(log(result["rop"], 2))
        if not math.isfinite(bits):
            raise ValueError("non-finite attack cost")
        return name, attack, {"bits": bits, "cost": str(result)}
    except Exception as exc:
        return name, attack, {"error": f"{type(exc).__name__}: {exc}"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parameter-selection", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--case", choices=CASES, action="append")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--target", type=float, default=128)
    args = parser.parse_args()
    if args.jobs < 1 or not math.isfinite(args.target) or args.target <= 0:
        parser.error("jobs and target must be positive")
    root = args.parameter_selection.resolve()
    names = args.case or list(CASES)
    report = {
        "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "parameter_selection_revision": revision(root),
        "estimator_revision": revision(root / "python" / "lattice-estimator"),
        "cost_model": "RC.BDGL16 (classical)",
        "samples": "unbounded; conservative generic LWE screen",
        "target_bits": args.target,
        "scope": "Four-attack LWE screen; not a protocol security certificate. "
                 "Known-support ring structure and circular/KDM security are not modeled.",
        "cases": {},
    }
    for name in names:
        n, log_q, sigma, alphabet = CASES[name]
        report["cases"][name] = dict(
            n=n, log_q=log_q, sigma=sigma, secret_alphabet=alphabet, attacks={}
        )
    tasks = [(str(root), name, attack) for name in names for attack in ATTACKS]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(estimate, task) for task in tasks]
        for future in concurrent.futures.as_completed(futures):
            name, attack, result = future.result()
            report["cases"][name]["attacks"][attack] = result
            print(name, attack, result.get("bits", result.get("error")), flush=True)
            args.output.write_text(json.dumps(report, indent=2) + "\n")
    complete = True
    for name, case in report["cases"].items():
        attacks = case["attacks"]
        case["complete"] = all("bits" in attacks[a] for a in ATTACKS)
        complete &= case["complete"]
        case["weakest_bits"] = min(
            (a["bits"] for a in attacks.values() if "bits" in a), default=None
        )
        case["passes_lwe_screen"] = case["complete"] and case["weakest_bits"] >= args.target
        status = ("INCOMPLETE" if not case["complete"] else
                  "PASS" if case["passes_lwe_screen"] else "FAIL")
        print(name, "weakest", case["weakest_bits"], status, flush=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    # Rejected candidates are an expected calibration result. An incomplete
    # estimator run is an execution failure and must not silently pass.
    return 0 if complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
