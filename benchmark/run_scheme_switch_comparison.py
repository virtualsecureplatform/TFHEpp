#!/usr/bin/env python3
"""Run scalar scheme comparisons sequentially, preserving reproducible evidence."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import subprocess


def git(repo, *args):
    return subprocess.check_output(["git", "-C", str(repo), *args], text=True).strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path,
                        help="new output directory; existing directories are never overwritten")
    parser.add_argument("--schemes", nargs="+", choices=("bfv", "clpx", "ckks"),
                        default=["bfv", "clpx", "ckks"])
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    executable = args.executable.resolve(strict=True)
    repo = Path(__file__).resolve().parents[1]
    args.output.mkdir(parents=True, exist_ok=False)
    cpu = platform.processor()
    if Path("/proc/cpuinfo").exists():
        match = re.search(r"model name\s*:\s*(.+)", Path("/proc/cpuinfo").read_text())
        if match:
            cpu = match[1]
    sources = ["test/clpx/scheme_switch_comparison.cpp", "include/bfv/scheme_switching.hpp",
               "include/bfv/bfv++.hpp", "include/tfhe/homdecomp.hpp", "include/tfhe/detwfa.hpp",
               "include/params/128bit.hpp", "include/clpx/bfv-clpx.hpp",
               "include/ckks/integer_multiplication.hpp"]
    report = {
        "utc": datetime.now(timezone.utc).isoformat(), "cpu": cpu,
        "threads": args.threads, "address_space_cap_bytes": 24 * 1024**3,
        "revision": git(repo, "rev-parse", "HEAD"), "status": git(repo, "status", "--short"),
        "source_sha256": {name: hashlib.sha256((repo / name).read_bytes()).hexdigest()
                          for name in sources},
        "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "method": "256 fixed 8-bit pairs; separate sequential processes; encryption, decryption and diagnostics excluded; key generation separate",
        "interfaces": {"bfv": "Boolean 8+8 inputs, 16 outputs", "clpx": "Boolean 8+8 inputs, 16 outputs",
                       "ckks": "one integer-phase TLWE per operand and product"},
        "runs": {}, "complete": False, "all_passed": False,
    }
    flags = executable.parents[2] / "test/CMakeFiles/scheme_switch_comparison.dir/flags.make"
    report["compiler_flags"] = flags.read_text() if flags.exists() else None
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=str(args.threads), OMP_DYNAMIC="FALSE",
               OMP_STACKSIZE="32M", MALLOC_ARENA_MAX="2")
    report["environment"] = {name: env[name] for name in
                             ("OMP_NUM_THREADS", "OMP_DYNAMIC", "OMP_STACKSIZE", "MALLOC_ARENA_MAX")}
    failed = False
    for scheme in args.schemes:
        log = args.output / f"{scheme}.log"
        resources = args.output / f"{scheme}.resources.txt"
        command = ["/usr/bin/time", "-v", "-o", str(resources), "prlimit",
                   f"--as={report['address_space_cap_bytes']}", "--", str(executable), scheme]
        print(f"Running {scheme}; log: {log}", flush=True)
        with log.open("w") as stream:
            completed = subprocess.run(command, env=env, stdout=stream, stderr=subprocess.STDOUT)
        output = log.read_text()
        fields = {name: float(value) for name, value in
                  re.findall(r"^([a-z_0-9]+)=([-+0-9.eE]+)$", output, re.MULTILINE)}
        rss = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", resources.read_text())
        passed = completed.returncode == 0 and "Products=256 Passed" in output
        report["runs"][scheme] = dict(fields, exit_status=completed.returncode, passed=passed,
                                     peak_rss_kib=int(rss[1]) if rss else None,
                                     command=command, log=log.name, resources=resources.name)
        report["complete"] = len(report["runs"]) == len(args.schemes)
        report["all_passed"] = report["complete"] and all(run["passed"] for run in report["runs"].values())
        (args.output / "comparison.json").write_text(json.dumps(report, indent=2) + "\n")
        failed |= not passed
        print(f"{scheme}: {'PASS' if passed else 'FAIL'}", flush=True)
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())
