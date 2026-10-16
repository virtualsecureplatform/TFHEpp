#!/usr/bin/env python3
"""Validate paired RAM receipts and report observed crossover, without fitting it."""
import argparse
import json
import math
import statistics
from pathlib import Path


def analyze(directory, bits):
    prefix = directory / f"crossover-{bits}x1"
    summary = json.loads(prefix.with_suffix(".json").read_text())
    records = [json.loads(line) for line in prefix.with_suffix(".jsonl").read_text().splitlines()]
    if len(records) != 10:
        raise ValueError(f"{bits}: incomplete trace")
    for position, record in enumerate(records):
        step, order = divmod(position, 2)
        backend = "ckks" if order ^ (step & 1) else "tfhe-serial"
        if (record["step"], record["order"], record["backend"], record["warmup"], record["bits"], record["correct"]) != (step, order, backend, step == 0, bits, True):
            raise ValueError(f"{bits}: inconsistent paired receipt at {position}")
        if not all(math.isfinite(record[k]) and record[k] >= 0 for k in ("seconds", "refresh_seconds", "max_phase_error")):
            raise ValueError(f"{bits}: invalid measurement")
    output = {"bits": bits, "warmups_per_backend": 1, "measured_accesses_per_backend": 4,
              "complete_operations_checked": 10, "cell_bits_checked": 10 * bits,
              "backends": {}}
    paired = {}
    for backend in summary["backends"]:
        name = backend["name"]
        actual = [r for r in records if r["backend"] == name and not r["warmup"]]
        for field in ("access_seconds", "refresh_seconds"):
            values = [r["seconds" if field == "access_seconds" else field] for r in actual]
            if values != backend[field]:
                raise ValueError(f"{bits}: summary and journal differ")
        paired[name] = backend["access_seconds"]
        output["backends"][name] = {"mean_access_seconds": statistics.mean(paired[name]),
                                    "min_access_seconds": min(paired[name]),
                                    "max_access_seconds": max(paired[name]),
                                    "mean_refresh_seconds": statistics.mean(backend["refresh_seconds"]),
                                    "max_phase_error": backend["max_phase_error"]}
    differences = [c - t for c, t in zip(paired["ckks"], paired["tfhe-serial"])]
    output["paired_ckks_minus_tfhe_seconds"] = differences
    output["ckks_faster_in_every_measured_pair"] = all(d < 0 for d in differences)
    output["tfhe_faster_in_every_measured_pair"] = all(d > 0 for d in differences)
    output["ckks_to_tfhe_mean_ratio"] = statistics.mean(paired["ckks"]) / statistics.mean(paired["tfhe-serial"])
    output["ckks_faster_on_mean"] = output["ckks_to_tfhe_mean_ratio"] < 1
    output["mean_ckks_stages_seconds"] = {k: v / 4 for k, v in summary["measured_ckks_stages"].items()}
    output["raw_summary"] = summary
    output["raw_records"] = records
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--bits", type=int, nargs="+", default=[256, 512, 1024, 8192])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    profiles = [analyze(args.directory, bits) for bits in sorted(set(args.bits))]
    brackets = [[left["bits"], right["bits"]] for left, right in zip(profiles, profiles[1:])
                if left["tfhe_faster_in_every_measured_pair"] and right["ckks_faster_in_every_measured_pair"]]
    result = {"scope": "Empirical end-to-end RAM crossover against the existing serial TFHE refresh; not a universal proof or comparison with parallel TFHE.",
              "security_validated": False, "profiles": profiles,
              "observed_all_pairs_crossover_brackets_bits": brackets,
              "consistent_crossover_observed": any(p["ckks_faster_in_every_measured_pair"] for p in profiles),
              "mean_crossover_observed": any(p["ckks_faster_on_mean"] for p in profiles),
              "absence_of_crossover_proved": False}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    for p in profiles:
        print(f"{p['bits']}: TFHE {p['backends']['tfhe-serial']['mean_access_seconds']:.6f}s, "
              f"CKKS {p['backends']['ckks']['mean_access_seconds']:.6f}s, ratio {p['ckks_to_tfhe_mean_ratio']:.4f}")


if __name__ == "__main__":
    main()
