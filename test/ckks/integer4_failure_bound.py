#!/usr/bin/env python3
"""Conservative coefficient-norm bound for the exact integer-only circuit.

The probability statement assumes independent ideal normal sampler outputs
rounded to integers as in TFHEpp. The conditional error bound itself uses
triangle inequalities, so it does not assume independent ciphertext errors
and remains valid when evaluation keys or input ciphertexts are reused.
"""

import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--multiplications", type=int, default=256)
    parser.add_argument("--input-bits", type=int, default=4)
    parser.add_argument("--cutoff-sigma", type=int, default=12)
    args = parser.parse_args()
    if args.multiplications <= 0:
        parser.error("multiplications must be positive")
    if not 1 <= args.input_bits <= 16 or args.cutoff_sigma <= 0:
        parser.error("input-bits must be 1..16 and cutoff-sigma positive")
    n, N, M = 2048, 16384, 2**(2*args.input_bits)
    # Keep extraction q=2^48, scale=2^40; output phase is M*u*v.
    postscale = 256*M
    delta = 2**64
    rows, digit_max, output_levels = 16, 128, 24
    # A one-unit allowance covers truncation or nearest-integer rounding.
    cutoff = args.cutoff_sigma
    ckks_error = math.ceil(cutoff * F(16, 5)) + 1
    tfhe_error = cutoff * 2**13 + 1
    eta = F(tfhe_error, 2**64)
    lift_round = F(1 + N, 2)
    lift_error = (n * ckks_error + lift_round) / delta
    # Actual unwrapped b-a*s has absolute value at most n+1 for ternary s.
    U = n + 1
    input_error = 2*U*eta + eta**2
    lift_product_error = 2*U*lift_error + N*lift_error**2
    # Three rescaled tensor components have weights 1, s, and s^2.
    # (1+N)^2/2 safely overbounds (1+N+N^2)/2.
    tensor_round = F((1 + N)**2, 2)
    relin_error = N * rows * digit_max * ckks_error
    arithmetic_error = (tensor_round + relin_error) / delta
    final_rescale_error = F(M, 256) * lift_round / 2**40
    direct_switch_error = F(M*N*output_levels*tfhe_error, 2**64)
    terms = {
        "input_tfhe_noise": M**2 * input_error,
        "ckks_phase_conversion": M**2 * lift_product_error,
        "tensor_rounding_and_relinearization": M**2 * arithmetic_error,
        "final_rescale_rounding": final_rescale_error,
        "direct_output_switching": direct_switch_error,
    }
    total = sum(terms.values(), F(0))
    draws = {
        "phase_key": n*N,
        "relinearization_key": rows*N,
        "output_switching_key": N*output_levels*2,
        "fresh_tfhe_inputs": 2*args.multiplications,
    }
    # P(|Z|>t) <= 2 exp(-t^2/2), then a union bound over all sampled errors.
    log2_probability = math.log2(sum(draws.values())) + 1 - cutoff**2/(2*math.log(2))
    # Product and postscale stay inside the centered active modulus before
    # the final intentional reduction. Check these using exact fractions.
    message_bound = (U + lift_error)**2 + arithmetic_error
    no_unintended_wrap = message_bound * delta * postscale < 2**127
    report = {
        "scope": "One key generation and at most the configured number of fresh-input "
                 f"multiplications. Unsigned {args.input_bits}-bit inputs encoded as x/{M}; "
                 f"exact {2*args.input_bits}-bit product. Widths other than 4 are analytical candidates.",
        "sampler_model": "Ideal normal samples with integer rounding allowance; "
                         "not a formal verification of std::normal_distribution or the PRNG.",
        "arithmetic": "Exact modular integer/NTT evaluation; no floating CKKS encoding or FFT.",
        "parameters": dict(ckks_n=N, tfhe_n=n, plaintext_modulus=M, log_delta=64,
                           ckks_input_log_q=256, phase_log_q=192, product_log_q=128,
                           extraction_log_q=48, extraction_log_delta=40),
        "input_bits": args.input_bits,
        "product_bits": 2*args.input_bits,
        "public_postscale": postscale,
        "multiplications": args.multiplications,
        "normal_cutoff_sigma": cutoff,
        "ckks_integer_noise_cutoff": ckks_error,
        "tfhe_integer_noise_cutoff": tfhe_error,
        "gaussian_draw_counts": draws,
        "integer_unit_error_terms": {k: float(v) for k,v in terms.items()},
        "integer_unit_error_bound": float(total),
        "exact_error_bound_fraction": str(total),
        "rounding_radius": 0.5,
        "no_unintended_wrap": no_unintended_wrap,
        "log2_failure_probability_bound": log2_probability,
        "target_log2_failure_probability": -40,
        "passes": total < F(1,2) and no_unintended_wrap and log2_probability <= -40,
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if report["passes"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
