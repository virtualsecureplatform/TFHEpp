#!/usr/bin/env python3
"""Reproduce subring algebra/mask checks and optional Parameter-Selection screens.

No private keys, retained ciphertexts, or evaluation-key payloads are read.
The mask test models the removed historical expander; it demonstrates its PRG
failure, not TFHE-secret recovery or a failure in the replacement BLAKE3 stream.
"""
import argparse
import concurrent.futures
import hashlib
import json
import math
from pathlib import Path
import random
import sys
import time

MASK = (1 << 64) - 1
STEP = 0x9e3779b97f4a7c15


def mix(x):
    x = ((x ^ (x >> 30)) * 0xbf58476d1ce4e5b9) & MASK
    x = ((x ^ (x >> 27)) * 0x94d049bb133111eb) & MASK
    return x ^ (x >> 31)


def unshift(x, shift):
    out = 0
    for offset in range(0, 64, shift):
        out ^= x >> offset
    return out


def unmix(x):
    x = unshift(x, 31)
    x = (x * pow(0x94d049bb133111eb, -1, 1 << 64)) & MASK
    x = unshift(x, 27)
    x = (x * pow(0xbf58476d1ce4e5b9, -1, 1 << 64)) & MASK
    return unshift(x, 30)


def fold(seed):
    state = 0x6a09e667f3bcc909
    for word in seed:
        state = mix(state ^ word)
    return state


def mul(a, b, q):
    out = [0] * len(a)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[(i+j) % len(a)] += x*y*(1 if i+j < len(a) else -1)
    return [v % q for v in out]


def algebra(repo):
    source = (repo/'include/ckks/ckks.hpp').read_text()
    rng = random.Random(102625)
    for _ in range(1000):
        x = rng.getrandbits(64)
        assert unmix(mix(x)) == x
    seed = [rng.getrandbits(64) for _ in range(4)]
    alternative = [rng.getrandbits(64) for _ in range(3)]
    state = 0x6a09e667f3bcc909
    for word in alternative:
        state = mix(state ^ word)
    alternative.append(state ^ unmix(fold(seed)))
    assert seed != alternative and fold(seed) == fold(alternative)
    # Enumerate only 2^12 high-bit suffixes of one synthetic q32 output.
    # A complete q32 distinguisher enumerates all 2^32 possible upper halves.
    high = rng.randrange(1 << 12)
    first = (high << 32) | rng.getrandbits(32)
    state = unmix(first)
    observed = [mix((state+i*STEP) & MASK) & 0xffffffff for i in range(4)]
    candidates = []
    for guess in range(1 << 12):
        recovered = unmix((guess << 32) | observed[0])
        if all((mix((recovered+i*STEP) & MASK) & 0xffffffff) == observed[i] for i in range(1, 4)):
            candidates.append(recovered)
    assert candidates == [state]
    # Independent integer check: deinterleaving commutes with multiplication
    # by an embedded small-ring secret. No primal/dual error scaling occurs.
    cases = 0
    for N, n in ((64, 8), (128, 16), (256, 16)):
        stride = N//n
        for _ in range(10):
            a = [rng.getrandbits(32) for _ in range(N)]
            t = [rng.randrange(-1, 2) for _ in range(n)]
            embedded = [t[i//stride] if i % stride == 0 else 0 for i in range(N)]
            product = mul(a, embedded, 1 << 32)
            for r in range(stride):
                assert product[r::stride] == mul(a[r::stride], t, 1 << 32)
                cases += 1
    sigma = 256
    rms = integer_rms(sigma)
    support = range(-12*sigma, 12*sigma+1)
    weights = [math.exp(-k*k/(2*rms*rms)) for k in support]
    normalizer = math.fsum(weights)
    def truncated_pmf(k):
        if k == 0:
            return math.erf(1/(sigma*math.sqrt(2)))
        k = abs(k)
        return .5*(math.erfc(k/(sigma*math.sqrt(2)))-math.erfc((k+1)/(sigma*math.sqrt(2))))
    tv = .5*math.fsum(abs(truncated_pmf(k)-w/normalizer) for k, w in zip(support, weights))
    return {'security_validated': False, 'source_sha256': hashlib.sha256(source.encode()).hexdigest(),
            'split_multiplication_checks': cases, 'distinct_seed_collision_verified': True,
            'mask_state_bits': 64, 'q32_mask_distinguisher_upper_bound_candidates': 1 << 32,
            'q32_restricted_suffix_demo_candidates': 1 << 12,
            'q32_restricted_suffix_demo_unique_recovery': True,
            'mask_attack_scope': 'Historical SplitMix model only; current production expansion is BLAKE3. Demonstrates a distribution failure, not encryption-secret recovery.',
            'sigma256_integer_rms': rms, 'single_coefficient_total_variation_from_variance_matched_discrete_gaussian': tv,
            'noise_scope': 'Ideal continuous-normal truncation PMF, numerical quadrature on +/-12 sigma; not a finite-precision sampler proof.'}


def integer_rms(sigma):
    return math.sqrt(math.fsum(k*k*(math.erfc(k/(sigma*math.sqrt(2)))-math.erfc((k+1)/(sigma*math.sqrt(2)))) for k in range(1, math.ceil(12*sigma)+1)))


def estimate(task):
    selection, case, sigma, samples, model, attack = task
    sys.path[:0] = [str(selection/'python/lattice-estimator'), str(selection/'python/estimates')]
    from estimator import LWE, ND, RC
    from narrow_ckks_security import configure_max_beta
    from sage.all import log, oo
    configure_max_beta(4096)
    p = LWE.Parameters(n=1024, q=1 << 32, Xs=ND.UniformMod(3),
                       Xe=ND.DiscreteGaussian(integer_rms(sigma)), m=oo if samples is None else samples)
    started = time.monotonic()
    try:
        if attack == 'mitm':
            result = LWE.mitm(p)
        elif attack in ('bdd_hybrid', 'bdd_mitm_hybrid'):
            mitm = attack == 'bdd_mitm_hybrid'
            result = LWE.primal_hybrid(p, mitm=mitm, babai=mitm, red_cost_model=getattr(RC, model))
        else:
            result = getattr(LWE, attack)(p, red_cost_model=getattr(RC, model))
        bits = float(log(result['rop'], 2))
        if not math.isfinite(bits):
            raise ValueError('nonfinite result does not establish security')
        output = {'bits': bits, 'cost': str(result)}
    except Exception as e:
        output = {'error': type(e).__name__+': '+str(e)}
    return case, model, attack, dict(output, seconds=time.monotonic()-started)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--mode', choices=('algebra', 'estimator'), default='algebra')
    p.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[3])
    p.add_argument('--selection', type=Path)
    args = p.parse_args()
    if args.mode == 'algebra':
        report = algebra(args.repo)
    else:
        if args.selection is None:
            p.error('--selection is required for estimator mode')
        original = ('primal_usvp', 'primal_bdd', 'dual', 'dual_hybrid')
        cases = {'sigma256_unbounded': (256, None, original+('mitm', 'bdd_hybrid', 'bdd_mitm_hybrid')),
                 'sigma128_unbounded': (128, None, original+('mitm', 'bdd_hybrid', 'bdd_mitm_hybrid')),
                 'sigma256_finite': (256, 262144, original)}
        tasks = [(args.selection, name, sigma, m, model, attack) for name, (sigma, m, attacks) in cases.items()
                 for model in ('BDGL16', 'MATZOV') for attack in attacks]
        report = {'security_validated': False, 'scope': 'Classical generic LWE proxies, conditional on uniform masks. Does not qualify existing SplitMix masks, ring-specific attacks, or cyclic key composition.',
                  'cases': {name: {'n': 1024, 'log_q': 32, 'sigma': s, 'integer_rms': integer_rms(s), 'samples': m, 'models': {}} for name, (s, m, _) in cases.items()}}
        with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
            for case, model, attack, result in pool.map(estimate, tasks):
                report['cases'][case]['models'].setdefault(model, {})[attack] = result
                args.output.write_text(json.dumps(report, indent=2)+'\n')
                print(case, model, attack, result.get('bits', result.get('error')), flush=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
