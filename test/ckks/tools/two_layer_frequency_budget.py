"""Exact integer frequency coverage only; this does not certify CKKS precision.

Run: python3 test/ckks/tools/two_layer_frequency_budget.py > plans.json
All seeds and intermediate frequencies are positive integers at most40.
"""
import itertools
import json


def signed(frequencies):
    return {0, *frequencies, *(-f for f in frequencies)}


def sums(frequencies):
    return {a + b for a in frequencies for b in frequencies}


def search():
    best_prefix = 0
    best_seeds = None
    candidates = []
    for seeds in itertools.combinations(range(1, 41), 3):
        level1 = sums(signed(seeds))
        level2 = sums(level1)
        prefix = 0
        while prefix + 1 in level2:
            prefix += 1
        if prefix > best_prefix:
            best_prefix, best_seeds = prefix, seeds
        if prefix >= 40:
            candidates.append((seeds, sorted(f for f in level1 if 0 < f <= 40)))
    plans = []
    for seeds, available in candidates:
        # n positive roots produce at most n+n^2 positive direct/sum/difference
        # frequencies, so n<=5 cannot cover all40 required frequencies.
        for count in range(6, len(available) + 1):
            hits = []
            for roots in itertools.combinations(available, count):
                covered = sums(signed(roots))
                if all(k in covered for k in range(1, 41)):
                    hits.append((len(set(roots) - set(seeds)), max(roots), roots))
            if hits:
                products, _, roots = min(hits)
                plans.append(dict(seeds=seeds, intermediate_count=count,
                                  intermediate_products=products, intermediate_roots=roots))
                break
    plans.sort(key=lambda p: (p['intermediate_count'], p['intermediate_products'], max(p['seeds'])))
    assert best_prefix == 46 and len(candidates) == 37
    assert plans[0]['seeds'] == (6, 9, 10)
    assert plans[0]['intermediate_roots'] == (3, 6, 9, 12, 15, 18, 19, 20)
    return dict(seed_range=[1, 40], seed_sets_checked=9880,
                maximum_contiguous_degree=best_prefix, maximum_degree_seeds=best_seeds,
                degree40_candidate_seed_sets=len(candidates), selected=plans[0], plans=plans,
                scope='Frequency identities only; retained encrypted37-bit profile fails the23-bit gate.')


if __name__ == '__main__':
    print(json.dumps(search(), indent=2))
