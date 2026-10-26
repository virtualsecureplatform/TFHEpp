#!/usr/bin/env python3
"""Validate complete ring-switch experiment receipts; report measured wins only."""
import argparse
import json
import math
import re
from pathlib import Path
import statistics

NAMES = ('tfhe-serial', 'direct', 'direct256', 'ring')
PERMUTATIONS = ((0, 1, 2), (2, 1, 0), (1, 0, 2), (2, 0, 1), (0, 2, 1), (1, 2, 0))
STAGES = ('packing', 'bootstrap', 'decoding', 'switching')


def finite(record, fields):
    if not all(math.isfinite(record[k]) and record[k] >= 0 for k in fields):
        raise ValueError('nonfinite or negative measurement')


def bridge_profile(record, expected):
    if (record['mode'] != expected or record['sigma'] != (128 if expected == 'direct' else 256)
            or not re.fullmatch('[0-9a-f]{64}', record['key_fingerprint'])
            or record['key_bytes'] != (2097152 if expected == 'ring' else 2149580800)
            or record['prepared_bytes'] != (4194304 if expected == 'ring' else 0)):
        raise ValueError('bridge profile mismatch')
    finite(record, ('setup_seconds', 'cumulative_kernel_seconds', 'cumulative_split_seconds',
                    'cumulative_extraction_seconds'))


def comparison(directory, bits):
    prefix = directory / f'ring-comparison-{bits}x1'
    records = [json.loads(l) for l in prefix.with_suffix('.jsonl').read_text().splitlines()]
    summary = json.loads(prefix.with_suffix('.json').read_text())
    if (len(records) != 20 or summary['bits'] != bits or not summary['all_cells_checked']
            or summary['warmups_per_backend'] != 1 or summary['measured_accesses_per_backend'] != 4
            or summary['security_validated'] is not False):
        raise ValueError('incomplete RAM comparison')
    for pos, record in enumerate(records):
        step, order = divmod(pos, 4)
        expected = NAMES[(order + step) % 4]
        if (record['bits'], record['backend'], record['step'], record['order'], record['warmup'], record['correct'], record['all_cells_checked']) != (bits, expected, step, order, step == 0, True, True):
            raise ValueError('RAM order/correctness mismatch')
        finite(record, ('seconds', 'refresh_seconds', 'max_phase_error', 'pre_switch_error'))
        if record['pre_switch_error'] > 1/16 or record['max_phase_error'] >= 1/8:
            raise ValueError('RAM phase budget exceeded')
        if expected != 'tfhe-serial':
            if set(record['stages']) != set(STAGES):
                raise ValueError('missing RAM stage')
            finite(record['stages'], STAGES)
    result = {'bits': bits, 'operations_checked': 20, 'cell_bits_checked': 20*bits, 'backends': {}}
    for backend in summary['backends']:
        name = backend['name']
        if name not in NAMES or name in result['backends']:
            raise ValueError('unexpected/duplicate backend')
        measured = [r for r in records if r['backend'] == name and not r['warmup']]
        values = [r['seconds'] for r in measured]
        if values != backend['access_seconds']:
            raise ValueError('RAM journal/summary mismatch')
        result['backends'][name] = {'mean_access_seconds': statistics.mean(values), 'access_seconds': values,
                                    'maximum_phase_error': max(r['max_phase_error'] for r in records if r['backend'] == name)}
        if name != 'tfhe-serial':
            bridge_profile(backend['bridge'], name)
            stages = backend['measured_stage_totals']
            if set(stages) != set(STAGES):
                raise ValueError('missing RAM stage totals')
            finite(stages, STAGES)
            for stage, value in stages.items():
                actual = sum(r['stages'][stage] for r in measured)
                if not math.isclose(value, actual, rel_tol=1e-9, abs_tol=1e-7):
                    raise ValueError('RAM stage sum mismatch')
            result['backends'][name]['mean_stages_seconds'] = {k: v/4 for k, v in stages.items()}
    if set(result['backends']) != set(NAMES):
        raise ValueError('missing RAM backend')
    ring = result['backends']['ring']
    result['ring_comparisons'] = {}
    for name in NAMES[:-1]:
        other = result['backends'][name]
        result['ring_comparisons'][name] = {'mean_speedup': other['mean_access_seconds']/ring['mean_access_seconds'],
            'ring_faster_in_every_measured_pair': all(r < d for r, d in zip(ring['access_seconds'], other['access_seconds']))}
    return result


def bridge(directory):
    records = [json.loads(l) for l in (directory / 'bridge.jsonl').read_text().splitlines()]
    summary = json.loads((directory / 'bridge.json').read_text())
    finite(summary, ('pre_switch_error',))
    if (len(records) != 84 or not summary['all_checks_passed'] or summary['pre_switch_error'] > 1/16
            or summary['security_validated'] is not False):
        raise ValueError('incomplete bridge benchmark')
    if len(summary['backends']) != 3:
        raise ValueError('missing bridge backend')
    for record, expected in zip(summary['backends'], NAMES[1:]):
        bridge_profile(record, expected)
    result = []
    for group, bits in enumerate((256, 512, 1024, 8192)):
        batch = records[group*21:(group+1)*21]
        for pos, record in enumerate(batch):
            trial = pos//3-1
            order = pos % 3
            backend = NAMES[1:][order if trial == -1 else PERMUTATIONS[trial][order]]
            if (record['bits'], record['trial'], record['order'], record['backend'], record['warmup'], record['correct']) != (bits, trial, order, backend, trial == -1, True):
                raise ValueError('bridge order/correctness mismatch')
            finite(record, ('seconds', 'maximum_added_phase_error', 'kernel_seconds', 'split_seconds', 'extraction_seconds'))
            if backend == 'ring' and record['maximum_added_phase_error'] > .025746934137548232:
                raise ValueError('bridge error budget exceeded')
        timings = {name: [r['seconds'] for r in batch if r['backend'] == name and not r['warmup']] for name in NAMES[1:]}
        means = {name: statistics.mean(values) for name, values in timings.items()}
        result.append({'bits': bits, 'mean_seconds': means,
                       'matched_noise_speedup': means['direct256']/means['ring'],
                       'ring_faster_in_every_matched_trial': all(r < d for r, d in zip(timings['ring'], timings['direct256']))})
    return result, summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory', type=Path)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    ram = [comparison(args.directory / 'comparison', bits) for bits in (256, 512, 1024, 8192)]
    micro, metadata = bridge(args.directory / 'bridge')
    capture = json.loads((args.directory / 'capture/capture.json').read_text())
    finite(capture, ('pre_switch_error',))
    if (not capture['correct'] or not capture['input_reused'] or capture['fresh_private_keys_generated']
            or capture['pre_switch_error'] > 1/16 or capture['output_bridge']['mode'] != 'ring'):
        raise ValueError('invalid full-batch capture')
    bridge_profile(capture['output_bridge'], 'ring')
    small = {}
    for name in NAMES[1:]:
        report = json.loads((args.directory / f'small-{name}/comparison.json').read_text())
        if report['operations'] != 16 or report['output_bridge']['mode'] != name or report['pre_switch_max_error'] > 1/16 or not all(b['correct'] and b['state_max_phase_error'] < 1/8 for b in report['backends']):
            raise ValueError('incomplete small RAM/feedback check')
        bridge_profile(report['output_bridge'], name)
        small[name] = report
    jobs = [json.loads(p.read_text()) for p in args.directory.glob('*/job.json')]
    names = {p.parent.name for p in args.directory.glob('*/job.json')}
    expected = {'pipeline', 'bridge'} if 'pipeline' in names else {'capture', 'small-direct', 'small-direct256', 'small-ring', 'comparison', 'bridge'}
    if names != expected or any(j.get('exit_code') != 0 for j in jobs):
        raise ValueError('incomplete/failed experiment job')
    integrity = [json.loads(p.read_text()) for p in args.directory.glob('integrity-*.json')]
    if not integrity or any(r['changed'] for r in integrity):
        raise ValueError('missing or failed retained-file integrity check')
    result = {'all_checks_passed': True, 'security_validated': False,
              'scope': 'Experimental ring bridge, independent persistent RAM states, same retained requests. Existing TFHE noise/KDM and whole-run failure qualification remain open.',
              'uncontended_timing_qualified': all(j['all_intervals_below_0_2_external_cores'] for j in jobs),
              'ram_profiles': ram, 'bridge_profiles': micro, 'bridge_metadata': metadata, 'capture': capture,
              'small_ram_feedback': small, 'jobs': jobs, 'integrity': integrity}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    for profile in ram:
        print(profile['bits'], 'ring vs matched direct:', profile['ring_comparisons']['direct256'])


if __name__ == '__main__':
    main()
