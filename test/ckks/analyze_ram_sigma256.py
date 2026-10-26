#!/usr/bin/env python3
"""Validate the complete opt-in sigma256 retained-fixture qualification run."""
import argparse
import json
import math
from pathlib import Path
import statistics


def finite(record, names):
    if any(not math.isfinite(record[n]) or record[n]<0 for n in names):
        raise ValueError('nonfinite or negative measurement')


def analyze(root):
    results=root/'results'
    capture=json.loads((results/'capture/capture.json').read_text())
    small=json.loads((results/'small-ring/comparison.json').read_text())
    summary=json.loads((results/'paired8192/crossover-8192x1.json').read_text())
    records=[json.loads(line) for line in (results/'paired8192/crossover-8192x1.jsonl').read_text().splitlines()]
    job=json.loads((results/'pipeline/job.json').read_text())
    integrity=json.loads((results/'integrity-sigma256.json').read_text())
    fixture=json.loads((root/'fixture-check.json').read_text())
    gates=json.loads((root/'synthetic-gates.json').read_text())
    if job.get('exit_code')!=0 or integrity['changed']:
        raise ValueError('failed job or fixture mutation')
    if not gates['correct'] or gates['lvl1_sigma']!=256 or gates['nand_checks']!=64 or gates['circuit_bootstrap_bit_checks']!=8192 or gates['ram_accesses']!=16 or not gates['feedback_checked']:
        raise ValueError('incomplete synthetic qualification')
    if not fixture['correct'] or fixture['lvl1_sigma']!=256 or not fixture['all_five_traces_checked']:
        raise ValueError('incomplete fixture qualification')
    if tuple(fixture[n] for n in ('lvl1_rows','lvl1_basebit','pks_rows','pks_basebit'))!=(4,5,12,2):
        raise ValueError('wrong decomposition profile')
    finite(fixture,('predicted_pre_switch_max_error',))
    if not fixture['large_noise_checked'] or fixture['predicted_pre_switch_max_error']>1/16:
        raise ValueError('large input noise check incomplete')
    for n in ('memory_noise_rms','private_switch_noise_rms'):
        finite(fixture,(n,))
        if not .95*256<fixture[n]<1.05*256:raise ValueError('fixture noise mismatch')
    if not capture['correct'] or capture['security_validated'] is not False:
        raise ValueError('capture failed')
    for report in (capture,small,summary):
        if report['output_bridge']['mode']!='ring' or report['output_bridge']['sigma']!=256:
            raise ValueError('wrong bridge profile')
        if report['security_validated'] is not False:raise ValueError('unexpected security claim')
    finite(capture,('pre_switch_error',))
    finite(small,('pre_switch_max_error',))
    if capture['pre_switch_error']>1/16 or small['pre_switch_max_error']>1/16:
        raise ValueError('input phase budget exceeded')
    if small['operations']!=16 or not small['read_before_write'] or len(small['backends'])!=2:
        raise ValueError('small RAM trace incomplete')
    if {b['name'] for b in small['backends']}!={'tfhe-cell','ckks-batch'}:
        raise ValueError('wrong small RAM backends')
    if any(not b['correct'] or not math.isfinite(b['state_max_phase_error']) or not 0<=b['state_max_phase_error']<1/8 for b in small['backends']):
        raise ValueError('small RAM state failed')
    if len(records)!=10 or summary['lvl1_sigma']!=256 or summary['address_bits']!=13 or summary['word_bits']!=1 or not summary['all_cells_checked'] or not summary['alternating_order'] or summary['warmups_per_backend']!=1 or summary['measured_accesses_per_backend']!=4:
        raise ValueError('paired full RAM trace incomplete')
    names=('tfhe-serial','ckks')
    for position,record in enumerate(records):
        step,order=divmod(position,2)
        expected=names[order^(step&1)]
        if (record['backend'],record['step'],record['order'],record['bits'],record['warmup'],record['correct'])!=(expected,step,order,8192,step==0,True):
            raise ValueError('paired order or correctness mismatch')
        finite(record,('seconds','refresh_seconds','max_phase_error','pre_switch_error'))
        if record['max_phase_error']>=1/8 or record['pre_switch_error']>1/16:
            raise ValueError('full RAM phase budget exceeded')
    by_name={b['name']:b for b in summary['backends']}
    if set(by_name)!=set(names) or len(summary['backends'])!=2:raise ValueError('wrong paired backends')
    for name in names:
        measured=[r for r in records if r['backend']==name and not r['warmup']]
        if [r['seconds'] for r in measured]!=by_name[name]['access_seconds'] or [r['refresh_seconds'] for r in measured]!=by_name[name]['refresh_seconds']:
            raise ValueError('paired journal/summary mismatch')
        if by_name[name]['max_phase_error']!=max(r['max_phase_error'] for r in records if r['backend']==name):
            raise ValueError('paired phase summary mismatch')
    if summary['pre_switch_max_error']!=max(r['pre_switch_error'] for r in records):
        raise ValueError('paired input summary mismatch')
    values={name:by_name[name]['access_seconds'] for name in names}
    means={name:statistics.mean(v) for name,v in values.items()}
    finite(summary['measured_ckks_stages'],('packing','bootstrap','decoding','switching'))
    return dict(all_checks_passed=True,security_validated=False,circular_security_assumed=True,
        diagnostic_owner_reused=True,original_sigma128_material_not_rehabilitated=True,
        mean_access_seconds=means,observed_ring_vs_serial_tfhe_speedup=means['tfhe-serial']/means['ckks'],
        ring_faster_in_each_measured_pair=all(a<b for a,b in zip(values['ckks'],values['tfhe-serial'])),
        mean_measured_ckks_stages={n:v/4 for n,v in summary['measured_ckks_stages'].items()},
        full_accesses_checked=10,full_cell_bits_checked=81920,
        uncontended_timing_qualified=job['all_intervals_below_0_2_external_cores'],
        capture=capture,small_feedback=small,paired_summary=summary,paired_records=records,
        synthetic_gates=gates,fixture_noise_audit=fixture,job=job,integrity=integrity)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    report=analyze(args.directory)
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(report['mean_access_seconds'])


if __name__=='__main__':main()
