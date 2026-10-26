#!/usr/bin/env python3
"""Assess the default RAM's TFHE noise families with Parameter-Selection.

Circular security is an explicit accepted assumption for this assessment.
Run under Sage. No keys or ciphertext payloads are read or generated here.
"""
import argparse
import concurrent.futures
import datetime
import hashlib
import math
from pathlib import Path
import sys
import time

import json

ORIGINAL = ('primal_usvp', 'primal_bdd', 'dual', 'dual_hybrid')
EXTENDED = ORIGINAL + ('mitm', 'bdd_hybrid', 'bdd_mitm_hybrid')


def integer_rms(sigma):
    return math.sqrt(math.fsum(k*k*(math.erfc(k/(sigma*math.sqrt(2))) -
                       math.erfc((k+1)/(sigma*math.sqrt(2))))
                       for k in range(1, math.ceil(12*sigma)+1)))


def cases():
    result = {}
    for sigma in (128, 132, 136, 144, 160, 192, 256):
        result[f'lvl1_sigma{sigma}'] = dict(n=1024, log_q=32, sigma=sigma,
             secret='ternary', m=None, attacks=EXTENDED)
    # One coefficient per independent raw PKS row, with known zero message.
    # (2048+1)*8*2^(3-1) = 65568; ring rotations are not counted as independent.
    for sigma in (128, 132, 192, 256):
        result[f'lvl1_sigma{sigma}_pks_finite'] = dict(n=1024, log_q=32, sigma=sigma,
             secret='ternary', m=65568, attacks=ORIGINAL)
    result['lvl1_sigma128_bk_finite'] = dict(n=1024, log_q=32, sigma=128,
             secret='ternary', m=630*3, attacks=ORIGINAL)
    result['lvl0_current'] = dict(n=630, log_q=16, sigma=65536*.0000925119974676756,
             secret='binary', m=None, attacks=EXTENDED)
    result['lvl2_current'] = dict(n=2048, log_q=64, sigma=8192,
             secret='ternary', m=None, attacks=EXTENDED)
    for value in result.values():
        value['integer_rms'] = integer_rms(value['sigma'])
    return result


def estimate(task):
    selection, name, profile, model, attack = task
    sys.path[:0] = [str(selection/'python/lattice-estimator'), str(selection/'python/estimates')]
    from estimator import LWE, ND, RC
    from narrow_ckks_security import configure_max_beta
    from sage.all import log, oo
    configure_max_beta(4096)
    p = LWE.Parameters(n=profile['n'], q=2**profile['log_q'],
            Xs=ND.UniformMod(3) if profile['secret']=='ternary' else ND.Binary,
            Xe=ND.DiscreteGaussian(profile['integer_rms']),
            m=oo if profile['m'] is None else profile['m'])
    before = time.monotonic()
    try:
        if attack=='mitm':
            cost = LWE.mitm(p)
        elif attack in ('bdd_hybrid', 'bdd_mitm_hybrid'):
            mitm = attack=='bdd_mitm_hybrid'
            cost = LWE.primal_hybrid(p, mitm=mitm, babai=mitm, red_cost_model=getattr(RC,model))
        else:
            cost = getattr(LWE,attack)(p, red_cost_model=getattr(RC,model))
        bits = float(log(cost['rop'],2))
        if not math.isfinite(bits):
            raise ValueError('nonfinite optimizer result does not establish security')
        output = dict(bits=bits, cost=str(cost))
    except Exception as e:
        output = dict(error=type(e).__name__+': '+str(e))
    return name, model, attack, dict(output, seconds=time.monotonic()-before)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--selection', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--jobs', type=int, default=4)
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error('--jobs must be positive')
    profiles = cases()
    report = dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        security_validated=False, circular_security_assumed=True,
        target_classical_bits=128, max_beta=4096,
        scope='Classical generic LWE cost estimates for variance-matched discrete-Gaussian proxies. Circular security accepted; no quantum or exact-sampler/ring-specific proof.',
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        cases={name:dict(profile, models={}) for name,profile in profiles.items()})
    tasks = [(args.selection,name,profile,model,attack)
             for name,profile in profiles.items() for model in ('BDGL16','MATZOV')
             for attack in profile['attacks']]
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(estimate,task) for task in tasks]
        for future in concurrent.futures.as_completed(futures):
            name,model,attack,result = future.result()
            report['cases'][name]['models'].setdefault(model,{})[attack] = result
            args.output.write_text(json.dumps(report,indent=2)+'\n')
            print(name,model,attack,result.get('bits',result.get('error')),flush=True)
    for profile in report['cases'].values():
        profile['summary'] = {}
        for model in ('BDGL16','MATZOV'):
            results = profile['models'].get(model,{})
            complete = all('bits' in results.get(a,{}) for a in profile['attacks'])
            weakest = min((v['bits'] for v in results.values() if 'bits' in v),default=None)
            profile['summary'][model] = dict(complete=complete, weakest_bits=weakest,
                  passes_128=complete and weakest>=128)
    args.output.write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    main()
