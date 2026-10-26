#!/usr/bin/env python3
"""Run retained-fixture ring bridge experiments and record load/provenance.

Requires --add-ring-switch-keys setup first. Never generates keys or inputs;
input32 is derived from the retained encrypted initial RAM, without encryption.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import time

REPO = Path(__file__).resolve().parents[3]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for data in iter(lambda: f.read(1 << 20), b''):
            h.update(data)
    return h.hexdigest()


def cpu_sample(pid):
    # Bracket child CPU reads with system reads, as in the existing monitors.
    def system():
        values = list(map(int, Path('/proc/stat').read_text().splitlines()[0].split()[1:]))
        return sum(values[:3]) + sum(values[5:8])
    first = system()
    text = Path(f'/proc/{pid}/stat').read_text().rsplit(')', 1)[1].split()
    ticks = int(text[11]) + int(text[12])
    last = system()
    return time.monotonic(), (first + last) / 2, ticks


def run(name, command, env, directory):
    directory.mkdir(parents=True, exist_ok=False)
    job = {'argv': command, 'environment': {k: env[k] for k in env if k.startswith(('OMP_', 'CKKS_RAM_'))},
           'started_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
           'external_cpu_monitor_scope': 'All host CPUs, subtracting child process CPU; observations, not proof of isolation.'}
    (directory / 'job.json').write_text(json.dumps(job, indent=2) + '\n')
    interval = []
    with (directory / 'run.log').open('w') as log:
        process = subprocess.Popen(command, env=env, cwd=REPO, stdout=log, stderr=subprocess.STDOUT)
        previous = cpu_sample(process.pid)
        while True:
            try:
                process.wait(timeout=2)
                break
            except subprocess.TimeoutExpired:
                try:
                    current = cpu_sample(process.pid)
                except FileNotFoundError:
                    continue
                dt = current[0] - previous[0]
                external = max(0., (current[1] - previous[1] - (current[2] - previous[2])) / os.sysconf('SC_CLK_TCK') / dt)
                interval.append(external)
                previous = current
    job.update(exit_code=process.returncode, external_cpu_cores_samples=interval,
               maximum_external_cpu_cores=max(interval, default=0.),
               all_intervals_below_0_2_external_cores=bool(interval) and all(v < .2 for v in interval),
               completed_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    (directory / 'job.json').write_text(json.dumps(job, indent=2) + '\n')
    print(name, 'exit', process.returncode, 'external CPU max', job['maximum_external_cpu_cores'], flush=True)
    if process.returncode:
        raise RuntimeError(f'{name} failed; see {directory / "run.log"}')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--bridge-benchmark', type=Path)
    p.add_argument('--fixture', type=Path, help='Override the RAM fixture from the baseline receipt')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--baseline-receipt', type=Path, default=REPO / 'docs/ckks-ram-crossover.json')
    p.add_argument('--phase', choices=['all', 'pipeline', 'capture', 'small', 'comparison', 'bridge', 'sigma256'], default='all')
    args = p.parse_args()
    if args.phase in ('all','bridge') and args.bridge_benchmark is None:
        p.error('--bridge-benchmark is required for all/bridge phases')
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    receipt = json.loads(args.baseline_receipt.read_text())
    argv = receipt['provenance']['run-job']['argv']
    binary_index = argv.index('0-15') + 1
    retained = argv[binary_index+1:]
    if args.fixture is not None:
        retained[3]=str(args.fixture.resolve())
    checkpoint, fixture = Path(retained[0]), Path(retained[3])
    env = dict(os.environ, OMP_NUM_THREADS='8', OMP_DYNAMIC='FALSE', OMP_PROC_BIND='close',
               OMP_PLACES='{0},{2},{4},{6},{8},{10},{12},{14}')
    env.pop('CKKS_RAM_CAPTURE_ONLY', None)
    resource.setrlimit(resource.RLIMIT_STACK, (134217728, resource.getrlimit(resource.RLIMIT_STACK)[1]))
    base_command = ['taskset', '-c', '0-15', str(args.binary.resolve()), *retained]
    tracked = [p for p in fixture.iterdir() if p.is_file()]
    manifest_path = output / 'retained-manifest.json'
    if not manifest_path.exists():
        manifest_path.write_text(json.dumps({str(p): digest(p) for p in tracked}, indent=2) + '\n')
    provenance = {'binary_sha256': digest(args.binary), 'bridge_benchmark_sha256': digest(args.bridge_benchmark) if args.bridge_benchmark else None,
                  'source_revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip(),
                  'source_diff': subprocess.check_output(['git', 'diff'], cwd=REPO, text=True),
                  'hardware': subprocess.check_output(['lscpu'], text=True),
                  'scope': 'Diagnostic prototype comparison under assumed circular security; complete security is not inferred from timing or correctness.'}
    changed_sources = subprocess.check_output(['git', 'ls-files', '-m', '-o', '--exclude-standard', 'include', 'test'], cwd=REPO, text=True).splitlines()
    provenance['source_sha256'] = {path: digest(REPO / path) for path in changed_sources if (REPO / path).is_file() and Path(path).suffix in ('.hpp', '.cpp', '.py', '.sh', '.txt')}
    cache = args.binary.resolve().parents[2] / 'CMakeCache.txt'
    if cache.exists():
        provenance['cmake_cache'] = cache.read_text()
    (output / f'provenance-{args.phase}.json').write_text(json.dumps(provenance, indent=2) + '\n')
    if args.phase == 'sigma256':
        directory=output/'pipeline'
        selected=dict(env, CKKS_RAM_PROFILE='sigma256-qualify', CKKS_RAM_REPORT_DIRECTORY=str(output))
        run('sigma256',base_command,selected,directory)
    if args.phase in ('all', 'pipeline'):
        directory = output / 'pipeline'
        selected = dict(env, CKKS_RAM_PROFILE='ring-prototype', CKKS_RAM_REPORT_DIRECTORY=str(output))
        run('pipeline', base_command, selected, directory)
    if args.phase == 'capture':
        directory = output / 'capture'
        selected = dict(env, CKKS_RAM_PROFILE='crossover-full', CKKS_RAM_OUTPUT_SWITCH='ring', CKKS_RAM_CAPTURE_ONLY='1', CKKS_RAM_REPORT_DIRECTORY=str(directory))
        run('capture', base_command, selected, directory)
    if args.phase == 'small':
        for mode in ('direct', 'direct256', 'ring'):
            directory = output / ('small-' + mode)
            selected = dict(env, CKKS_RAM_PROFILE='4x2', CKKS_RAM_OUTPUT_SWITCH=mode, CKKS_RAM_REPORT_DIRECTORY=str(directory))
            run('small-' + mode, base_command, selected, directory)
    if args.phase == 'comparison':
        directory = output / 'comparison'
        selected = dict(env, CKKS_RAM_PROFILE='ring-compare', CKKS_RAM_REPORT_DIRECTORY=str(directory))
        run('comparison', base_command, selected, directory)
    if args.phase in ('all', 'bridge'):
        directory = output / 'bridge'
        command = ['taskset', '-c', '0-15', str(args.bridge_benchmark.resolve()), str(checkpoint), str(fixture), str(output / 'capture/output-input32.bin'), str(directory)]
        run('bridge', command, env, directory)
    before = json.loads(manifest_path.read_text())
    changed = [path for path, sha in before.items() if digest(path) != sha]
    (output / f'integrity-{args.phase}.json').write_text(json.dumps({'retained_files_checked': len(before), 'changed': changed}, indent=2) + '\n')
    if changed:
        raise RuntimeError('retained fixture mutation')


if __name__ == '__main__':
    main()
