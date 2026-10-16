"""Generate one warm-up and five GeLU latency trials sharing retained caches.

Usage: python3 test/ckks/tools/make_gelu_latency_replay.py EMPTY_DIRECTORY
Compile replay.cpp with normal retained-replay flags and test/ckks include path.
Use the original-degree grouped GeLU profile with eight OpenMP threads.
Only the replay loop and diagnostic metadata change; no keys/input are generated.
"""
import hashlib
import json
from pathlib import Path
import re
import sys


def replace_once(text, old, new):
    assert text.count(old) == 1, old
    return text.replace(old, new, 1)


if len(sys.argv) != 2:
    raise SystemExit(__doc__)
repo = Path(__file__).resolve().parents[3]
out = Path(sys.argv[1]).resolve()
out.mkdir(parents=True, exist_ok=True)
if any(out.iterdir()):
    raise SystemExit('Use an empty experiment directory')
source = repo / 'test/ckks/ckks_narrow_ll14_basis_ntt_fbs.cpp'
main = source.read_text()
pair_line = next(line for line in main.splitlines()
                 if 'paired replay controls must be separate' in line)
pairs = re.findall(r'int\((\w+)\)', pair_line)
assert len(pairs) >= 10
needle = ' const auto start=std::chrono::steady_clock::now();'
guard = (' if(requested!="gelu"||sparse_tanh||reduced_degree||'
         '!basis_accumulate||!basis_group_tensors||!native_auto_batch||'
         '!dd_windowed||!direct_rns||narrow_benchmark::benchmarkThreads()!=8||'
         + '||'.join(pairs)
         + ')throw std::invalid_argument("latency replay requires fixed original-degree grouped GeLU with eight threads");\n')
main = replace_once(main, needle, guard + needle)
needle = ' const auto sweep_tree_evaluator=key->narrow.product_tree_evaluator;'
main = replace_once(main, needle, ' functions.assign(6,requested);\n' + needle)
lines = main.splitlines()
for i, line in enumerate(lines):
    if line.startswith(' const auto artifact=function+'):
        assert line.endswith(';')
        lines[i] = line[:-1] + '+"-latency"+std::to_string(trial_index);'
main = '\n'.join(lines) + '\n'
needle = r'<<",\"factor_stage_seconds\":"'
main = replace_once(main, needle,
                    r'<<",\"latency_trial\":"<<trial_index'
                    r'<<",\"latency_warmup\":"<<(trial_index==0?"true":"false")' + needle)
(out / 'replay.cpp').write_text(main)
files = [source, repo / 'include/ckks/narrow_factor_rns32.hpp',
         repo / 'include/ckks/narrow_bootstrapping.hpp', Path(__file__).resolve()]
(out / 'source-manifest.json').write_text(json.dumps({
    str(p.relative_to(repo)): hashlib.sha256(p.read_bytes()).hexdigest()
    for p in files}, indent=2) + '\n')
