"""Generate a retained-key comparison of direct and reused root preparation.

Usage: python3 test/ckks/tools/make_table_roots_replay.py EMPTY_DIRECTORY
Compile replay.cpp with normal retained-replay flags and test/ckks include path.
Warm direct/reused, then measure direct/reused/reused/direct twice. Both policies
share all prepared keys and the retained input. No keys generated or encrypted.
"""
import hashlib
import json
from pathlib import Path
import re
import sys

if len(sys.argv)!=2:
    raise SystemExit(__doc__)
repo=Path(__file__).resolve().parents[3]
out=Path(sys.argv[1]).resolve()
out.mkdir(parents=True,exist_ok=True)
if any(out.iterdir()):
    raise SystemExit('Use an empty experiment directory')
source=repo/'test/ckks/ckks_narrow_ll14_basis_ntt_fbs.cpp'
main=source.read_text()
pair_line=next(line for line in main.splitlines() if 'paired replay controls must be separate' in line)
pairs=re.findall(r'int\((\w+)\)',pair_line)
assert len(pairs)>=10
needle=' const auto start=std::chrono::steady_clock::now();'
assert main.count(needle)==1
main=main.replace(needle,' if(requested!="tanh"||!sparse_tanh||!basis_accumulate||!basis_group_tensors||!native_auto_batch||!dd_windowed||!direct_rns||'+ '||'.join(pairs)+')throw std::invalid_argument("root reuse experiment requires fixed grouped sparse Tanh");\n'+needle,1)
needle=' const auto sweep_tree_evaluator=key->narrow.product_tree_evaluator;'
assert main.count(needle)==1
main=main.replace(needle,' functions.assign(10,requested);\n'+needle,1)
needle=' for(const auto& function:functions) {'
assert main.count(needle)==1
main=main.replace(needle,needle+'''
 constexpr std::array<bool,10> root_modes{false,true,false,true,true,false,false,true,true,false};
 const auto reuse_table_roots=root_modes.at(trial_index);
 key->narrow.reuse_table_roots=reuse_table_roots;''',1)
lines=main.splitlines()
for i,line in enumerate(lines):
    if line.startswith(' const auto artifact=function+'):
        assert line.endswith(';')
        lines[i]=line[:-1]+'+"-roots"+std::to_string(trial_index);'
main='\n'.join(lines)+'\n'
main=main.replace('  root_stage_times.leaf_factors_and_inline_alignment_seconds+=t.leaf_factors_and_inline_alignment_seconds;', '  root_stage_times.input_and_factor_cache_seconds+=t.input_and_factor_cache_seconds;\n  root_stage_times.leaf_factors_and_inline_alignment_seconds+=t.leaf_factors_and_inline_alignment_seconds;',1)
needle=r'<<",\"factor_stage_seconds\":"'
assert main.count(needle)==1
main=main.replace(needle,'<<",\\"root_input_stage_seconds\\":"<<root_stage_times.input_and_factor_cache_seconds'+r'<<",\"table_roots_reused\":"<<(reuse_table_roots?"true":"false")<<",\"table_roots_warmup\":"<<(trial_index<2?"true":"false")'+needle,1)
(out/'replay.cpp').write_text(main)
files=[source,repo/'include/ckks/narrow_bootstrapping.hpp',Path(__file__).resolve()]
(out/'source-manifest.json').write_text(json.dumps({str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},indent=2)+'\n')
