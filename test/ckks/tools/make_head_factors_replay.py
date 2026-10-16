"""Generate a retained-key comparison of head-leaf construction policies.

Usage: python3 test/ckks/tools/make_head_factors_replay.py EMPTY_DIRECTORY
Compile replay.cpp with normal retained-replay flags and test/ckks include path.
Warm preinitialized/appended, then measure preinitialized/appended/appended/
preinitialized twice. All prepared keys, root reuse and retained input are shared.
"""
import hashlib
import json
from pathlib import Path
import re
import sys

def replace_once(text,old,new):
    assert text.count(old)==1,old
    return text.replace(old,new,1)

if len(sys.argv)!=2:
    raise SystemExit(__doc__)
repo=Path(__file__).resolve().parents[3]
out=Path(sys.argv[1]).resolve()
out.mkdir(parents=True,exist_ok=True)
if any(out.iterdir()):
    raise SystemExit('Use an empty experiment directory')
source=repo/'test/ckks/ckks_narrow_ll14_basis_ntt_fbs.cpp'
main='#include <sys/resource.h>\n'+source.read_text()
pair_line=next(line for line in main.splitlines() if 'paired replay controls must be separate' in line)
pairs=re.findall(r'int\((\w+)\)',pair_line)
assert len(pairs)>=10
needle=' const auto start=std::chrono::steady_clock::now();'
guard=' if(requested!="tanh"||!sparse_tanh||!basis_accumulate||!basis_group_tensors||!native_auto_batch||!dd_windowed||!direct_rns||'+'||'.join(pairs)+')throw std::invalid_argument("head append experiment requires fixed grouped sparse Tanh");\n'
main=replace_once(main,needle,guard+needle)
needle=' const auto sweep_tree_evaluator=key->narrow.product_tree_evaluator;'
main=replace_once(main,needle,' functions.assign(10,requested);\n'+needle)
needle=' for(const auto& function:functions) {'
main=replace_once(main,needle,needle+'''
 constexpr std::array<bool,10> head_modes{false,true,false,true,true,false,false,true,true,false};
 const auto append_head_factors=head_modes.at(trial_index);
 key->narrow.append_head_factors=append_head_factors;
 key->narrow.reuse_table_roots=true;''')
root_line=next(line for line in main.splitlines() if line.startswith(' std::vector<Narrow::Output> roots;CKKSNarrowBlindRotateBatch'))
main=replace_once(main,root_line,' rusage root_usage_before{},root_usage_after{};getrusage(RUSAGE_SELF,&root_usage_before);\n'+root_line+'\n getrusage(RUSAGE_SELF,&root_usage_after);')
lines=main.splitlines()
for i,line in enumerate(lines):
    if line.startswith(' const auto artifact=function+'):
        assert line.endswith(';')
        lines[i]=line[:-1]+'+"-head"+std::to_string(trial_index);'
main='\n'.join(lines)+'\n'
needle='  root_stage_times.leaf_factors_and_inline_alignment_seconds+=t.leaf_factors_and_inline_alignment_seconds;'
main=replace_once(main,needle,'  root_stage_times.input_and_factor_cache_seconds+=t.input_and_factor_cache_seconds;\n'+needle)
needle=r'<<",\"factor_stage_seconds\":"'
fields=(r'<<",\"root_input_stage_seconds\":"<<root_stage_times.input_and_factor_cache_seconds'
        r'<<",\"root_minor_faults\":"<<(root_usage_after.ru_minflt-root_usage_before.ru_minflt)'
        r'<<",\"head_factors_appended\":"<<(append_head_factors?"true":"false")'
        r'<<",\"head_factors_warmup\":"<<(trial_index<2?"true":"false")')
main=replace_once(main,needle,fields+needle)
assert main.index('root_minor_faults')>main.index('getrusage(RUSAGE_SELF,&root_usage_after);')
assert main.count('root_minor_faults')==1
(out/'replay.cpp').write_text(main)
files=[source,repo/'include/ckks/narrow_bootstrapping.hpp',Path(__file__).resolve()]
(out/'source-manifest.json').write_text(json.dumps({str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},indent=2)+'\n')
