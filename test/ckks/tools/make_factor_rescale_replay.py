"""Generate a shared-cache serial/parallel factor-rescale replay.

Usage: python3 test/ckks/tools/make_factor_rescale_replay.py EMPTY_DIRECTORY
Compile replay.cpp with retained replay flags and test/ckks on the include path.
Run the fixed grouped sparse-Tanh profile with native batches and windowed DD.
Warm both policies, then measure serial/parallel/parallel/serial. No keys generated.
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
frozen=Path(__file__).with_name('factor_rescale_source_candidate.hpp')
main='#include "factor_rescale_source_candidate.hpp"\n'+source.read_text()
main=main.replace(' if(direct_rns){',' decltype(key->narrow.factor_source) serial_projection_source,parallel_projection_source;\n if(direct_rns){',1)
needle='  using Source=CKKSNarrowDeferredRescaleFactorSourceRNS32<Narrow,TableSource,TableSource>;'
assert main.count(needle)==1
main=main.replace(needle,needle+'''
  serial_projection_source=std::make_shared<factor_rescale_experiment::SerialProjectionFactorSourceRNS32<Narrow,TableSource,TableSource>>(raw,retained,matrix_tiled,parallel_conjugation,matrix_dot8,factor_short_ntt,true,matrix_retained_rows,matrix_joint,parallel_table_ntt,matrix_fused_components,table_batch_encode);
  parallel_projection_source=std::make_shared<Source>(raw,retained,matrix_tiled,parallel_conjugation,matrix_dot8,factor_short_ntt,true,matrix_retained_rows,matrix_joint,parallel_table_ntt,matrix_fused_components,table_batch_encode);''',1)
pair_line=next(line for line in main.splitlines() if 'paired replay controls must be separate' in line)
pairs=re.findall(r'int\((\w+)\)',pair_line)
assert len(pairs)>=10
needle=' const auto start=std::chrono::steady_clock::now();'
assert main.count(needle)==1
main=main.replace(needle,' if(requested!="tanh"||!sparse_tanh||!basis_accumulate||!basis_group_tensors||!native_auto_batch||!dd_windowed||!direct_rns||!parallel_conjugation||'+ '||'.join(pairs)+')throw std::invalid_argument("factor rescale experiment requires fixed grouped sparse Tanh and parallel factors");\n'+needle,1)
needle=' const auto sweep_tree_evaluator=key->narrow.product_tree_evaluator;'
assert main.count(needle)==1
main=main.replace(needle,' functions.assign(6,requested);\n'+needle,1)
needle=' for(const auto& function:functions) {'
assert main.count(needle)==1
main=main.replace(needle,needle+'''
 constexpr std::array<bool,6> projection_modes{false,true,false,true,true,false};
 const auto parallel_projection=projection_modes.at(trial_index);
 key->narrow.factor_source=parallel_projection?parallel_projection_source:serial_projection_source;''',1)
lines=main.splitlines()
for i,line in enumerate(lines):
    if line.startswith(' const auto artifact=function+'):
        assert line.endswith(';')
        lines[i]=line[:-1]+'+"-rescale"+std::to_string(trial_index);'
main='\n'.join(lines)+'\n'
needle=r'<<",\"factor_stage_seconds\":"'
assert main.count(needle)==1
main=main.replace(needle,r'<<",\"factor_rescale_parallel\":"<<(parallel_projection?"true":"false")<<",\"factor_rescale_warmup\":"<<(trial_index<2?"true":"false")'+needle,1)
(out/'replay.cpp').write_text(main)
(out/frozen.name).write_text(frozen.read_text())
files=[source,frozen,repo/'include/ckks/narrow_factor_rns32.hpp',Path(__file__).resolve()]
(out/'source-manifest.json').write_text(json.dumps({str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},indent=2)+'\n')
