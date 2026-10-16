"""Build an isolated independent/allocated/reused DD comparison source.

Usage: python3 test/ckks/tools/make_dd_scratch_replay.py EMPTY_DIRECTORY [--gadget-workspace]
Use --gadget-workspace to compare fresh/reused gadget buffers instead of outer
DD buffers. Mode1 is the selected frozen control and mode2 is production.
Compile replay.cpp using the normal retained-replay flags and add test/ckks to
its include path. Run the fixed grouped sparse-Tanh/native-batch profile with
--dd-windowed 1. Warm0/1/2, then measure0/1/2/2/1/0 with one prepared key chain.
"""
import hashlib
import json
from pathlib import Path
import re
import sys

if len(sys.argv) not in (2,3) or (len(sys.argv)==3 and sys.argv[2]!='--gadget-workspace'):
    raise SystemExit(__doc__)
repo = Path(__file__).resolve().parents[3]
out = Path(sys.argv[1]).resolve()
out.mkdir(parents=True, exist_ok=True)
if any(out.iterdir()):
    raise SystemExit('Use an empty experiment directory')
source = repo / 'test/ckks/ckks_narrow_ll14_basis_ntt_fbs.cpp'
gadget = len(sys.argv)==3
old = Path(__file__).with_name('dd_gadget_candidate.hpp' if gadget else 'dd_scratch_candidate.hpp')
main = source.read_text()
main = f'#include "{old.name}"\n#include <sys/resource.h>\n' + main
adapter = '''template<class Base>struct ReallocatedDDChain {
 std::shared_ptr<Base> base;
 template<std::size_t Level>decltype(auto)get()const{return base->template get<Level>();}
 template<std::size_t Level,class Outputs,class Inputs>
 void multiplyBatch(Outputs&next,const Inputs&factors,const std::vector<std::pair<std::size_t,std::size_t>>&pairs)const
  requires requires{std::declval<const Base&>().template get<Level>().key;}
 {
  const auto relin=base->template get<Level>();
  constexpr auto Q=std::remove_cvref_t<decltype(relin)>::log_q;
  std::vector<CKKSRNS32Product<P,Q,Narrow::log_delta>> work;
  work.reserve(pairs.size());
  for(const auto&[tree,pair]:pairs)work.push_back({&next[tree][pair],&factors[tree][2*pair],&factors[tree][2*pair+1]});
  dd_scratch_experiment::multiplyReallocated<P,Q,Narrow::log_delta>(work,relin.key);
 }
};
'''
if gadget:
    adapter=adapter.replace('dd_scratch_experiment::multiplyReallocated','dd_gadget_experiment::multiplyFreshGadgetWorkspace')
assert main.count('int main(') == 1
main = main.replace('int main(', adapter + 'int main(', 1)
main = main.replace(' auto windowed_dd_calls=', ' long scratch_minor_faults=0;\n auto windowed_dd_calls=', 1)
main = main.replace('independent_dd_batch_tree,windowed_dd_batch_tree;',
                    'independent_dd_batch_tree,windowed_dd_batch_tree,reallocated_dd_batch_tree;', 1)
needle = '    windowed_dd_calls->fetch_add(1,std::memory_order_relaxed);\n   };\n'
assert main.count(needle) == 1
main = main.replace(needle, needle + '''   using Chain=typename decltype(prepared)::element_type;
   auto old=std::make_shared<ReallocatedDDChain<Chain>>(ReallocatedDDChain<Chain>{prepared});
   reallocated_dd_batch_tree=[old,windowed_dd_calls](std::vector<Narrow::TreeResult>&out,std::vector<std::vector<Narrow::Factor>> factors,const Narrow::RelinKeys&raw){
    narrow_detail::productTreeBatch<Narrow,0,ReallocatedDDChain<Chain>,true>(out,std::move(factors),raw,old.get(),true);
    windowed_dd_calls->fetch_add(1,std::memory_order_relaxed);
   };
''', 1)
pair_line = next(line for line in main.splitlines() if 'paired replay controls must be separate' in line)
pair_flags = re.findall(r'int\((\w+)\)', pair_line)
assert len(pair_flags) >= 10
needle = ' const auto start=std::chrono::steady_clock::now();'
assert main.count(needle) == 1
main = main.replace(needle, ' if(requested!="tanh"||!sparse_tanh||!basis_accumulate||!basis_group_tensors||!native_auto_batch||!dd_windowed||' +
                    '||'.join(pair_flags) + ')throw std::invalid_argument("scratch experiment requires a fixed grouped sparse Tanh profile and windowed DD");\n' + needle, 1)
needle = ' const auto sweep_tree_evaluator=key->narrow.product_tree_evaluator;'
assert main.count(needle) == 1
main = main.replace(needle, ' functions.assign(9,requested);\n' + needle, 1)
lines = main.splitlines()
for i,line in enumerate(lines):
    if line.startswith(' const bool use_dd_windowed='):
        lines[i] = ''' constexpr std::array<unsigned,9> scratch_modes{0,1,2,0,1,2,2,1,0};
 const auto scratch_mode=scratch_modes.at(trial_index);
 const bool use_dd_windowed=scratch_mode!=0;'''
    if line.startswith(' if(dd_windowed||dd_windowed_pair)key->narrow.product_tree_batch_evaluator='):
        lines[i] += '''
 if(scratch_mode==1)key->narrow.product_tree_batch_evaluator=reallocated_dd_batch_tree;
 scratch_minor_faults=0;
 const auto scratch_backend=key->narrow.product_tree_batch_evaluator;
 key->narrow.product_tree_batch_evaluator=[scratch_backend,&scratch_minor_faults](auto&out,auto factors,const auto&raw){
  rusage before{},after{};getrusage(RUSAGE_SELF,&before);
  scratch_backend(out,std::move(factors),raw);
  getrusage(RUSAGE_SELF,&after);scratch_minor_faults+=after.ru_minflt-before.ru_minflt;
 };'''
    if line.startswith(' const auto artifact=function+'):
        assert line.endswith(';')
        lines[i] = line[:-1] + '+"-scratch"+std::to_string(trial_index);'
main = '\n'.join(lines) + '\n'
needle = r'<<",\"dd_windowed\":"'
assert main.count(needle) == 1
main = main.replace(needle, r'<<",\"dd_scratch_minor_faults\":"<<scratch_minor_faults<<",\"dd_scratch_mode\":"<<scratch_mode<<",\"dd_scratch_warmup\":"<<(trial_index<3?"true":"false")' + needle, 1)
(out / 'replay.cpp').write_text(main)
(out / old.name).write_text(old.read_text())
files = [source,old,repo/'include/ckks/prepared_dd_batch_rns32.hpp',Path(__file__).resolve()]
(out / 'source-manifest.json').write_text(json.dumps({str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},indent=2)+'\n')
