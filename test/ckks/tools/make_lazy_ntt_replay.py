"""Generate an isolated nine-trial NTT experiment; never edit production kernels.

Usage: python3 test/ckks/tools/make_lazy_ntt_replay.py OUTPUT_DIRECTORY
Compile the generated replay.cpp with include/ before the repository include
path and add the repository test/ckks include path. Run only the documented
retained sparse-Tanh profile. Modes warm0/1/2 then measure0/1/2/2/1/0.
"""
import hashlib
import json
from pathlib import Path
import re
import sys

repo = Path(__file__).resolve().parents[3]
out = Path(sys.argv[1]).resolve()
out.mkdir(parents=True, exist_ok=True)
if any(out.iterdir()):
    raise SystemExit('Use an empty experiment directory')
include = out / 'include'
include.mkdir()
header_path = repo / 'include/modular_ntt32.hpp'
rns_path = repo / 'include/modular_rns32.hpp'
candidate_path = Path(__file__).with_name('lazy_ntt_candidate.hpp')
main_path = repo / 'test/ckks/ckks_narrow_ll14_basis_ntt_fbs.cpp'
header = header_path.read_text()
needle = '\nstruct LazyNTTAccess;\n'
assert header.count(needle) == 1
header = header.replace('#include <algorithm>', '#include <algorithm>\n#include <atomic>', 1)
header = header.replace(needle, needle + '''
class NegacyclicNTTPlan;
inline std::atomic<unsigned> experimental_ntt_mode{0};
void experimentalNTTForward(const NegacyclicNTTPlan&,std::span<std::uint32_t>);
void experimentalNTTInverse(const NegacyclicNTTPlan&,std::span<std::uint32_t>);
''', 1)
for method, condition, function in [
    ('forwardBitReversed', '>0', 'experimentalNTTForward'),
    ('inverseBitReversed', '>1', 'experimentalNTTInverse')]:
    pattern = rf'(void {method}\(std::span<std::uint32_t> values\) const\s*\{{)'
    header, count = re.subn(pattern, r'\1\n        if(experimental_ntt_mode.load(std::memory_order_relaxed)' + condition +
                           '){' + function + '(*this,values);return;}', header)
    assert count == 1
candidate = '\n'.join(line for line in candidate_path.read_text().splitlines()
                      if not line.startswith('#include') and line != '#pragma once')
header += '\n' + candidate + '''
namespace TFHEpp::modular_ntt32 {
inline void experimentalNTTForward(const NegacyclicNTTPlan&p,std::span<std::uint32_t>x)
{LazyNTTAccess::forward(p,x,true);}
inline void experimentalNTTInverse(const NegacyclicNTTPlan&p,std::span<std::uint32_t>x)
{LazyNTTAccess::inverse(p,x,true);}
}
'''
(include / 'modular_ntt32.hpp').write_text(header)
(include / 'modular_rns32.hpp').write_text(rns_path.read_text())
main = main_path.read_text()
needle = ' const auto sweep_tree_evaluator=key->narrow.product_tree_evaluator;'
assert main.count(needle) == 1
main = main.replace(needle, ''' if(requested!="tanh"||!sparse_tanh||!basis_accumulate||!basis_group_tensors||!native_auto_batch||dd_windowed||dd_windowed_pair)
  throw std::invalid_argument("NTT experiment requires fixed grouped sparse Tanh, native auto batches and independent DD");
 functions.assign(9,requested);
''' + needle, 1)
lines = main.splitlines()
for i, line in enumerate(lines):
    if line.startswith(' const auto artifact=function+'):
        assert line.endswith(';')
        lines[i] = ''' constexpr std::array<unsigned,9> ntt_modes{0,1,2,0,1,2,2,1,0};
 const auto ntt_mode=ntt_modes.at(trial_index);
 modular_ntt32::experimental_ntt_mode.store(ntt_mode,std::memory_order_relaxed);
''' + line[:-1] + '+"-ntt"+std::to_string(trial_index);'
        break
else:
    raise AssertionError('missing artifact selection')
main = '\n'.join(lines) + '\n'
needle = r'<<",\"basis_group_tensors\":"'
assert main.count(needle) == 1
main = main.replace(needle, r'<<",\"ntt_experiment_mode\":"<<ntt_mode<<",\"ntt_experiment_warmup\":"<<(trial_index<3?"true":"false")' + needle, 1)
(out / 'replay.cpp').write_text(main)
(out / 'source-manifest.json').write_text(json.dumps({str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest()
    for p in [header_path,rns_path,candidate_path,main_path]},indent=2)+'\n')
