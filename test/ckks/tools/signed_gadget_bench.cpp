// Deterministic public-residue fixtures; no key generation or encryption.
// Build with the same defines, include paths and libraries as the configured
// ckks_dd_gadget_rns32 target. The validation receipt records the exact command.
// No argument: one apply with its internal team. Any argument:32 independent
// applies with one outer team, matching the batched product-tree context.
#include "../ckks_narrow_basis_benchmark.hpp"
#include <ckks/prepared_gadget_rns32.hpp>
#include <iostream>
#include <iomanip>
using namespace TFHEpp;
template<class P,uint32_t Q,uint32_t Bits>struct PublicRows {
    static constexpr size_t rows=(Q+Bits-1)/Bits;
    std::array<TRLWE<P>,rows> ciphertexts;
};
template<uint32_t Q,uint32_t Bits,uint32_t Chunk,uint32_t Lift,bool Batch>void benchmark()
{
    using P=narrow_benchmark::BenchmarkParam<14,7>;
    auto key=std::make_unique<PublicRows<P,Q,Bits>>();std::mt19937_64 random(20261017+Q+Bits);
    for(auto&r:key->ciphertexts)for(auto&c:r)for(auto&v:c){for(auto&w:v.limb)w=random();v=ckks_detail::reduceToLevel<P,Q>(v);}
    CKKSPreparedGadgetRNS32<P,Q,Bits,Chunk> prepared(*key);
    auto input=std::make_unique<Polynomial<P>>();auto reference=std::make_unique<TRLWE<P>>();
    for(auto&v:*input){for(auto&w:v.limb)w=random();v=ckks_detail::reduceToLevel<P,Q>(v<<Lift);}
    std::vector<TRLWE<P>> outputs(Batch?32:1);prepared.apply(*reference,*input,true,false);
    const auto run=[&](bool batch_encode){
        const auto start=std::chrono::steady_clock::now();
        if constexpr(Batch){
#pragma omp parallel for schedule(dynamic,1)
            for(size_t j=0;j<outputs.size();++j)prepared.apply(outputs[j],*input,true,batch_encode);
        }else prepared.apply(outputs[0],*input,true,batch_encode);
        const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        for(const auto&out:outputs)if(out!=*reference)throw std::runtime_error("signed batch encoding changed gadget output");
        return seconds;
    };
    run(false);run(true);
    for(bool batch_encode:{false,true,true,false})for(int trial=0;trial<5;++trial)
        std::cout<<"gadget q="<<Q<<" bits="<<Bits<<" chunk="<<Chunk<<" lift="<<Lift<<" batch="<<Batch<<" encode="<<batch_encode<<" trial="<<trial<<" seconds="<<run(batch_encode)<<" identical\n";
}
template<bool Batch>void all(){
    benchmark<433,40,145,48,Batch>();benchmark<389,40,80,48,Batch>();
    benchmark<345,40,115,48,Batch>();benchmark<301,40,80,48,Batch>();
    benchmark<257,40,145,48,Batch>();benchmark<435,44,112,37,Batch>();
    benchmark<435,44,112,270,Batch>();
}
int main(int argc,char**){std::cout<<std::setprecision(17);if(argc>1)all<true>();else all<false>();}
