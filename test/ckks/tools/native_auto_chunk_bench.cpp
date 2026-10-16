// Deterministic noncryptographic rows, reused across layouts. No encryption.
// Build with the configured ckks_dd_gadget_rns32 flags and libraries.
#include "../ckks_narrow_basis_benchmark.hpp"
#include <ckks/public_auto_ntt.hpp>
#include <iomanip>
#include <iostream>
using namespace TFHEpp;
using P=narrow_benchmark::BenchmarkParam<14,7>;
template<std::uint32_t Chunk>using Cache=PreparedPublicAutoNTT<P,435,44,20,true,Chunk>;
template<std::uint32_t DataQ>void benchmark(const Cache<80>&c80,const Cache<112>&c112,
    const Cache<145>&c145,const Cache<218>&c218,std::size_t batch)
{
    auto input=std::make_unique<TRLWE<P>>();std::mt19937_64 random(20261025+DataQ);
    for(auto&component:*input)for(auto&v:component)for(auto&w:v.limb)w=random();
    (*input)[0][0]={};(*input)[0][1]=P::T{}-P::T{1};
    (*input)[0][2]=P::T{1}<<(DataQ-1);(*input)[0][3]=(P::T{1}<<DataQ)+P::T{3};
    auto reference=std::make_unique<TRLWE<P>>();c112.template apply<DataQ>(*reference,*input);
    std::vector<TRLWE<P>> outputs(batch);
    const auto run=[&](const auto&cache){
        const auto start=std::chrono::steady_clock::now();
        if(batch==1)cache.template apply<DataQ>(outputs[0],*input);
        else {
#pragma omp parallel for schedule(static) num_threads(batch)
            for(std::size_t j=0;j<batch;++j)cache.template apply<DataQ>(outputs[j],*input);
        }
        const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        for(const auto&out:outputs)if(out!=*reference)throw std::runtime_error("native automorphism chunks changed output");
        return seconds;
    };
    run(c80);run(c112);run(c145);run(c218);
    for(int chunk:{112,80,145,218,218,145,80,112})for(int trial=0;trial<3;++trial){
        const auto seconds=chunk==80?run(c80):chunk==112?run(c112):chunk==145?run(c145):run(c218);
        std::cout<<"data_q="<<DataQ<<" batch="<<batch<<" chunk="<<chunk<<" trial="<<trial<<" seconds="<<seconds<<" identical\n";
    }
}
int main(int argc,char**argv)
{
    const auto batch=argc>1?std::stoul(argv[1]):1;
    if(batch!=1&&batch!=4&&batch!=8)throw std::invalid_argument("batch must be1,4 or8");
    auto raw=std::make_unique<PublicAutoGadgetRows<P,435,44>>();raw->exponent=2*P::n-1;
    std::mt19937_64 random(20261025);
    for(auto&row:raw->ciphertexts)for(auto&component:row)for(auto&v:component){
        for(auto&w:v.limb)w=random();
        v=ckks_detail::reduceToLevel<P,435>(v);
    }
    Cache<80>c80(*raw);Cache<112>c112(*raw);Cache<145>c145(*raw);Cache<218>c218(*raw);
    std::cout<<std::setprecision(17)<<"bytes80="<<c80.byteSize()<<" bytes112="<<c112.byteSize()<<" bytes145="<<c145.byteSize()<<" bytes218="<<c218.byteSize()<<"\n";
    benchmark<398>(c80,c112,c145,c218,batch);benchmark<385>(c80,c112,c145,c218,batch);
    benchmark<165>(c80,c112,c145,c218,batch);benchmark<93>(c80,c112,c145,c218,batch);
}
