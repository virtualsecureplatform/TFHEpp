// Deterministic coefficient fixtures only; no encryption or key generation.
#include "../ckks_narrow_basis_benchmark.hpp"
#include "dd_batch_candidate.hpp"
#include <iomanip>
#include <iostream>
using namespace TFHEpp;
template<class P,std::uint32_t Q,std::uint32_t Bits>struct Rows {
    static constexpr std::size_t rows=(Q+Bits-1)/Bits;
    std::array<TRLWE<P>,rows> ciphertexts;
};
template<std::uint32_t Nbit,std::uint32_t Q,std::uint32_t Chunk>
void run(std::size_t size,bool timed){
    using P=narrow_benchmark::BenchmarkParam<Nbit,7>;
    using In=CKKSCiphertext<P,Q,44>;using Out=CKKSCiphertext<P,Q-44,44>;
    using Key=CKKSAuxiliaryDDRelinKeyNTT<P,Q,48,40,20,true,Chunk>;
    auto rows=std::make_unique<Rows<P,Q+48,40>>();std::mt19937_64 random(20261030+Q);
    const auto fill=[&](auto&poly,auto bits){for(auto&v:poly){for(auto&w:v.limb)w=random();v=bits(v);}};
    for(auto&r:rows->ciphertexts)for(auto&c:r)fill(c,[](auto v){return ckks_detail::reduceToLevel<P,Q+48>(v);});
    Key key(std::make_shared<typename Key::Cache>(*rows));
    CKKSRelinKeyTensorRNS32<P,Key> original{key};
    std::vector<In> left(size),right(size);std::vector<Out> expected(size),actual(size);
    for(std::size_t j=0;j<size;++j)for(std::size_t c=0;c<2;++c){
        fill(left[j].ct[c],[](auto v){return ckks_detail::reduceToLevel<P,Q>(v);});
        fill(right[j].ct[c],[](auto v){return ckks_detail::reduceToLevel<P,Q>(v);});
    }
    if(!timed&&size>=4){
        left[0].ct={};right[1].ct={};
        for(auto&c:left[2].ct)c.fill((typename P::T{1}<<(Q-1))+typename P::T{1});
        for(auto&c:right[3].ct)c.fill(typename P::T{}-typename P::T{1});
    }
    const auto normal_threads=std::min<int>(size,omp_get_max_threads());
#pragma omp parallel for schedule(dynamic,1) num_threads(normal_threads)
    for(std::size_t j=0;j<size;++j)CKKSMultRelinThenRescale<P>(expected[j],left[j],right[j],original);
    const auto measure=[&](int mode){
        auto start=std::chrono::steady_clock::now();
        if(mode<=1){
            const auto threads=mode?omp_get_max_threads():normal_threads;
#pragma omp parallel for schedule(dynamic,1) num_threads(threads)
            for(std::size_t j=0;j<size;++j)CKKSMultRelinThenRescale<P>(actual[j],left[j],right[j],original);
        }else dd_batch_experiment::multiply<P,Q,44,48,40,Chunk>(actual,left,right,key,mode);
        const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        for(std::size_t j=0;j<size;++j)if(actual[j].ct!=expected[j].ct)throw std::runtime_error("DD batch changed ciphertext coefficients");
        return seconds;
    };
    if(!timed){
        for(int window:{3,8,16,32})measure(window);
        std::vector<Out> empty;std::vector<In> no_input;
        dd_batch_experiment::multiply<P,Q,44,48,40,Chunk>(empty,no_input,no_input,key);
        bool rejected=false;try{dd_batch_experiment::multiply<P,Q,44,48,40,Chunk>(empty,left,right,key);}catch(const std::invalid_argument&){rejected=true;}
        if(!rejected)throw std::runtime_error("DD batch shape accepted");
        rejected=false;try{dd_batch_experiment::multiply<P,Q,44,48,40,Chunk>(actual,left,right,key,0);}catch(const std::invalid_argument&){rejected=true;}
        if(!rejected)throw std::runtime_error("DD zero window accepted");
        std::cout<<"check N="<<P::n<<" q="<<Q<<" batch="<<size<<" identical\n";
    }else{
        for(int mode:{0,1,8,16,32})measure(mode);
        for(int mode:{0,1,8,16,32,32,16,8,1,0})for(int trial=0;trial<3;++trial)
            std::cout<<"q="<<Q<<" batch="<<size<<" mode="<<mode<<" trial="<<trial<<" seconds="<<measure(mode)<<" identical\n";
    }
}
int main(int argc,char**argv)try{
    const bool check=argc>1&&std::string(argv[1])=="--check";std::cout<<std::setprecision(10);
    if(check){
        run<4,385,145>(7,false);run<4,341,80>(7,false);run<4,297,115>(7,false);run<4,253,80>(7,false);run<4,209,145>(7,false);
        run<9,385,145>(7,false);run<14,385,145>(7,false);run<14,341,80>(7,false);run<14,297,115>(7,false);run<14,253,80>(7,false);run<14,209,145>(7,false);
    }else for(auto size:{32u,16u,8u,4u,2u}){
        run<14,385,145>(size,true);run<14,341,80>(size,true);run<14,297,115>(size,true);run<14,253,80>(size,true);run<14,209,145>(size,true);
    }
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
