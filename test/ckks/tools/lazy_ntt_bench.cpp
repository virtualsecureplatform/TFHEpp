// Exact finite-field fixtures only; no encryption or key generation.
#include "lazy_ntt_candidate.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <omp.h>
using namespace TFHEpp;
using Access=modular_ntt32::LazyNTTAccess;
std::vector<modular_ntt::PrimeModulus> primes(){
    std::vector<modular_ntt::PrimeModulus> p(modular_rns32::primes.begin(),modular_rns32::primes.end());
    for(auto x:modular_ntt32::primes)if(std::none_of(p.begin(),p.end(),[&](auto y){return x.value==y.value;}))p.push_back(x);
    return p;
}
std::size_t reverse(std::size_t i,std::size_t n){std::size_t v=0;for(;n>1;n>>=1){v=2*v+(i&1);i>>=1;}return v;}
void check(){
    std::mt19937_64 random(20261030);std::size_t cases=0;
    for(auto prime:primes())for(bool shorts:{false,true})for(std::size_t n=2;n<=65536;n*=2){
        modular_ntt32::NegacyclicNTTPlan plan(n,prime,shorts);
        std::vector<std::uint32_t> in(n),a,b;
        for(std::size_t i=0;i<n;++i)in[i]=i%8==0?0:i%8==1?prime.value-1:random()%prime.value;
        for(bool bit:{false,true}){
            a=b=in;
            if(bit)plan.forwardBitReversed(a);else plan.forward(a);
            Access::forward(plan,b,bit,true);
            if(a!=b)throw std::runtime_error("lazy forward differs from original spectrum");
            if(n<=16){
                const auto psi=modular_ntt::power(prime.primitive_root,(prime.value-1)/(2*n),prime.value);
                for(std::size_t k=0;k<n;++k){
                    const auto root=modular_ntt::power(psi,2*k+1,prime.value);std::uint64_t sum=0,power=1;
                    for(auto x:in){sum=(sum+x*power)%prime.value;power=power*root%prime.value;}
                    if(b[bit?reverse(k,n):k]!=sum)throw std::runtime_error("lazy forward differs from independent evaluation");
                }
            }
            Access::inverse(plan,b,bit,true);if(b!=in)throw std::runtime_error("lazy round trip differs");
            a=b=in;
            if(bit)plan.inverseBitReversed(a);else plan.inverse(a);
            Access::inverse(plan,b,bit,true);if(a!=b)throw std::runtime_error("lazy inverse differs on arbitrary spectrum");
            a=in;b=in;Access::forward(plan,a,bit,true);Access::forward(plan,b,bit,true);
            for(std::size_t i=0;i<n;++i)a[i]=std::uint64_t(a[i])*b[i]%prime.value;
            Access::inverse(plan,a,bit,true);
            if(n<=16)for(std::size_t k=0;k<n;++k){
                std::uint64_t expected=0;
                for(std::size_t i=0;i<n;++i){
                    const auto j=(k+n-i)%n;const auto product=std::uint64_t(in[i])*in[j]%prime.value;
                    expected=(expected+(i+j<n?product:prime.value-product))%prime.value;
                }
                if(a[k]!=expected)throw std::runtime_error("lazy convolution differs from schoolbook");
            }
            bool rejected=false;try{std::vector<std::uint32_t>bad(n+1);Access::forward(plan,bad,bit);}catch(const std::invalid_argument&){rejected=true;}
            if(!rejected)throw std::runtime_error("bad lazy NTT shape accepted");++cases;
        }
    }
    std::cout<<"checked "<<cases<<" transform cases; every stage range, canonical output, original spectrum/inverse, independent small-ring evaluation and schoolbook convolution pass\n";
}
void benchmark(std::size_t batch,int threads,bool bit,bool shorts){
    constexpr std::size_t n=16384;
    const auto ps=primes();std::vector<modular_ntt32::NegacyclicNTTPlan> plans;
    for(auto p:ps)plans.emplace_back(n,p,shorts);
    std::vector<std::vector<std::uint32_t>> in(batch),out(batch),expected(batch);
    std::mt19937_64 random(20261030);
    for(std::size_t j=0;j<batch;++j){in[j].resize(n);for(auto&v:in[j])v=random()%ps[j%ps.size()].value;expected[j]=in[j];if(bit)plans[j%ps.size()].forwardBitReversed(expected[j]);else plans[j%ps.size()].forward(expected[j]);}
    const auto run=[&](int mode){
        for(std::size_t j=0;j<batch;++j)out[j]=in[j];
        auto start=std::chrono::steady_clock::now();
#pragma omp parallel for schedule(dynamic,1) num_threads(threads)
        for(std::size_t j=0;j<batch;++j){const auto&p=plans[j%plans.size()];
            if(mode)Access::forward(p,out[j],bit);else if(bit)p.forwardBitReversed(out[j]);else p.forward(out[j]);}
        const auto forward=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        if(out!=expected)throw std::runtime_error("timed forward mismatch");
        start=std::chrono::steady_clock::now();
#pragma omp parallel for schedule(dynamic,1) num_threads(threads)
        for(std::size_t j=0;j<batch;++j){const auto&p=plans[j%plans.size()];
            if(mode)Access::inverse(p,out[j],bit);else if(bit)p.inverseBitReversed(out[j]);else p.inverse(out[j]);}
        const auto inverse=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        if(out!=in)throw std::runtime_error("timed inverse mismatch");return std::pair{forward,inverse};
    };
    run(0);run(1);
    for(int mode:{0,1,1,0})for(int trial=0;trial<8;++trial){auto [f,i]=run(mode);std::cout<<"batch="<<batch<<" threads="<<threads<<" bit="<<bit<<" short="<<shorts<<" mode="<<mode<<" trial="<<trial<<" forward="<<f<<" inverse="<<i<<" identical\n";}
}
int main(int argc,char**argv)try{
    std::cout<<std::setprecision(17);
    if(argc>1&&std::string(argv[1])=="--check"){check();return 0;}
    if(argc>1&&std::string(argv[1])=="--bit-reversed-only"){
        for(bool shorts:{true,false}){benchmark(17,1,true,shorts);benchmark(17*8,24,true,shorts);}return 0;
    }
    for(bool bit:{true,false})for(bool shorts:{true,false}){benchmark(17,1,bit,shorts);benchmark(17*8,24,bit,shorts);}
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
