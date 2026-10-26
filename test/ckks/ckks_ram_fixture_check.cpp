#include "ckks_ram_fixture.hpp"
#include <iomanip>
#include <numbers>
#ifdef __linux__
#include <sys/resource.h>
#endif
#ifdef USE_LVL1_SIGMA256
#if defined(__GNUC__)
[[gnu::noinline]]
#endif
double largeNoiseCheck(const std::filesystem::path&root,const TFHEpp::Key<ckks_ram::L>&secret,const TFHEpp::EvalKey&ek)
{
    using namespace TFHEpp;using namespace ckks_ram;
    auto fixture=std::make_unique<Fixture<13,1>>();load(root/"trace-8192x1.bin",*fixture);
    auto ram=std::make_unique<Fixture<13,1>::RAM>();ram->memory=fixture->initial;
    double maximum=0;
    const auto out=ram->access(fixture->addresses[0],fixture->enables[0],fixture->data[0],ek,[&](const auto&updated){
        for(std::size_t i=0;i<updated.size();++i){const double phase=std::ldexp(double(std::bit_cast<std::int32_t>(tlweSymPhase<L>(updated[i],secret))),-32);const bool bit=i==0||(i&1);const double restored=.125*std::sin(4*std::numbers::pi*phase);maximum=std::max(maximum,std::abs(restored-(bit?.125:-.125)));}
        return updated;
    });
    if(tlweSymDecrypt<L>(out[0],secret))throw std::runtime_error("warmup read mismatch");
    if(maximum>1./16)throw std::runtime_error("ideal-sine pre-switch noise check failed: "+std::to_string(maximum));
    return maximum;
}
#endif
template<std::size_t A,std::size_t W>
double check(const std::filesystem::path &root,const TFHEpp::Key<ckks_ram::L> &secret,const char *name)
{
    using namespace TFHEpp;
    auto fixture=std::make_unique<ckks_ram::Fixture<A,W>>();ckks_ram::load(root/name,*fixture);
    const auto plain=ckks_ram::trace<A,W>();
    if(fixture->addresses.size()!=plain.size()||fixture->enables.size()!=plain.size()||fixture->data.size()!=plain.size())throw std::runtime_error("trace length mismatch");
    for(std::size_t i=0;i<plain.size();++i){
        for(std::size_t a=0;a<A;++a)if(tlweSymDecrypt<ckks_ram::L>(fixture->addresses[i][a],secret)!=bool((plain[i].address>>a)&1))throw std::runtime_error("address mismatch");
        for(std::size_t w=0;w<W;++w)if(tlweSymDecrypt<ckks_ram::L>(fixture->data[i][w],secret)!=bool((plain[i].data>>w)&1))throw std::runtime_error("data mismatch");
        if(tlweSymDecrypt<ckks_ram::L>(fixture->enables[i],secret)!=plain[i].enable)throw std::runtime_error("enable mismatch");
    }
    long double square=0;
    for(std::size_t w=0;w<W;++w)for(std::size_t i=0;i<(std::size_t{1}<<A);++i){
        TLWE<ckks_ram::L> cell;SampleExtractIndex<ckks_ram::L>(cell,fixture->initial[w][i],0);
        if(tlweSymDecrypt<ckks_ram::L>(cell,secret)!=bool((i>>w)&1))throw std::runtime_error("initial memory mismatch");
        SampleExtractIndex<ckks_ram::L>(cell,fixture->initial[w][i],1);
        const double noise=std::bit_cast<std::int32_t>(tlweSymPhase<ckks_ram::L>(cell,secret));square+=noise*noise;
    }
    return std::sqrt(double(square/(W*(std::size_t{1}<<A))));
}
int main(int argc,char **argv)try{
    if(argc!=2 && !(argc==3&&std::string(argv[2])=="--ram-noise"))throw std::runtime_error("usage: ckks_ram_fixture_check FIXTURE_DIRECTORY [--ram-noise]");
#ifndef USE_LVL1_SIGMA256
    if(argc==3)throw std::runtime_error("--ram-noise requires the sigma256 profile");
#endif
    auto owner=std::make_unique<ckks_ram::SavedOwner>();ckks_ram::load(std::filesystem::path(argv[1])/"owner.bin",*owner);
    check<2,2>(argv[1],owner->secret(),"trace-4x2.bin");check<8,1>(argv[1],owner->secret(),"trace-256x1.bin");
#ifdef USE_LVL1_SIGMA256
    check<9,1>(argv[1],owner->secret(),"trace-512x1.bin");check<10,1>(argv[1],owner->secret(),"trace-1024x1.bin");
    const auto memory_rms=check<13,1>(argv[1],owner->secret(),"trace-8192x1.bin");
    auto ek=std::make_unique<TFHEpp::EvalKey>();ckks_ram::load(std::filesystem::path(argv[1])/"tfhe-eval.bin",*ek);
    using PKS=ckks_ram::PKS;using L=ckks_ram::L;
    const auto &key=ek->getprivksk<PKS>("privksk4cb_1");
    constexpr std::size_t digits=1<<(PKS::basebit-1),total=(PKS::domainP::n+1)*PKS::t*digits;
    long double square=0;
    for(std::size_t sample=0;sample<4096;++sample){const auto index=sample*(total-1)/4095;TFHEpp::TLWE<L> value;TFHEpp::SampleExtractIndex<L>(value,key[index/(PKS::t*digits)][(index/digits)%PKS::t][index%digits],1);const double noise=std::bit_cast<std::int32_t>(TFHEpp::tlweSymPhase<L>(value,owner->secret()));square+=noise*noise;}
    const double private_rms=std::sqrt(double(square/4096));
    if(memory_rms<.95*256||memory_rms>1.05*256||private_rms<.95*256||private_rms>1.05*256)throw std::runtime_error("retuned fixture noise audit failed");
    double predicted_error=0;
    if(argc==3){
#ifdef __linux__
        rlimit limit{};constexpr rlim_t needed=128*1024*1024;
        if(getrlimit(RLIMIT_STACK,&limit))throw std::runtime_error("cannot read stack limit");
        if(limit.rlim_cur<needed){if(limit.rlim_max<needed)throw std::runtime_error("large RAM diagnostic requires a 128 MiB stack limit");limit.rlim_cur=needed;if(setrlimit(RLIMIT_STACK,&limit))throw std::runtime_error("cannot raise RAM diagnostic stack limit");}
#endif
        predicted_error=largeNoiseCheck(argv[1],owner->secret(),*ek);
    }
    std::cout<<std::setprecision(17)<<"{\"lvl1_sigma\":256,\"lvl1_rows\":"<<L::l<<",\"lvl1_basebit\":"<<L::Bgbit<<",\"pks_rows\":"<<PKS::t<<",\"pks_basebit\":"<<PKS::basebit<<",\"all_five_traces_checked\":true,\"memory_noise_rms\":"<<memory_rms<<",\"private_switch_noise_rms\":"<<private_rms<<",\"predicted_pre_switch_max_error\":"<<predicted_error<<",\"large_noise_checked\":"<<(argc==3?"true":"false")<<",\"correct\":true}\n";
#else
    std::cout<<"Constructor-free owner archive and both saved traces passed"<<std::endl;
#endif
}catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
