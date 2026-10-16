#include "ckks_ram_fixture.hpp"
template<std::size_t A,std::size_t W>
void check(const std::filesystem::path &root,const TFHEpp::Key<ckks_ram::L> &secret,const char *name)
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
    for(std::size_t w=0;w<W;++w)for(std::size_t i=0;i<(std::size_t{1}<<A);++i){
        TLWE<ckks_ram::L> cell;SampleExtractIndex<ckks_ram::L>(cell,fixture->initial[w][i],0);
        if(tlweSymDecrypt<ckks_ram::L>(cell,secret)!=bool((i>>w)&1))throw std::runtime_error("initial memory mismatch");
    }
}
int main(int argc,char **argv)try{
    if(argc!=2)throw std::runtime_error("usage: ckks_ram_fixture_check FIXTURE_DIRECTORY");
    auto owner=std::make_unique<ckks_ram::SavedOwner>();ckks_ram::load(std::filesystem::path(argv[1])/"owner.bin",*owner);
    check<2,2>(argv[1],owner->secret(),"trace-4x2.bin");check<8,1>(argv[1],owner->secret(),"trace-256x1.bin");
    std::cout<<"Constructor-free owner archive and both saved traces passed"<<std::endl;
}catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
