#include "ckks_ram_fixture.hpp"
#include <tfhe++.hpp>
#include <iomanip>

int main()try{
    using namespace TFHEpp;using namespace ckks_ram;
    if(std::ldexp(L::α,32)!=256 || PKS::targetP::α!=L::α)throw std::runtime_error("sigma256 target not selected");
    auto owner=std::make_unique<SecretKey>();const auto target=owner->key.getIndependent<L>();
    auto ek=std::make_unique<EvalKey>();
    ek->emplaceiksk<KS>(*owner);ek->emplacebkfft<CB>(*owner);ek->emplacebkfft<BR>(*owner);ek->emplaceprivksk4cb<PKS>(*owner);
    const auto &pks=ek->getprivksk<PKS>("privksk4cb_1");
    constexpr std::size_t digits=1<<(PKS::basebit-1),total=(PKS::domainP::n+1)*PKS::t*digits;
    long double square=0;
    for(std::size_t sample=0;sample<4096;++sample){
        const auto index=sample*(total-1)/4095,i=index/(PKS::t*digits),j=(index/digits)%PKS::t,u=index%digits;
        TLWE<L> value;SampleExtractIndex<L>(value,pks[i][j][u],1);
        const double noise=std::bit_cast<std::int32_t>(tlweSymPhase<L>(value,target));square+=noise*noise;
    }
    const double rms=std::sqrt(double(square/4096));
    if(rms<.95*256||rms>1.05*256)throw std::runtime_error("private-switch noise was not retuned");
    for(int repeat=0;repeat<16;++repeat)for(int a=0;a<2;++a)for(int b=0;b<2;++b){
        TLWE<L> left,right,out;tlweSymEncrypt<L>(left,a?L::μ:-L::μ,target);tlweSymEncrypt<L>(right,b?L::μ:-L::μ,target);
        HomNAND<KS,BR>(out,left,right,*ek);
        if(tlweSymDecrypt<L>(out,target)!=!(a&&b))throw std::runtime_error("sigma256 NAND mismatch");
    }
    for(int repeat=0;repeat<4;++repeat)for(int bit=0;bit<2;++bit){
        TLWE<L> input;tlweSymEncrypt<L>(input,bit?L::μ:-L::μ,target);
        TRGSWFFT<L> select,inverse;CircuitBootstrappingWithInv<KS,CB,PKS>(select,inverse,input,*ek);
        Polynomial<L> zero,one;zero.fill(-L::μ);one.fill(L::μ);
        TRLWE<L> c0,c1,out;trlweSymEncrypt<L>(c0,zero,target);trlweSymEncrypt<L>(c1,one,target);
        CMUXFFT<L>(out,select,c1,c0);
        const auto decoded=trlweSymDecrypt<L>(out,target);
        for(bool value:decoded)if(value!=bool(bit))throw std::runtime_error("sigma256 circuit bootstrap mismatch");
    }
    using RAM=Fixture<2,2>::RAM;auto ram=std::make_unique<RAM>();
    for(std::size_t w=0;w<2;++w)for(std::size_t i=0;i<4;++i){Polynomial<L> message{};message[0]=(i>>w)&1?L::μ:-L::μ;trlweSymEncrypt<L>(ram->memory[w][i],message,target);}
    const auto refresh=[&](const std::vector<TLWE<L>>&in){std::vector<TLWE<L>> out(in.size());for(std::size_t i=0;i<in.size();++i){TLWE<KS::targetP> low;TRLWE<L> cell;IdentityKeySwitch<KS>(low,in[i],ek->getiksk<KS>());BlindRotate<BR>(cell,low,ek->getbkfft<BR>(),μpolygen<L,L::μ>());SampleExtractIndex<L>(out[i],cell,0);}return out;};
    std::array<std::uint32_t,4> expected{0,1,2,3};RAM::Word previous{};double maximum=0;
    const auto requests=trace<2,2>();
    for(std::size_t step=0;step<requests.size();++step){const auto request=requests[step];RAM::Address address;RAM::Word data;TLWE<L> enable;
        for(std::size_t i=0;i<2;++i){tlweSymEncrypt<L>(address[i],(request.address>>i)&1?L::μ:-L::μ,target);tlweSymEncrypt<L>(data[i],(request.data>>i)&1?L::μ:-L::μ,target);}
        tlweSymEncrypt<L>(enable,request.enable?L::μ:-L::μ,target);
        const auto output=ram->access(address,enable,step==9?previous:data,*ek,refresh);
        for(std::size_t w=0;w<2;++w)if(tlweSymDecrypt<L>(output[w],target)!=bool((expected[request.address]>>w)&1))throw std::runtime_error("sigma256 RAM read mismatch");
        previous=output;if(request.enable)expected[request.address]=request.data;
        for(std::size_t w=0;w<2;++w)for(std::size_t i=0;i<4;++i){TLWE<L> cell;SampleExtractIndex<L>(cell,ram->memory[w][i],0);const bool bit=(expected[i]>>w)&1;if(tlweSymDecrypt<L>(cell,target)!=bit)throw std::runtime_error("sigma256 RAM cell mismatch");const double phase=std::ldexp(double(std::bit_cast<std::int32_t>(tlweSymPhase<L>(cell,target))),-32);maximum=std::max(maximum,std::abs(phase-(bit?.125:-.125)));}
    }
    std::cout<<std::setprecision(17)<<"{\"lvl1_sigma\":256,\"private_switch_noise_rms\":"<<rms<<",\"noise_samples\":4096,\"nand_checks\":64,\"circuit_bootstrap_bit_checks\":8192,\"ram_accesses\":16,\"feedback_checked\":true,\"all_cells_checked\":true,\"maximum_ram_phase_error\":"<<maximum<<",\"correct\":true}\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
