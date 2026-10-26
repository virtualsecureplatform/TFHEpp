#include <tfhe/lweParams.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/array.hpp>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace TFHEpp;
template<class F> void rejects(F action,const char*expected){try{action();}catch(const std::exception&e){if(std::string(e.what()).find(expected)!=std::string::npos)return;throw;}throw std::runtime_error("incompatible noise profile accepted");}
int main(int argc,char**argv)try{
    lweParams parameters;
    if(parameters.lvl1.α!=TFHEParameterNoiseValue<lvl1param>())throw std::runtime_error("level1 noise is uninitialized");
    if(argc==3){
        if(std::string(argv[1])=="--write"){
            std::ofstream file(argv[2],std::ios::binary);cereal::PortableBinaryOutputArchive ar(file);ar(parameters);return 0;
        }
        if(std::string(argv[1])=="--reject"){
            rejects([&]{std::ifstream file(argv[2],std::ios::binary);cereal::PortableBinaryInputArchive ar(file);ar(parameters);},"noise profile mismatch");return 0;
        }
        throw std::runtime_error("unknown noise-profile command");
    }
    if(argc!=1)throw std::runtime_error("invalid noise-profile arguments");
    std::stringstream saved;{cereal::PortableBinaryOutputArchive ar(saved);ar(parameters);}
    lweParams restored;{cereal::PortableBinaryInputArchive ar(saved);ar(restored);}
    if(!(restored==parameters)||!restored.noise_profile_bound)throw std::runtime_error("profile roundtrip mismatch");
    auto changed=parameters;changed.lvl1.α*=2;
    if(changed==parameters)throw std::runtime_error("noise omitted from equality");
    rejects([&]{std::stringstream f;cereal::PortableBinaryOutputArchive ar(f);ar(changed);},"cannot relabel");
    std::stringstream wrong;{cereal::PortableBinaryOutputArchive ar(wrong);auto profile=TFHEArchiveNoiseProfile();profile[18]^=1;ar(std::uint32_t{0x54465032},profile);}
    rejects([&]{cereal::PortableBinaryInputArchive ar(wrong);ar(restored);},"noise profile mismatch");
    std::stringstream payload;{cereal::PortableBinaryOutputArchive ar(payload);ar(std::uint32_t{0x54465032},TFHEArchiveNoiseProfile(),changed.lvl0.n);changed.archiveTail(ar);}
    rejects([&]{cereal::PortableBinaryInputArchive ar(payload);ar(restored);},"parameter payload mismatch");
    std::stringstream legacy;{cereal::PortableBinaryOutputArchive ar(legacy);ar(parameters.lvl0.n);parameters.archiveTail(ar);}
#ifdef USE_LVL1_SIGMA256
    rejects([&]{cereal::PortableBinaryInputArchive ar(legacy);ar(restored);},"sigma256 requires");
#else
    {cereal::PortableBinaryInputArchive ar(legacy);ar(restored);}
    if(restored.noise_profile_bound)throw std::runtime_error("legacy archive silently upgraded");
    rejects([&]{std::stringstream f;cereal::PortableBinaryOutputArchive ar(f);ar(restored);},"cannot relabel");
#endif
    std::cout<<"TFHE noise initialization, binding, payload and legacy checks passed\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
