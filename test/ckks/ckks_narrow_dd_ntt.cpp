#include <ckks/prepared_dd_ntt.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>

namespace {
struct P {
    static constexpr std::int32_t key_value_max = 1, key_value_min = -1;
    static constexpr std::uint32_t nbit = 14, n = 1 << nbit, k = 1;
    static constexpr std::uint32_t l = 1, l_a = 1, lₐ = 1;
    static constexpr std::uint32_t Bgbit = 4, Bg_abit = 4, Bgₐbit = 4;
    using T = TFHEpp::MultiLimbUInt<7>;
    static constexpr T Bg = T{1} << Bgbit, Bg_a = Bg, Bgₐ = Bg;
    static constexpr TFHEpp::ErrorDistribution errordist = TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double α = std::ldexp(1.0, -448);
    static constexpr double ckks_min_noise_stddev = 3.2;
    static constexpr T μ = T{1} << 253;
    static constexpr std::uint32_t plain_modulusbit = 20;
    static constexpr T plain_modulus = T{786433};
    static constexpr double Δ = 0.0;
    static constexpr std::uint32_t l̅ = 112, l̅ₐ = l̅, B̅gbit = 4, B̅gₐbit = 4;
};

using K=TFHEpp::CKKSCompactDDRelinKey<P,433,64,16>;
using Chain=TFHEpp::CKKSCompactDDRelinKeyChain<P,376,43,5>;
std::filesystem::path wide_path, baseline_path, factors_path;

int run(const std::filesystem::path&root){using namespace TFHEpp;using Clock=std::chrono::steady_clock;
 const std::string origin="f70c9f110231eac0b4fffe300cccd9ba89a519809e1cc89c45e8814f8f7e53af";
 const auto path=wide_path;auto wide=std::make_unique<K>();const bool reused=true;
 if(std::filesystem::file_size(path)!=1+8+64+56+sizeof(wide->data))throw std::runtime_error("wide key length mismatch");
 {std::ifstream f(path,std::ios::binary);cereal::PortableBinaryInputArchive ar(f);uint64_t length;ar(length);if(length!=64)throw std::runtime_error("wide key origin length mismatch");std::string stored(64,' ');ar(cereal::binary_data(stored.data(),stored.size()),*wide);if(stored!=origin||f.peek()!=std::char_traits<char>::eof())throw std::runtime_error("wide key mismatch");}
 std::cout<<"public DD433 primary64 saved; reused="<<reused<<std::endl;
 auto wide32=makeCKKSCompactDDRelinKeyWithSecondaryBase<32>(*wide);auto before=Clock::now();CKKSPreparedDDRelinKeyNTT<P,433,64,32> ntt(*wide32);const double preparation=std::chrono::duration<double>(Clock::now()-before).count();
 auto chain=std::make_unique<Chain>();{std::ifstream f(baseline_path,std::ios::binary);cereal::PortableBinaryInputArchive ar(f);std::array<uint64_t,6> meta;std::string stored;auto post=std::make_unique<CKKSCompactDDRelinKey<P,101>>();ar(meta,stored,*chain,*post);if(meta!=std::array<uint64_t,6>{0x4444434841494e31ULL,P::n,448,376,43,5}||stored!=origin||f.peek()!=std::char_traits<char>::eof())throw std::runtime_error("baseline DD mismatch");}
 auto baseline16=makeCKKSCompactDDRelinKeyWithSecondaryBase<16>(chain->get<0>());auto fft=makeCKKSDDRelinKeyFFT(*baseline16);
 auto a=std::make_unique<CKKSCiphertext<P,333,43>>(),b=std::make_unique<CKKSCiphertext<P,333,43>>();{std::ifstream f(factors_path,std::ios::binary);cereal::PortableBinaryInputArchive ar(f);std::array<uint64_t,4> meta;auto prepared=std::make_unique<CKKSCiphertext<P,48,44>>();ar(meta,*prepared,*a,*b);if(meta!=std::array<uint64_t,4>{0x524f4f544e4f4931ULL,P::n,448,36})throw std::runtime_error("factor fixture mismatch");}
 auto tensor=std::make_unique<TRLWE3<P>>();CKKSTensorProductRescale<P,333,333,0>(*tensor,a->ct,b->ct);
 auto lifted=std::make_unique<TRLWE3<P>>();for(size_t c=0;c<3;++c)for(size_t i=0;i<P::n;++i)(*lifted)[c][i]=ckks_detail::reduceToLevel<P,433>((*tensor)[c][i]<<100);
 auto switched=std::make_unique<TRLWE<P>>(),baseline_switched=std::make_unique<TRLWE<P>>();
 ntt.apply(*switched,(*lifted)[2]);CKKSRelinKeySwitch<P>(*baseline_switched,(*tensor)[2],*fft);
 before=Clock::now();for(size_t j=0;j<8;++j)ntt.apply(*switched,(*lifted)[2]);const double ntt_seconds=std::chrono::duration<double>(Clock::now()-before).count()/8;
 before=Clock::now();for(size_t j=0;j<8;++j)CKKSRelinKeySwitch<P>(*baseline_switched,(*tensor)[2],*fft);const double fft_seconds=std::chrono::duration<double>(Clock::now()-before).count()/8;
 auto relin=std::make_unique<CKKSCiphertext<P,433,186>>();for(size_t c=0;c<2;++c)for(size_t i=0;i<P::n;++i)relin->ct[c][i]=ckks_detail::reduceToLevel<P,433>((*lifted)[c][i]+(*switched)[c][i]);
 auto unscaled=std::make_unique<CKKSCiphertext<P,333,86>>();CKKSScaleDown<P,433,186,100>(*unscaled,*relin);auto output=std::make_unique<CKKSCiphertext<P,290,43>>();CKKSScaleDown<P,333,86,43>(*output,*unscaled);
 auto shared=makeCKKSDDRelinKeyNTT(*wide);CKKSAuxiliaryDDRelinKeyNTT<P,333,100> auxiliary(shared);auto api_result=std::make_unique<CKKSCiphertext<P,290,43>>();CKKSMultRelinThenRescale<P>(*api_result,*a,*b,auxiliary);if(api_result->ct!=output->ct)throw std::runtime_error("generic multiplication differs from qualified auxiliary sequence");std::cout<<"generic multiplication API byte-identical to manual auxiliary sequence"<<std::endl;
 auto old=std::make_unique<CKKSCiphertext<P,333,86>>();for(size_t c=0;c<2;++c)for(size_t i=0;i<P::n;++i)old->ct[c][i]=ckks_detail::reduceToLevel<P,333>((*tensor)[c][i]+(*baseline_switched)[c][i]);auto oldscaled=std::make_unique<CKKSCiphertext<P,290,43>>();CKKSScaleDown<P,333,86,43>(*oldscaled,*old);
 {std::ofstream f(root/"outputs.bin",std::ios::binary);cereal::PortableBinaryOutputArchive ar(f);ar(std::array<uint64_t,8>{0x5744444352595031ULL,P::n,448,333,43,433,100,64},*a,*b,*tensor,*relin,*unscaled,*output,*oldscaled);}
 std::ofstream report(root/"report.json");report<<std::setprecision(17)<<"{\"n\":16384,\"private_keys_generated\":0,\"new_public_dd_keys_generated\":"<<(reused?0:1)<<",\"public_dd_archive_reused\":"<<(reused?"true":"false")<<",\"wide_key_log_q\":433,\"primary_bits\":64,\"secondary_bits\":32,\"auxiliary_log_scale\":100,\"ntt_cache_bytes\":"<<ntt.spectra.size()*sizeof(modular_ntt::ShoupMultiplier)<<",\"ntt_preparation_seconds\":"<<preparation<<",\"ntt_switch_seconds_per_call\":"<<ntt_seconds<<",\"baseline_fft_switch_seconds_per_call\":"<<fft_seconds<<",\"switch_speedup\":"<<fft_seconds/ntt_seconds<<"}";
 std::cout<<"wide DD switch seconds="<<ntt_seconds<<" baseline="<<fft_seconds<<" speedup="<<fft_seconds/ntt_seconds<<std::endl;return 0;
}
}
int main(int argc,char**argv)try{
 if(argc!=5)throw std::runtime_error("usage: SAVED_DD433 SAVED_BASELINE_DD_CHAIN SAVED_ROOT_FACTORS NEW_OUTPUT_DIRECTORY");
 wide_path=argv[1];baseline_path=argv[2];factors_path=argv[3];
 if(!std::filesystem::create_directory(argv[4]))throw std::runtime_error("output directory must be new");
 return run(argv[4]);
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}
