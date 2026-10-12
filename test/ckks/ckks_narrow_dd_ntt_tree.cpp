#include "ckks_narrow_basis_benchmark.hpp"
#include <ckks/prepared_dd_ntt_chain.hpp>
#include <iomanip>
using namespace TFHEpp;
using P=narrow_benchmark::BenchmarkParam<14,7>;
struct Schedule : CKKSNarrowBootstrapSchedule<P,48,373,43,40,31,true,12,true> {
 static constexpr uint32_t selector_matrix_digit_bits=56;
};
int main(int argc,char **argv)try {
 if(argc!=5)throw std::runtime_error("usage: SAVED_DD433 SAVED_DD_CHAIN SAVED_FACTORS NEW_OUTPUT_DIRECTORY");
 const std::filesystem::path root=argv[4];
 if(!std::filesystem::create_directory(root))throw std::runtime_error("output directory must be new");
 const std::string origin="f70c9f110231eac0b4fffe300cccd9ba89a519809e1cc89c45e8814f8f7e53af";
 using Master=CKKSCompactDDRelinKey<P,433,64,16>;
 auto master=std::make_unique<Master>();
 if(std::filesystem::file_size(argv[1])!=1+8+64+56+sizeof(master->data))throw std::runtime_error("master key length mismatch");
 {std::ifstream f(argv[1],std::ios::binary);cereal::PortableBinaryInputArchive ar(f);uint64_t length;ar(length);if(length!=64)throw std::runtime_error("master origin length mismatch");std::string stored(64,' ');ar(cereal::binary_data(stored.data(),stored.size()),*master);if(stored!=origin || f.peek()!=std::char_traits<char>::eof())throw std::runtime_error("master origin/length mismatch");}
 auto keys=std::make_unique<Schedule::RelinKeys>();
 {std::ifstream f(argv[2],std::ios::binary);cereal::PortableBinaryInputArchive ar(f);std::array<uint64_t,6> meta;std::string stored;auto post=std::make_unique<CKKSCompactDDRelinKey<P,101>>();ar(meta,stored,*keys,*post);if(meta!=std::array<uint64_t,6>{0x4444434841494e31ULL,P::n,448,376,43,5} || stored!=origin || f.peek()!=std::char_traits<char>::eof())throw std::runtime_error("chain profile mismatch");}
 std::vector<Schedule::Factor> factors(Schedule::blocks+1);
 {std::ifstream f(argv[3],std::ios::binary);cereal::PortableBinaryInputArchive ar(f);std::array<uint64_t,4> meta;auto prepared=std::make_unique<CKKSCiphertext<P,48,44>>();ar(meta,*prepared);if(meta!=std::array<uint64_t,4>{0x524f4f544e4f4931ULL,P::n,448,36})throw std::runtime_error("factor profile mismatch");for(size_t i=0;i<Schedule::blocks;++i)ar(factors[i]);}
 auto one=std::make_unique<CKKSSlotVector<P>>();one->fill({1,0});ckksSlotEncode<P,333,43>(factors.back().ct[P::k],*one);
 const auto start=std::chrono::steady_clock::now();
 CKKSPreparedAuxiliaryDDRelinKeyChainNTT<P,333,43,5> ntt(*master);
 const double setup=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 auto fft=CKKSNarrowMakePreparedRelinearization<Schedule>(*keys);
 auto result=std::make_unique<Schedule::TreeResult>(),baseline=std::make_unique<Schedule::TreeResult>();
 auto t=std::chrono::steady_clock::now();narrow_detail::productTree<Schedule>(*result,factors,*keys,&ntt);
 const double ntt_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();
 auto parallel=std::make_unique<Schedule::TreeResult>();
 t=std::chrono::steady_clock::now();
 narrow_detail::productTree<Schedule,0,decltype(ntt),true>(*parallel,factors,*keys,&ntt);
 const double parallel_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();
 if(parallel->ct!=result->ct)throw std::runtime_error("parallel NTT tree differs from serial NTT tree");
 std::cout<<"parallel tree byte-identical; serial="<<ntt_seconds<<" parallel="<<parallel_seconds<<std::endl;
 t=std::chrono::steady_clock::now();narrow_detail::productTree<Schedule>(*baseline,factors,*keys,fft.get());
 const double fft_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();
 {std::ofstream f(root/"outputs.bin",std::ios::binary);cereal::PortableBinaryOutputArchive ar(f);ar(std::array<uint64_t,8>{0x44444e5454524531ULL,P::n,448,333,43,118,5,factors.size()});for(const auto &factor:factors)ar(factor);ar(*result,*baseline);}
 {std::ofstream f(root/"report.json");f<<std::setprecision(17)<<"{\"new_private_keys_generated\":0,\"new_encrypted_evaluation_keys_generated\":0,\"ntt_chain_cache_bytes\":"<<ntt.byteSize()<<",\"ntt_chain_preparation_seconds\":"<<setup<<",\"ntt_tree_seconds\":"<<ntt_seconds<<",\"fft_tree_seconds\":"<<fft_seconds<<",\"speedup\":"<<fft_seconds/ntt_seconds<<"}";}
 std::cout<<"five-level tree NTT="<<ntt_seconds<<" FFT="<<fft_seconds<<" speedup="<<fft_seconds/ntt_seconds<<std::endl;
}catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
