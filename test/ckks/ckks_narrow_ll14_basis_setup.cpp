#include "ckks_narrow_basis_benchmark.hpp"
#include <ckks/prepared_dd_ntt.hpp>
#include <iomanip>
using namespace TFHEpp;
using P=narrow_benchmark::BenchmarkParam<14,7>;
struct EncryptionP:P { static constexpr uint32_t B̅gbit=16,l̅=28; };
using Direct=CKKSNarrowBootstrapSchedule<P,48,440,57,42,31,false,0,true>;
int main(int argc,char**argv)try {
 if(argc!=4)throw std::runtime_error("CHECKPOINT NEW_DIRECTORY pilot|full");
 const std::filesystem::path root=argv[2];const bool full=std::string(argv[3])=="full";
 if(!full && std::string(argv[3])!="pilot")throw std::runtime_error("invalid setup mode");
 if(!std::filesystem::create_directory(root))throw std::runtime_error("setup directory must be new");
 auto dense=std::make_unique<Key<P>>(),sparse=std::make_unique<Key<P>>();
 const std::string origin="f70c9f110231eac0b4fffe300cccd9ba89a519809e1cc89c45e8814f8f7e53af";
 {std::ifstream f(std::filesystem::path(argv[1])/"checkpoint.bin",std::ios::binary);if(!f)throw std::runtime_error("checkpoint missing");
  cereal::PortableBinaryInputArchive ar(f);std::array<uint64_t,9> meta;std::string function;ar(meta,function);
  if(meta!=std::array<uint64_t,9>{0x4e415252504b3031ULL,P::n,85,41,379,49,34,100,448}||function!="tanh")throw std::runtime_error("checkpoint mismatch");
  ar(*dense,*sparse);}
 const auto positions=narrow_detail::structuredPositions<Direct>(*sparse);
 const auto started=std::chrono::steady_clock::now();double generation_seconds=0,write_seconds=0;
 for(bool tail:{false,true})for(size_t group=0;group<(full?(tail?Direct::tail_packs:Direct::blocks):(tail?0:1));++group) {
  auto chunk=std::make_unique<std::vector<Direct::Selector>>(Direct::block_width);
  const auto before=std::chrono::steady_clock::now();
#pragma omp parallel for schedule(static)
  for(size_t y=0;y<Direct::block_width;++y) {
   auto slots=std::make_unique<CKKSSlotVector<P>>();bool nonzero=false;
   const auto indicator=[&](size_t x) {
    if(y>=Direct::block_size(x)||positions[x]!=int(y))return std::complex<double>{};
    nonzero=true;return (*sparse)[x*Direct::block_width+y]==P::T{1}?std::complex<double>{1,0}:std::complex<double>{0,1};
   };
   if(!tail)std::fill_n(slots->begin(),Direct::head_slots,indicator(group));
   else for(size_t ell=0;ell<Direct::tail_blocks_per_pack;++ell) {
    const auto x=group*Direct::tail_blocks_per_pack+ell;if(x>=Direct::blocks)break;
    std::fill_n(slots->begin()+(2*ell+1)*Direct::tail_width,Direct::tail_width,indicator(x));
   }
   auto polynomial=std::make_unique<Polynomial<P>>();
   if(nonzero)ckksSlotEncode<P,440,57>(*polynomial,*slots);
   ckks_detail::encryptPolynomialAtLevel<EncryptionP,440>((*chunk)[y].ct,*polynomial,*dense,{std::ldexp(102.4,-440),0});
  }
  const auto generated=std::chrono::steady_clock::now();generation_seconds+=std::chrono::duration<double>(generated-before).count();
  {std::ofstream f(root/((tail?"tail-":"head-")+std::to_string(group)+".bin"),std::ios::binary);cereal::PortableBinaryOutputArchive ar(f);
   std::array<uint64_t,9> meta{0x4e4152524f573031ULL,P::n,440,57,42,Direct::block_width,448,tail,group};ar(meta);
   ar(cereal::binary_data(reinterpret_cast<const uint64_t*>(chunk->data()),chunk->size()*sizeof(Direct::Selector)));}
  write_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-generated).count();
  std::cout<<"saved "<<(tail?"tail":"head")<<" group="<<group<<" generation_seconds="<<generation_seconds<<" write_seconds="<<write_seconds<<std::endl;
  if(!full) {
   auto slots=std::make_unique<CKKSSlotVector<P>>();double error=0;
   for(size_t y:{size_t{0},size_t(positions[0]),Direct::block_width-1}) {
    ckksSlotDecrypt<P>(*slots,(*chunk)[y],*dense);
    const auto value=positions[0]==int(y)?((*sparse)[y]==P::T{1}?std::complex<double>{1,0}:std::complex<double>{0,1}):std::complex<double>{};
    for(size_t i=0;i<P::n/2;++i)error=std::max(error,std::abs((*slots)[i]-(i<Direct::head_slots?value:std::complex<double>{})));
   }
   std::cout<<"pilot maximum selector error="<<std::setprecision(17)<<error<<std::endl;
   if(error>=std::ldexp(1.0,-23))throw std::runtime_error("pilot selector accuracy failed");
  }
 }
 auto dd=std::make_unique<CKKSCompactDDRelinKey<P,433,20,20>>();
 CKKSCompactDDRelinKeyGenForNTT(*dd,*dense);
 {std::ofstream f(root/"dd433-primary20.bin",std::ios::binary);cereal::PortableBinaryOutputArchive ar(f);ar(origin,*dd);}
 std::ofstream report(root/"setup.json");report<<std::setprecision(17)<<"{\"full\":"<<(full?"true":"false")<<",\"private_keys_generated\":0,\"private_keys_reused\":true,\"selector_noise_stddev\":102.4,\"selector_log_q\":440,\"selector_log_delta\":57,\"dd_log_q\":433,\"dd_primary_bits\":20,\"dd_secondary_bits\":20,\"generation_seconds\":"<<generation_seconds<<",\"write_seconds\":"<<write_seconds<<",\"total_setup_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<"}";
 std::cout<<"setup complete"<<std::endl;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
