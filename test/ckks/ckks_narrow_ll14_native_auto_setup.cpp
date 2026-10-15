#include "ckks_narrow_basis_benchmark.hpp"
#include <ckks/native_auto_file.hpp>
#include <iomanip>
#include <set>
using namespace TFHEpp;
using P=narrow_benchmark::BenchmarkParam<14,7>;
using Schedule=CKKSNarrowBootstrapSchedule<P,48,427,44,42,31,true,13,true>;
int main(int argc,char**argv)try{
 if(argc!=3)throw std::runtime_error("usage: checkpoint-directory native-auto-directory");
 const std::filesystem::path directory=argv[2];std::filesystem::create_directories(directory);
 const std::string origin="f70c9f110231eac0b4fffe300cccd9ba89a519809e1cc89c45e8814f8f7e53af";
 const double sigma=3.2;auto dense=std::make_unique<Key<P>>();
 {std::ifstream file(std::filesystem::path(argv[1])/"checkpoint.bin",std::ios::binary);
  if(!file)throw std::runtime_error("retained checkpoint missing");
  cereal::PortableBinaryInputArchive archive(file);std::array<uint64_t,9> metadata{};std::string function;
  archive(metadata,function);
  if(metadata!=std::array<uint64_t,9>{0x4e415252504b3031ULL,P::n,85,41,379,49,34,100,448}||function!="tanh")throw std::runtime_error("retained checkpoint profile mismatch");
  archive(*dense);
 }
 std::set<uint32_t> exponents{2*P::n-1};
 for(std::size_t ell=0;ell<Schedule::tail_blocks_per_pack;++ell)exponents.insert(ckks_detail::rotationAutomorphismForSteps<P>(2*(ell+1)*Schedule::tail_width));
 const auto started=std::chrono::steady_clock::now();std::size_t generated=0,reused=0;
 for(auto exponent:exponents){
  const auto target=directory/(std::to_string(exponent)+".bin");
  if(std::filesystem::exists(target)){
   std::ifstream file(target,std::ios::binary);CKKSReadNativeAutoRows<P,435,44>(file,exponent,origin,sigma);++reused;
  }else{
   const auto pending=std::filesystem::path(target.string()+".pending");
   if(!std::filesystem::create_directory(pending))throw std::runtime_error("native auto staging directory already exists");
   struct Cleanup{std::filesystem::path path;~Cleanup(){std::error_code error;std::filesystem::remove_all(path,error);}} cleanup{pending};
   std::filesystem::permissions(pending,std::filesystem::perms::owner_all,std::filesystem::perm_options::replace);
   auto rows=std::make_unique<PublicAutoGadgetRows<P,435,44>>();CKKSGenerateNativeAutoGadgetRows(*rows,exponent,*dense);
   {std::ofstream file(pending/"payload.bin",std::ios::binary);CKKSWriteNativeAutoRows(file,*rows,origin,sigma);}
   if(std::filesystem::exists(target))throw std::runtime_error("native auto archive appeared during setup");
   std::filesystem::rename(pending/"payload.bin",target);++generated;
  }
  std::cout<<"native auto exponent="<<exponent<<" generated="<<generated<<" reused="<<reused<<std::endl;
 }
 std::ofstream report(directory/"setup.json");report<<std::setprecision(17)<<"{\"private_keys_generated\":0,\"private_key_reused\":true,\"evaluation_key_sets_generated\":"<<generated<<",\"encrypted_gadget_rows_generated\":"<<generated*PublicAutoGadgetRows<P,435,44>::rows<<",\"evaluation_key_sets_reused\":"<<reused<<",\"log_q\":435,\"primary_bits\":44,\"noise_stddev\":3.2,\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<",\"security_validated\":false}";
}catch(const std::exception&error){std::cerr<<error.what()<<'\n';return 1;}
