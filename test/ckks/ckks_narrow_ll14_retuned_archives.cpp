#include "ckks_narrow_basis_benchmark.hpp"
#include <ckks/narrow_ll14_noise_profile.hpp>
#include <ckks/native_auto_file.hpp>
#include <sstream>
using namespace TFHEpp;
using P=narrow_benchmark::BenchmarkParam<9,7>;
constexpr auto native_bits=CKKSNarrowLL14RetunedNoise::native_auto_primary_bits;
using Legacy=CKKSNarrowBootstrapSchedule<P,48,440,57,42,31,false,0,true>;
struct Retuned:Legacy {static constexpr uint64_t selector_archive_magic=CKKSNarrowLL14RetunedNoise::selector_magic;};
template<class F> void rejects(F action) {
 bool rejected=false;try{action();}catch(const std::runtime_error&){rejected=true;}
 if(!rejected)throw std::runtime_error("incompatible noise archive accepted");
}
template<class S,class Other> void checkSelectors(const std::filesystem::path&path) {
 using File=CKKSNarrowSelectorFile<S>;File owner(path);
 std::vector<typename S::Selector> values(S::block_width);values[0].ct[1][0]=42;
 for(bool tail:{false,true})for(size_t g=0;g<(tail?S::tail_packs:S::blocks);++g)owner.write(tail,g,values);
 owner.finish();File reopened(path,typename File::Existing{});
 std::vector<typename S::Selector> loaded;reopened.loadChunk(false,0,0,1,loaded);
 if(loaded.size()!=1||loaded[0].ct[1][0]!=P::T{42})throw std::runtime_error("selector round-trip failed");
 rejects([&]{CKKSNarrowSelectorFile<Other> wrong(path,typename CKKSNarrowSelectorFile<Other>::Existing{});});
}
int main()try {
 const auto root=std::filesystem::temp_directory_path()/("ckks-noise-archives-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
 if(!std::filesystem::create_directory(root))throw std::runtime_error("test directory exists");
 struct Cleanup{std::filesystem::path path;~Cleanup(){std::error_code e;std::filesystem::remove_all(path,e);}} cleanup{root};
 checkSelectors<Legacy,Retuned>(root/"legacy");checkSelectors<Retuned,Legacy>(root/"retuned");
 auto rows=std::make_unique<PublicAutoGadgetRows<P,435,native_bits>>();rows->exponent=2*P::n-1;rows->ciphertexts[0][1][0]=123;
 std::stringstream stream(std::ios::in|std::ios::out|std::ios::binary);const std::string origin(64,'0');
 CKKSWriteNativeAutoRows(stream,*rows,origin,4096.0);stream.seekg(0);
 auto loaded=CKKSReadNativeAutoRows<P,435,native_bits>(stream,rows->exponent,origin,4096.0);
 if(loaded->ciphertexts[0][1][0]!=P::T{123})throw std::runtime_error("native archive round-trip failed");
 stream.clear();stream.seekg(0);rejects([&]{CKKSReadNativeAutoRows<P,435,native_bits>(stream,rows->exponent,origin,3.2);});
 std::cout<<"legacy/retuned selector separation, chunk round-trip, native noise mismatch: passed; no key generation\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
