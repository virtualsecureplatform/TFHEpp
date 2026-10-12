#include <ckks/narrow_refresh.hpp>
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
    static constexpr std::uint32_t nbit = 9, n = 1 << nbit, k = 1;
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
struct FullP : P {
    static constexpr std::uint32_t nbit = 14, n = 1 << nbit;
};
using Schedule = TFHEpp::CKKSNarrowBootstrapSchedule<P, 48, 379, 49, 34, 31>;
using Packed = TFHEpp::CKKSNarrowPackedKey<Schedule>;
constexpr std::array<std::uint64_t, 9> metadata{
    0x5041434b43463031ULL, P::n, 379, 49, 34, 31, 448, 4, 4};
static_assert(Packed::log_q == 422);
using Clock = std::chrono::steady_clock;
double elapsed(Clock::time_point start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}


void run(const std::filesystem::path &root) {
 using namespace TFHEpp;
 auto packed=std::make_unique<Packed>();
 std::ifstream f(root/"packed.bin",std::ios::binary);
 {cereal::PortableBinaryInputArchive a(f);std::array<uint64_t,9> m{};a(m,*packed);if(m!=metadata)throw std::runtime_error("profile mismatch");}
 CKKSNarrowPreparePackedSubstitutions<Schedule>(*packed,P::nbit);
 std::ifstream reference(root/"terms.bin",std::ios::binary);cereal::PortableBinaryInputArchive ar(reference);
 struct Record {bool tail;uint64_t index;Schedule::Selector term;};
 std::vector<Record> records(4*P::n);
 for(auto &r:records)ar(r.tail,r.index,r.term);
 if(reference.peek()!=std::char_traits<char>::eof())throw std::runtime_error("reference trailing data");
 size_t checked=0;
 for(const auto [begin,end]:std::array<std::pair<size_t,size_t>,5>{{{0,1},{17,33},{P::n-1,P::n},{99,99},{0,P::n}}}) {
  std::vector<const Record*> expected;
  for(size_t j=0;j<records.size();++j){size_t i=j%P::n,x=0;for(size_t bit=0;bit<P::nbit;++bit)x=(x<<1)|((i>>bit)&1);if(x>=begin&&x<end)expected.push_back(&records[j]);}
  size_t count=0;auto start=Clock::now();
  CKKSNarrowExtractSelectorTerms<Schedule>(*packed,[&](bool tail,size_t index,const auto &term){
   if(count>=expected.size())throw std::runtime_error("unexpected term");
   const auto &r=*expected[count++];if(r.tail!=tail||r.index!=index||r.term.ct!=term.ct)throw std::runtime_error("bounded extraction differs from full reference");
  },begin,end);
  if(count!=expected.size())throw std::runtime_error("missing terms");checked+=count;
  std::cout<<"range "<<begin<<' '<<end<<" terms "<<count<<" seconds "<<elapsed(start)<<std::endl;
 }
 size_t rejected=0;
 for(const auto [begin,end]:std::array<std::pair<size_t,size_t>,2>{{{9,8},{0,P::n+1}}}){
  try{CKKSNarrowExtractSelectorTerms<Schedule>(*packed,[](bool,size_t,const auto&){throw std::runtime_error("invalid range evaluated");},begin,end);}
  catch(const std::invalid_argument&){++rejected;}
 }
 if(rejected!=2)throw std::runtime_error("range validation failed");
 std::cout<<"{\"n\":512,\"fresh_keys_generated\":false,\"terms_byte_identical\":true,\"terms_checked\":"<<checked<<",\"invalid_ranges_rejected\":"<<rejected<<"}"<<std::endl;
}
}
int main(int argc,char **argv)try{if(argc!=2)throw std::runtime_error("usage: packed-range FIXTURE");run(argv[1]);}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}
