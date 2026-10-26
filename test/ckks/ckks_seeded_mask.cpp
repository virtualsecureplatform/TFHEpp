#include "ckks_narrow_basis_benchmark.hpp"
#include <ckks/scheme_switching.hpp>
#include <bfv/compact-cover-bgv.hpp>
#include <sstream>
#include <type_traits>

namespace {
void require(bool ok,const char *message){if(!ok)throw std::runtime_error(message);}
template<class F> void rejects(F action){bool threw=false;try{action();}catch(const std::exception&){threw=true;}require(threw,"incompatible mask archive accepted");}
#ifdef USE_BLAKE3
using P=narrow_benchmark::BenchmarkParam<6,2>;
using Row=TFHEpp::CKKSSeededKeySwitchRow<P,32>;
void stream(){
    static_assert(!std::is_default_constructible_v<BLAKE3PRNG::BLAKE3PRNG<std::uint64_t>>);
    static_assert(!std::is_constructible_v<BLAKE3PRNG::BLAKE3PRNG<std::uint64_t>,std::uint64_t>);
    const TFHEpp::CKKSSeed seed{1,2,3,4};
    const std::string_view domain="TFHEpp.CKKS.seeded-mask.v2";
    TFHEpp::ckks_detail::CKKSSeedStream stream(seed);
    std::array<std::uint8_t,32> bytes{};
    for(std::size_t i=0;i<4;++i)for(std::size_t j=0;j<8;++j)bytes[8*i+j]=std::uint8_t(seed[i]>>(8*j));
    blake3_hasher hash;blake3_hasher_init_keyed(&hash,bytes.data());
    blake3_hasher_update(&hash,domain.data(),domain.size());
    std::vector<std::uint8_t> expected(4*4096);
    blake3_hasher_finalize_seek(&hash,0,expected.data(),expected.size());
    for(std::size_t i=0;i<expected.size();i+=8){
        std::uint64_t word=0;for(std::size_t j=0;j<8;++j)word|=std::uint64_t(expected[i+j])<<(8*j);
        require(stream.next()==word,"seeded mask differs from reference XOF across refill");
    }
    TFHEpp::SeededMaskStream ckks(seed,domain),bfv(seed,"TFHEpp.BFV.cover-mask.v3");
    require(ckks.next()!=bfv.next(),"mask domains are not separated");
    // These distinct seeds collided in the removed 64-bit-state expander.
    TFHEpp::ckks_detail::CKKSSeedStream first(seed),second({5,6,7,UINT64_C(3599532394160262945)});
    require(first.next()!=second.next(),"historical seed collision persists");
    for(std::size_t i=0;i<4;++i){auto other=seed;other[i]^=UINT64_C(1)<<63;TFHEpp::ckks_detail::CKKSSeedStream a(seed),b(other);require(a.next()!=b.next(),"seed word not incorporated");}
}
void archives(){
    using namespace TFHEpp;
    auto row=std::make_unique<Row>();Key<P> secret{};Polynomial<P> message{},expanded{},phase{};
    for(std::size_t i=0;i<P::n;++i){secret[i]=P::T::from_signed_i64(int(i%3)-1);message[i]=P::T{std::uint32_t(i*8191)};}
    ckks_scheme_switch_detail::exactEncryptPolynomialAtLevel<P,32>(*row,message,secret,{0,0});
    ckks_detail::fillSeededUniformPolynomialAtLevel<P,32>(expanded,row->mask_seeds[0]);
    ckks_detail::exactActiveTorusMulBySmall<P,32>(phase,expanded,secret);
    for(std::size_t i=0;i<P::n;++i)require(ckks_detail::reduceToLevel<P,32>(row->body[i]-phase[i])==message[i],"new seeded encryption phase mismatch");
    std::stringstream saved;{cereal::PortableBinaryOutputArchive ar(saved);ar(*row);}
    Row loaded;{cereal::PortableBinaryInputArchive ar(saved);ar(loaded);}
    require(loaded.mask_seeds==row->mask_seeds&&loaded.body==row->body,"seeded archive roundtrip");
    std::stringstream legacy;{cereal::PortableBinaryOutputArchive ar(legacy);ar(row->mask_seeds,row->body);}
    rejects([&]{cereal::PortableBinaryInputArchive ar(legacy);ar(loaded);});
    std::stringstream wrong;{cereal::PortableBinaryOutputArchive ar(wrong);std::array<std::uint64_t,5> header{0x434b534d41534b32ULL,P::n,1,31,1};ar(header,row->mask_seeds,row->body);}
    rejects([&]{cereal::PortableBinaryInputArchive ar(wrong);ar(loaded);});
    namespace cc=compact_cover_bgv;
    const auto path=std::filesystem::temp_directory_path()/("legacy-seeded-mask-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    cc::TransitionKeyHeader old;old.reserved=1;old.width=1;old.limbs=1;old.rows=1;
    {std::ofstream file(path,std::ios::binary);file.write(reinterpret_cast<const char*>(&old),sizeof(old));}
    try{rejects([&]{cc::SeededTransitionKeyProvider reader(path);});}
    catch(...){std::filesystem::remove(path);throw;}
    std::filesystem::remove(path);
}
#endif
}
int main()try{
#ifdef USE_BLAKE3
    stream();archives();
#else
    rejects([]{TFHEpp::SeededMaskStream stream({1,2,3,4},"test");});
#endif
    std::cout<<"Cryptographic mask streams and legacy archive rejection passed\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
