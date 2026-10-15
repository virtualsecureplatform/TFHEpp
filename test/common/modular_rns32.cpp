#include <modular_rns32.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <chrono>
#include <iostream>
#include <random>

using namespace TFHEpp;
using Integer=boost::multiprecision::cpp_int;
using T=MultiLimbUInt<7>;
template<std::size_t L> Integer integer(const MultiLimbUInt<L>&x){Integer v=0;for(std::size_t j=L;j-->0;){v<<=64;v+=x.limb[j];}return v;}
template<std::size_t Count>void checkPrefix(){
    using Small=MultiLimbUInt<2>;
    constexpr std::size_t n=521;
    std::mt19937_64 random(20261015+Count);Integer modulus=1;
    for(std::size_t p=0;p<Count;++p)modulus*=modular_rns32::primes[p].value;
    const Integer half=modulus/2,mask=(Integer{1}<<128)-1;
    std::array<std::vector<std::uint32_t>,Count> rows;
    for(auto &row:rows)row.resize(n);
    std::vector<Integer> expected(n);
    for(std::size_t k=0;k<n;++k){
        Integer v=0;for(int j=0;j<4;++j){v<<=64;v+=random();}v%=modulus;if(v>half)v-=modulus;
        if(k<7)v=std::array<Integer,7>{0,1,-1,half,half-1,-half,-half+1}[k];
        expected[k]=v&mask;
        for(std::size_t p=0;p<Count;++p){Integer r=v%modular_rns32::primes[p].value;if(r<0)r+=modular_rns32::primes[p].value;rows[p][k]=r.convert_to<std::uint32_t>();}
    }
    std::array<std::span<std::uint32_t>,Count> spans;for(std::size_t p=0;p<Count;++p)spans[p]=rows[p];
    std::vector<Small> output(n);modular_rns32::CenteredCRT<2,Count> crt;crt.reconstruct(spans,output);
    for(std::size_t k=0;k<n;++k)if(integer(output[k])!=expected[k])throw std::runtime_error("prefix centered CRT mismatch");
    std::cout<<"prefix primes="<<Count<<" centered CRT values="<<n<<" passed\n";
}
template<std::size_t Limbs,std::size_t Count>void checkDirect(){
    using Output=MultiLimbUInt<Limbs>;
    constexpr std::size_t n=521;
    std::mt19937_64 random(20261018+Limbs*101+Count);
    Integer modulus=1;for(std::size_t p=0;p<Count;++p)modulus*=modular_rns32::primes[p].value;
    const Integer half=modulus/2,mask=(Integer{1}<<(64*Limbs))-1;
    std::array<std::vector<std::uint32_t>,Count> rows;
    for(auto&row:rows)row.resize(n);
    std::vector<Integer> expected(n);
    for(std::size_t k=0;k<n;++k){
        Integer v=0;for(int j=0;j<8;++j){v<<=64;v+=random();}v%=modulus;if(v>half)v-=modulus;
        if(k<7)v=std::array<Integer,7>{0,1,-1,half,half-1,-half,-half+1}[k];
        if(k>=7&&k<47){const auto d=(k-7)/2;if(k&1)v=half-d;else v=-half+d;}
        if(k>=100&&k<100+Count)v=-modulus/modular_rns32::primes[k-100].value;
        if(k==255)v=half;
        if(k==256)v=-half;
        expected[k]=v&mask;
        for(std::size_t p=0;p<Count;++p){Integer r=v%modular_rns32::primes[p].value;if(r<0)r+=modular_rns32::primes[p].value;rows[p][k]=r.convert_to<std::uint32_t>();}
    }
    std::array<std::span<std::uint32_t>,Count> spans;for(std::size_t p=0;p<Count;++p)spans[p]=rows[p];
    std::vector<Output> output(n);modular_rns32::DirectCenteredCRT<Limbs,Count> crt;
    const auto fallbacks=crt.reconstruct(spans,output);
    for(std::size_t k=0;k<n;++k)if(integer(output[k])!=expected[k])throw std::runtime_error("direct quotient CRT differs from arbitrary precision oracle");
    if constexpr(Count>1)if(!fallbacks)throw std::runtime_error("centering endpoints did not exercise exact fallback");
    std::array<std::span<std::uint32_t>,Count> empty;
    if(crt.reconstruct(empty,{})!=0)throw std::runtime_error("empty CRT reconstruction required fallback");
    bool rejected=false;spans[0]=spans[0].first(n-1);
    try{crt.reconstruct(spans,output);}catch(const std::invalid_argument&){rejected=true;}
    if(!rejected)throw std::runtime_error("direct CRT accepted inconsistent residue lengths");
    std::cout<<"direct limbs="<<Limbs<<" primes="<<Count<<" values="<<n<<" fallbacks="<<fallbacks<<" passed\n";
}
template<std::size_t Count>void checkDirectWidths(){
    checkDirect<1,Count>();checkDirect<3,Count>();checkDirect<7,Count>();checkDirect<9,Count>();
}
int main(){
    checkPrefix<5>();checkPrefix<6>();
    []<std::size_t... I>(std::index_sequence<I...>){(checkDirect<2,I+1>(),...);}(std::make_index_sequence<17>{});
    checkDirectWidths<4>();checkDirectWidths<6>();checkDirectWidths<7>();checkDirectWidths<9>();
    checkDirectWidths<12>();checkDirectWidths<14>();checkDirectWidths<17>();
    std::mt19937_64 random(20261015);
    Integer modulus=1;std::size_t checks=0;
    modular_rns32::CanonicalEncoder<7,440> encoder;
    for(std::size_t i=0;i<modular_rns32::prime_count;++i){
        const auto p=modular_rns32::primes[i].value;modulus*=p;
        for(std::uint64_t d=2;d*d<=p;++d)if(p%d==0)throw std::runtime_error("composite RNS modulus");
        const auto root=modular_ntt::power(modular_rns32::primes[i].primitive_root,(p-1)/131072,p);
        if(modular_ntt::power(root,65536,p)!=p-1)throw std::runtime_error("RNS root order mismatch");
        for(std::size_t j=0;j<1024;++j){
            const auto a=std::uint32_t(random()),b=std::uint32_t(random()%p);
            if(modular_rns32::Multiplier(b,p).apply(a,p)!=(std::uint64_t(a)*b)%p)
                throw std::runtime_error("32-bit Shoup mismatch");
            T input;for(auto &word:input.limb)word=random();
            const Integer active=integer(input)&((Integer{1}<<440)-1);
            if(encoder.residue(input,i)!=(active%p).convert_to<std::uint32_t>())
                throw std::runtime_error("full selector RNS encode mismatch");
            checks+=2;
        }
    }
    constexpr std::size_t n=16389;
    std::array<std::vector<std::uint32_t>,17> rows;
    for(auto&row:rows)row.resize(n);
    std::vector<Integer> expected(n);
    const Integer half=modulus/2,mask=(Integer{1}<<448)-1;
    for(std::size_t k=0;k<n;++k){
        Integer value=0;for(int j=0;j<8;++j){value<<=64;value+=random();}
        value%=modulus;if(value>half)value-=modulus;
        if(k<7)value=std::array<Integer,7>{0,1,-1,half,half-1,-half,-half+1}[k];
        expected[k]=value&mask;
        for(std::size_t i=0;i<17;++i){const auto p=modular_rns32::primes[i].value;Integer r=value%p;if(r<0)r+=p;rows[i][k]=r.convert_to<std::uint32_t>();}
    }
    std::array<std::span<std::uint32_t>,17> spans;
    for(std::size_t i=0;i<17;++i)spans[i]=rows[i];
    std::vector<T> output(n);modular_rns32::CenteredCRT<7> crt;
    const auto saved=rows;
    const auto start=std::chrono::steady_clock::now();crt.reconstruct(spans,output);
    const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    for(std::size_t k=0;k<n;++k)if(integer(output[k])!=expected[k])throw std::runtime_error("centered wide CRT mismatch");
    modular_rns32::DirectCenteredCRT<7> direct;
    const auto expected_output=output;
    for(bool fast:{false,true,true,false})for(int trial=0;trial<5;++trial){
        rows=saved;
        for(std::size_t p=0;p<17;++p)spans[p]=rows[p];
        const auto begin=std::chrono::steady_clock::now();
        std::size_t fallback=0;
        if(fast)fallback=direct.reconstruct(spans,output);else crt.reconstruct(spans,output);
        const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
        if(output!=expected_output)throw std::runtime_error("timed direct CRT mismatch");
        std::cout<<"CRT variant="<<(fast?"direct":"garner")<<" trial="<<trial<<" seconds="<<seconds<<" fallbacks="<<fallback<<" identical\n";
    }
    std::cout<<"prime/root/encode/CRT checks="<<checks+n<<" CRT seconds="<<elapsed<<" all passed\n";
}
