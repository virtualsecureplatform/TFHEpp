#include <modular_rns32.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <chrono>
#include <iostream>
#include <random>

using namespace TFHEpp;
using Integer=boost::multiprecision::cpp_int;
using T=MultiLimbUInt<7>;
Integer integer(const T&x){Integer v=0;for(std::size_t j=7;j-->0;){v<<=64;v+=x.limb[j];}return v;}
int main(){
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
    const auto start=std::chrono::steady_clock::now();crt.reconstruct(spans,output);
    const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    for(std::size_t k=0;k<n;++k)if(integer(output[k])!=expected[k])throw std::runtime_error("centered wide CRT mismatch");
    std::cout<<"prime/root/encode/CRT checks="<<checks+n<<" CRT seconds="<<elapsed<<" all passed\n";
}
