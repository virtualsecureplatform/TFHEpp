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
template<std::size_t Bits,std::size_t SourceQ=448,std::size_t Shift=0>
void checkBatchEncoding(){
    constexpr std::size_t n=521,active=std::min(Bits,SourceQ-Shift);
    std::mt19937_64 random(20261019+Bits*1009+Shift);
    std::vector<T> input(n+2);
    for(auto&v:input)for(auto&w:v.limb)w=random();
    input[1]={};input[2]=T{0}-T{1};input[3]=T{1}<<447;
    modular_rns32::CanonicalEncoder<7,Bits> encoder;
    const Integer mask=(Integer{1}<<active)-1;
    for(std::size_t p=0;p<17;++p){
        const auto modulus=modular_rns32::primes[p].value;
        modular_ntt32::Montgomery arithmetic(modulus);
        std::vector<std::uint32_t> expected(n);
        for(std::size_t i=0;i<n;++i)expected[i]=(((integer(input[i+1])>>Shift)&mask)%modulus).convert_to<std::uint32_t>();
        for(std::size_t length:{0,1,7,8,9,255,256,257,521}){
            std::vector<std::uint32_t> ordinary(length+2,0xdeadbeef),montgomery=ordinary;
            encoder.template residues<SourceQ,Shift>(std::span(ordinary).subspan(1,length),std::span<const T>(input).subspan(1,length),p);
            encoder.template residues<SourceQ,Shift,true>(std::span(montgomery).subspan(1,length),std::span<const T>(input).subspan(1,length),p);
            for(std::size_t i=0;i<length;++i)
                if(ordinary[i+1]!=expected[i] || montgomery[i+1]!=arithmetic.encode(expected[i]))
                    throw std::runtime_error("batched shifted RNS encoding differs from integer oracle");
            if(ordinary.front()!=0xdeadbeef || ordinary.back()!=0xdeadbeef || montgomery.front()!=0xdeadbeef || montgomery.back()!=0xdeadbeef)
                throw std::runtime_error("batched RNS encoding overwrote guard coefficients");
        }
    }
    std::vector<std::uint32_t> output(n);bool rejected=false;
    try{encoder.template residues<SourceQ,Shift>(output,std::span<const T>(input).first(n-1),0);}catch(const std::invalid_argument&){rejected=true;}
    if(!rejected)throw std::runtime_error("batched RNS encoding accepted inconsistent lengths");
    rejected=false;
    try{encoder.template residues<SourceQ,Shift>(output,std::span<const T>(input).first(n),17);}catch(const std::out_of_range&){rejected=true;}
    if(!rejected)throw std::runtime_error("batched RNS encoding accepted invalid prime");
    std::cout<<"batch bits="<<Bits<<" source="<<SourceQ<<" shift="<<Shift<<" all primes/domains/tails passed\n";
}
template<std::size_t Count>void checkSignedEncoding()
{
    constexpr std::size_t n=521;
    std::mt19937_64 random(20261021+Count);
    modular_rns32::SignedEncoder<Count> encoder;
    std::vector<std::int64_t> input(n+2);
    for(unsigned bits:{20,40,44,64}){
        for(auto&v:input){
            auto word=random();
            if(bits<64)word=(word&((std::uint64_t{1}<<bits)-1))-(std::uint64_t{1}<<(bits-1));
            v=std::bit_cast<std::int64_t>(word);
        }
        for(std::size_t p=0;p<Count;++p){
            const auto prime=std::int64_t(modular_rns32::primes[p].value);
            const std::int64_t edges[]={0,1,-1,INT64_MAX,INT64_MIN,prime,-prime,2*prime,-2*prime,(prime<<32)-1,-((prime<<32)-1)};
            for(std::size_t i=0;i<128;++i)input[i+1]=edges[i%std::size(edges)];
            for(std::size_t length:{0,1,7,8,9,255,256,257,521}){
                std::vector<std::uint32_t> expected(length+2,0xdeadbeef);
                for(std::size_t i=0;i<length;++i){Integer value=input[i+1];Integer r=value%prime;if(r<0)r+=prime;expected[i+1]=r.convert_to<std::uint32_t>();}
                for(bool batch:{false,true}){
                    auto actual=expected;std::fill(actual.begin()+1,actual.end()-1,0);
                    encoder.residues(std::span(actual).subspan(1,length),std::span<const std::int64_t>(input).subspan(1,length),p,batch);
                    if(actual!=expected)throw std::runtime_error("signed RNS residues differ from arbitrary precision oracle or overwrite guards");
                }
            }
        }
    }
    bool rejected=false;std::vector<std::uint32_t> output(n);
    try{encoder.residues(output,std::span<const std::int64_t>(input).first(n-1),0);}catch(const std::invalid_argument&){rejected=true;}
    if(!rejected)throw std::runtime_error("signed RNS batch accepted inconsistent lengths");
    rejected=false;
    try{encoder.residues(output,std::span<const std::int64_t>(input).first(n),Count);}catch(const std::out_of_range&){rejected=true;}
    if(!rejected)throw std::runtime_error("signed RNS batch accepted invalid prime");
    std::cout<<"signed RNS primes="<<Count<<" full int64 domains, unaligned buffers, tails and guards passed\n";
}

int main(){
    checkSignedEncoding<1>();checkSignedEncoding<17>();
    checkPrefix<5>();checkPrefix<6>();
    checkBatchEncoding<1>();checkBatchEncoding<20>();checkBatchEncoding<32>();checkBatchEncoding<44>();
    checkBatchEncoding<64>();checkBatchEncoding<193>();checkBatchEncoding<440>();checkBatchEncoding<448>();
    checkBatchEncoding<83,165,83>();checkBatchEncoding<105,209,105>();checkBatchEncoding<127,253,127>();
    checkBatchEncoding<149,297,149>();checkBatchEncoding<171,341,171>();checkBatchEncoding<193,385,193>();
    checkBatchEncoding<1,1,1>();checkBatchEncoding<20,448,440>();
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
