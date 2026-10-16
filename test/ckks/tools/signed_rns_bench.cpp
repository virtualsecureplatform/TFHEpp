// Public signed-integer conversion benchmark; no keys or ciphertexts.
// c++ -O3 -march=native -funroll-loops -std=c++20 -Iinclude \
//   -Ithirdparties/cereal/include test/ckks/tools/signed_rns_bench.cpp -o /tmp/signed-rns-bench
#include <modular_rns32.hpp>
#include <random>
#include <chrono>
#include <iostream>
#include <iomanip>
using namespace TFHEpp;
int main()
{
    std::mt19937_64 random(20261017);
    modular_rns32::SignedEncoder<> encoder;
    std::vector<std::int64_t> input(16389);
    std::vector<std::uint32_t> output(input.size());
    std::array<std::vector<std::uint32_t>,17> expected;
    std::cout<<std::setprecision(17);
    for(unsigned bits:{40,44,64}){
        for(auto&v:input){
            auto word=random();
            if(bits<64)word=(word&((std::uint64_t{1}<<bits)-1))-(std::uint64_t{1}<<(bits-1));
            v=std::bit_cast<std::int64_t>(word);
        }
        for(std::size_t p=0;p<17;++p){expected[p].resize(input.size());encoder.residues(expected[p],input,p,false);}
        for(bool batch:{false,true,true,false})for(int trial=0;trial<5;++trial){
            double seconds=0;
            for(int repeat=0;repeat<16;++repeat)for(std::size_t p=0;p<17;++p){
                const auto start=std::chrono::steady_clock::now();
                encoder.residues(output,input,p,batch);
                seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
                // Consume every output outside each timed call.
                if(output!=expected[p])throw std::runtime_error("timed signed encoding changed residues");
            }
            std::cout<<"bits="<<bits<<" batch="<<batch<<" trial="<<trial<<" seconds="<<seconds<<" identical\n";
        }
    }
}
