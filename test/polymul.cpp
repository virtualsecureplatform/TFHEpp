#include<tfhe++.hpp>
#include<cassert>
#include<random>
#include<algorithm>
#include<iostream>

using namespace std;

int main(){
    const uint32_t num_test = 100;
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> Bgdist(0, DEF_Bg);
    uniform_int_distribution<uint32_t> Torusdist(0, UINT32_MAX);

    cout<<"Start LVL1 test."<<endl;
    for(int i; i<num_test;i++){
        array<uint32_t,DEF_N> a;
        for(int i = 0;i<DEF_N;i++)a[i] = Torusdist(engine);
        array<double,DEF_N> resfft;
        TFHEpp::TwistIFFTlvl1(resfft,a);
        array<uint32_t,DEF_N> res;
        TFHEpp::TwistFFTlvl1(res,resfft);
        for(int i = 0;i<DEF_N;i++)assert(abs(static_cast<int32_t>(a[i] - res[i]))<=1);
    }
    cout<<"FFT Passed"<<endl;

    for(int i; i<num_test;i++){
        array<uint32_t,DEF_N> a;
        for(int i = 0;i<DEF_N;i++)a[i] = Bgdist(engine) - DEF_Bg/2;
        array<uint32_t,DEF_N> b;
        for(int i = 0;i<DEF_N;i++)b[i] = Torusdist(engine);

        array<uint32_t,DEF_N> polymul;
        TFHEpp::PolyMullvl1(polymul, a,b);
        array<uint32_t,DEF_N> naieve;
        for(int i = 0;i<DEF_N; i++){
            uint32_t ri = 0;
            for(int j = 0;j<=i; j++) ri += static_cast<int32_t>(a[j])*b[i-j];
            for(int j = i+1;j<DEF_N; j++) ri -= static_cast<int32_t>(a[j])*b[DEF_N+i-j];
            naieve[i] = ri;
        }
    for(int i = 0;i<DEF_N; i++) assert(abs(static_cast<int>(naieve[i] - polymul[i])) <= 1);
    }
    cout<<"PoluMul Passed"<<endl;

    return 0;
}