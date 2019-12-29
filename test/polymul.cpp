#include<tfhe++.hpp>
#include<cassert>
#include<random>
#include<algorithm>
#include<iostream>

using namespace std;

int main(){
    const uint32_t num_test = 1000;
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> Bgdist(0, DEF_Bg);
    uniform_int_distribution<uint32_t> Torus32dist(0, UINT32_MAX);

    cout<<"Start LVL1 test."<<endl;
    for(int i; i<num_test;i++){
        array<uint32_t,DEF_N> a;
        for(int i = 0;i<DEF_N;i++)a[i] = Torus32dist(engine);
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
        for(int i = 0;i<DEF_N;i++)b[i] = Torus32dist(engine);

        array<uint32_t,DEF_N> polymul;
        TFHEpp::PolyMullvl1(polymul, a,b);
        array<uint32_t,DEF_N> naieve;
        for(int i = 0;i<DEF_N; i++){
            uint32_t ri = 0;
            for(int j = 0;j<=i; j++) ri += static_cast<int32_t>(a[j])*b[i-j];
            for(int j = i+1;j<DEF_N; j++) ri -= static_cast<int32_t>(a[j])*b[DEF_N+i-j];
            naieve[i] = ri;
        }
    for(int i = 0;i<DEF_N; i++) assert(abs(static_cast<int32_t>(naieve[i] - polymul[i])) <= 1);
    }
    cout<<"PoluMul Passed"<<endl;



    uniform_int_distribution<uint64_t> Bgbardist(0, DEF_Bgbar);
    uniform_int_distribution<uint64_t> Torus64dist(0, UINT64_MAX);

    cout<<"Start LVL2 test."<<endl;
    for(int i; i<num_test;i++){
        array<uint64_t,DEF_nbar> a;
        for(int i = 0;i<DEF_nbar;i++)a[i] = Torus64dist(engine);
        array<double,DEF_nbar> resfft;
        TFHEpp::TwistIFFTlvl2(resfft,a);
        array<uint64_t,DEF_nbar> res;
        TFHEpp::TwistFFTlvl2(res,resfft);
        for(int i = 0;i<DEF_N;i++)assert(abs(static_cast<int64_t>(a[i] - res[i]))<=(1<<14));
    }
    cout<<"FFT Passed"<<endl;

    for(int i; i<num_test;i++){
        array<uint64_t,DEF_nbar> a;
        for(int i = 0;i<DEF_nbar;i++)a[i] = Bgbardist(engine) - DEF_Bgbar/2;
        array<uint64_t,DEF_nbar> b;
        for(int i = 0;i<DEF_nbar;i++)b[i] = Torus64dist(engine);

        array<uint64_t,DEF_nbar> polymul;
        TFHEpp::PolyMullvl2(polymul, a,b);
        array<uint64_t,DEF_nbar> naieve;
        for(int i = 0;i<DEF_nbar; i++){
            uint64_t ri = 0;
            for(int j = 0;j<=i; j++) ri += static_cast<int64_t>(a[j])*b[i-j];
            for(int j = i+1;j<DEF_nbar; j++) ri -= static_cast<int64_t>(a[j])*b[DEF_nbar+i-j];
            naieve[i] = ri;
        }
    for(int i = 0;i<DEF_nbar; i++) assert(abs(static_cast<int64_t>(naieve[i] - polymul[i])) <= (1<<27));
    }
    cout<<"PoluMul Passed"<<endl;

    return 0;
}