#include<tfhe++.hpp>
#include<assert.h>
#include<random>
#include<algorithm>
#include<iostream>

using namespace std;

int main(){
    const uint32_t num_test = 100;
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> Bgdist(-DEF_Bg, DEF_Bg);
    uniform_int_distribution<uint32_t> Torusdist(0, UINT32_MAX);

    cout<<"Start LVL1 test."<<endl;
    for(int i; i<num_test;i++){
        cout<<i<<endl;
        array<uint32_t,DEF_N> a;
        for(int i = 0;i<DEF_N;i++)a[i] = Bgdist(engine);
        array<uint32_t,DEF_N> b;
        for(int i = 0;i<DEF_N;i++)b[i] = Torusdist(engine);

        array<uint32_t,DEF_N> polymul = TFHEpp::PolyMullvl1(a,b);
        array<uint32_t,DEF_N> naieve;
        for(int j = 0;j<DEF_N; j++) naieve[j] = a[j]*b[0];
        for(int i = 1;i<DEF_N; i++){
            rotate(a.begin(),a.end(),a.end());
            a[0] *= -1;
            for(int j = 0;j<DEF_N; j++) naieve[j] += a[j]*b[i];
        }
    for(int i = 0;i<DEF_N; i++) assert(naieve[i] != polymul[i]);
    }

    return 0;
}