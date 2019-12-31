#pragma once

#include<array>
#include<algorithm>
#include<cmath>
#include<random>
#include<array>

#include<params.hpp>

namespace TFHEpp{
    using namespace std;
    static random_device generator;

    inline uint32_t dtot32(double d) {
        return int32_t(int64_t((d - int64_t(d))*(1L<<32)));
    }

    inline uint32_t gaussian32(uint32_t message, double α){
        normal_distribution<double> distribution(0.,α); //TODO: can we create a global distrib of param 1 and multiply by sigma?
        double err = distribution(generator);
        return message + dtot32(err);
    }
    
    template<uint32_t N>
    inline void MulInFD(array<double,N> &res,const array<double,N> &a,const array<double,N> &b){
        for(int i = 0;i < N/2; i++){
            double aimbim = a[i+N/2] * b[i+N/2];
            double arebim = a[i] * b[i+N/2];
            res[i] = a[i] * b[i] - aimbim;
            res[i+N/2] = a[i+N/2] * b[i] + arebim;
        }
    }

    template<uint32_t N>
    inline void FMAInFD(array<double,N> &res,const array<double,N> &a,const array<double,N> &b){
        for(int i = 0;i < N/2; i++){
            res[i] = a[i+N/2] * b[i+N/2] - res[i];
            res[i] = a[i] * b[i] - res[i];
            res[i+N/2] += a[i] * b[i+N/2];
            res[i+N/2] += a[i+N/2] * b[i];
        }
    }
}