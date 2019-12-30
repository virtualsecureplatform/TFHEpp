#pragma once

#include<array>
#include<algorithm>
#include<cmath>
#include<random>

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
}