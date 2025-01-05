#include<fnt.hpp>
#include<random>
#include<array>
#include<iostream>
#include<cassert>

int main(){
    constexpr uint32_t num_test = 1000;
    constexpr unsigned int Nbit = 6;
    constexpr unsigned int N = 1u << Nbit;

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<int64_t> Pdist(0, FNTpp::P);

    std::cout<< "Start ModLShift Test"<< std::endl;
    for(int test = 0; test < num_test; test++){
        const int64_t a = Pdist(engine);
        const uint shift = std::uniform_int_distribution<uint>(0, 63)(engine);
        const int64_t res = FNTpp::ModLshift(a, shift);
        if(res != (static_cast<__int128_t>(a) << shift) % FNTpp::P)
            std::cout << "a: " << a << " shift: " << shift << " res: " << res << " expected: " << static_cast<int64_t>((static_cast<__int128_t>(a) << shift) % FNTpp::P) << std::endl;
        assert(res == (static_cast<__int128_t>(a) << shift) % FNTpp::P);
    }
    std::cout<< "Passed ModLShift"<< std::endl;

    std::cout << "invN Test" << std::endl;
    assert(1 == FNTpp::ModLshift(N, 2*FNTpp::K-Nbit));
    std::cout << "Passed invN" << std::endl;

    std::cout << "Start univariable FNT only test." << std::endl;
    for(int test = 0; test < num_test; test++){
        std::array<int64_t, N> a;
        for(int i = 0; i < N; i++) a[i] = Pdist(engine);
        std::array<int64_t, N> res;
        res = a;
        FNTpp::FNT<Nbit>(res);
        FNTpp::IFNT<Nbit>(res);
        FNTpp::MulInvN<Nbit>(res);
        for(int i = 0; i < N; i++)
            if(a[i] != res[i])   std::cout << "i: "<< i << " a: " << a[i] << " res: " << res[i] << std::endl;
        for(int i = 0; i < N; i++) assert(a[i] == res[i]);
    }
    std::cout << "Univariable FNT only test Passed" << std::endl;

    return 0;
}