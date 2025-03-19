#include <array>
#include <cassert>
#include <fnt.hpp>
#include <iostream>
#include <random>

int main()
{
    constexpr uint32_t num_test = 1000;

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<int64_t> Pdist(0, FNTpp::P);

    {
        std::cout << "Start ModLShift Test" << std::endl;

        constexpr unsigned int Nbit = FNTpp::Kbit + 1;
        constexpr unsigned int N = 1u << Nbit;
        for (int test = 0; test < num_test; test++) {
            const int64_t a = Pdist(engine);
            const uint shift =
                std::uniform_int_distribution<uint>(0, 63)(engine);
            const int64_t res = FNTpp::ModLshift(a, shift);
            if (res != (static_cast<__int128_t>(a) << shift) % FNTpp::P)
                std::cout << "a: " << a << " shift: " << shift
                          << " res: " << res << " expected: "
                          << static_cast<int64_t>(
                                 (static_cast<__int128_t>(a) << shift) %
                                 FNTpp::P)
                          << std::endl;
            assert(res == (static_cast<__int128_t>(a) << shift) % FNTpp::P);
        }
        std::cout << "Passed ModLShift" << std::endl;

        std::cout << "invN Test" << std::endl;
        assert(1 == FNTpp::ModLshift(N, 2 * FNTpp::K - Nbit));
        std::cout << "Passed invN" << std::endl;

        std::cout << "Start univariable FNT only test." << std::endl;
        for (int test = 0; test < num_test; test++) {
            std::array<int64_t, N> a;
            for (int i = 0; i < N; i++) a[i] = Pdist(engine);
            std::array<int64_t, N> res;
            res = a;
            FNTpp::FNT<Nbit>(res);
            FNTpp::IFNT<Nbit>(res);
            FNTpp::MulInvN<Nbit>(res);
            for (int i = 0; i < N; i++)
                if (a[i] != res[i])
                    std::cout << "i: " << i << " a: " << a[i]
                              << " res: " << res[i] << std::endl;
            for (int i = 0; i < N; i++) assert(a[i] == res[i]);
        }
        std::cout << "Univariable FNT only test Passed" << std::endl;
    }

    {
        std::cout << "Start univariable TwistFNT only test." << std::endl;
        constexpr unsigned int Nbit = FNTpp::Kbit;
        constexpr unsigned int N = 1u << Nbit;
        for (int test = 0; test < num_test; test++) {
            std::array<int64_t, N> a;
            for (int i = 0; i < N; i++) a[i] = Pdist(engine);
            std::array<int64_t, N> res, temp;
            res = a;
            FNTpp::TwistFNT<Nbit>(temp, res);
            FNTpp::TwistIFNT<Nbit>(res, temp);
            for (int i = 0; i < N; i++)
                if (a[i] != res[i])
                    std::cout << "i: " << i << " a: " << a[i]
                              << " res: " << res[i] << std::endl;
            for (int i = 0; i < N; i++) assert(a[i] == res[i]);
        }
        std::cout << "Univariable TwistFNT only test Passed" << std::endl;

        std::cout << "Start univariable TwistFNT Mul test." << std::endl;
        for (int test = 0; test < num_test; test++) {
            std::array<int64_t, N> a, b;
            for (int i = 0; i < N; i++) a[i] = Pdist(engine);
            for (int i = 0; i < N; i++) b[i] = Pdist(engine);
            std::array<int64_t, N> res, naieve, fnta, fntb;
            naieve = {};
            for (int i = 0; i < N; i++) {
                for (int j = 0; j <= i; j++)
                    naieve[i] =
                        (naieve[i] + static_cast<__int128_t>(a[j]) * b[i - j]) %
                        FNTpp::P;
                for (int j = i + 1; j < N; j++)
                    naieve[i] = (naieve[i] -
                                 static_cast<__int128_t>(a[j]) * b[N + i - j]) %
                                FNTpp::P;
            }
            for (int i = 0; i < N; i++)
                naieve[i] = (naieve[i] + FNTpp::P) % FNTpp::P;
            FNTpp::TwistFNT<Nbit>(fnta, a);
            FNTpp::TwistFNT<Nbit>(fntb, b);
            for (int i = 0; i < N; i++)
                fnta[i] =
                    (static_cast<__int128_t>(fnta[i]) * fntb[i]) % FNTpp::P;
            FNTpp::TwistIFNT<Nbit>(res, fnta);
            for (int i = 0; i < N; i++)
                if (naieve[i] != res[i])
                    std::cout << "i: " << i << " naieve: " << naieve[i]
                              << " res: " << res[i] << std::endl;
            for (int i = 0; i < N; i++) assert(naieve[i] == res[i]);
        }
        std::cout << "Univariable TwistFNT Mul test Passed" << std::endl;
    }

    std::cout << "Start multivariable TwistFNT only test." << std::endl;
    constexpr unsigned int Nbit = 10;
    constexpr unsigned int N = 1u << Nbit;
    constexpr unsigned int NinFNT = 1u << (Nbit + (Nbit / (FNTpp::Kbit + 1)));
    for (int test = 0; test < num_test; test++) {
        std::array<int64_t, N> a;
        for (int i = 0; i < N; i++) a[i] = Pdist(engine);
        std::array<int64_t, N> res;
        std::array<int64_t, NinFNT> temp;
        res = a;
        FNTpp::TwistFNT<Nbit>(temp, res);
        FNTpp::TwistIFNT<Nbit>(res, temp);
        for (int i = 0; i < N; i++)
            if (a[i] != res[i])
                std::cout << "i: " << i << " a: " << a[i] << " res: " << res[i]
                          << std::endl;
        for (int i = 0; i < N; i++) assert(a[i] == res[i]);
    }
    std::cout << "Multivariable TwistFNT only test Passed" << std::endl;

    std::cout << "Start multivariable TwistFNT Mul test." << std::endl;
    for (int test = 0; test < num_test; test++) {
        std::array<int64_t, N> a, b;
        for (int i = 0; i < N; i++) a[i] = Pdist(engine);
        for (int i = 0; i < N; i++) b[i] = Pdist(engine);
        std::array<int64_t, N> res, naieve;
        naieve = {};
        for (int i = 0; i < N; i++) {
            for (int j = 0; j <= i; j++)
                naieve[i] =
                    (naieve[i] + static_cast<__int128_t>(a[j]) * b[i - j]) %
                    FNTpp::P;
            for (int j = i + 1; j < N; j++)
                naieve[i] =
                    (naieve[i] - static_cast<__int128_t>(a[j]) * b[N + i - j]) %
                    FNTpp::P;
        }
        for (int i = 0; i < N; i++)
            naieve[i] = (naieve[i] + FNTpp::P) % FNTpp::P;

        std::array<int64_t, NinFNT> fnta, fntb;
        FNTpp::TwistFNT<Nbit>(fnta, a);
        FNTpp::TwistFNT<Nbit>(fntb, b);
        for (int i = 0; i < NinFNT; i++)
            fnta[i] = (static_cast<__int128_t>(fnta[i]) * fntb[i]) % FNTpp::P;
        FNTpp::TwistIFNT<Nbit>(res, fnta);
        for (int i = 0; i < N; i++)
            if (naieve[i] != res[i])
                std::cout << "i: " << i << " naieve: " << naieve[i]
                          << " res: " << res[i] << std::endl;
        for (int i = 0; i < N; i++) assert(naieve[i] == res[i]);
    }
    std::cout << "Multivariable TwistFNT Mul test Passed" << std::endl;

    return 0;
}