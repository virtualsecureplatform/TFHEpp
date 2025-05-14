#include <cassert>
#include <chrono>
#include <iostream>
#include <memory>
#include <random>
#include <tfhe++.hpp>

int main()
{
    using brP = TFHEpp::lvl01param;
    using iksP = TFHEpp::lvl10param;
    using cbiksP = TFHEpp::lvl20param;
    using cbbrP = TFHEpp::lvl02param;
    // using brP = TFHEpp::lvl02param;
    // using iksP = TFHEpp::lvl20param;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);
    std::uniform_int_distribution<uint8_t> bytegen(0, 255);
    std::unique_ptr<TFHEpp::SecretKey> sk(new TFHEpp::SecretKey());
    constexpr uint num_test = 1;
    std::vector<std::array<TFHEpp::TLWE<typename iksP::domainP>, 128>> cin(
        num_test);
    std::vector<std::array<uint8_t,16>> pin(num_test);
    std::vector<std::array<std::array<TFHEpp::TLWE<typename iksP::domainP>, 128>,11>> cexpandedkey(
        num_test);
    std::vector<std::array<uint8_t,16>> key(num_test);

    for (int i = 0; i < num_test; i++) {
        for (int j = 0; j < 16; j++)
            key[i][j] = bytegen(engine);
        std::array<uint8_t, 4*TFHEpp::Nb*(TFHEpp::Nr+1)> w;
        TFHEpp::KeyExpansion(w, key[i]);
        for (int j = 0; j < 11; j++)
            for(int k = 0; k < 16; k++)
                for(int l = 0; l < 8; l++)
                    cexpandedkey[i][j][k*8+l] =
                        TFHEpp::tlweSymEncrypt<typename iksP::domainP>((w[j*16+k]>>l)&0x1?1ULL << (std::numeric_limits<typename iksP::domainP::T>::digits - 2):-(1ULL << (std::numeric_limits<typename iksP::domainP::T>::digits - 2)), *sk);

        for (int j = 0; j < 16; j++)
            pin[i][j] = bytegen(engine);
        for(int j = 0; j < 16; j++)
            for(int k = 0; k < 8; k++)
                cin[i][j*8+k] =
                    TFHEpp::tlweSymEncrypt<typename iksP::domainP>((pin[i][j]>>k)&0x1?1ULL << (std::numeric_limits<typename iksP::domainP::T>::digits - 2):-(1ULL << (std::numeric_limits<typename iksP::domainP::T>::digits - 2)), *sk);
    }
    std::vector<std::array<TFHEpp::TLWE<typename brP::targetP>, 128>> cres(
        num_test);
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<brP>(*sk);
    ek.emplacebkfft<cbbrP>(*sk);
    ek.emplaceiksk<iksP>(*sk);
    ek.emplaceiksk<cbiksP>(*sk);
    ek.emplaceahk<typename cbbrP::targetP>(*sk);
    ek.emplacecbsk<typename cbbrP::targetP>(*sk);

    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();
    for (int test = 0; test < num_test; test++) {
        std::cout << "test: " << test << std::endl;
        TFHEpp::AESEnc<iksP, brP, cbiksP, cbbrP>(cres[test], cin[test], cexpandedkey[test], ek);
        // TFHEpp::AESEnc<iksP, brP>(cres[test], cin[test], cexpandedkey[test], ek);
    }

    end = std::chrono::system_clock::now();
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::cout << elapsed / num_test << "ms" << std::endl;
    for (int i = 0; i < num_test; i++) {
        AES aes(AESKeyLength::AES_128);
        auto c = aes.EncryptECB(pin[i].data(), 128, key[i].data());
        for(int j = 0; j < 16; j++) {
            uint8_t pres = 0;
            for (int k = 0; k < 8; k++) {
                if(TFHEpp::tlweSymDecrypt<typename brP::targetP>(
                     cres[i][j*8+k], *sk)) 
                    pres |= (1 << k);
            }
            // std::cout << "j: " << j << ",c: " << (int)c[j] << " pres: " << (int)pres << std::endl;
            assert(pres == c[j]);
        }
    }
    std::cout << "Passed" << std::endl;
}