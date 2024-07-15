#include <cassert>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

int main()
{
    // using CBbsP = lvl02param;
    // using CBprivksP = lvl22param;
    // using ksP = lvl20param;
    using CBbsP = lvl02param;
    using CBprivksP = lvl21param;
    using ksP = lvl10param;

    constexpr uint32_t address_bit = 8;  // Address by words.
    constexpr uint32_t words_bit = 5;
    constexpr uint32_t word = 1 << words_bit;
    constexpr uint32_t width_bit =
        ksP::domainP::nbit -
        words_bit;  // log_2 of how many words are in one TRLWE message.
    static_assert(address_bit >= width_bit);
    // static_assert(address_bit > width_bit+2);
    // constexpr uint32_t width = 1 << width_bit;
    static_assert(address_bit > width_bit);
    constexpr uint32_t num_trlwe = 1 << (address_bit - width_bit);
    random_device seeder;
    default_random_engine engine(seeder());
    uniform_int_distribution<uint8_t> binary(0, 1);

    SecretKey *sk = new SecretKey;
    TFHEpp::EvalKey ek;
    ek.emplaceiksk<ksP>(*sk);
    ek.emplacebkfft<CBbsP>(*sk);
    ek.emplaceprivksk4cb<CBprivksP>(*sk);
    vector<array<uint8_t, ksP::domainP::n>> pmemory(num_trlwe);
    vector<array<typename ksP::domainP::T, ksP::domainP::n>> pmu(num_trlwe);
    vector<uint8_t> address(address_bit);
    vector<uint8_t> pres(word);

    for (array<uint8_t, ksP::domainP::n> &i : pmemory)
        for (uint8_t &p : i) p = binary(engine);
    for (int i = 0; i < num_trlwe; i++)
        for (int j = 0; j < ksP::domainP::n; j++)
            pmu[i][j] = pmemory[i][j]
                            ? 2 * ksP::domainP::μ
                            : -2 * ksP::domainP::μ;  // This will increase noise
                                                     // torellance.
    for (uint8_t &p : address) p = binary(engine);

    alignas(64) array<TRGSWFFT<typename ksP::domainP>, address_bit> bootedTGSW;
    vector<TLWE<typename ksP::domainP>> encaddress(address_bit);
    array<TRLWE<typename ksP::domainP>, num_trlwe> encmemory;
    vector<TLWE<typename ksP::domainP>> encres(word);

    encaddress = bootsSymEncrypt(address, *sk);
    for (int i = 0; i < num_trlwe; i++)
        encmemory[i] = trlweSymEncrypt<typename ksP::domainP>(
            pmu[i], (*sk).key.get<typename ksP::domainP>());

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    for (int i = 0; i < width_bit; i++)
        CircuitBootstrappingFFT<ksP, CBbsP, CBprivksP>(bootedTGSW[i],
                                                       encaddress[i], ek);
    for (int i = width_bit; i < address_bit; i++)
        CircuitBootstrappingFFTInv<ksP, CBbsP, CBprivksP>(bootedTGSW[i],
                                                          encaddress[i], ek);
    TRLWE<typename ksP::domainP> encumemory;

    UROMUX<typename ksP::domainP, address_bit, width_bit>(
        encumemory, bootedTGSW, encmemory);
    LROMUX<typename ksP::domainP, address_bit, width_bit>(encres, bootedTGSW,
                                                          encumemory);
    end = chrono::system_clock::now();

    pres = bootsSymDecrypt(encres, *sk);
    uint32_t uaddress = 0;
    uint32_t laddress = 0;
    for (int i = 0; i < (address_bit - width_bit); i++)
        uaddress += address[i + width_bit] << i;
    array<bool, ksP::domainP::n> umemory;
    umemory = trlweSymDecrypt<typename ksP::domainP>(
        encumemory, (*sk).key.get<typename ksP::domainP>());

    for (int i = 0; i < width_bit; i++)
        laddress += static_cast<uint32_t>(address[i]) << (i + words_bit);
    for (uint32_t i = 0; i < word; i++)
        assert(static_cast<int>(pres[i]) ==
               static_cast<int>(pmemory[uaddress][laddress + i]));
    cout << "Passed" << endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed << "ms" << endl;
}