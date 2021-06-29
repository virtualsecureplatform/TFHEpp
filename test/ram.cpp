#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

int main()
{
    using CBbsP = lvl02param;
    using CBprivksP = lvl22param;
    using ksP = lvl20param;

    constexpr uint32_t address_bit = 8;  // Address by bytes.
    constexpr uint32_t memsize = 1 << address_bit;
    random_device seeder;
    default_random_engine engine(seeder());
    uniform_int_distribution<uint8_t> binary(0, 1);

    SecretKey *sk = new SecretKey;
    CloudKey<CBbsP, CBprivksP, ksP> *ck =
        new CloudKey<CBbsP, CBprivksP, ksP>(*sk);
    vector<uint8_t> pmemory(memsize);
    vector<array<ksP::domainP::T, ksP::domainP::n>> pmu(memsize);
    vector<uint8_t> address(address_bit);
    uint8_t pres;
    uint8_t wrflag;
    uint8_t writep;

    for (uint8_t &p : pmemory) p = binary(engine);

    for (int i = 0; i < memsize; i++) {
        pmu[i] = {};
        pmu[i][0] = pmemory[i] ? ksP::domainP::μ : -ksP::domainP::μ;
    }
    for (uint8_t &p : address) p = binary(engine);
    uint32_t addressint = 0;
    for (int i = 0; i < address_bit; i++)
        addressint += static_cast<uint32_t>(address[i]) << i;

    wrflag = binary(engine);
    writep = pmemory[addressint] > 0 ? 0 : 1;

    array<array<TRGSWFFT<typename ksP::domainP>, address_bit>, 2> *bootedTGSW =
        new array<array<TRGSWFFT<typename ksP::domainP>, address_bit>, 2>;
    vector<TLWE<typename ksP::targetP>> encaddress(address_bit);
    array<TRLWE<typename ksP::domainP>, memsize> *encmemory =
        new array<TRLWE<typename ksP::domainP>, memsize>;
    TLWE<typename ksP::domainP> encreadreshigh;
    TLWE<typename ksP::targetP> encreadres;
    TRLWE<typename ksP::domainP> encumemory;
    TLWE<typename ksP::targetP> cs;
    TLWE<typename ksP::targetP> c1;
    TRLWE<typename ksP::domainP> writed;

    encaddress = bootsSymEncrypt(address, *sk);
    for (int i = 0; i < memsize; i++)
        (*encmemory)[i] = trlweSymEncrypt<typename ksP::domainP>(
            pmu[i], ksP::domainP::α, (*sk).key.get<typename ksP::domainP>());
    cs = tlweSymEncrypt<typename ksP::targetP>(
        wrflag ? ksP::targetP::μ : -ksP::targetP::μ, ksP::targetP::α,
        (*sk).key.get<typename ksP::targetP>());
    c1 = tlweSymEncrypt<typename ksP::targetP>(
        writep ? ksP::targetP::μ : -ksP::targetP::μ, ksP::targetP::α,
        (*sk).key.get<typename ksP::targetP>());

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    // Addres CB
    for (int i = 0; i < address_bit; i++) {
        CircuitBootstrappingFFTwithInv<CBbsP, CBprivksP>(
            (*bootedTGSW)[1][i], (*bootedTGSW)[0][i], encaddress[i], (*ck).ck);
    }

    // Read
    RAMUX<typename ksP::domainP, address_bit>(encumemory, (*bootedTGSW)[0],
                                              *encmemory);
    SampleExtractIndex<typename ksP::domainP>(encreadreshigh, encumemory, 0);
    IdentityKeySwitch<ksP>(encreadres, encreadreshigh, (*ck).ksk);

    // Write
    HomMUXwoSE<CBbsP>(writed, cs, c1, encreadres, (*ck).ck.bkfft);
    for (int i = 0; i < memsize; i++) {
        TRLWE<typename ksP::domainP> temp;
        TFHEpp::RAMwriteBar<typename ksP::domainP, address_bit>(
            temp, writed, (*encmemory)[i], i, *bootedTGSW);
        TLWE<typename ksP::domainP> temp2;
        SampleExtractIndex<typename ksP::domainP>(temp2, temp, 0);
        TLWE<typename ksP::targetP> temp3;
        IdentityKeySwitch<ksP>(temp3, temp2, (*ck).ksk);
        GateBootstrappingTLWE2TRLWEFFT<CBbsP>((*encmemory)[i], temp3,
                                              (*ck).ck.bkfft);
    }

    end = chrono::system_clock::now();
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed << "ms" << endl;
    pres = tlweSymDecrypt<typename ksP::targetP>(
        encreadres, (*sk).key.get<typename ksP::targetP>());

    assert(static_cast<int>(pres) == static_cast<int>(pmemory[addressint]));

    array<bool, ksP::domainP::n> pwriteres =
        trlweSymDecrypt<typename ksP::domainP>(
            (*encmemory)[addressint], (*sk).key.get<typename ksP::domainP>());
    assert(static_cast<int>(pwriteres[0]) ==
           static_cast<int>((wrflag > 0) ? writep : pmemory[addressint]));

    cout << "Passed" << endl;
}