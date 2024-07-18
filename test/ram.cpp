#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

int main()
{
    using CBbsP = TFHEpp::lvl02param;
    using CBprivksP = TFHEpp::lvl21param;
    using ksP = TFHEpp::lvl10param;
    using brP = TFHEpp::lvl01param;
    using CBiksP = TFHEpp::lvl10param;

    constexpr uint32_t address_bit = 8;  // Address by bytes.
    constexpr uint32_t memsize = 1 << address_bit;
    random_device seeder;
    default_random_engine engine(seeder());
    uniform_int_distribution<uint8_t> binary(0, 1);

    SecretKey *sk = new SecretKey;
    TFHEpp::KeySwitchingKey<CBiksP> *iksk =
        new TFHEpp::KeySwitchingKey<CBiksP>();
    TFHEpp::ikskgen<CBiksP>(*iksk, *sk);
    TFHEpp::EvalKey ek;
    ek.emplaceiksk<ksP>(*sk);
    ek.emplacebkfft<CBbsP>(*sk);
    ek.emplacebkfft<brP>(*sk);
    ek.emplaceprivksk4cb<CBprivksP>(*sk);
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

    array<array<TRGSWFFT<typename CBprivksP::targetP>, address_bit>, 2>
        *bootedTGSW = new (std::align_val_t(64))
            array<array<TRGSWFFT<typename CBprivksP::targetP>, address_bit>, 2>;
    vector<TLWE<typename CBprivksP::targetP>> encaddress(address_bit);
    array<TRLWE<typename CBprivksP::targetP>, memsize> *encmemory =
        new (std::align_val_t(64))
            array<TRLWE<typename CBprivksP::targetP>, memsize>;
    TLWE<typename CBprivksP::targetP> encreadres;
    TRLWE<typename CBiksP::domainP> encumemory;
    TLWE<typename ksP::domainP> cs;
    TLWE<typename ksP::domainP> c1;
    TRLWE<typename CBiksP::domainP> writed;

    encaddress = bootsSymEncrypt(address, *sk);
    for (int i = 0; i < memsize; i++)
        (*encmemory)[i] = trlweSymEncrypt<typename ksP::domainP>(
            pmu[i], (*sk).key.get<typename ksP::domainP>());
    cs = tlweSymEncrypt<typename ksP::domainP>(
        wrflag ? ksP::domainP::μ : -ksP::domainP::μ,
        (*sk).key.get<typename ksP::domainP>());
    c1 = tlweSymEncrypt<typename ksP::domainP>(
        writep ? ksP::domainP::μ : -ksP::domainP::μ,
        (*sk).key.get<typename ksP::domainP>());

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    // Addres CB
    for (int i = 0; i < address_bit; i++) {
        CircuitBootstrappingFFTwithInv<ksP, CBbsP, CBprivksP>(
            (*bootedTGSW)[1][i], (*bootedTGSW)[0][i], encaddress[i], ek);
    }

    // Read
    RAMUX<typename CBprivksP::targetP, address_bit>(
        encumemory, (*bootedTGSW)[0], *encmemory);
    SampleExtractIndex<typename CBprivksP::targetP>(encreadres, encumemory, 0);

    // Write
    HomMUXwoSE<CBiksP, brP>(writed, cs, c1, encreadres, ek);
    for (int i = 0; i < memsize; i++) {
        TRLWE<typename CBiksP::domainP> temp;
        TFHEpp::RAMwriteBar<typename CBiksP::domainP, address_bit>(
            temp, writed, (*encmemory)[i], i, *bootedTGSW);
        TLWE<typename CBiksP::domainP> temp2;
        SampleExtractIndex<typename CBiksP::domainP>(temp2, temp, 0);
        TLWE<typename CBiksP::targetP> temp3;
        IdentityKeySwitch<CBiksP>(temp3, temp2, *iksk);
        BlindRotate<brP>((*encmemory)[i], temp3, ek.getbkfft<brP>(),
                         μpolygen<typename brP::targetP, brP::targetP::μ>());
    }

    end = chrono::system_clock::now();
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed << "ms" << endl;
    pres = tlweSymDecrypt<typename CBiksP::domainP>(
        encreadres, (*sk).key.get<typename CBiksP::domainP>());

    assert(static_cast<int>(pres) == static_cast<int>(pmemory[addressint]));

    array<bool, CBprivksP::targetP::n> pwriteres =
        trlweSymDecrypt<typename CBprivksP::targetP>(
            (*encmemory)[addressint],
            (*sk).key.get<typename CBprivksP::targetP>());
    assert(static_cast<int>(pwriteres[0]) ==
           static_cast<int>((wrflag > 0) ? writep : pmemory[addressint]));

    cout << "Passed" << endl;
}