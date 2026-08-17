#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

int main()
{
    using CBbsP = TFHEpp::cblvl02param;
    using CBprivksP = TFHEpp::cblvl21param;
    using dataP = typename CBprivksP::targetP;
    using ksP = TFHEpp::independentlvl10param;
    using brP = TFHEpp::independentlvl01param;

    constexpr uint32_t address_bit = 8;  // Address by bytes.
    constexpr uint32_t memsize = 1 << address_bit;
    random_device seeder;
    default_random_engine engine(seeder());
    uniform_int_distribution<uint8_t> binary(0, 1);

    SecretKey* sk = new SecretKey;
    TFHEpp::EvalKey ek;
    ek.emplaceiksk<ksP>(*sk);
    ek.emplacebkfft<CBbsP>(*sk);
    ek.emplacebkfft<brP>(*sk);
    ek.emplaceprivksk4cb<CBprivksP>(*sk);
    vector<uint8_t> pmemory(memsize);
    vector<array<dataP::T, dataP::n>> pmu(memsize);
    vector<uint8_t> address(address_bit);
    uint8_t pres;
    uint8_t wrflag;
    uint8_t writep;

    for (uint8_t& p : pmemory) p = binary(engine);

    for (int i = 0; i < memsize; i++) {
        pmu[i] = {};
        pmu[i][0] = pmemory[i] ? dataP::μ : -dataP::μ;
    }
    for (uint8_t& p : address) p = binary(engine);
    uint32_t addressint = 0;
    for (int i = 0; i < address_bit; i++)
        addressint += static_cast<uint32_t>(address[i]) << i;

    wrflag = binary(engine);
    writep = pmemory[addressint] > 0 ? 0 : 1;

    array<array<TRGSWFFT<dataP>, address_bit>, 2>*
        bootedTGSW = new (std::align_val_t(64))
            array<array<TRGSWFFT<dataP>, address_bit>, 2>;
    vector<TLWE<dataP>> encaddress(address_bit);
    array<TRLWE<dataP>, memsize>* encmemory =
        new (std::align_val_t(64)) array<TRLWE<dataP>, memsize>;
    TLWE<dataP> encreadres;
    TRLWE<dataP> encumemory;
    TLWE<dataP> cs;
    TLWE<dataP> c1;
    TRLWE<dataP> writed;

    bootsSymEncrypt<dataP>(encaddress, address,
                           sk->key.getIndependent<dataP>());
    for (int i = 0; i < memsize; i++)
        trlweSymEncrypt<dataP>((*encmemory)[i], pmu[i],
                               sk->key.getIndependent<dataP>());
    tlweSymEncrypt<dataP>(cs, wrflag ? dataP::μ : -dataP::μ,
                          sk->key.getIndependent<dataP>());
    tlweSymEncrypt<dataP>(c1, writep ? dataP::μ : -dataP::μ,
                          sk->key.getIndependent<dataP>());

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    // Addres CB
    for (int i = 0; i < address_bit; i++) {
        CircuitBootstrappingWithInv<ksP, CBbsP, CBprivksP>(
            (*bootedTGSW)[1][i], (*bootedTGSW)[0][i], encaddress[i], ek);
    }

    // Read
    RAMUX<dataP, address_bit>(encumemory, (*bootedTGSW)[0], *encmemory);
    SampleExtractIndex<dataP>(encreadres, encumemory, 0);

    // Write
    HomMUXwoSE<ksP, brP>(writed, cs, c1, encreadres, ek);
    for (int i = 0; i < memsize; i++) {
        TRLWE<dataP> temp;
        TFHEpp::RAMwriteBar<dataP, address_bit>(
            temp, writed, (*encmemory)[i], i, *bootedTGSW);
        TLWE<dataP> temp2;
        SampleExtractIndex<dataP>(temp2, temp, 0);
        TLWE<typename ksP::targetP> temp3;
        IdentityKeySwitch<ksP>(temp3, temp2, ek.getiksk<ksP>());
        BlindRotate<brP>((*encmemory)[i], temp3, ek.getbkfft<brP>(),
                         μpolygen<typename brP::targetP, brP::targetP::μ>());
    }

    end = chrono::system_clock::now();
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed << "ms" << endl;
    pres = tlweSymDecrypt<dataP>(encreadres,
                                 sk->key.getIndependent<dataP>());

    assert(static_cast<int>(pres) == static_cast<int>(pmemory[addressint]));

    array<bool, dataP::n> pwriteres = trlweSymDecrypt<dataP>(
        (*encmemory)[addressint], sk->key.getIndependent<dataP>());
    assert(static_cast<int>(pwriteres[0]) ==
           static_cast<int>((wrflag > 0) ? writep : pmemory[addressint]));

    cout << "Passed" << endl;
}
