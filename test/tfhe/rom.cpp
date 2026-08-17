#include <bit>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>
#include <sstream>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

int main()
{
    // using CBbsP = lvl02param;
    // using CBprivksP = lvl22param;
    // using ksP = lvl20param;
    using CBbsP = cblvl02param;
    using CBprivksP = cblvl21param;
    using dataP = typename CBprivksP::targetP;
    using ksP = lvl10param;

    constexpr uint32_t address_bit = 8;  // Address by words.
    constexpr uint32_t words_bit = 5;
    constexpr uint32_t word = 1 << words_bit;
    constexpr uint32_t width_bit =
        dataP::nbit -
        words_bit;  // log_2 of how many words are in one TRLWE message.
    static_assert(address_bit >= width_bit);
    // static_assert(address_bit > width_bit+2);
    // constexpr uint32_t width = 1 << width_bit;
    static_assert(address_bit > width_bit);
    constexpr uint32_t num_trlwe = 1 << (address_bit - width_bit);
    random_device seeder;
    default_random_engine engine(seeder());
    uniform_int_distribution<uint8_t> binary(0, 1);

    SecretKey* sk = new SecretKey;
    TFHEpp::EvalKey ek;
    ek.emplaceiksk<ksP>(*sk);
    ek.emplacebk<lvl01param>(*sk);
    ek.emplacesubiksk<ksP>(*sk);
    ek.emplacebkfft<CBbsP>(*sk);
    ek.emplacebkfft<lvl01param>(*sk);
    ek.emplaceprivksk4cb<CBprivksP>(*sk);
    std::stringstream serialized_key{ios::binary | ios::out | ios::in};
    {
        cereal::PortableBinaryOutputArchive archive{serialized_key};
        archive(ek);
    }
    TFHEpp::EvalKey deserialized_ek;
    {
        cereal::PortableBinaryInputArchive archive{serialized_key};
        archive(deserialized_ek);
    }
    auto& test_ek = deserialized_ek;
    vector<array<uint8_t, dataP::n>> pmemory(num_trlwe);
    vector<array<typename dataP::T, dataP::n>> pmu(num_trlwe);
    vector<uint8_t> address(address_bit);
    vector<uint8_t> pres(word);

    for (array<uint8_t, dataP::n>& i : pmemory)
        for (uint8_t& p : i) p = binary(engine);
    for (int i = 0; i < num_trlwe; i++)
        for (int j = 0; j < dataP::n; j++)
            pmu[i][j] = pmemory[i][j]
                            ? 2 * dataP::μ
                            : -2 * dataP::μ;  // This will increase noise
                                                     // torellance.
    for (uint8_t& p : address) p = binary(engine);

    alignas(64) array<TRGSWFFT<dataP>, address_bit> bootedTGSW;
    vector<TLWE<dataP>> encaddress(address_bit);
    vector<TLWE<typename CBbsP::domainP>> directEncaddress(address_bit);
    array<TRLWE<dataP>, num_trlwe> encmemory;
    vector<TLWE<dataP>> encres(word);

    bootsSymEncrypt<dataP>(encaddress, address,
                           sk->key.getIndependent<dataP>());
    bootsSymEncrypt<typename CBbsP::domainP>(directEncaddress, address, *sk);
    for (int i = 0; i < num_trlwe; i++)
        trlweSymEncrypt<dataP>(encmemory[i], pmu[i],
                               (*sk).key.getIndependent<dataP>());

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    for (int i = 0; i < width_bit; i++)
        CircuitBootstrapping<CBbsP, CBprivksP>(bootedTGSW[i],
                                               directEncaddress[i], test_ek);
    for (int i = width_bit; i < address_bit; i++)
        CircuitBootstrappingInv<CBbsP, CBprivksP>(bootedTGSW[i],
                                                  directEncaddress[i], test_ek);
    TRLWE<dataP> encumemory;

    UROMUX<dataP, address_bit, width_bit>(encumemory, bootedTGSW, encmemory);
    LROMUX<dataP, address_bit, width_bit>(encres, bootedTGSW, encumemory);
    end = chrono::system_clock::now();

    pres = bootsSymDecrypt<dataP>(encres, sk->key.getIndependent<dataP>());
    uint32_t uaddress = 0;
    uint32_t laddress = 0;
    for (int i = 0; i < (address_bit - width_bit); i++)
        uaddress += address[i + width_bit] << i;
    array<bool, dataP::n> umemory;
    umemory = trlweSymDecrypt<dataP>(encumemory,
                                     (*sk).key.getIndependent<dataP>());

    for (int i = 0; i < width_bit; i++)
        laddress += static_cast<uint32_t>(address[i]) << (i + words_bit);
    for (uint32_t i = 0; i < word; i++)
        assert(static_cast<int>(pres[i]) ==
               static_cast<int>(pmemory[uaddress][laddress + i]));

    for (uint32_t i = 0; i < word; i++) {
        TLWE<typename ksP::targetP> switched;
        IdentityKeySwitch<ksP>(switched, encres[i], test_ek.getiksk<ksP>());
        assert(tlweSymDecrypt<typename ksP::targetP>(switched, *sk) ==
               pmemory[uaddress][laddress + i]);
    }

    // Iyokan allows a ROM address width smaller than the number of words in a
    // TRLWE. Exercise its horizontal mux schedule and final key switching.
    constexpr uint32_t iyokan_address_bit = 4;
    constexpr uint32_t iyokan_output_bit = 8;
    constexpr uint32_t iyokan_log_words =
        dataP::nbit - std::countr_zero(iyokan_output_bit);
    array<uint8_t, 128> iyokan_memory{};
    for (uint32_t i = 0; i < iyokan_memory.size(); i++)
        iyokan_memory[i] = (i / 8) & 1;
    array<typename dataP::T, dataP::n> iyokan_mu{};
    for (uint32_t i = 0; i < iyokan_memory.size(); i++)
        iyokan_mu[i] = iyokan_memory[i] ? dataP::μ : -dataP::μ;
    TRLWE<dataP> iyokan_data;
    trlweSymEncrypt<dataP>(iyokan_data, iyokan_mu,
                           sk->key.getIndependent<dataP>());
    array<uint8_t, iyokan_address_bit> iyokan_address{1, 1, 1, 0};
    array<TLWE<typename CBbsP::domainP>, iyokan_address_bit> iyokan_encaddress;
    for (uint32_t i = 0; i < iyokan_address_bit; i++)
        tlweSymEncrypt<typename CBbsP::domainP>(
            iyokan_encaddress[i],
            iyokan_address[i] ? CBbsP::domainP::μ : -CBbsP::domainP::μ,
            sk->key.getSubset<typename CBbsP::domainP>());
    array<TRGSWFFT<dataP>, iyokan_address_bit> iyokan_tgsw;
    for (uint32_t i = 0; i < iyokan_address_bit; i++)
        CircuitBootstrapping<CBbsP, CBprivksP>(iyokan_tgsw[i],
                                               iyokan_encaddress[i], test_ek);
    TRLWE<dataP> iyokan_acc = iyokan_data, iyokan_temp;
    for (uint32_t bit = 1; bit <= iyokan_log_words; bit++) {
        if (iyokan_log_words - bit >= iyokan_address_bit) continue;
        for (int component = 0; component <= dataP::k; component++)
            PolynomialMulByXaiMinusOne<dataP>(
                iyokan_temp[component], iyokan_acc[component],
                2 * dataP::n - (dataP::n >> bit));
        ExternalProduct<dataP>(iyokan_temp, iyokan_temp,
                               iyokan_tgsw[iyokan_log_words - bit]);
        for (uint32_t i = 0; i < dataP::n; i++)
            for (int component = 0; component <= dataP::k; component++)
                iyokan_acc[component][i] += iyokan_temp[component][i];
    }
    for (uint32_t i = 0; i < iyokan_output_bit; i++) {
        TLWE<dataP> extracted;
        TLWE<typename ksP::targetP> switched;
        SampleExtractIndex<dataP>(extracted, iyokan_acc, i);
        IdentityKeySwitch<ksP>(switched, extracted, test_ek.getiksk<ksP>());
        assert(tlweSymDecrypt<typename ksP::targetP>(switched, *sk) ==
               iyokan_memory[7 * iyokan_output_bit + i]);
    }
    cout << "Passed" << endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed << "ms" << endl;
}
