#include <bitset>
#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

template <class P, uint32_t address_bit, uint32_t words_bit>
void combUROMUX(TRLWE<P>& res,
                const array<TRGSWFFT<P>, address_bit>& invaddress,
                const std::vector<TRLWE<P>>& data)
{
    constexpr uint32_t width_bit =
        P::nbit -
        words_bit;  // log_2 of how many words are in one TRLWElvl1 message.
    constexpr uint32_t Ubit = address_bit - width_bit;
    constexpr uint32_t num_trlwe = 1 << (Ubit);
    array<TRLWE<P>, num_trlwe / 2> temp;

    for (uint32_t index = 0; index < num_trlwe / 2; index++) {
        CMUXFFT<P>(temp[index], invaddress[width_bit], data[2 * index],
                   data[2 * index + 1]);
    }

    for (uint32_t bit = 0; bit < (Ubit - 2); bit++) {
        const uint32_t stride = 1 << bit;
        for (uint32_t index = 0; index < (num_trlwe >> (bit + 2)); index++) {
            CMUXFFT<P>(
                temp[(2 * index) * stride], invaddress[width_bit + bit + 1],
                temp[(2 * index) * stride], temp[(2 * index + 1) * stride]);
        }
    }

    constexpr uint32_t stride = 1 << (Ubit - 2);
    CMUXFFT<P>(res, invaddress[address_bit - 1], temp[0], temp[stride]);
}

template <class P, uint32_t address_bit, uint32_t words_bit>
void combRAMUX(
    std::vector<TRLWE<P>>& res,
    const array<TRGSWFFT<P>, address_bit>& invaddress,
    const std::vector<array<TRLWE<P>, 1 << address_bit>>& data)
{
    constexpr uint32_t words = 1U << words_bit;
    constexpr uint32_t num_trlwe = 1 << address_bit;
    for (int i = 0; i < words; i++) {
        array<TRLWE<P>, num_trlwe / 2> temp;

        for (uint32_t index = 0; index < num_trlwe / 2; index++) {
            CMUXFFT<P>(temp[index], invaddress[0], data[i][2 * index],
                       data[i][2 * index + 1]);
        }

        for (uint32_t bit = 0; bit < (address_bit - 2); bit++) {
            const uint32_t stride = 1 << bit;
            for (uint32_t index = 0; index < (num_trlwe >> (bit + 2));
                 index++) {
                CMUXFFT<P>(
                    temp[(2 * index) * stride], invaddress[bit + 1],
                    temp[(2 * index) * stride], temp[(2 * index + 1) * stride]);
            }
        }
        constexpr uint32_t stride = 1 << (address_bit - 2);
        CMUXFFT<P>(res[i], invaddress[address_bit - 1], temp[0], temp[stride]);
    }
}

template <class P, class iksP, class bkP, uint32_t address_bit,
          uint32_t words_bit>
void combWRAM(std::vector<array<TRLWE<P>, 1U << address_bit>>& encram,
              const array<array<TRGSWFFT<P>, address_bit>, 2>& address,
              const std::vector<TRLWE<P>>& encwritep, const EvalKey& ek)
{
    constexpr uint32_t memsize = 1U << address_bit;
    constexpr uint32_t words = 1U << words_bit;

#pragma omp parallel for
    for (int i = 0; i < memsize; i++) {
        const bitset<address_bit> addressbitset(i);
        for (int j = 0; j < words; j++) {
            TRLWE<P> temp = encwritep[j];
            for (int k = 0; k < address_bit; k++)
                CMUXFFT<P>(temp, address[addressbitset[k]][k], temp,
                           encram[j][i]);
            TLWE<P> temp2;
            SampleExtractIndex<P>(temp2, temp, 0);
            TLWE<lvl0param> temp3;
            IdentityKeySwitch<iksP>(temp3, temp2, ek.getiksk<iksP>());
            BlindRotate<bkP>(encram[j][i], temp3, ek.getbkfft<bkP>(),
                             μpolygen<P, P::μ>());
        }
    }
}

int main()
{
    using CBbrP = cblvl02param;
    using CBprivP = cblvl21param;
    using dataP = typename CBprivP::targetP;
    using gateiksP = independentlvl10param;
    using gatebkP = independentlvl01param;

#ifdef USE_CONCRETE
    // The concrete parameter set is selected for seven RAM address levels.
    // The final bit below selects RAM versus ROM. A 16-bit word leaves the two
    // upper-address bits required by combUROMUX while staying at that depth.
    constexpr uint32_t address_bit = 8;
    constexpr uint32_t words_bit = 4;
#else
    constexpr uint32_t address_bit = 10;
    constexpr uint32_t words_bit = 3;
#endif
    constexpr uint32_t words = 1U << words_bit;
    constexpr uint32_t memsize = 1 << address_bit;
    constexpr uint32_t numromtrlwe =
        1U << (address_bit - 1 + words_bit - dataP::nbit);
    constexpr uint32_t numramtrlwe = 1U << (address_bit - 1);
    constexpr uint32_t width_bit =
        dataP::nbit -
        words_bit;  // log_2 of how many words are in one TRLWE message.
    random_device seeder;
    default_random_engine engine(seeder());
    uniform_int_distribution<uint8_t> binary(0, 1);

    SecretKey* sk = new SecretKey;
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<gatebkP>(*sk);
    ek.emplacebkfft<CBbrP>(*sk);
    ek.emplaceiksk<gateiksP>(*sk);
    ek.emplaceprivksk4cb<CBprivP>(*sk);
    vector<uint8_t> ramp(memsize / 2 * words);  // unit of memsize is byte(8bit)
    vector<uint8_t> romp(memsize / 2 * words);
    vector<array<array<typename dataP::T, dataP::n>, numramtrlwe>>
        ramu(words);
    vector<array<typename dataP::T, dataP::n>> romu(numromtrlwe);
    vector<uint8_t> address(address_bit);
    array<uint8_t, words> pres;
    array<uint8_t, words> writep;

    for (uint8_t wrflag = 0; wrflag <= 1; wrflag++) {
        for (uint8_t msbaddress = 0; msbaddress <= 1; msbaddress++) {
            for (uint8_t& p : address) p = binary(engine);
            address[address_bit - 1] = msbaddress;
            for (uint8_t& p : ramp) p = binary(engine);
            for (uint8_t& p : romp) p = binary(engine);

            for (int i = 0; i < words; i++) {
                for (int j = 0; j < numramtrlwe; j++) {
                    ramu[i][j] = {};
                    ramu[i][j][0] =
                        ramp[j * words + i] ? dataP::μ : -dataP::μ;
                }
            }
            for (int i = 0; i < numromtrlwe; i++) {
                for (int j = 0; j < dataP::n; j++) {
                    romu[i][j] =
                        romp[i * dataP::n + j] ? dataP::μ : -dataP::μ;
                }
            }

            uint32_t addressint = 0;
            for (int i = 0; i < address_bit - 1; i++)
                addressint += static_cast<uint32_t>(address[i]) << i;

            for (int i = 0; i < words; i++)
                writep[i] = ramp[addressint * words + i] > 0 ? 0 : 1;

            std::array<array<TRGSWFFT<dataP>, address_bit - 1>, 2>*
                bootedTGSW = new (std::align_val_t(64))
                    array<array<TRGSWFFT<dataP>, address_bit - 1>,
                          2>;  // MSB of address is evaluated by HomMUX, not
                               // CMUX.
            vector<TLWE<dataP>> encaddress(address_bit);
            std::vector<TRLWE<dataP>> encrom(numromtrlwe);
            std::vector<array<TRLWE<dataP>, numramtrlwe>> encram(words);

            std::vector<TRLWE<dataP>> encramread(words);
            std::vector<TLWE<dataP>> encramreadres(words);

            TRLWE<dataP> encumemory;
            std::vector<TLWE<dataP>> encromreadres(words);

            std::vector<TLWE<dataP>> encreadres(words);

            TLWE<dataP> encwrflag;
            std::vector<TLWE<dataP>> encwritep(words);
            std::vector<TRLWE<dataP>> writed(words);

            bootsSymEncrypt<dataP>(encaddress, address,
                                   sk->key.getIndependent<dataP>());
            for (int i = 0; i < numromtrlwe; i++)
                trlweSymEncrypt<dataP>(encrom[i], romu[i],
                                       sk->key.getIndependent<dataP>());
            for (int i = 0; i < words; i++)
                for (int j = 0; j < numramtrlwe; j++)
                    trlweSymEncrypt<dataP>(encram[i][j], ramu[i][j],
                                           sk->key.getIndependent<dataP>());

            tlweSymEncrypt<dataP>(
                encwrflag, (wrflag > 0) ? dataP::μ : -dataP::μ,
                sk->key.getIndependent<dataP>());

            for (int i = 0; i < words; i++)
                tlweSymEncrypt<dataP>(
                    encwritep[i], writep[i] ? dataP::μ : -dataP::μ,
                    sk->key.getIndependent<dataP>());

            chrono::system_clock::time_point start, end;
            start = chrono::system_clock::now();
            // Addres CB
            for (int i = 0; i < address_bit - 1; i++) {
                CircuitBootstrappingWithInv<gateiksP, CBbrP, CBprivP>(
                    (*bootedTGSW)[1][i], (*bootedTGSW)[0][i], encaddress[i],
                    ek);
            }

            // Read
            combRAMUX<dataP, address_bit - 1, words_bit>(
                encramread, (*bootedTGSW)[0], encram);
            for (int i = 0; i < words; i++) {
                SampleExtractIndex<dataP>(encramreadres[i], encramread[i], 0);
            }

            combUROMUX<dataP, address_bit - 1, words_bit>(
                encumemory, (*bootedTGSW)[0], encrom);

            LROMUX<dataP, address_bit - 1, width_bit>(
                encromreadres, (*bootedTGSW)[1], encumemory);

            for (int i = 0; i < words; i++)
                HomMUX<dataP>(encreadres[i], encaddress[address_bit - 1],
                              encramreadres[i], encromreadres[i], ek);

            // Controll
            TLWE<dataP> cs;
            HomAND<gateiksP, gatebkP>(cs, encwrflag,
                                      encaddress[address_bit - 1], ek);
            for (int i = 0; i < words; i++)
                HomMUXwoSE<gateiksP, gatebkP>(writed[i], cs, encwritep[i],
                                              encramreadres[i], ek);

            // Write
            combWRAM<dataP, gateiksP, gatebkP, address_bit - 1, words_bit>(
                encram, *bootedTGSW, writed, ek);

            end = chrono::system_clock::now();
            double elapsed =
                std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                      start)
                    .count();
            cout << elapsed << "ms" << endl;

            // test
            for (int i = 0; i < words; i++)
                pres[i] = tlweSymDecrypt<dataP>(
                    encreadres[i], sk->key.getIndependent<dataP>());

            for (int i = 0; i < words; i++)
                assert(static_cast<int>(ramp[addressint * words + i]) ==
                       static_cast<int>(tlweSymDecrypt<dataP>(
                           encramreadres[i], sk->key.getIndependent<dataP>())));
            for (int i = 0; i < words; i++)
                assert(static_cast<int>(romp[addressint * words + i]) ==
                       static_cast<int>(tlweSymDecrypt<dataP>(
                           encromreadres[i], sk->key.getIndependent<dataP>())));

            for (int i = 0; i < words; i++)
                assert(static_cast<int>(pres[i]) ==
                       static_cast<int>(address[address_bit - 1] > 0
                                            ? ramp[addressint * words + i]
                                            : romp[addressint * words + i]));

            array<array<bool, dataP::n>, words> pwriteres;
            for (int i = 0; i < words; i++)
                pwriteres[i] = trlweSymDecrypt<dataP>(
                    encram[i][addressint], sk->key.getIndependent<dataP>());

            cout << static_cast<int>(wrflag > 0) << ":"
                 << static_cast<int>(tlweSymDecrypt<dataP>(
                        encwrflag, sk->key.getIndependent<dataP>()))
                 << endl;
            cout << static_cast<int>(address[address_bit - 1] > 0) << ":"
                 << static_cast<int>(
                        tlweSymDecrypt<dataP>(encaddress[address_bit - 1],
                                              sk->key.getIndependent<dataP>()))
                 << endl;
            bool csp = ((wrflag > 0) & (address[address_bit - 1] > 0));
            cout << static_cast<int>(csp) << ":"
                 << static_cast<int>(tlweSymDecrypt<dataP>(
                        cs, sk->key.getIndependent<dataP>()))
                 << endl;
            assert(static_cast<int>(tlweSymDecrypt<dataP>(
                       cs, sk->key.getIndependent<dataP>())) ==
                   static_cast<int>(csp));
            array<array<bool, dataP::n>, words> writedp;
            for (int i = 0; i < words; i++)
                writedp[i] = trlweSymDecrypt<dataP>(
                    writed[i], sk->key.getIndependent<dataP>());

            for (int i = 0; i < words; i++)
                assert(static_cast<int>(writep[i]) ==
                       static_cast<int>(tlweSymDecrypt<dataP>(
                           encwritep[i], sk->key.getIndependent<dataP>())));
            for (int i = 0; i < words; i++)
                assert(static_cast<int>(writedp[i][0]) ==
                       static_cast<int>(csp ? writep[i]
                                            : ramp[addressint * words + i]));
            for (int i = 0; i < words; i++)
                assert(static_cast<int>(pwriteres[i][0]) ==
                       static_cast<int>(csp ? writep[i]
                                            : ramp[addressint * words + i]));
        }
    }
    cout << "Passed" << endl;
}
