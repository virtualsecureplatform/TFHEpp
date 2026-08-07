#include <cassert>
#include <iostream>
#include <tfhe++.hpp>

int main()
{
    {
        using iksP = TFHEpp::lvl10param;
        using bkP = TFHEpp::lvl02param;
        using sskP = TFHEpp::lvl22param;
        constexpr int validbit = 8;

        TFHEpp::SecretKey sk;
        TFHEpp::EvalKey ek;
        ek.emplaceiksk<iksP>(sk);
        ek.emplacebkfft<bkP>(sk);

        TFHEpp::AnnihilateKey<typename bkP::targetP> ahk;
        TFHEpp::annihilatekeygen<typename bkP::targetP>(ahk, sk);

        const auto digits = TFHEpp::EncodeHatEncoderInt8<typename sskP::targetP, validbit>(13);
        std::vector<TFHEpp::TLWE<typename iksP::domainP>> tlwes;
        TFHEpp::bootsSymEncrypt<typename iksP::domainP>(tlwes, digits, sk);
        assert(tlwes.size() == digits.size());

        TFHEpp::TRLWE<typename sskP::targetP> switched;
        TFHEpp::TLWES2CLPXIKS<iksP, bkP, sskP, 4, 0>(switched, tlwes, ahk, ek);

        bool nonzero = false;
        for (const auto &poly : switched)
            for (const auto value : poly) nonzero = nonzero || (value != 0);
        assert(nonzero);
    }

    {
        using iksP10 = TFHEpp::lvl1hparam;
        using iksP21 = TFHEpp::lvl21param;
        using bkP01 = TFHEpp::lvlh1param;
        using bkP02 = TFHEpp::lvlh2param;
        using iksP20 = TFHEpp::lvl2hparam;
        using bigP = TFHEpp::SS2CLPXlvl2param;
        constexpr uint32_t validbit = 16;

        TFHEpp::SecretKey sk;
        TFHEpp::EvalKey ek;
        ek.emplaceiksk<iksP10>(sk);
        ek.emplaceiksk<iksP21>(sk);
        ek.emplaceiksk<iksP20>(sk);
        ek.emplacebkfft<bkP01>(sk);
        ek.emplacebkfft<bkP02>(sk);

        for (const uint32_t plaintext : {0x2029U, 0x6029U, 0xa029U,
                                         0xe029U}) {
            const auto encoded = TFHEpp::EncodeHatEncoderP<bigP>(plaintext);
            const auto big =
                TFHEpp::clpxSymIntEncrypt<bigP>(encoded, sk.key.get<bigP>());

            std::vector<TFHEpp::TLWE<typename iksP10::domainP>> out(validbit);
            TFHEpp::CLPX2TLWESIKSanybit<iksP10, iksP21, bkP01, bkP02,
                                        iksP20, 9, 2>(out, big, ek, sk);

            uint32_t decoded = 0;
            for (uint32_t bit = 0; bit < validbit; bit++) {
                const bool actual =
                    TFHEpp::tlweSymDecrypt<typename bkP01::targetP>(
                        out[bit], sk.key.get<typename bkP01::targetP>());
                decoded |= static_cast<uint32_t>(actual) << bit;
            }
            if (decoded != plaintext) {
                std::cerr << "CLPX2TFHE decoded " << decoded << ", expected "
                          << plaintext << std::endl;
                return 1;
            }
        }
    }

    std::cout << "Passed" << std::endl;
}
