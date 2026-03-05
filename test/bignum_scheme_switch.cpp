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
        TFHEpp::TLWES2BigNumIKSezM<iksP, bkP, sskP, 4, 0>(switched, tlwes, ahk, ek,
                                                          sk);

        bool nonzero = false;
        for (const auto &poly : switched)
            for (const auto value : poly) nonzero = nonzero || (value != 0);
        assert(nonzero);
    }

    {
        using iksP10 = TFHEpp::lvl10param;
        using iksP21 = TFHEpp::lvl21param;
        using bkP01 = TFHEpp::lvl01param;
        using bkP02 = TFHEpp::lvl02param;
        using iksP20 = TFHEpp::lvl22param;
        using bigP = typename iksP20::domainP;
        constexpr uint32_t validbit = 8;

        TFHEpp::SecretKey sk;
        TFHEpp::EvalKey ek;
        ek.emplaceiksk<iksP10>(sk);
        ek.emplaceiksk<iksP21>(sk);
        ek.emplaceiksk<iksP20>(sk);
        ek.emplacebkfft<bkP01>(sk);
        ek.emplacebkfft<bkP02>(sk);

        const auto encoded = TFHEpp::EncodeHatEncoderP<bigP>(9);
        const auto big = TFHEpp::bigNumSymIntEncrypt<bigP>(encoded, sk.key.get<bigP>());

        std::vector<TFHEpp::TLWE<typename bkP01::targetP>> out(validbit);
        TFHEpp::BIGNUM2TLWESIKSanybit<iksP10, iksP21, bkP01, bkP02, iksP20, 4, 2>(
            out, big, ek, sk);

        bool nonzero = false;
        for (const auto &tlwe : out)
            for (const auto value : tlwe) nonzero = nonzero || (value != 0);
        assert(nonzero);
    }

    std::cout << "Passed" << std::endl;
}
