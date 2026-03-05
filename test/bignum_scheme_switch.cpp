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

    std::cout << "Passed" << std::endl;
}
