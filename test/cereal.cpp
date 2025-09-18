#include <cassert>
#include <fstream>
#include <iostream>
#include <tfhe++.hpp>

int main()
{
    TFHEpp::SecretKey sk;
    {
        {
            std::ofstream ofs{"./secret.key", std::ios::binary};
            cereal::PortableBinaryOutputArchive ar(ofs);
            sk.serialize(ar);
        }
        TFHEpp::SecretKey ski;
        {
            std::ifstream ifs{"./secret.key", std::ios::binary};
            cereal::PortableBinaryInputArchive ar(ifs);
            ski.serialize(ar);
        }
        for (int i = 0; i < TFHEpp::lvl0param::n; i++)
            assert(sk.key.get<TFHEpp::lvl0param>()[i] ==
                   ski.key.get<TFHEpp::lvl0param>()[i]);
        for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            assert(sk.key.get<TFHEpp::lvl1param>()[i] ==
                   ski.key.get<TFHEpp::lvl1param>()[i]);
        for (int i = 0; i < TFHEpp::lvl2param::n; i++)
            assert(sk.key.get<TFHEpp::lvl2param>()[i] ==
                   ski.key.get<TFHEpp::lvl2param>()[i]);
        std::cout << "n" << ":" << ski.params.lvl0.n << std::endl;
        assert(sk.params == ski.params);
    }

    {
        TFHEpp::IOpacket iopacket;
        iopacket.tlwelvl0ports["test"].resize(1);
        iopacket.tlwelvl0ports["test"][0] =
            TFHEpp::tlweSymEncrypt<TFHEpp::lvl0param>(
                TFHEpp::lvl0param::μ, TFHEpp::lvl0param::α,
                sk.key.get<TFHEpp::lvl0param>());
        {
            std::ofstream ofs{"./iopacket.data", std::ios::binary};
            cereal::PortableBinaryOutputArchive ar(ofs);
            iopacket.serialize(ar);
        }
    }

    {
        TFHEpp::EvalKey ek(sk);
        ek.emplacebkfft<TFHEpp::lvl01param>(sk);
        ek.emplacebkntt<TFHEpp::lvl01param>(sk);
        ek.emplaceiksk<TFHEpp::lvl10param>(sk);
        {
            std::ofstream ofs{"./gatekey.key", std::ios::binary};
            cereal::PortableBinaryOutputArchive ar(ofs);
            ek.serialize(ar);
        }
    }
    {
        TFHEpp::EvalKey ek;
        {
            std::ifstream ifs{"./gatekey.key", std::ios::binary};
            cereal::PortableBinaryInputArchive ar(ifs);
            ek.serialize(ar);
            std::cout << (ek.params.lvl0.n) << std::endl;
        }
    }
}