#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/vector.hpp>
#include <fstream>
#include <memory>
#include <random>
#include <tfhe++.hpp>
#include <vector>

using iksP = TFHEpp::lvl10param;
using brP = TFHEpp::lvl01param;

void HalfAdder(TFHEpp::TLWE<TFHEpp::lvl1param>& sum,
               TFHEpp::TLWE<TFHEpp::lvl1param>& carry,
               const TFHEpp::TLWE<TFHEpp::lvl1param>& ca,
               const TFHEpp::TLWE<TFHEpp::lvl1param>& cb,
               const TFHEpp::EvalKey& ek)
{
    HomXOR<iksP, brP, brP::targetP::μ>(sum, ca, cb, ek);
    HomAND<iksP, brP, brP::targetP::μ>(carry, ca, cb, ek);
}

// elementary full comparator gate that is used to compare the i-th bit:
//   input: cai and cbi the i-th bit of a and b
//          cc: carry of lower bit
//   algo: sum and carry
void FullAdder(TFHEpp::TLWE<TFHEpp::lvl1param>& sum,
               TFHEpp::TLWE<TFHEpp::lvl1param>& carry,
               const TFHEpp::TLWE<TFHEpp::lvl1param>& cai,
               const TFHEpp::TLWE<TFHEpp::lvl1param>& cbi,
               const TFHEpp::TLWE<TFHEpp::lvl1param>& cci,
               const TFHEpp::EvalKey& ek)
{
    TFHEpp::TLWE<TFHEpp::lvl1param> hsum, hcarryab, hcarryabc;
    HalfAdder(hsum, hcarryab, cai, cbi, ek);
    HalfAdder(sum, hcarryabc, hsum, cci, ek);
    TFHEpp::HomOR<iksP, brP, brP::targetP::μ>(carry, hcarryabc, hcarryab, ek);
}

int main()
{
    TFHEpp::EvalKey ek;

    // reads the cloud key from file
    {
        const std::string path = "./cloud.key";
        std::ifstream ifs("./cloud.key", std::ios::binary);
        cereal::PortableBinaryInputArchive ar(ifs);
        ek.serialize(ar);
    }
    // import input
    std::vector<TFHEpp::TLWE<TFHEpp::lvl1param>> clientcipher;
    {
        std::ifstream ifs{"cloud.data", std::ios::binary};
        cereal::PortableBinaryInputArchive ar(ifs);
        ar(clientcipher);
    }

    // get cloutd input
    uint16_t cloud_input;
    std::cout << "Type cloud input (16bit unsigned interger)" << std::endl;
    std::cin >> cloud_input;

    // encrypt the input (trivial cypher)
    std::vector<TFHEpp::TLWE<TFHEpp::lvl1param>> cloudcipher(16);
    for (int i = 0; i < 16; i++) {
        if ((cloud_input >> i) & 1)
            TFHEpp::HomCONSTANTONE<TFHEpp::lvl1param>(cloudcipher[i]);
        else
            TFHEpp::HomCONSTANTZERO<TFHEpp::lvl1param>(cloudcipher[i]);
    }

    TFHEpp::TLWE<TFHEpp::lvl1param> flagequiv, flagclient;
    TFHEpp::HomCONSTANTONE<TFHEpp::lvl1param>(flagequiv);
    TFHEpp::HomCONSTANTONE<TFHEpp::lvl1param>(flagclient);
    for (int i = 0; i < 16; i++) {
        TFHEpp::TLWE<TFHEpp::lvl1param> tmpsum;
        cloudcipher[i][TFHEpp::lvl1param::k * TFHEpp::lvl1param::n] *=
            -1;  // NOT
        FullAdder(tmpsum, flagclient, clientcipher[i], cloudcipher[i],
                  flagclient, ek);
        TFHEpp::HomANDYN<iksP, brP, brP::targetP::μ>(flagequiv, flagequiv,
                                                     tmpsum, ek);
    }

    // export the result ciphertexts to a file
    std::vector<TFHEpp::TLWE<TFHEpp::lvl1param>> rescipher(2);
    rescipher[0] = flagequiv;
    rescipher[1] = flagclient;
    {
        std::ofstream ofs{"result.data", std::ios::binary};
        cereal::PortableBinaryOutputArchive ar(ofs);
        ar(rescipher);
    };
}