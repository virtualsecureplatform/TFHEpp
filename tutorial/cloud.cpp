#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/vector.hpp>
#include <fstream>
#include <memory>
#include <random>
#include <tfhe++.hpp>
#include <vector>

void HalfAdder(TFHEpp::TLWElvl0& sum, TFHEpp::TLWElvl0& carry,
               const TFHEpp::TLWElvl0& ca, const TFHEpp::TLWElvl0& cb,
               const TFHEpp::GateKey& gk)
{
    HomXOR(sum, ca, cb, gk);
    HomAND(carry, ca, cb, gk);
}

// elementary full comparator gate that is used to compare the i-th bit:
//   input: cai and cbi the i-th bit of a and b
//          cc: carry of lower bit
//   algo: sum and carry
void FullAdder(TFHEpp::TLWElvl0& sum, TFHEpp::TLWElvl0& carry,
               const TFHEpp::TLWElvl0& cai, const TFHEpp::TLWElvl0& cbi,
               const TFHEpp::TLWElvl0& cci, const TFHEpp::GateKey& gk)
{
    TFHEpp::TLWElvl0 hsum, hcarryab, hcarryabc;
    HalfAdder(hsum, hcarryab, cai, cbi, gk);
    HalfAdder(sum, hcarryabc, hsum, cci, gk);
    TFHEpp::HomOR(carry, hcarryabc, hcarryab, gk);
}

int main()
{
    // To avoid stack limitation, GateKey should be allocated like this
    std::unique_ptr<TFHEpp::GateKey> gk(new TFHEpp::GateKey());

    // reads the cloud key from file
    {
        const std::string path = "./cloud.key";
        std::ifstream ifs("./cloud.key", std::ios::binary);
        cereal::PortableBinaryInputArchive ar(ifs);
        gk->serialize(ar);
    }
    // import input
    std::vector<TFHEpp::TLWElvl0> clientcipher;
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
    std::vector<TFHEpp::TLWElvl0> cloudcipher(16);
    for (int i; i < 16; i++) {
        if ((cloud_input >> i) & 1)
            TFHEpp::HomCONSTANTONE(cloudcipher[i]);
        else
            TFHEpp::HomCONSTANTZERO(cloudcipher[i]);
    }

    TFHEpp::TLWElvl0 flagequiv, flagcloud;
    TFHEpp::HomCONSTANTONE(flagequiv);
    TFHEpp::HomCONSTANTONE(flagcloud);
    for (int i; i < 16; i++) {
        TFHEpp::TLWElvl0 tmpsum;
        cloudcipher[i][TFHEpp::DEF_n] *= -1;  // NOT
        FullAdder(tmpsum, flagcloud, clientcipher[i], cloudcipher[i], flagcloud,
                  *gk);
        TFHEpp::HomANDYN(flagequiv, flagequiv, tmpsum, *gk);
    }

    // export the result ciphertexts to a file
    std::vector<TFHEpp::TLWElvl0> rescipher(2);
    rescipher[0] = flagequiv;
    rescipher[1] = flagcloud;
    {
        std::ofstream ofs{"result.data", std::ios::binary};
        cereal::PortableBinaryOutputArchive ar(ofs);
        ar(rescipher);
    };
}