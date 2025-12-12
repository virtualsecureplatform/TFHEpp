# Serialization

TFHEpp uses the cereal library to serialize keys and ciphertexts. Both
`SecretKey` and `EvalKey` are cereal‑serializable.

## Saving keys

```cpp
#include <fstream>
#include <cereal/archives/portable_binary.hpp>
#include <tfhe++.hpp>

TFHEpp::SecretKey sk;
TFHEpp::EvalKey ek;
ek.emplacebkfft<TFHEpp::lvl01param>(sk);
ek.emplaceiksk<TFHEpp::lvl10param>(sk);

std::ofstream os("keys.bin", std::ios::binary);
cereal::PortableBinaryOutputArchive ar(os);
ar(sk, ek);
```

## Loading keys

```cpp
TFHEpp::SecretKey sk;
TFHEpp::EvalKey ek;

std::ifstream is("keys.bin", std::ios::binary);
cereal::PortableBinaryInputArchive ar(is);
ar(sk, ek);
```

Ciphertext types (`TLWE`, `TRLWE`, `TRGSW`) are plain `std::array`s and can be
archived directly as long as you include the appropriate cereal headers.

