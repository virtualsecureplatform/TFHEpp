# Getting Started

This page shows how to build TFHEpp and run a minimal Boolean‑gate example.

## Prerequisites

- **Compiler:** C++20 with UTF‑8 identifiers and `extern template` support.
  TFHEpp is primarily tested with GCC ≥ 11 and Clang ≥ 16.
- **CPU:** AVX2 is required for the default SPQLIOS FFT backend.
  If your platform lacks AVX2, use an alternative backend (see `Build.md`).
- **CMake:** ≥ 3.16
- **Submodules:** third‑party code lives under `thirdparties/` and is already
  vendored in this repo.

## Build

```bash
git clone <this repo>
cd TFHEpp
mkdir build && cd build
cmake .. -DENABLE_TEST=ON
make -j
```

This produces a static (default) or shared library target named `tfhe++`,
plus optional tests/tutorials depending on CMake flags.

## Minimal Example (Homomorphic NAND)

This is a compact version of `test/nand.cpp`.

```cpp
#include <iostream>
#include <tfhe++.hpp>

using namespace TFHEpp;

int main() {
    SecretKey sk;
    EvalKey ek;
    ek.emplacebkfft<lvl01param>(sk);  // bootstrapping key (lvl0 -> lvl1)
    ek.emplaceiksk<lvl10param>(sk);   // identity key switch (lvl1 -> lvl0)

    std::vector<uint8_t> pa = {1};
    std::vector<uint8_t> pb = {0};

    auto ca = bootsSymEncrypt(pa, sk);  // TLWE<lvl1param>
    auto cb = bootsSymEncrypt(pb, sk);

    TLWE<lvl1param> cres;
    HomNAND(cres, ca[0], cb[0], ek);

    auto pres = bootsSymDecrypt(std::vector<TLWE<lvl1param>>{cres}, sk);
    std::cout << "NAND result = " << int(pres[0]) << "\n";
}
```

### Build the example with CMake

In your own project:

```cmake
add_subdirectory(path/to/TFHEpp)
add_executable(my_example main.cpp)
target_link_libraries(my_example PRIVATE tfhe++)
```

Or compile against an installed TFHEpp:

```bash
g++ -std=c++20 main.cpp -ltfhe++ -I<install-prefix>/include -L<install-prefix>/lib
```

## Next Steps

- Read `API.md` for the core workflow and gate set.
- If you need different security/performance tradeoffs, see `Parameters.md`.
- For non‑Boolean workloads or advanced TFHE features, see `Advanced.md`.

