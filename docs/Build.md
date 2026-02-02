# Build and Integration

TFHEpp uses CMake and builds a library target named `tfhe++`.

## Standard build

```bash
mkdir build && cd build
cmake .. <options>
make -j
```

### Useful CMake options

Security / parameters:

- `-DUSE_80BIT_SECURITY=ON` — legacy 80‑bit parameters (fast, not recommended).
- `-DUSE_CGGI19=ON` — CGGI19 parameter set.
- `-DUSE_CONCRETE=ON` — Zama CONCRETE parameter set (fast 128‑bit option).
- `-DUSE_TFHE_RS=ON` — TFHE‑rs‑style parameter set.
- `-DUSE_COMPRESS=ON` — compressed ciphertext parameters.
- `-DUSE_TERNARY=ON` — ternary secret keys.
- `-DUSE_SUBSET_KEY=ON` — subset/partial key switching support.
- `-DUSE_TERNARY_CMUX=ON` — ternary CMUX (implies `USE_TERNARY`).
- `-DUSE_KEY_BUNDLE=ON` — key‑bundle algorithm for bootstrapping.

FFT / backends:

- default: SPQLIOS (AVX2) from `thirdparties/spqlios`
- `-DUSE_AVX512=ON` — AVX512 SPQLIOS variant.
- `-DUSE_FFTW3=ON` — FFTW3 backend (requires system `libfftw3-dev`).
  **Caution:** linking FFTW3 implies GPLv3+ for your binary.
- `-DUSE_MKL=ON` — Intel MKL backend (requires MKL installed and `MKLROOT`).
- `-DUSE_CONCRETE_FFT=ON` — Rust concrete‑fft backend.
- `-DUSE_SPQLIOX_AARCH64=ON` — AArch64 backend (needs `xbyak_aarch64`).
- `-DUSE_SPQLIOS_ARITHMETIC=ON` — spqlios‑arithmetic backend.
- `-DUSE_HEXL=ON` — Intel HEXL NTT backend (AVX512).

Build products:

- `-DENABLE_TEST=ON` — build tests under `test/`.
- `-DENABLE_BENCHMARK=ON` — build benchmarks under `benchmark/`.
- `-DENABLE_TUTORIAL=ON` — build tutorial programs under `tutorial/`.
- `-DENABLE_SHARED=ON` — build `tfhe++` as a shared library.
- `-DDEBUG=ON` — disable `-O3`, enable debug flags.

Toolchain / CPU flags:

- `-DUSE_MARCH_NATIVE=ON` (default) — enable `-march=native` for CPU‑tuned builds.
- `-DUSE_MARCH_NATIVE=OFF` — build portable binaries (useful in CI or mixed CPUs).

Randomness:

- `-DUSE_BLAKE3=ON` (default) — BLAKE3‑XOF CSPRNG.
- `-DUSE_RANDEN=ON` — legacy Randen CSPRNG (deprecated).

## Integration patterns

### Add TFHEpp as a subproject

```cmake
add_subdirectory(external/TFHEpp)
target_link_libraries(my_target PRIVATE tfhe++)
```

TFHEpp exports include paths and links any selected backend automatically.

### Install and use via find_package

```bash
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/prefix
make install
```

Then:

```cmake
find_library(TFHEPP_LIB tfhe++ REQUIRED)
target_link_libraries(my_target PRIVATE ${TFHEPP_LIB})
target_include_directories(my_target PRIVATE /path/to/prefix/include)
```

## Platform notes

- By default TFHEpp assumes AVX2. If you need portability, select FFTW3 or the
  AArch64 backend.
- OpenMP is enabled when available; it accelerates batch encryption/decryption
  and some bootstrapping routines.
