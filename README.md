[![Test](https://github.com/virtualsecureplatform/TFHEpp/actions/workflows/test.yml/badge.svg)](https://github.com/virtualsecureplatform/TFHEpp/actions/workflows/test.yml)
# TFHEpp
TFHEpp is a full scratched pure C++ version of TFHE. TFHEpp is slightly(about 10%) faster than the original [TFHE implementation](https://github.com/tfhe/tfhe). In addition to that, THFEpp supports [Circuit Bootstrapping](https://eprint.iacr.org/2018/421), [Programable Boootstrapping many LUT](https://eprint.iacr.org/2021/729), Partial Key([Klemsa's version](https://eprint.iacr.org/2021/634), [Zama's version](https://eprint.iacr.org/2023/979) . Named as subset key in the code), [Modified Chen's Packing](https://eprint.iacr.org/2024/1318) (We call it as annihilate key switching in our code), [mean compensation in rounding noises](https://eprint.iacr.org/2025/809), and [Double Decomposition](https://eprint.iacr.org/2023/771) (bivariate polynomial representation) with 128-bit Torus support. 
We also include partial support for B/FV, written in include/bfv++.hpp. For the implementation of the Modified Chen's Packing, we also used the idea of [dividing TRLWE at each recursion](https://eprint.iacr.org/2025/1088), though I developed that independently. 
TFHEpp depends on AVX2 because we use SPQLIOS FMA. If you want to run TFHEpp without AVX2, see the spqlios++ branch. It includes a pure C++ implementation of SPQLIOS as a header-only library, but it is slow. 

# Supported Compiler

This code includes UTF-8 identifiers like α, using `extern template`, and std::make_unique_for_overwrite. Therefore, GCC11 or later and Clang16 or later are primarily supported compilers. 

# Parameter
The default parameter (128bit.hpp) is 128-bit security. Please add -DUSE_80BIT_SECURITY=ON to use a faster but less secure parameter.
The fastest parameter with 128-bit security is concrete.hpp. Please add -DUSE_CONCRETE=ON to use a faster parameter. 
If you need your own parameter, please modify the files under include/params/ or add your own hpp and add change the macro definition in include/param.hpp. 

## Ternary Key
As you will see in the files under include/params, we use a ternary key for lvl1 or higher-level parameters. 

# FFTW3 Support
Some environments that do not support AVX2 cannot use spqlios. Instead of spqlios, TFHEpp can use fftw3.
To use fftw3,  install `libfftw3-dev` and add `-DUSE_FFTW3=ON` to the compile option.

# Third-party libraries
Codes under thirdparties directory contain third-party libraries, Randen, BLAKE3, Cereal, and SPQLIOS. See the corresponding directory to check the licenses.

## Randen
Previous (till Version 9 in release tags) TFHEpp uses this as a Cryptographically Secure Pseudo-Random Number Generator (CSPRNG). Original repository is [here](https://github.com/google/randen).
I just removed some unnecessary codes, with no modification. 
This is now deprecated. To use this, set `-DUSE_BLAKE3=OFF -DUSE_RANDEN=ON` explicitly. 

## BLAKE3
This is used as another CSPRNG implementation using its eXtended Output Function (XOF) mode. Because Randen is not peer-reviewed algorithm and BLAKE3 is now bit faster than RANDEN, we are now using this as a default CSPRNG. 
FYI, CRYSTAL-Kyber is using SHA-3's XOF as a CSPRNG generator. We are currently not implemening this because it is slow. 
Micorosft SEAL implements BLAKE2 as one of the supported CSPRNGs. 

## Cereal
cereal is a header-only C++11 serialization library. TFHEpp uses this to export ciphertexts and keys. Cereal is treated by the git submodule.

## SPQLIOS
SPQLIOS is the FFT library using AVX2 dedicated to the ring R\[X\]/(X^N+1) for N a power of 2. These codes come from [experimental-tfhe](https://github.com/tfhe/experimental-tfhe/tree/master/circuit-bootstrapping/src/spqlios). We just renamed instances to adapt to our codes.

## SPQLIOS-AVX512
This AVX512 version of SPQLIOS developed in [MOSFHET](https://github.com/antoniocgj/MOSFHET). I confirmed this is faster than SPQLIOS on Intel i5-11400 and AMD Ryzen 7700X. 

## SPQlios-Arithmetic
[SPQlios Arithmetic](https://github.com/tfhe/spqlios-arithmetic) seems to be the successor of SPQLIOS. Currently, only FFT is supported. 

## FFTW3
[FFTW](https://www.fftw.org/) is one of the most famous FFT libraries. 

**CAUTION**: REDISTRIBUTING BINARY WHICH DEPENDS ON FFTW3 MEANS YOU AGREED WITH GPLv3 OR LATER.

## MKL
Intel MKL is the library provided by Intel and including FFTW compatible interface for FFT.
We assume to install MKL by [this procedure](https://www.intel.com/content/www/us/en/developer/articles/guide/installing-free-libraries-and-python-apt-repo.html) and already ran `source /opt/intel/mkl/bin/mklvars.sh`.

Add `-DUSE_MKL` to the CMake option to use MKL

### FFTW3 API
To use this, you have to also add `-DUSE_FFTW3`.

### Native API
Instead of FFTW3 API, I also added native API version. This will be enabled if `-DUSE_FFTW3` is not specified with `-DUSE_MKL`.

## concrete-fft
concrete-fft is the pure Rust FFT library developed by Zama.ai. This can be enabled by `-DUSE_CONCRETE_FFT`.

## spqliox_aarch64
spqliox_aarch64 is the FFT library for aarch64 forked from SPQLIOS.
This is slightly faster than FFTW3(average 1ms).
This library requires [xbyak_aarch64](https://github.com/fujitsu/xbyak_aarch64), and
to use this library, add `-DUSE_SPQLIOX_AARCH64=on` to the CMake option.

<center>

| FFTW3    | spqliox_aarch64 |
| :------: | :-------------: |
| 15.801ms | 14.368ms        |

</center>

## HEXL

[HEXL](https://github.com/intel/hexl.git) is the NTT library optimized for AVX512. 
To use this library, add `-DUSE_HEXL=on` to the CMake option.

# Speed Test

Following Code measure how many time homomorphic NAND takes on your computer with TFHEpp. 
```
git clone https://github.com/virtualsecureplatform/TFHEpp
cd TFHEpp
mkdir build
cd build
cmake .. -DENABLE_TEST=ON
make
ulimit -s unlimited
./test/nand 
```

If you want to run semantically equivalent test on original TFHE, run below code.
```
git clone https://github.com/tfhe/tfhe.git --recursive
cd tfhe
mkdir build
cd build
cmake ../src -DENABLE_TESTS=on -DENABLE_NAYUKI_PORTABLE=off -DENABLE_NAYUKI_AVX=off -DENABLE_SPQLIOS_AVX=off -DENABLE_SPQLIOS_FMA=on -DCMAKE_BUILD_TYPE=optim
make
./test/test-gate-bootstrapping-spqlios-fma
```

If you have Docker on your system, this will do above on docker.

```
git clone https://github.com/virtualsecureplatform/TFHEpp
cd TFHEpp
docker build -t tfheppbench .
```

This is for TFHE-10ms. Because TFHE-10ms only supports 80-bit security parameter, this is not included in Dockerfile
```
git clone https://github.com/virtualsecureplatform/tfhe-10ms.git --recursive
cd tfhe-10ms
mkdir build
cd build
cmake ../src -DENABLE_TESTS=on -DENABLE_NAYUKI_PORTABLE=off -DENABLE_NAYUKI_AVX=off -DENABLE_SPQLIOS_AVX=off -DENABLE_SPQLIOS_FMA=on -DCMAKE_BUILD_TYPE=optim
make
./test/test-bootstrapping-fft-spqlios-fma 
```

# Theory

Here is the slides (in japanese).
https://nindanaoto.github.io/

# Citation

For the people who want to cite this library directly (may be in addition to [VSP paper](https://www.usenix.org/conference/usenixsecurity21/presentation/matsuoka)), I give a below example of bibtex citation.

```
@misc{TFHEpp,
	author = {Kotaro Matsuoka},
	title = {{TFHEpp: pure C++ implementation of TFHE cryptosystem}},
  	year = {2020},
	howpublished = {\url{https://github.com/virtualsecureplatform/TFHEpp}}
}
```
