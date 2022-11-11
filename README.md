# oveus-tfhe
This repository contains the TFHE library for oveus based on [TFHEpp](https://github.com/virtualsecureplatform/TFHEpp).
Some performance improvements described in [this paper](https://doi.org/10.1145/3474366.3486927) and some UX improvements are applied.
The files in this repository are licensed under AGPLv3. For third parties' programs, see corresponding directories under /thirdparties.
If you need the commercial license for our AGPLv3 parts, please contact to oveus@axell.co.jp

**CAUTION:** TRADEMARKS INCLUDING THE NAME "oveus" ARE NOT LICENSED.

## Patents
This repository contains some codes implementing the contents of AXELL CORPORATION's patents.

Below is the README for TFHEpp.

# TFHEpp
TFHEpp is full Scracthed pure C++ Ver. of TFHE. TFHEpp is slightly(about 10%) faster than original [TFHE implementation](https://github.com/tfhe/tfhe). In addition to that, THFEpp supports Circuit Bootstrapping and [Private Boootstrapping many LUT](https://eprint.iacr.org/2021/729).
TFHEpp depends on AVX2 because we use SPQLIOS FMA. If you want run TFHEpp without AVX2, see spqlios++ branch. It include pure C++ implementation of SPQLIOS as header only library, but slow.

# Supported Compiler

This code includes utf-8 identifiers like α. Therefore, Clang and GCC10 or later are primarily supported compilers. GCC9 is not supported.

# Parameter
The default parameter is 128 bit security. Please add -DUSE_80BIT_SECURITY=ON to use faster but less secure parameter.

# FFTW3 Support
Some environments which do not support AVX2 cannot use spqlios. Instead of spqlios, TFHEpp can use fftw3.
To use fftw3,  install `libfftw3-dev` and add `-DUSE_FFTW3=ON` to the compile option.

# Third party libraries
Codes under thirdparties directory contain third-party libraries, Randen, Cereal, and SPQLIOS. See the corresponding directory to check the licenses.

## Randen
TFHEpp uses this as a Cryptographically Secure Pseudo-Random Number Generator (CSPRNG). Original repository is [here](https://github.com/google/randen).
I just removed some unnecessary codes, no modification.

## Cereal
cereal is a header-only C++11 serialization library. TFHEpp uses this to export ciphertexts and keys. Cereal is treated by git submodule.

## SPQLIOS
SPQLIOS is the FFT library using AVX2 that is dedicated to the ring R\[X\]/(X^N+1) for N a power of 2. These codes come from [experimental-tfhe](https://github.com/tfhe/experimental-tfhe/tree/master/circuit-bootstrapping/src/spqlios). We just renamed instances to adapt with our codes.

## SPWLIOS-AVX512
This is the AVX512 version of SPQLIOS developed in [MOSFHET](https://github.com/antoniocgj/MOSFHET). I confirmed that this is faster than SPQLIOS on Intel i5-11400.

## FFTW3
[FFTW](https://www.fftw.org/) is one of the most famous FFT libraries. 

**CAUTION**: REDISTRIBUTING BINARY WHICH DEPENDS ON FFTW3 MEANS YOU AGREED WITH GPLv3 OR LATER.

## MKL
Intel MKL is the library provided by Intel and including FFTW compatible interface for FFT.
We assume to install MKL by [this procedure](https://www.intel.com/content/www/us/en/developer/articles/guide/installing-free-libraries-and-python-apt-repo.html) and already ran `source /opt/intel/mkl/bin/mklvars.sh`.

Add `-DUSE_MKL` to the compile option to use MKL

## spqliox_aarch64
spqliox_aarch64 is the FFT library for aarch64 forked from SPQLIOS.
This is slightly faster than FFTW3(average 1ms).
This library requires [xbyak_aarch64](https://github.com/fujitsu/xbyak_aarch64), and
to use this library, add `-DUSE_SPQLIOX_AARCH64=on` to the compile option.

<center>

| FFTW3    | spqliox_aarch64 |
| :------: | :-------------: |
| 15.801ms | 14.368ms        |

</center>

# Speed Test

Following Code measure how many time homomorphic NAND takes on your computer with TFHEpp. 
```
git clone https://github.com/virtualsecureplatform/TFHEpp
cd TFHEpp
mkdir build
cd build
cmake .. -DENABLE_TEST=ON
make
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
