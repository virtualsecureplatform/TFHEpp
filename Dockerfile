FROM ubuntu:20.04

LABEL maintainer="nindanaoto <dirty.knife603@gmail.com>"

RUN apt-get update && apt-get upgrade -y && apt-get install -y build-essential clang libomp-dev cmake git && git clone --recursive --depth 1 https://github.com/virtualsecureplatform/TFHEpp && git clone --recursive --depth 1 https://github.com/tfhe/tfhe.git && git clone --recursive --depth 1 https://github.com/virtualsecureplatform/tfhe-10ms.git && mkdir TFHEpp/build && mkdir tfhe/build && mkdir tfhe-10ms/build

WORKDIR TFHEpp/build

RUN cmake ..&& make

WORKDIR /tfhe/build

RUN cmake ../src -DENABLE_TESTS=on -DENABLE_NAYUKI_PORTABLE=off -DENABLE_NAYUKI_AVX=off -DENABLE_SPQLIOS_AVX=off -DENABLE_SPQLIOS_FMA=on -DCMAKE_BUILD_TYPE=optim && make

# WORKDIR /tfhe-10ms/build

# RUN cmake ../src -DENABLE_TESTS=on -DENABLE_NAYUKI_PORTABLE=off -DENABLE_NAYUKI_AVX=off -DENABLE_SPQLIOS_AVX=off -DENABLE_SPQLIOS_FMA=on -DCMAKE_BUILD_TYPE=optim && make

WORKDIR /

RUN echo "TFHEpp" &&./TFHEpp/build/test/nand && echo "THFE" && ./tfhe/build/test/test-gate-bootstrapping-spqlios-fma && echo "TFHE-10ms" && tfhe-10ms/build/test/test-bootstrapping-fft-spqlios-fma 