# Types in TFHEpp
On this page, some types used in TFHEpp are introduced to describe the interfaces of TFHEpp.
Most of the types are defined in [params.hpp](../include/params.hpp).

## Philosophy
Most of TFHEpp's types are defined by `std::array`. This is because I hope to let the compiler detect type mismatch for homomorphic operations and do the compile time optimizations like loop unrolling.
This also makes it easier to export TFHEpp's data to other libraries and vice versa.
Carrying some additional information like current noise can be useful for some applications, but I decided to keep the data as simple as possible to reduce the cost of understanding whole code structures.

## Ciphertext Types

There are three ciphertext types in TFHE: TLWE(Torus Learning With Errors), TRLWE(Torus Ring LWE), and TRGSW(Torus Ring Gentry-Sahai-Walters).

### TLWE

This is the most fundamental type in TFHE. In the actual implementation, Torus is discretized into integers, so TLWE is Integer LWE with a power of two modoli. 
In TFHEpp, TLWE is represented by an array of unsigned integers with length n+1, like `std::array<uint32_t,n+1>`.
This TLWE ciphertext is representing ($\mathbf{a},b$). 

### TRLWE

This is the ring version of TLWE. Though I use the term TRLWE in this library because of historical reasons, this is a TMRLWE(Torus Module Ring LWE) ciphertext which means the dimension of the ciphertext can be larger than two.
The form of TRLWE is $(\mathbf{a}[X],b[X])$, so the $k+1$-dimensional vector of $N$ degree Torus polynomials. 
This is represented by the type like `std::array<std::array<uint32_t,N>k+1>`.

### TRGSW

This is the ring version of the Gentry-Sahai-Walters homomorphic encryption scheme.