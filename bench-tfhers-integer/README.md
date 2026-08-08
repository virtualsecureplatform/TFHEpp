# tfhe-rs integer multiplication benchmark

This benchmark compares the TFHE2CLPX `64 x 64 -> 128` pipeline with a
single-threaded TFHE-rs integer multiplication having the same input and output
widths. TFHE-rs `FheUint64` multiplication truncates modulo `2^64`, so the
encrypted inputs are widened to `FheUint128` before multiplication. Both the
widening cost and the 128-bit multiplication cost are reported.

The dependency is pinned to TFHE-rs 1.6.1 (`99ded5fc5222b0ed24f4fc5420e0a04c56b5c88a`).
The Gaussian `2M64` two-message-bit/two-carry-bit parameters are selected to
match the paper's decryption-failure target more closely than TFHE-rs's `2M128`
default. Rayon is restricted to one worker thread.

Run on the same pinned physical core used for TFHEpp:

```sh
git clone --depth 1 --branch tfhe-rs-1.6.1 \
  https://github.com/zama-ai/tfhe-rs.git tfhe-rs
RAYON_NUM_THREADS=1 taskset -c 0 cargo run --release \
  --manifest-path bench-tfhers-integer/Cargo.toml
```

On an AMD Ryzen 9 7950X3D, five runs produced a mean of 34,085.921 ms
(47.187 ms sample standard deviation). The two widening casts account for
0.358 ms of that total. Key generation, encryption, warm-up, and decryption are
outside the timed region.
